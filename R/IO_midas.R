
#' Read a MIDAS formated file into a tidy data frame
#'
#' A MIDAS (Minimum Information for Data Analysis in Systems Biology) is a
#' specially formated csv file that almost nobody uses. However, we use it some
#' time particularly when there is a need to use \emph{DataRail} or \emph{CellNOpt}.
#'
#' @param MIDASfile A path to a correctly formated MIDAS file.
#'
#' @details For more information about the MIDAS format you can read
#'     \href{http://www.cellnopt.org/doc/cnodocs/midas.html}{here}.
#'
#'     Up to this point, I have never seen a MIDAS file with prefix
#'     other than \emph{TR, DA, DV}, so any other prefix will be removed with a warning.
#'
#' @return A tidy data frame (tibble to be exact) where all the "DA" fields have
#'     been "gathered" into a "Time" variable and all the "DV" field into a "Value" variable.
#'     Column prefices (TR, DA, DV) have been removed.
#'
#' @author Teo Sakel
#' @export
#'
#' @examples
read_MIDAS <- function(MIDASfile) {
    #  based on https://github.com/Bioconductor-mirror/CellNOptR/blob/master/R/readMIDAS.R

    ## A few checks not to fill like engineers!
    if (!file.exists(MIDASfile)) stop(paste("No file:", MIDASfile))
    ncols = count.fields(MIDASfile, sep = ",")
    if (ncols[1] < 3) stop(paste("Too few columns in", MIDASfile))
    if (any(ncols[1] != ncols)) stop(paste("Cannot read:", MIDASfile, "rows have unequal length"))

    midas <- read.csv(MIDASfile, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")

    # Check column names: Prefix
    prefix <- grepl("^TR:|^DA:|^DV:", colnames(midas))
    if (!all(prefix)) {
        warn_msg <- paste("Unknown column prefices:",
                          paste(colnames(midas)[!prefix], collapse = ","),
                          "They will be removed.")
        warning(warn_msg)
        midas <- midas[prefix]
    }

    # Check column names: Fields
    colon_count <- stringr::str_count(colnames(midas), ":")
    if (any(colon_count) > 2 || any(colon_count) < 1) {
        stop("Ill-formated header. Only options are: \"XX:Specy and XX:userword:Specy\"")
    }

    TRcol <- grep(pattern = "^TR:", colnames(midas), value = TRUE)
    DAcol <- grep(pattern = "^DA:", colnames(midas), value = TRUE)
    DVcol <- grep(pattern = "^DV:", colnames(midas), value = TRUE)

    # The MIDAS may contain the DA:ALL special time
    if (length(DAcol) == 1 && DAcol == "DA:ALL") {
        # Easy case - Just tidy analyte values
        colnames(midas) <- gsub("^DA:ALL$", "Time", colnames(midas))
        midas <- tidyr::gather_(midas, "Analyte", "Value", DVcol)
        midas$Analyte <- gsub("^DV:", "", midas$Analyte)
        col_order <- c(TRcol, "Time", "Analyte", "Value")
        midas <- midas[col_order]
    } else {
        # Check if Analytes match between DA and DV
        analytes_dv <- vapply(stringr::str_split(DVcol, ":"), "[[", "", 2)
        analytes_da <- vapply(stringr::str_split(DAcol, ":"), "[[", "", 2)
        if (!setequal(analytes_dv, analytes_da)) {
            mismatch <- setdiff(union(analytes_dv, analytes_da), intersect(analytes_dv, analytes_da))
            stop(paste0("There is a mismatch in analytes between DA and DV.\n",
                        "The following analytes don't match:\n",
                        paste(mismatch, collapse = ", ")))
        }

        # Tidy!
        pair_values <- function(x) {
            data.frame(Analyte = x,
                       Time = midas[[paste0("DA:", x)]],
                       Value = midas[[paste0("DV:", x)]], stringsAsFactors = FALSE)
        }
        DAV <- dplyr::bind_rows(lapply(analytes_da, pair_values))
        midas <- dplyr::bind_cols(midas[rep(1:nrow(midas), length(analytes_da)), TRcol], DAV)
    }

    # Some final cleaning
    if (any(stringr::str_count(midas[["Analyte"]], ":") > 0)) {
        midas <- tidyr::separate_(midas, "Analyte", c("Analyte", "Group"), sep = ":")
    }
    colnames(midas) <- gsub("^TR:", "", colnames(midas))
    colnames(midas) <- sub(":", ".", colnames(midas))
    midas <- midas[!is.na(midas$Value), , drop = FALSE]

    return(tibble::as_tibble(midas))
}

#' Write a tidy data frame to a MIDAS-formated csv file
#'
#' A wrapper function, it calls \code{make_MIDAS} using all the arguments but \code{filename}
#' and then uses \code{write.csv} to write the result into \code{filename}.
#'
#' @param df a tidy data frame
#' @param filename a path naming a file or connection (it will be passed to write.csv)
#' @param ... arguments passed to \code{make_MIDAS}
#'
#' @details Column-field delimiters are substituted from dots to semicolons to confront with the original protocol.
#'
#' @seealso The \code{\link{make_MIDAS}} for details about the conversion.
#' @export
#'
#' @examples
write_MIDAS <- function(df, filename, ...) {
    df <- make_MIDAS(df, ...)
    colnames(df) <- sub(".", ":", colnames(df), fixed = TRUE)
    write.csv(df, file = filename, row.names = FALSE)
}

#' Convert a tidy data frame into MIDAS-formated data frame
#'
#' A MIDAS (Minimum Information for Data Analysis in Systems Biology) is a
#' specially formated csv file that almost nobody uses. However, we use it some
#' time particularly when there is a need to use \emph{DataRail} or \emph{CellNOpt}.
#'
#' @param df a tidy data frame to be converted
#' @param analyte_col the column of analytes,
#'   it will be used to name the DA- and DV-columns
#' @param time_col the column that encodes the time,
#'   it will be "spread" into the DA-columns
#' @param value_col the column that encodes the expression values,
#'   it will be spread into the DV-columns
#' @param well_id a column specifying which rows belong to the same treatment.
#'   To avoid conflicts when "spreading". If \code{NULL} each group of treatments
#'   is assumed be unique (no replicates)
#'
#' @return A MIDAS formated data frame. The columns of this data frame are grouped
#'     into 3 categories and each category is tagged with a two letter prefix:
#'     \itemize{
#'         \item \code{TR} for treatments: describe the condition of the experiment
#'         \item \code{DA} for data acquisition: describe the timing of measurements
#'         \item \code{DV} for data values: describe the results of the expreriment
#'     }
#'     DA- and DV-columns are suffixed by the analyte they refer to.
#'     TR-columns may have an extra field specifying the category of the treatment (eg stimulus, inhibitor)
#'
#'     Unlike the original MIDAS the different fields are delimited by a dot "\code{.}"
#'     and not semicolons "\code{:}" as described by the protocol to be more \code{R} friendly.
#'
#'     All the columns not-specified are considered TR-columns. They must be either "binary"
#'     columns or be factor-coersible otherewise an error is thrown.
#'     Factor-like columns will be expanded with the following format: \code{TR.<colname>.<value>}.
#'
#' @seealso For more details about the MIDAS format you can read
#'     \href{http://www.cellnopt.org/doc/cnodocs/midas.html}{here}
#' @export
#'
#' @examples
make_MIDAS <- function(df, analyte_col="Analyte", time_col="Time",
                       value_col="Value", well_id=NULL) {

    TRcol <- setdiff(colnames(df), c(time_col, analyte_col, value_col, well_id))
    analytes <- unique(df[[analyte_col]])
    if (is.null(well_id)) {
        # assume that everything is unique
        df['__uid__'] <- dplyr::group_indices_(df, TRcol)  # pythonic convention ;P
        well_id <- '__uid__'
    }

    # First fix the treatments
    # condense = not binary = have multiple values
    is_condensed <- function(x) {
        dfx <- df[[x]]
        !(is.numeric(dfx) && all(unique(dfx) %in% c(0, 1)))
    }
    uncondense <- function(x) {
        xval <- df[[x]]
        if (is.numeric(xval)) {
            stop(paste0("Only factor-like columns can be used as treatment. ", x, " is numeric"))
        } else {
            xval <- paste0(x, ".", xval)  # results in TR.x.x_val
        }
        as.data.frame(
            vapply(unique(xval), "==", rep(FALSE, length(xval)), xval) + 0
        )
    }
    TRcondensed <- TRcol[vapply(TRcol, is_condensed, FALSE)]
    if (length(TRcondensed) > 0) {
        df_uncondensed <- dplyr::bind_cols(lapply(TRcondensed, uncondense))
        df <- dplyr::bind_cols(df_uncondensed, df[setdiff(colnames(df), TRcondensed)])
        TRcol <- setdiff(colnames(df), c(time_col, analyte_col, value_col, well_id))
    }

    # Fix Time & Values (there is probably a better way...)
    df <- tidyr::unite_(df, "__pairs__", c(time_col, value_col), sep = "_")
    df <- tidyr::spread_(df, analyte_col, "__pairs__", fill = NA)
    df[well_id] <- NULL
    for (analyte in analytes) {
        df <- tidyr::separate_(df, analyte, paste0(c("DA.", "DV."), analyte),
                               sep = "_", convert = TRUE)
    }

    colnames(df)[match(TRcol, colnames(df))] <- paste0("TR.", TRcol)
    corder <- unlist(lapply(c("^TR.", "^DA.", "^DV."), grep, x = colnames(df)))
    df[corder]
}
