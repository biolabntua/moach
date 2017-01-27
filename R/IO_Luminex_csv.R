
#' Read a Luminex-Xponent generated csv file
#'
#' @param path The path to a Luminx csv file
#'
#' @return A list of lists with all the information from the file.
#'     At the top level there are two lists
#'     \itemize{
#'         \item \code{BatchHeader} that stores all the information before the \emph{Results}
#'         section as well as the CRC and
#'         \item \code{AssayData} that stores all the information in the \emph{Results} section
#'     }
#'
#'     \code{BatchHeader} stores data as key-value pairs with keys taken from the first column of the csv file and
#'     values from the rest. Values are mostly character-variables and there are few character vectors.
#'     The only exception is \code{CALInfo} with is a list of dataframes containing information about:
#'     \itemize{
#'         \item Last Calibration in \code{LastCalibration}
#'         \item Classification Calibrator in \code{ClassificationCalibratorExtended}
#'         \item Reporter Calibratorin (you guessed it!) \code{ReporterCalibrator}
#'     }
#'
#'     \code{AssayData} stores data as dataframes, except for \code{Run} data which is a list of dataframes
#'     (Audit Logs and Warning/Errors). There are 4 dataframes in \code{AssayData}:
#'     \itemize{
#'         \item \strong{Exprs} Expression data, ie \emph{raw} data, including \strong{Medians, Counts, Net MFI}, etc
#'         \item \strong{Average} Averages for replicate wells
#'         \item \strong{Well} Information about each well (currently only the \emph{Dilution Factor})
#'         \item \strong{Bead} Information about the beads, eg BeadID
#'     }
#'     more information about the different data generated can be found
#'     \href{https://www.luminexcorp.com/blog/2012/06/19/its-all-about-the-stats/}{here}.
#'
#' @details The format details (eg comma or semicolon delimiter) are infered by the function.
#'     There is a known edge-case, for some semicolon-separated files,
#'     where some calibration information is misplaced. This is the **only** case where you should
#'     mess with the csv file, you should **not mess** with it.
#'
#' @author Teo Sakel
#' @export
#'
#' @examples
#' ## Read a csv file for a 384 plate
#' path <- system.file("extdata", "test0_384.csv", package="moach")
#' lumcsv <- read_Xponent_csv(path)
read_Xponent_csv <- function(path) {

    # Read and Clean Lines
    Lines <- readLines(path)
    delim <- ifelse(stringr::str_count(Lines[1], ",") > 0, ",", ";")
    delim_regex <- paste0(delim, "(?=([^\"]*\"[^\"]*\")*[^\"]*$)")  # csv regex
    Lines <- stringr::str_replace(Lines, paste0(delim, "*$"), "") # Remove trailing delimiters
    Lines <- stringr::str_split(Lines, delim_regex)
    Lines <- lapply(Lines, stringr::str_replace_all, "^\"|\"$", "")
    Lines <- vapply(Lines, paste, "", collapse = "\t")

    # Organize Lines into Blocks
    empty_lines <- which(Lines == "") + 1
    Lines <- .splitAt(Lines, empty_lines)
    Lines[vapply(Lines, function(x) all(x == ""), FALSE)] <- NULL  # remove empty blocks
    header_end <- grep("^Results", vapply(Lines, "[", "", 1))
    header_lines <- Lines[1:(header_end - 1)]
    results_lines <- Lines[(header_end + 1):length(Lines)]


    Xcsv <- list(
        BatchHeader = .parse_Xcsv_header(header_lines),
        AssayData = .parse_Xcsv_results(results_lines)
    )
    # special case
    Xcsv$BatchHeader$CRC <- as.character(tryCatch(results_lines[[length(results_lines)]][2],
                                                  error = function(e) NA))

    class(Xcsv) <- "XponentCSV"
    return(Xcsv)
}


# SubMethods ----------------------------------------------------------------------------------

.parse_Xcsv_header <- function(Lines) {
    ## Comments:
        # Should we pass the Dates in the header?

    # Parameters
    delim <- "\t"
    fields <- vapply(Lines, "[", "", 1)

    # Parse Batch Info
    batch_info <- c(Lines[[1]], Lines[[2]])  # first two fields are Batch Info
    batch_info <- .parse_Xcsv_BatchInfo(batch_info)

    # Parse Calibration Info
    lastcal <- grep("^Most Recent Calibration and Verification Results:", fields)
    calinfo <- grep("^CALInfo:", fields)
    CALInfo <- .parse_Xcsv_CALInfo(list(unlist(Lines[lastcal]), unlist(Lines[calinfo])))

    # Parse Sample Field
    sample_field <- Lines[[length(Lines)]]
    sample_field <- stringr::str_split(sample_field[grep("^Samples", sample_field)], "\t")[[1]]
    Samples <- as.integer(sample_field[2])
    MinEvents <- paste(sample_field[4:length(sample_field)], collapse = " ")

    return(list(
        BatchInfo = batch_info,
        CALInfo = CALInfo,
        Samples = Samples,
        MinEvents = MinEvents
    ))
}

.parse_Xcsv_BatchInfo <- function(Lines) {

    # Parameters
    list_fields <- c("ProtocolPlate", "ProtocolMicrosphere")
    parse_list_fields <- function(x) {
        y <- setNames(x[seq(2, length(x), 2)], x[seq(1, length(x), 2)])
        y[y == ""] <- NA
        y
    }
    parse_simple_fields <- function(x) {
        as.character(ifelse(any(stringr::str_detect(x, "\\S+")), paste(x, collapse = " "),  NA))
    }

    # Clean up
    Lines <- stringr::str_split(Lines, "\t")  # split into list of vectors
    fields <- vapply(Lines, "[", "", 1)
    fields <- stringr::str_replace_all(fields, "\\s+", "")
    Lines <- Lines[fields != ""]; fields <- fields[fields != ""]  # remove empty lines
    simple_fields <- setdiff(fields, list_fields)
    # Final Form!
    BatchInfo <- setNames(lapply(Lines, function(x) x[2:max(2, length(x))]), fields)
    BatchInfo[simple_fields] <- lapply(BatchInfo[simple_fields], parse_simple_fields)
    BatchInfo[list_fields] <- lapply(BatchInfo[list_fields], parse_list_fields)

    return(BatchInfo)
}

.parse_Xcsv_CALInfo <- function(Lines) {
    # I am not sure about the format. Added a lot of tryCatch to avoid trouble...

    LastCAL <- NA
    if (!is.null(Lines[[1]])) {
        LastCAL <- .chr_vector_to_df(Lines[[1]][-1], sep = "\t", header = FALSE,
                                     col.names = c("Description", "Results"),
                                     colClasses = "character")
    }

    CALInfo <- NA
    if (!is.null(Lines[[2]])) {
        CALInfo <- Lines[[2]][-1]
        cce <- which(CALInfo == "Classification Calibrator Extended")
        rcal <- which(CALInfo == "Reporter Calibrator")
        Classification_Calibrator_Extended <- tryCatch(
            .chr_vector_to_df(CALInfo[(cce+1):(rcal -1)], sep = "\t", header=TRUE),
            error = function(e) NA
        )

        Reporter_Calibrator <- tryCatch({
            rcal <- CALInfo[(rcal+1):length(CALInfo)]
            lots <- grep("^Lot", rcal)
            rcal <- .splitAt(rcal, lots)
            rcal <- lapply(rcal, .chr_vector_to_df, sep = "\t", header=TRUE)
            dplyr::bind_rows(rcal)
        }, error = function(e) NA)

    }

    return(list(
        LastCalibration = LastCAL,
        ClassificationCalibratorExtended = Classification_Calibrator_Extended,
        ReporterCalibrator = Reporter_Calibrator
    ))

}

.parse_Xcsv_results <- function(Lines) {
    ## Comments:
        # Break into smaller functions?
        # Export as tibble?

    # Parameters
    id_vars <- c("Location", "Sample", "Analyte")
    exprs_var <- c("Median", "Count", "Net MFI", "Result", "Range", "% Recovery", "Comments")
    aver_vars <- c("Avg Net MFI", "Avg Result", "Avg Range",
                   "%CV Replicates", "Standard Expected Concentration")
    well_vars <- c("Dilution Factor")
    bead_vars <- c("Units", "Per Bead Count", "Analysis Types")
    run_vars <- c("Audit Logs", "Warnings/Errors")

    # Prepare Data
    fields <- vapply(Lines, "[", "", 1)
    fields <- stringr::str_split(fields, "\t")
    fields <- vapply(fields, "[", "", 2)
    names(Lines) <- fields

    # Parse Well Variables
    exprs_var <- intersect(exprs_var, fields)
    Exprs <- lapply(exprs_var,
        function(type) {
            df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t")
            df$Total.Events <- NULL
            df$Location <- stringr::str_extract(df$Location, "[A-Z]+[0-9]+")
            analytes <- setdiff(colnames(df), id_vars[-3])
            df <- tidyr::gather_(df, id_vars[3], type, analytes)
            colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "_")
            df
        }
    )
    Exprs <- Reduce(
        function(...) merge.data.frame(by = id_vars, sort = FALSE, ...),
        Exprs
    )

    # Parse Average Data
    aver_vars <- intersect(aver_vars, fields)
    Average <- NA
    if (length(aver_vars) > 0) {
        Average <- lapply(aver_vars,
            function(type) {
                df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t")
                analytes <- setdiff(colnames(df), id_vars[-3])
                df <- tidyr::gather_(df, id_vars[3], type, analytes)
                colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "_")
                df
            }
        )
        Average <- Reduce(
            function(...) merge.data.frame(by = id_vars[-1], sort = FALSE, ...),
            Average
        )
    }

    # Parse Sample Data
    well_vars <- intersect(well_vars, fields)
    Wells <- NA
    if (length(well_vars) > 0) {
        Wells <- lapply(well_vars,
            function(type) {
                df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t")
                df$Location <- stringr::str_extract(df$Location, "[A-Z]+[0-9]+")
                colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "_")
                df
            }
        )
        Wells <- Reduce(
            function(...) merge.data.frame(by = id_vars[-3], sort = FALSE, ...),
            Wells
        )
    }

    # Parse Bead Data
    bead_vars <- intersect(bead_vars, fields)
    Beads <- lapply(bead_vars,
        function(type) {
            cols <- stringr::str_split(Lines[[type]][-1], "\t")[1:3]
            col_names <- stringr::str_replace_all(vapply(cols, "[", "", 1), ":", "")
            df <- setNames(lapply(cols, function(x) x[2:max(2, length(x))]), col_names)
            if ("BeadID" %in% col_names) df[["BeadID"]] <- as.integer(df[["BeadID"]])
            if (type == "Per Bead Count") df[[3]] <- as.integer(df[[3]])
#            df[[type]] <- c(df[[type]], rep(NA, length(df$BeadID) - length(df[[type]])))
            df <- data.frame(df, stringsAsFactors = FALSE)
            colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "_")
            df
        })
    Beads <- Reduce(
        function(...) merge.data.frame(sort = FALSE, ...),
        Beads
    )

    # RunData
    run_vars <- intersect(run_vars, fields)
    Run <- lapply(run_vars,
        function(type) {
            df <- .chr_vector_to_df(Lines[[type]][-1], header = TRUE, sep = "\t")
            if ("Location" %in% colnames(df)) df$Location <- stringr::str_extract(df$Location, "[A-Z]+[0-9]+")
            colnames(df) <- stringr::str_replace_all(colnames(df), "\\s+", "_")
            df
        }
    )
    names(Run) <- run_vars

    return(list(
        Exprs = Exprs,
        Average = Average,
        Wells = Wells,
        Beads = Beads,
        Run = Run
    ))
}

