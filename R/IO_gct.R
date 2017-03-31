
read_gct <- function(path, na.strings = c("", "na", "NA", "-666", "-666.0"),
                     tidy = FALSE) {
    header <- readLines(path, n = 2)
    dims <- as.integer(strsplit(header[2], "\t")[[1]])
    if (header[1] == "#1.2") {
       dims <- c(dims, 1, 0)
    } else if (header[1] != "#1.3") {
        stop(paste("Unknown verstion of gct", header[1]))
    }
    Ncols <- dims[2] + dims[3] + 1  # main + meta + cids

    ### Meta Data for Columns
    meta_cols <- readLines(path, n = 3 + dims[4])[3:(3 + dims[4])]
    meta_cols <- strsplit(meta_cols, "\t")
    names(meta_cols) <- vapply(meta_cols, "[[", "", 1)
    names(meta_cols)[1] <- "cid"
    meta_rows_header <- c("rid", meta_cols$cid[2:(dims[3] + 1)])  # for later use
    meta_cols <- lapply(meta_cols, "[", (dims[3] + 2):Ncols)  # drop dead-area
    meta_cols <- lapply(meta_cols, type.convert, na.strings = na.strings, as.is = TRUE)
    meta_cols <- as.data.frame(meta_cols, stringsAsFactors = FALSE)

    ### Meta Data for Rows
    main_rows <- read.table(path, sep = "\t", na.strings = na.strings,
                            skip = 3 + dims[4], stringsAsFactors = FALSE)
    meta_rows <- main_rows[, 1:(dims[3] + 1)]
    colnames(meta_rows) <- meta_rows_header

    ### Main matrix
    mat <- unname(as.matrix.data.frame(main_rows[, (dims[3] + 2):Ncols]))
    dimnames(mat) <- list(rid = meta_rows$rid, cid = meta_cols$cid)

    ### Return
    gct <- list(matrix = mat, meta = list(rows = meta_rows, cols = meta_cols))
    if (tidy) return(tidy_gct(gct))
    return(gct)
}
