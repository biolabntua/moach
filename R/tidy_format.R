
tidy_gct <- function(gct, value_col = "value") {
    df <- tibble::as_data_frame(gct$mat)
    df <- tibble::rownames_to_column(df, "rid")
    df <- tidyr::gather_(df, "cid", value_col, colnames(df)[2:ncol(df)])
    df <- dplyr::left_join(df, gct$meta$cols, "cid")
    df <- dplyr::left_join(df, gct$meta$rows, "rid")
    df <- dplyr::select_(df, c(colnames(gct$meta$rows), colnames(gct$meta$cols), value_col))
    return(df)
}
