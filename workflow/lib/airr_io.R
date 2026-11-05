# Serialise AIRR tables with T/F booleans, empty nulls and unquoted TSV fields.
# Column order is preserved. These helpers use base R.

#' Convert logical columns to T/F and missing or blank cells to empty strings.
#'
#' Returns the frame in input column order. Literal NA identifiers stay
#' unchanged, and the caller's numeric-formatting options are restored on exit.
normalise_airr_frame <- function(df) {
  previous <- options(scipen = TABLE_SCIPEN)
  on.exit(options(previous), add = TRUE)
  for (col in names(df)) {
    v <- df[[col]]
    if (is.logical(v)) {
      df[[col]] <- ifelse(is.na(v), "", ifelse(v, "T", "F"))
    } else {
      cv <- as.character(v)
      # Literal strings such as NA and None can be identifiers. Only actual
      # missing values and blank cells become AIRR nulls.
      cv[is_airr_null(v)] <- ""
      df[[col]] <- cv
    }
  }
  df
}

#' Write an AIRR Rearrangement TSV and return its path invisibly.
#'
#' Uses explicit T/F and null serialisation before the shared atomic writer.
#' airr::write_rearrangement() applies package-dependent type coercion, so this
#' writer keeps the output tokens stable across package versions.
write_airr_table <- function(df, path) {
  write_delimited_table(normalise_airr_frame(df), path)
}
