# Read pipeline TSVs as strings and write deterministic tables atomically.
# Shared parsers distinguish blank AIRR nulls, literal NA identifiers, booleans
# and non-negative read support. AIRR-specific writers remain in airr_io.R.
# These helpers use base R.

#' Return TRUE for actual missing values and cells blank after trimming.
#'
#' Literal NA remains a string: AIRR nulls use empty fields.
is_airr_null <- function(x) is.na(x) | trimws(as.character(x)) == ""

#' Parse each AIRR boolean as true, false, missing or invalid.
#'
#' Accepts logical values and case-insensitive T/TRUE/1 and F/FALSE/0 tokens.
#' AIRR nulls are missing; all other tokens, including literal NA, are invalid.
parse_airr_boolean <- function(x) {
  if (is.logical(x)) {
    return(ifelse(is.na(x), "missing", ifelse(x, "true", "false")))
  }
  up <- toupper(trimws(as.character(x)))
  out <- rep("invalid", length(up))
  out[is_airr_null(x)] <- "missing"
  out[up %in% c("T", "TRUE", "1")] <- "true"
  out[up %in% c("F", "FALSE", "0")] <- "false"
  out
}

#' Convert a consensus_count field to numeric TRUST4 read support.
#'
#' Read support may be zero, but missing, non-numeric, infinite or negative
#' values are invalid because every downstream consensus_count metric reads this
#' same quantity.
as_read_support <- function(x, field = "read support") {
  value <- suppressWarnings(as.numeric(x))
  invalid <- !is.finite(value) | value < 0
  if (any(invalid)) {
    stop(sprintf("%s contains %d invalid value(s); expected finite non-negative counts",
                 field, sum(invalid)))
  }
  value
}

# Writers use empty fields for nulls and favour fixed numeric notation.
TABLE_SCIPEN <- 12L

# Create the parent directory of every path in a possibly nested list of
# outputs. Side effect only; returns nothing.
ensure_output_dirs <- function(outputs) {
  for (path in unname(unlist(outputs, use.names = FALSE))) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
  invisible(NULL)
}

#' Stop if df lacks required_columns, naming context in the error.
#'
#' Returns nothing. A zero-row frame is valid if its columns are present.
require_columns <- function(df, required_columns, context) {
  missing_columns <- setdiff(required_columns, colnames(df))
  if (length(missing_columns)) {
    stop(context, " is missing required columns: ",
         paste(missing_columns, collapse = ", "))
  }
  invisible(NULL)
}

#' Read a TSV as character columns, preserving literal NA and identifier zeros.
#'
#' Checks required_columns and uses context (the path by default) in errors.
#' Blank fields stay empty strings; a header-only input returns zero rows.
read_pipeline_table <- function(path, required_columns = character(0), context = path) {
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                          colClasses = "character", na.strings = character(0))
  require_columns(df, required_columns, context)
  df
}

#' Build a frame with the given columns and no rows.
#'
#' `columns` are the column names, all typed as character. Used where a stage has
#' nothing to report but still owes its consumers a table with the released
#' header.
make_empty_frame <- function(columns) {
  data.frame(stats::setNames(rep(list(character(0)), length(columns)), columns),
             stringsAsFactors = FALSE, check.names = FALSE)
}

#' Combine frames that may not share the same columns.
#'
#' `frames` is a list of data frames; empty ones are dropped. `columns` fixes the
#' output columns, which matters when a frame that contributed no rows is still
#' the only one declaring a column; leave it NULL to use the union of whatever
#' the non-empty frames carry. Columns a frame lacks are filled with NA.
#'
#' Returns one frame, or a zero-row frame with `columns` when nothing was left.
bind_frames_on_column_union <- function(frames, columns = NULL) {
  frames <- frames[vapply(frames, nrow, integer(1)) > 0L]
  if (!length(frames)) {
    return(if (is.null(columns)) data.frame() else make_empty_frame(columns))
  }
  if (is.null(columns)) {
    columns <- Reduce(union, lapply(frames, names))
  }
  frames <- lapply(frames, function(frame) {
    for (column in setdiff(columns, names(frame))) frame[[column]] <- NA
    frame[, columns, drop = FALSE]
  })
  do.call(rbind, frames)
}

#' Read a column, or a constant where the frame does not have it.
#'
#' Returns a vector as long as the frame has rows, so the result can be assigned
#' straight back into it.
column_or_default <- function(df, column, default = NA) {
  if (column %in% names(df)) df[[column]] else rep(default, nrow(df))
}

#' Write a metric table.
#'
#' columns selects the columns and their output order; NULL keeps the input order.
#' Writes the file and returns its path invisibly.
write_pipeline_table <- function(df, path, columns = NULL) {
  write_delimited_table(df, path, columns = columns)
}

#' Write a TSV through a temporary file in the destination directory.
#'
#' columns selects and orders required columns when supplied. Missing values become
#' empty fields. Returns the path invisibly after an atomic rename, so existing
#' hard links retain the previous contents.
#'
#' TABLE_SCIPEN favours fixed notation for integer readers while allowing
#' scientific notation for small values. Restores the caller's options on exit.
write_delimited_table <- function(df, path, columns = NULL) {
  if (!is.null(columns)) {
    require_columns(df, columns, path)
    df <- df[, columns, drop = FALSE]
  }
  previous <- options(scipen = TABLE_SCIPEN)
  on.exit(options(previous), add = TRUE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "-"),
                        tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.table(df, temporary, sep = "\t", row.names = FALSE, quote = FALSE,
                     na = "")
  if (!file.rename(temporary, path)) stop("Could not publish table: ", path)
  invisible(path)
}
