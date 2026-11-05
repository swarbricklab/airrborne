# Namespace library-local sequence IDs as library_id~sequence_id for donor pools.
# Use the same key for joins to library annotations and scope membership.
# These helpers use base R.

# `~` sorts after every character a library or sequence ID can contain, and
# pRESTO and Change-O do not use it as an annotation delimiter.
POOLED_SEQUENCE_ID_SEPARATOR <- "~"

#' The pooled identifier of one library-local record.
#'
#' `library_id` names the library and `sequence_id` is the record's library-local
#' AIRR identifier; both are vectorised. Returns the `<library_id>~<sequence_id>`
#' string that pool_scope_donor_igh.R writes into the pooled table.
make_pooled_sequence_id <- function(library_id, sequence_id) {
  paste(library_id, sequence_id, sep = POOLED_SEQUENCE_ID_SEPARATOR)
}

#' Attach per-library IGH evidence to pooled clone members in membership order.
#'
#' features uses library-local sequence_id; membership uses library_id~sequence_id.
#' Each member must match one feature row. Unused common-QC feature rows are valid.
#' region_levels names the IMGT V regions whose rates, mutation counts and coverage
#' counts are attached alongside whole-V evidence and qualified_c_call.
#'
#' Numeric evidence accepts blank text or an actual R NA as missing. The literal
#' text "NA" and other malformed numeric text stop the join with the column and
#' pooled record ID. Empty membership receives
#' typed empty columns, including when features has no rows or columns.
attach_igh_features_to_membership <- function(membership, features,
                                             region_levels) {
  numeric_columns <- c("shm_frequency", "informative_v",
                       paste0("shm_frequency_", region_levels),
                       paste0("mutations_", region_levels),
                       paste0("informative_", region_levels),
                       paste0("informative_supported_", region_levels),
                       paste0("germline_occupied_", region_levels),
                       paste0("median_majority_votes_", region_levels))
  if (nrow(features) > 0L) {
    require_columns(features, c("library_id", "sequence_id", "qualified_c_call",
                                numeric_columns), "IGH feature table")
    features$sequence_id <- make_pooled_sequence_id(features$library_id,
                                                      features$sequence_id)
    if (anyDuplicated(features$sequence_id)) {
      stop("sequence_id is not unique across IGH feature tables")
    }
  }
  feature_index <- match(membership$sequence_id, features$sequence_id)
  if (anyNA(feature_index)) {
    stop("Lineage membership contains ", sum(is.na(feature_index)),
         " sequence(s) missing from IGH features")
  }
  membership$qualified_c_call <- as.character(features$qualified_c_call[feature_index])
  for (column in numeric_columns) {
    values <- features[[column]][feature_index]
    numbers <- suppressWarnings(as.numeric(values))
    invalid <- !is_airr_null(values) & is.na(numbers) & !is.nan(numbers)
    if (any(invalid)) {
      i <- which(invalid)[1L]
      stop("IGH feature column ", column, " must contain numeric evidence or missing values; ",
           "record ", membership$sequence_id[i], " has value '", values[i], "'")
    }
    membership[[column]] <- numbers
  }
  membership
}
