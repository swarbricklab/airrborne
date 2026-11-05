# Read string identifiers first; convert only named measurements to numbers.
# Run from the checkout with Rscript examples/read_release.R RELEASE [METADATA_CSV].
read_rows <- function(path, sep = "\t") {
  read.table(path, header = TRUE, sep = sep, quote = if (sep == ",") "\"" else "", comment.char = "",
             colClasses = "character", na.strings = "", check.names = FALSE)
}

number <- function(x) as.numeric(x)

row_keys <- function(rows, keys) {
  if (anyNA(rows[, keys, drop = FALSE])) stop("Missing join key")
  do.call(paste, c(rows[, keys, drop = FALSE], sep = "\r"))
}

unique_keys <- function(rows, keys) {
  value <- row_keys(rows, keys)
  if (anyDuplicated(value)) stop("Duplicate lookup key")
  value
}

add_metadata <- function(headline, metadata, key = "tissue_id") {
  index <- match(headline[[key]], unique_keys(metadata, key))
  columns <- setdiff(names(metadata), key)
  overlap <- intersect(columns, names(headline))
  for (column in overlap) {
    present <- !is.na(index)
    if (anyNA(metadata[index[present], column]) ||
        any(headline[present, column] != metadata[index[present], column])) {
      stop("Metadata conflicts with scientific columns")
    }
  }
  columns <- setdiff(columns, overlap)
  joined <- cbind(headline, metadata[index, columns, drop = FALSE])
  stopifnot(nrow(joined) == nrow(headline), identical(joined$library_id, headline$library_id))
  joined
}

join_member_airr <- function(members, airr) {
  representatives <- airr[airr$is_member_representative == "T", , drop = FALSE]
  unique_keys(members, c("scope_id", "sequence_id"))
  index <- match(row_keys(members, c("library_id", "sequence_id")),
                 unique_keys(representatives, c("library_id", "sequence_id")))
  stopifnot(!anyNA(index))
  matched <- representatives[index, , drop = FALSE]
  stopifnot(all(number(matched$member_consensus_count) == number(members$member_consensus_count)))
  joined <- cbind(members, sequence = matched$sequence,
                   source_consensus_count = matched$consensus_count)
  stopifnot(nrow(joined) == nrow(members),
            sum(number(joined$member_consensus_count)) == sum(number(members$member_consensus_count)))
  joined
}

main <- function(args) {
  if (!length(args)) stop("Supply the release directory and optionally a metadata CSV")
  release <- args[1]
  headline <- read_rows(file.path(release, "metrics/headline.tsv"))
  unique_keys(headline, c("scope_id", "library_id"))
  if (length(args) > 1) headline <- add_metadata(headline, read_rows(args[2], sep = ","))
  headline$input_read_pairs <- number(headline$input_read_pairs)
  print(headline[, c("library_id", "input_read_pairs", "igh_grouped_consensus_count")])
  shm <- read_rows(file.path(release, "metrics/shm.tsv"))
  shm <- shm[shm$isotype_group == "all_IGH" & shm$v_region == "whole_v", , drop = FALSE]
  stopifnot(!anyDuplicated(row_keys(shm, c("scope_id", "library_id"))))
  members <- read_rows(file.path(release, "metrics/detail/igh/clone_members.tsv"))
  for (library in unique(members$library_id)) {
    selected <- members[members$library_id == library, , drop = FALSE]
    airr <- read_rows(gzfile(file.path(release, "airr", paste0(library, "_rearrangements.tsv.gz"))))
    joined <- join_member_airr(selected, airr)
    cat(library, nrow(joined), "IGH members; transferred support", sum(number(joined$member_consensus_count)), "\n")
  }
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
