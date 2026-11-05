# Construct the shared exact-junction and resolved V/J key using base R.
# The separator cannot occur in a nucleotide sequence or gene token. ASCII `~`
# sorts after every character of those parts in the C collation, so keys order
# as their parts do. pRESTO and Change-O do not use it as a delimiter.
.GROUP_SEP <- "~"

#' Build the exact rearrangement-group key (vectorised).
#'
#' Join the uppercased A/C/G/T junction and canonical v_gene and j_gene tokens
#' with `~`. Return one key per record, or NA when the junction is invalid or
#' either entity is missing. Empty vectors return character(0). The caller
#' resolves complete V/J calls before constructing this library-independent key.
#' build_rearrangement_group_table() in rearrangement_groups.R applies non-IGH
#' eligibility before grouping on this key. build_clone_membership_frame() in
#' igh_clones.R also writes it for sequence matching, not IGH clone assignment.
#'
#' An amino-acid grouping change needs more than replacing junction here: review
#' non-IGH eligibility, representative annotations and support in
#' build_rearrangement_group_table(), downstream key comparisons and tests.
#' Different nucleotide junctions can encode the same amino-acid sequence.
exact_rearrangement_group_key <- function(junction, v_token, j_token) {
  junction <- toupper(trimws(as.character(junction)))
  v_token <- as.character(v_token)
  j_token <- as.character(j_token)
  key <- paste(junction, v_token, j_token, sep = .GROUP_SEP)
  key[is.na(junction) | !grepl("^[ACGT]+$", junction) |
        is_airr_null(v_token) | is_airr_null(j_token)] <- NA_character_
  key
}

