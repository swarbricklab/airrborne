# Sum TRUST4-reported read support before downstream QC over the seven canonical
# loci, including productive and non-productive reported observations and
# records whose CDR3 TRUST4 imputed. Unknown loci are outside these totals.
# TRUST4's own report selection remains embedded in the input. The sums measure
# recovered read support, not transcripts or cells. Base R only.

SUPPORT_BCR_LOCI <- c("IGH", "IGK", "IGL")
SUPPORT_TCR_LOCI <- c("TRA", "TRB", "TRD", "TRG")
SUPPORT_LOCI <- c(SUPPORT_BCR_LOCI, SUPPORT_TCR_LOCI)

#' Build one library's TRUST4-reported support totals before pipeline QC.
#'
#' support_by_locus contains summed recovered read support for every SUPPORT_LOCI
#' entry, already validated with as_read_support(). Return one row with library_id,
#' total, BCR and TCR support and one <locus>_trust4_consensus_count per locus.
#' A library with no assemblies supplies a zero-filled vector and still gets a row.
build_repertoire_support_row <- function(sample_label, support_by_locus) {
  tcr <- sum(support_by_locus[SUPPORT_TCR_LOCI])
  bcr <- sum(support_by_locus[SUPPORT_BCR_LOCI])
  result <- data.frame(
    library_id = sample_label,
    trust4_consensus_count = tcr + bcr,
    bcr_trust4_consensus_count = bcr,
    tcr_trust4_consensus_count = tcr,
    stringsAsFactors = FALSE
  )
  for (locus in SUPPORT_LOCI) {
    result[[paste0(tolower(locus), "_trust4_consensus_count")]] <- support_by_locus[[locus]]
  }
  result
}

#' Join input depth to the support row and add yield per million input read pairs.
#'
#' depth is the library's one-row input depth table. Return library_id, then
#' n_fastq_file_pairs and input_read_pairs as the depth table spells them,
#' input_bases from its total_bases, yield_status, the total, BCR and TCR support,
#' then for each locus <locus>_trust4_consensus_count and
#' <locus>_yield_per_million_read_pairs. TRUST4-reported support retains TRUST4's
#' own report selection. The denominator counts read pairs, so both mates of one
#' fragment are one unit; dividing by total reads instead would halve every rate.
#' A zero input depth has unavailable yield; a positive depth with no recovered
#' support has zero yield.
add_repertoire_yield <- function(support, depth) {
  pairs <- as_read_support(depth$input_read_pairs, "input_read_pairs")
  if (length(pairs) != 1L || pairs != floor(pairs)) {
    stop("input_read_pairs must be one non-negative integer")
  }
  result <- data.frame(
    library_id = support$library_id,
    n_fastq_file_pairs = depth$n_fastq_file_pairs,
    input_read_pairs = pairs,
    input_bases = depth$total_bases,
    yield_status = if (pairs > 0) "available" else "no_input_read_pairs",
    trust4_consensus_count = support$trust4_consensus_count,
    bcr_trust4_consensus_count = support$bcr_trust4_consensus_count,
    tcr_trust4_consensus_count = support$tcr_trust4_consensus_count,
    stringsAsFactors = FALSE
  )
  for (locus in tolower(SUPPORT_LOCI)) {
    count <- paste0(locus, "_trust4_consensus_count")
    result[[count]] <- support[[count]]
    result[[paste0(locus, "_yield_per_million_read_pairs")]] <- if (pairs > 0) {
      support[[count]] / pairs * 1e6
    } else NA_real_
  }
  result
}
