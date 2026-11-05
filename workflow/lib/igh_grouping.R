# Prepare IGH grouping once per library beside unchanged common-QC rows.
# Fitting, assignment, audit and export read the same representative links.
# These helpers use base R.

#' Select the sequence rows available for IGH clone fitting and assignment.
#'
#' This preparation does not assign clone IDs. IGH clone assignment links similar
#' nucleotide junctions at the chosen threshold; a chain of links can join
#' sequences farther apart than that threshold. fit_donor_threshold.R
#' estimates donor thresholds, scope_clone_threshold.R selects the scope threshold
#' and assign_scope_donor_clones.R calls SCOPer hierarchicalClones for each donor.
#'
#' `records` are the library's common-QC IGH rows and `max_n` is
#' `params.lineage.max_n`, the most N or ? symbols a grouping-eligible junction
#' may carry. Eligibility requires a valid junction and one accepted V and J
#' assignment, either a single gene or an identical-reference group of genes.
#' clone_junction_eligible() in igh_clones.R checks the nucleotide alphabet and
#' missing-symbol allowance. vj_tie_kind() in gene_calls.R supplies the V/J statuses.
#' Present C calls must be locus-compatible and have complete candidate lists,
#' checked by grouping_constant_status() in gene_calls.R. A missing C call passes.
#' The dominant common-QC backbone variant is chosen before eligibility, without
#' substituting a minor eligible variant. Status records semicolon-separated
#' failure reasons or ineligible_backbone_representative for blocked siblings.
#' Retained representatives then merge by exact full sequence. Return the mapping
#' from prepare_assembly_mapping(), including a typed empty map for no records.
prepare_igh_grouping <- function(records, library_id, max_n) {
  require_columns(records, c("junction", "v_gene_status", "j_gene_status"),
                  "Common-QC IGH")
  if (any(records$locus != "IGH")) stop("IGH preparation requires only IGH records")
  junction_failed <- !clone_junction_eligible(records$junction, max_n)
  v_failed <- !records$v_gene_status %in% c("resolved", "assigned_group")
  j_failed <- !records$j_gene_status %in% c("resolved", "assigned_group")
  constant_status <- grouping_constant_status(records)
  constant_failed <- constant_status != "eligible"
  # V and J failures name the segment and its status, as non-IGH grouping does.
  reasons <- lapply(seq_len(nrow(records)), function(i) {
    failures <- c(if (junction_failed[i]) "ineligible_junction",
                  if (v_failed[i]) paste0("v_call_", records$v_gene_status[i]),
                  if (j_failed[i]) paste0("j_call_", records$j_gene_status[i]),
                  if (constant_failed[i]) constant_status[i])
    if (length(failures)) paste(failures, collapse = ";") else "eligible"
  })
  status <- unlist(reasons, use.names = FALSE)
  eligible <- !junction_failed & !v_failed & !j_failed & !constant_failed
  prepare_assembly_mapping(records, library_id, eligible, status, collapse = TRUE)
}

#' Return one representative per distinct retained full sequence in each library.
#'
#' mapping must describe every common-QC IGH record. inference_support selects
#' transferred support or the representative's own support; neither distToNearest
#' nor hierarchicalClones uses it. Membership applies the transfer once for release.
#' No representatives returns an empty records frame with the same columns.
prepared_igh_representatives <- function(records, mapping, inference_support) {
  mapping <- validate_grouping_mapping(records, mapping)
  if (any(records$locus != "IGH")) stop("Prepared donor pool must contain only IGH")
  keep <- mapping$is_inference_representative
  out <- records[keep, , drop = FALSE]
  if (inference_support) out$consensus_count <- mapping$inference_consensus_count[keep]
  out
}
