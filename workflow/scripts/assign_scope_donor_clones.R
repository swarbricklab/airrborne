#!/usr/bin/env Rscript
# Group one donor's prepared IGH sequences into clones within one scope.
# SCOPer 1.5.0 hierarchicalClones performs the assignment (Gupta et al. 2017,
# doi:10.4049/jimmunol.1601850). The saved scope cutoff controls how different
# nucleotide junctions may be; donor sequences never join across donor boundaries.
# Preparation in igh_grouping.R has already selected eligible representatives.
# This script owns the SCOPer call, then writes membership and matching records
# with eligible source read support transferred once by igh_clones.R. Empty or
# wholly ineligible pools still write both tables with their declared columns.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "igh_clones",
  "igh_grouping", "qc_audit", "rearrangement_keys", "scope_clone_threshold"))
stage_log <- start_stage_log("assign_scope_donor_clones")
suppressPackageStartupMessages(library(optparse))
opts <- parse_script_arguments(list(
  make_option("--input_file", type = "character", help = "Full donor common-QC IGH pool. Required."),
  make_option("--grouping_file", type = "character", help = "Prepared donor grouping. Required."),
  make_option("--threshold_file", type = "character", help = "Operative scope threshold. Required."),
  make_option("--scope_id", type = "character", help = "Scope identifier. Required."),
  make_option("--donor_id", type = "character", help = "Donor identifier. Required."),
  make_option("--max_n", type = "integer", help = "Configured SCOPer missing-base allowance. Required."),
  make_option("--membership_file", type = "character", help = "Clone membership output. Required."),
  make_option("--records_file", type = "character", help = "Matching sequence-record output. Required."),
  make_option("--qc_audit_file", type = "character", help = "Assigned population counts. Required.")
))
threshold <- read_scope_clone_threshold(opts$threshold_file, opts$scope_id)
max_n <- opts$max_n
if (!is.finite(max_n) || max_n < 0L) stop("max_n must be a non-negative integer")
source <- read_pipeline_table(opts$input_file)
mapping <- validate_grouping_mapping(source, read_pipeline_table(opts$grouping_file))
# SCOPer receives one representative per distinct retained full sequence in
# each library. The membership builder adds contributing source support once.
cluster_db <- prepared_igh_representatives(source, mapping, inference_support = FALSE)
db <- source[mapping$grouping_eligible, , drop = FALSE]
eligible_mapping <- mapping[mapping$grouping_eligible, , drop = FALSE]
representative <- match(eligible_mapping$inference_representative_sequence_id, db$sequence_id)
out <- make_empty_frame(union(CLONE_MEMBER_COLUMNS, CLONE_RECORD_COLUMNS))
if (nrow(db)) {
  suppressPackageStartupMessages(library(scoper))
  cluster_db <- igh_inference_calls(cluster_db)
  # method = "nt" compares nucleotide junctions; normalize = "len" divides the
  # mismatch count by junction length. SCOPer compares equal-length junctions
  # within matching V/J groups. igh_inference_calls() gives each resolved V/J
  # reference-sequence group one temporary call, also used during SHazaM fitting.
  # Single linkage joins sequences through a chain of within-cutoff neighbours;
  # two members of a clone need not be within the cutoff of each other.
  # Change distance choices here and in fit_donor_distance_threshold() together.
  # cell_id stays NULL because bulk RNA-seq supplies no cell identities. cdr3 =
  # FALSE keeps the whole supplied junction; mod3 = FALSE adds no frame filter
  # after common QC. max_n is the same allowance used during preparation.
  res <- hierarchicalClones(
    cluster_db,
    threshold = threshold,
    method = "nt",
    linkage = "single",
    normalize = "len",
    junction = "junction",
    v_call = "v_inference_call",
    j_call = "j_inference_call",
    clone = "clone_id",
    fields = NULL,
    cell_id = NULL,
    first = FALSE,
    cdr3 = FALSE,
    mod3 = FALSE,
    max_n = max_n,
    summarize_clones = FALSE,
    verbose = FALSE,
    nproc = 1
  )

  if (!is.data.frame(res)) res <- as.data.frame(res)
  if (!("clone_id" %in% names(res))) {
    stop("hierarchicalClones did not return a clone_id column")
  }

  db$clone_id <- propagate_clone_id_to_siblings(db, representative, res)
  out <- build_clone_membership_frame(db, mapping, opts$scope_id, opts$donor_id)
}
write_pipeline_table(out, opts$membership_file, CLONE_MEMBER_COLUMNS)
write_pipeline_table(out, opts$records_file, CLONE_RECORD_COLUMNS)
audit <- qc_audit_rows(
  NA_character_, "donor_clones", length(unique(out$clone_id)),
  support_in = sum(as_read_support(out$consensus_count)), analysis = "grouping",
  stage_index = 35, input_unit = "donor_clone", grain = "donor",
  scope_id = opts$scope_id, donor_id = opts$donor_id)
for (library in sort(unique(source$library_id))) {
  local <- out[out$library_id == library, , drop = FALSE]
  audit <- rbind(audit, qc_audit_rows(
    library, "clone_observations", length(unique(local$clone_id)),
    support_in = sum(as_read_support(local$consensus_count)), analysis = "grouping",
    stage_index = 34, input_unit = "library_clone_observation",
    scope_id = opts$scope_id, donor_id = opts$donor_id))
}
write_qc_audit(audit, opts$qc_audit_file)
log_info("Assigned %d representatives to %d donor clones", nrow(out), length(unique(out$clone_id)))
finish_stage_log(stage_log)
