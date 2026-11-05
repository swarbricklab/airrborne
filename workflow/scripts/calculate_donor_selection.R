#!/usr/bin/env Rscript
# Estimate BASELINe selection by library and qualified isotype group within a
# donor/scope. Record input counts before clone collapse and estimate availability
# afterwards. The pinned Immcantation image supplies SHazaM 1.3.2; the method and
# its attribution are in workflow/lib/baseline_selection.R. That module owns
# selection_record_eligible(), clone consensus, incomplete-codon masking and
# BASELINe estimation. Selection reads alignments and qualified calls directly;
# regional SHM coverage does not decide which records enter selection.

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse)
  library(shazam)
})
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------ args ----
source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "baseline_selection", "gene_calls", "isotype_calls",
  "isotype_resolution", "pooled_records", "qc_audit"))
stage_log <- start_stage_log("calculate_donor_selection")

script_options <- parse_script_arguments(list(
  make_option("--pooled_file", type = "character",
              help = "Donor-pooled QC-pass IGH table for one scope. Required."),
  make_option("--membership_file", type = "character",
              help = paste("Clone membership table for the same scope and",
                           "donor. Required.")),
  make_option("--manifest_file", type = "character",
              help = paste("Scope manifest naming the donor's libraries, with",
                           "their tissue. Required.")),
  make_option("--qc_audit_file", type = "character", help = "Selection input counts. Required."),
  make_option("--out_file", type = "character",
              help = paste("Path to write the library x isotype BASELINe",
                           "shard to. Required.")),
  make_option("--scope_id", type = "character",
              help = "Analysis scope every output row is stamped with. Required."),
  make_option("--donor_id", type = "character",
              help = "Donor whose lineages are being summarised. Required."),
  make_option("--threads", type = "integer",
              help = paste("Cores handed to shazam's collapseClones,",
                           "calcBaseline, groupBaseline and summarizeBaseline;",
                           "a positive integer. Required.")),
  make_option("--collapse_min_freq", type = "double",
              help = paste("Frequency a base must reach to enter a clone's",
                           "effective sequence, in (0, 1]. Required.")),
  make_option("--isotype_call_files", type = "character",
              help = paste("Comma-separated qualified IGH call tables, one per",
                           "library of this donor. Required."))
))
validate_selection_collapse_frequency(script_options$collapse_min_freq)
pooled_file <- script_options$pooled_file
membership_file <- script_options$membership_file
manifest_file <- script_options$manifest_file
out_file <- script_options$out_file
scope_id <- script_options$scope_id
donor_id <- script_options$donor_id
threads <- script_options$threads
collapse_min_freq <- script_options$collapse_min_freq
isotype_call_files <- split_path_list(script_options$isotype_call_files)

if (!is.finite(threads) || threads < 1L) stop("threads must be a positive integer")

log_info("Scope %s, donor %s, nproc %d, collapse_min_freq %s",
         scope_id, donor_id, threads, collapse_min_freq)

# =========================================================== sample metadata ==
manifest_required <- c("library_id", "donor_id", "tissue_id")
manifest <- read_pipeline_table(manifest_file,
                                required_columns = manifest_required,
                                context = "Scope manifest")
samples <- manifest[
  as.character(manifest$donor_id) == donor_id,
  manifest_required,
  drop = FALSE
]
if (!nrow(samples)) stop("Scope manifest contains no libraries for donor ", donor_id)
if (anyDuplicated(samples$library_id)) stop("Scope manifest library_id is not unique")
samples <- samples[order(samples$library_id), , drop = FALSE]

# ================================================================ read data ==
db <- read_pipeline_table(
  pooled_file,
  required_columns = c("sequence_id", "library_id", "sequence_alignment",
                       "germline_alignment_d_mask", "v_call", "j_call", "consensus_count"),
  context = "Pooled IGH"
)
if (anyDuplicated(db$sequence_id)) stop("Pooled sequence_id is not unique")
outside_manifest <- setdiff(unique(as.character(db$library_id)), as.character(samples$library_id))
if (length(outside_manifest)) {
  stop("Pooled IGH contains a library outside this scope/donor: ",
       paste(outside_manifest, collapse = ", "))
}

mem <- read_pipeline_table(
  membership_file,
  required_columns = c("library_id", "sequence_id", "clone_id", "consensus_count"),
  context = "Lineage membership"
)
if (anyDuplicated(mem$sequence_id)) stop("Lineage membership sequence_id is not unique")
log_info("Read %d pooled QC-pass IGH record(s) from %s and %d lineage member(s) from %s",
         nrow(db), pooled_file, nrow(mem), membership_file)

# Join the qualification sidecars directly. Selection does not consume SHM
# measurements. Namespace library-local keys to match the unchanged donor pool.
calls <- do.call(rbind, lapply(isotype_call_files, read_pipeline_table))
require_columns(calls, c("library_id", "sequence_id", "qualified_c_call"), "IGH isotype calls")
calls$sequence_id <- make_pooled_sequence_id(calls$library_id, calls$sequence_id)
db <- join_qualified_isotype_calls(db, calls)

source_index <- match(as.character(mem$sequence_id), as.character(db$sequence_id))
if (anyNA(source_index)) stop("Lineage membership contains a sequence outside the common-QC pool")
if (any(as.character(db$library_id[source_index]) != as.character(mem$library_id))) {
  stop("Pooled IGH and lineage membership disagree on library_id")
}

# Membership already applies grouping and collapse. Selection eligibility starts
# from these members, retaining both original and transferred read support.
seqs <- db[source_index, , drop = FALSE]
seqs$clone_id <- as.character(mem$clone_id)
seqs$source_consensus_count <- as_read_support(seqs$consensus_count)
seqs$consensus_count <- as_read_support(mem$consensus_count)
seqs$isotype_group <- resolve_analysis_igh_isotype(
  as.character(seqs$qualified_c_call)
)

# ========================================== eligibility and isotype grouping ==
eligible_records <- seqs[
  selection_record_eligible(seqs$clone_id, seqs$sequence_alignment,
                            seqs$germline_alignment_d_mask),
  , drop = FALSE
]

audit <- list()
for (i in seq_len(nrow(samples))) {
  library <- samples$library_id[i]
  for (group in ANALYSIS_IGH_ISOTYPE_GROUPS) {
    before <- seqs[seqs$library_id == library &
      (group == "all_IGH" | (!is.na(seqs$isotype_group) & seqs$isotype_group == group)), , drop = FALSE]
    after <- eligible_records[eligible_records$library_id == library &
      (group == "all_IGH" | (!is.na(eligible_records$isotype_group) &
                            eligible_records$isotype_group == group)), , drop = FALSE]
    for (support_column in c("consensus_count", "source_consensus_count")) {
      stage <- if (support_column == "consensus_count") "selection_input" else
        "selection_input_source_support"
      audit[[length(audit) + 1L]] <- qc_audit_rows(
        library, stage, nrow(before), nrow(after),
        sum(before[[support_column]]), sum(after[[support_column]]),
        analysis = "selection", stage_index = 70, kind = "eligibility",
        input_unit = "inference_representative", output_unit = "selection_input_record",
        scope_id = scope_id, donor_id = donor_id, tissue_id = samples$tissue_id[i],
        isotype_group = group, support_conserved = FALSE)
    }
  }
}
write_qc_audit(do.call(rbind, audit), script_options$qc_audit_file)

# A donor can hold no record a density could be built from: none of its pooled
# IGH joined a clone, or none carried a D-masked germline. An empty repertoire
# is a valid result, so the donor's libraries are reported as contributing no
# clone to any isotype group and the estimator is skipped.
if (nrow(eligible_records) == 0L) {
  selection_rows <- build_empty_selection_rows(samples, scope_id)
  write_pipeline_table(selection_rows, out_file)
  log_info(paste0("Clone members %d, none with an alignment and a D-masked ",
                  "germline. Wrote %d row(s) reporting no eligible lineage to %s"),
           nrow(seqs), nrow(selection_rows), out_file)
  finish_stage_log(stage_log)
  quit(save = "no", status = 0)
}

# Each record enters all-IGH and, if its C call resolved, its own isotype group.
# Grouping repeats records only for these separate estimates.
eligible <- build_isotype_grouped_selection_input(eligible_records)
log_info("Clone members %d, with an alignment and a D-masked germline %d, grouped %d",
         nrow(seqs), nrow(eligible_records), nrow(eligible))

# =============================================== collapse, estimate, report ==
log_info("Building one effective sequence per library, isotype and lineage from %d grouped record(s)",
         nrow(eligible))
collapsed <- collapse_clones_to_effective_sequences(
  eligible, collapse_min_freq, threads)
log_info("Effective sequences %d, with gaps or incomplete codons masked in %d",
         nrow(collapsed$effective_sequences), collapsed$n_masked)
representatives <- collapsed$effective_sequences

baseline_results <- setNames(lapply(
  ANALYSIS_IGH_ISOTYPE_GROUPS,
  function(group) estimate_baseline_selection(representatives, group, threads)
), ANALYSIS_IGH_ISOTYPE_GROUPS)

selection_rows <- build_selection_result_rows(baseline_results, representatives,
                                              samples, scope_id)

write_pipeline_table(selection_rows, out_file)
log_info("Wrote %d library x isotype group rows to %s", nrow(selection_rows), out_file)

# One line per released row, so the sigma estimates and the clone counts
# behind them are readable in the log without opening the table.
for (index in seq_len(nrow(selection_rows))) {
  row <- selection_rows[index, , drop = FALSE]
  log_info(
    paste0("%s %s: n_clone_consensus_sequences=%d used_cdr=%d used_fwr=%d sigma_cdr=%.4f ",
           "sigma_fwr=%.4f sigma_cdr_status=%s sigma_fwr_status=%s"),
    row$library_id, row$isotype_group, row$n_clone_consensus_sequences,
    row$n_clone_consensus_sequences_used_cdr, row$n_clone_consensus_sequences_used_fwr,
    row$sigma_cdr, row$sigma_fwr, row$sigma_cdr_status, row$sigma_fwr_status)
}

finish_stage_log(stage_log)
