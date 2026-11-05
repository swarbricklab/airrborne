#!/usr/bin/env Rscript
# Build exact non-IGH clonotypes from common-QC AIRR and write their prepared
# grouping map and audit. build_rearrangement_group_table() in
# rearrangement_groups.R owns eligibility and the nucleotide group definition.
# prepare_assembly_mapping() in assembly_backbone.R transfers support;
# exact_rearrangement_group_key() in rearrangement_keys.R builds comparison keys.

suppressPackageStartupMessages(library(optparse))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys",
  "rearrangement_groups", "qc_audit", "pooled_records"))
stage_log <- start_stage_log("build_clonotypes")

NONIGH_LOCI <- c("IGK", "IGL", "TRA", "TRB", "TRD", "TRG")

script_options <- parse_script_arguments(c(
  list(
    make_option("--input_file", type = "character",
                help = "Common QC-pass AIRR table for one library. Required."),
    make_option("--library_id", type = "character",
                help = "Library the records belong to. Required."),
    make_option("--membership_file", type = "character",
                help = "Path to write representative-to-clonotype membership. Required."),
    make_option("--grouping_file", type = "character",
                help = "Path to write the non-IGH grouping map to. Required."),
    make_option("--qc_audit_file", type = "character",
                help = "Path to write the grouping population audit to. Required.")),
  lapply(NONIGH_LOCI, function(locus) {
    make_option(paste0("--", tolower(locus), "_groups_file"), type = "character",
                help = sprintf("Path to write the %s clonotype table to. Required.", locus))
  })
))
library_id <- script_options$library_id
input_path <- script_options$input_file
records <- read_pipeline_table(
  input_path,
  required_columns = c("sequence_id", "locus", "junction", "v_call", "j_call",
                       "v_gene_status", "j_gene_status",
                       "v_gene", "j_gene",
                       "consensus_count"),
  context = "common QC-pass AIRR table")
log_info("Loaded %d common-QC records from %s", nrow(records), input_path)

mappings <- list()
memberships <- list()
audit <- qc_audit_empty()
for (locus in NONIGH_LOCI) {
  output_path <- script_options[[paste0(tolower(locus), "_groups_file")]]
  locus_records <- records[as.character(records$locus) == locus, , drop = FALSE]
  built <- build_rearrangement_group_table(locus_records, library_id, locus)

  mappings[[locus]] <- built$mapping
  memberships[[locus]] <- built$membership
  audit <- rbind(audit, qc_audit_grouping_rows(library_id, locus_records, built$mapping, locus),
    qc_audit_rows(library_id, "exact_groups", nrow(built$groups),
      support_in = sum(built$groups$consensus_count), locus = locus,
      analysis = "grouping", stage_index = 34, input_unit = "clonotype"))
  if (built$n_records == 0L) {
    log_info("[%s] 0 records -> 0 clonotypes", locus)
  } else {
    log_info("[%s] Excluded %d of %d records from grouping after eligibility checks and backbone selection",
             locus, built$n_dropped, built$n_records)
    log_info("[%s] Of those, %d named several V or J germline entities",
             locus, built$n_dropped_ambiguous_vj)
    if (built$n_eligible == 0L) {
      log_info("[%s] 0 eligible records -> 0 clonotypes", locus)
    } else {
      log_info("[%s] %d eligible AIRR records -> %d nucleotide clonotypes",
               locus, built$n_eligible, nrow(built$groups))
    }
  }

  # The representative is a member, so its released ID takes the same library
  # qualification as clonotype_members.sequence_id and the AIRR export.
  groups <- built$groups
  groups$representative_sequence_id <- make_pooled_sequence_id(
    groups$library_id, groups$representative_sequence_id)
  write_pipeline_table(groups, output_path,
                       columns = REARRANGEMENT_GROUP_COLUMNS)
  log_info("[%s] Wrote %d clonotypes to %s", locus,
           nrow(built$groups), output_path)
}

mapping <- do.call(rbind, mappings)
source_ids <- records$sequence_id[records$locus != "IGH"]
mapping <- mapping[match(source_ids, mapping$sequence_id), , drop = FALSE]
write_pipeline_table(mapping, script_options$grouping_file, GROUPING_COLUMNS)
membership <- do.call(rbind, memberships)
membership$sequence_id <- if (nrow(membership)) {
  make_pooled_sequence_id(membership$library_id, membership$sequence_id)
} else character()
write_pipeline_table(membership, script_options$membership_file, CLONOTYPE_MEMBER_COLUMNS)
write_qc_audit(audit, script_options$qc_audit_file)
finish_stage_log(stage_log)
