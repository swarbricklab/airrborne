#!/usr/bin/env Rscript
# Prepare one grouping-map row per common-QC IGH record in a library.
# prepare_igh_grouping() in igh_grouping.R checks junction and gene assignments;
# params.lineage.max_n controls the junction's N or ? allowance.
# prepare_assembly_mapping() in assembly_backbone.R chooses the dominant
# backbone variant before eligibility, then merges identical retained full DNA.
# The mapping preserves source links and assigns each eligible source's read
# support once. No IGH records produces an empty mapping with the same columns.
# Common-QC AIRR stays unchanged. assign_scope_donor_clones.R assigns clones
# later, using the scope threshold and these prepared representative links.

suppressPackageStartupMessages(library(optparse))

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "igh_clones",
  "igh_grouping", "qc_audit"))
stage_log <- start_stage_log("prepare_igh_grouping")

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "Common QC-pass AIRR table for one library. Required."),
  make_option("--library_id", type = "character",
              help = "Library the records belong to. Required."),
  make_option("--max_n", type = "integer",
              help = "Most N or ? symbols a grouping-eligible junction may carry. Required."),
  make_option("--grouping_file", type = "character",
              help = "Path to write the IGH grouping map to. Required."),
  make_option("--qc_audit_file", type = "character",
              help = "Path to write the grouping population audit to. Required.")
))
library_id <- script_options$library_id

records <- read_pipeline_table(script_options$input_file)
records <- records[records$locus == "IGH", , drop = FALSE]
mapping <- prepare_igh_grouping(records, library_id, script_options$max_n)
write_pipeline_table(mapping, script_options$grouping_file, GROUPING_COLUMNS)
write_qc_audit(qc_audit_grouping_rows(library_id, records, mapping, "IGH"),
               script_options$qc_audit_file)
log_info("Prepared %d IGH records, %d eligible, %d representatives with %s support",
         nrow(mapping), sum(mapping$grouping_eligible),
         sum(mapping$is_inference_representative), sum(mapping$inference_consensus_count))
finish_stage_log(stage_log)
