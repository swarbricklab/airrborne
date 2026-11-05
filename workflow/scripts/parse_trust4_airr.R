#!/usr/bin/env Rscript
# Preserve native TRUST4 observations after removing imputed junctions.
# Native means the rows TRUST4 reports before this pipeline reannotates them.
# One assembly can report alternative CDR3 sequences; their IDs and support
# remain separate here, including when full sequences are identical. Later
# grouping selects representatives through prepare_assembly_mapping() in
# workflow/lib/assembly_backbone.R.
#
# This script reads one library's AIRR table and cdr3.out, writes the retained
# AIRR rows and counts native yield, imputation losses and locus routing in the
# QC audit. flag_imputed_cdr3_records() in workflow/lib/trust4_parsing.R owns
# the fixed TRUST4 score rule. IgBLAST reannotation and common QC run afterwards.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "airr_io", "qc_audit", "trust4_parsing"))
stage_log <- start_stage_log("parse_trust4_airr")

suppressPackageStartupMessages(library(optparse))

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "TRUST4's own AIRR table for one library. Required."),
  make_option("--output_file", type = "character",
              help = "Path to write the normalised AIRR table to. Required."),
  make_option("--cdr3_file", type = "character",
              help = paste("TRUST4's cdr3.out for the same library; its column",
                           "10 carries the per-CDR3 score. Required.")),
  make_option("--qc_audit_file", type = "character",
              help = "Path to write the parser QC audit to. Required."),
  make_option("--library_id", type = "character",
              help = "Library identifier the qc_audit rows are keyed on. Required.")
))
input_file <- script_options$input_file
output_file <- script_options$output_file
cdr3_file <- script_options$cdr3_file
qc_audit_file <- script_options$qc_audit_file
library_id <- script_options$library_id

# The audit keeps the empty library's zero counts as well as record losses.
# Invalid input stops the script before the outputs are written.
qc_audit <- qc_audit_empty()

native <- read_pipeline_table(
  input_file,
  required_columns = c("sequence_id", "sequence", "locus", "consensus_count"),
  context = "TRUST4 AIRR"
)
if (any(is_airr_null(native$sequence_id))) stop("TRUST4 sequence_id is missing")
if (anyDuplicated(native$sequence_id)) stop("TRUST4 sequence_id is not unique")
if (any(is_airr_null(native$sequence))) stop("TRUST4 sequence is missing")
if (any(is_airr_null(native$locus))) stop("TRUST4 locus is missing")

# TRUST4 can complete a CDR3 with inferred bases. Remove those records before
# IgBLAST could treat the completed junction as an observed sequence.
cdr3 <- read_trust4_cdr3_scores(cdr3_file)
imputed <- flag_imputed_cdr3_records(native$sequence_id, cdr3$cdr3_key, cdr3$cdr3_score)
if (any(imputed)) {
  by_locus <- table(native$locus[imputed])
  log_info(
    "Dropping %d imputed-CDR3 record(s) (TRUST4 score %s): %s",
    sum(imputed), TRUST4_IMPUTED_CDR3_SCORE,
    paste(sprintf("%s=%d", names(by_locus), as.integer(by_locus)), collapse = ", "))
}

# Read support is coerced before anything is removed, because the qc_audit needs
# the pre-drop denominators.
support <- as_read_support(native$consensus_count, "consensus_count")

# Head of the QC audit. The two assembly counts come from cdr3.out because it is
# the assembly-level file this stage reads. Contigs TRUST4 assembled without any
# CDR3 call do not appear in it, and rows in cdr3.out that TRUST4's own simplerep
# stage did not report as clonotypes are why the AIRR row count can be below the
# call count.
qc_audit <- rbind(
  qc_audit,
  qc_audit_rows(library_id, "trust4_contigs_with_cdr3_call",
                length(unique(cdr3$consensus_id)), locus = "all",
                analysis = "assembly", stage_index = 1, input_unit = "contig",
                locus_basis = "not_applicable", support_conserved = NA),
  qc_audit_rows(library_id, "trust4_cdr3_calls", nrow(cdr3), locus = "all",
                analysis = "assembly", stage_index = 2, input_unit = "cdr3_call",
                locus_basis = "not_applicable", support_conserved = NA),
  # A level, not a gate: the mask keeps every record, so the counts in and out
  # agree. The locus sums cover the seven canonical loci, so a native record on
  # any other locus first appears in the routing row below.
  qc_audit_gate_rows(library_id, "trust4_reported_airr",
                   rep(TRUE, nrow(native)), support, native$locus,
                   analysis = "common", stage_index = 3, locus_basis = "native"),
  qc_audit_gate_rows(library_id, "imputed_cdr3_filter", !imputed, support, native$locus,
                   analysis = "common", stage_index = 4, locus_basis = "native",
                   kind = "eligibility", support_conserved = FALSE)
)

native <- native[!imputed, , drop = FALSE]
support <- support[!imputed]

# Identical sequences must carry one TRUST4 locus. Records are not merged here.
first_sequence <- match(native$sequence, native$sequence)
if (any(native$locus != native$locus[first_sequence])) {
  stop("Identical TRUST4 sequence has conflicting locus assignments")
}

# Keep assembly IDs and original support for later representative selection and
# position-vote measurements. This parser neither chooses a dominant CDR3
# alternative nor transfers support between records.
parsed <- native
parsed_support <- support
write_airr_table(parsed, output_file)

# The routing row counts the parsed records against the subset IgBLAST is given, because
# run_igblast_all_loci.sh writes a FASTA only for the seven canonical loci. The
# difference between the two is what IgBLAST never sees.
qc_audit <- rbind(
  qc_audit,
  qc_audit_rows(library_id, "igblast_locus_routing", nrow(parsed),
                sum(parsed$locus %in% QC_AUDIT_LOCI), sum(parsed_support),
                sum(parsed_support[parsed$locus %in% QC_AUDIT_LOCI]), locus = "all",
                analysis = "common", stage_index = 6, locus_basis = "native",
                kind = "eligibility", support_conserved = FALSE)
)
write_qc_audit(qc_audit, qc_audit_file)

log_info(
  "TRUST4 parsing: preserved %d non-imputed source rows in %s",
  nrow(parsed), output_file
)

finish_stage_log(stage_log)
