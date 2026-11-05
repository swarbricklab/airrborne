#!/usr/bin/env Rscript
# Combine IgBLAST/MakeDb annotations into one table for common record QC.
# run_igblast_all_loci.sh supplies the IG and TR tables, reconstructed IG
# germlines and C-region alignments. TRUST4 rows supply IDs and read support.
# Missing germline or C-alignment evidence stays missing; it does not remove
# an otherwise annotated row here.
#
# workflow/lib/reannotated_airr.R owns shared TRA/TRD locus resolution, CIGAR
# conversion and candidate-list completeness. Its build_reannotated_airr_frame()
# keeps raw V/J allele lists beside assignments from gene_calls.R. Multiple
# named genes resolve together only when their named alleles have identical,
# complete, unambiguous reference sequences. Equally best matches that reach the
# alignment cap prevent resolution, even if the visible call names one gene. The
# QC audit records annotation losses; filter_airr.R applies common QC next.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "airr_io", "gene_calls", "reannotated_airr",
  "qc_audit"))
stage_log <- start_stage_log("build_reannotated_airr")

suppressPackageStartupMessages({
  library(airr)
  library(optparse)
})

# The five IgBLAST/Change-O products are each given as the empty string when the
# step that would have written them had nothing to run on, which is why they are
# required arguments that accept an empty value rather than optional ones.
script_options <- parse_script_arguments(list(
  make_option("--ig_db_file", type = "character",
              help = paste("MakeDb-pass table for the IG loci, or the empty",
                           "string when IgBLAST annotated no IG record.",
                           "Required.")),
  make_option("--tr_db_file", type = "character",
              help = paste("MakeDb-pass table for the TR loci, or the empty",
                           "string when IgBLAST annotated no TR record.",
                           "Required.")),
  make_option("--ig_germ_file", type = "character",
              help = paste("CreateGermlines germline-pass table for the IG",
                           "loci, or the empty string when no germline was",
                           "reconstructed. Required.")),
  make_option("--ig_fmt7_file", type = "character",
              help = paste("IgBLAST fmt7 hit table for the IG loci, read for",
                           "the first best-ranked C alignment, or the empty",
                           "string when IgBLAST was not run. Required.")),
  make_option("--tr_fmt7_file", type = "character",
              help = paste("IgBLAST fmt7 hit table for the TR loci, read for",
                           "the first best-ranked C alignment, or the empty",
                           "string when IgBLAST was not run. Required.")),
  make_option("--trust4_file", type = "character",
              help = paste("Normalised TRUST4 AIRR table, which supplies the",
                           "sequence identifiers and read support. Required.")),
  make_option("--reference_metadata_file", type = "character",
              help = "Original reference allele metadata and sequence classes. Required."),
  make_option("--num_alignments_v", type = "integer", help = "Configured V alignment cap. Required."),
  make_option("--num_alignments_d", type = "integer", help = "Configured D alignment cap. Required."),
  make_option("--num_alignments_j", type = "integer", help = "Configured J alignment cap. Required."),
  make_option("--num_alignments_c", type = "integer", help = "Configured C alignment cap. Required."),
  make_option("--output_file", type = "character",
              help = "Path to write the combined all-locus AIRR table to. Required."),
  make_option("--library_id", type = "character",
              help = "Library identifier stamped on every output row. Required."),
  make_option("--qc_audit_in_file", type = "character",
              help = paste("QC audit shard parse_trust4_airr.R wrote, which",
                           "this stage extends. Required.")),
  make_option("--qc_audit_file", type = "character",
              help = "Path to write the extended QC audit shard to. Required.")
))
ig_db_file <- script_options$ig_db_file
tr_db_file <- script_options$tr_db_file
ig_germ_file <- script_options$ig_germ_file
ig_fmt7_file <- script_options$ig_fmt7_file
tr_fmt7_file <- script_options$tr_fmt7_file
trust4_file <- script_options$trust4_file
reference_metadata_file <- script_options$reference_metadata_file
output_file <- script_options$output_file
library_id <- script_options$library_id
qc_audit_in_file <- script_options$qc_audit_in_file
qc_audit_file <- script_options$qc_audit_file


native <- read_pipeline_table(
  trust4_file, required_columns = c("sequence_id", "consensus_count", "locus"),
  context = "TRUST4 AIRR"
)
if (anyDuplicated(native$sequence_id)) stop("TRUST4 sequence_id is not unique")
support <- as_read_support(native$consensus_count, "consensus_count")

db <- bind_frames_on_column_union(list(
  read_optional_airr_table(ig_db_file, "IG MakeDb-pass"),
  read_optional_airr_table(tr_db_file, "TR MakeDb-pass")
))

# Which routed sequences IgBLAST and MakeDb could annotate, keyed on the TRUST4
# locus they were routed under. Records whose locus is outside the canonical
# seven were never handed to IgBLAST and are already accounted for by the
# routing row, so keying on TRUST4's locus here excludes them automatically.
annotated <- as.character(native$sequence_id) %in% as.character(db$sequence_id)
qc_audit <- rbind(
  read_qc_audit(qc_audit_in_file),
  qc_audit_gate_rows(library_id, "igblast_makedb",
                   annotated, support, native$locus,
                   analysis = "common", stage_index = 7, locus_basis = "native",
                   kind = "eligibility", support_conserved = FALSE)
)

if (!nrow(db)) {
  qc_audit <- rbind(qc_audit, qc_audit_gate_rows(
    library_id, "igblast_locus_agreement", logical(), numeric(), character(),
    analysis = "common", stage_index = 8, kind = "diagnostic", support_conserved = FALSE))
  write_airr_table(make_empty_frame(REANNOTATED_COLUMNS), output_file)
  log_info("Wrote 0 all-locus MakeDb-pass records to %s", output_file)
  write_qc_audit(qc_audit, qc_audit_file)
  finish_stage_log(stage_log)
  quit(save = "no", status = 0)
}
if (anyDuplicated(db$sequence_id)) stop("MakeDb-pass sequence_id is not unique")

c_evidence <- bind_frames_on_column_union(list(
  read_optional_c_evidence(ig_fmt7_file, "IG IgBLAST fmt7"),
  read_optional_c_evidence(tr_fmt7_file, "TR IgBLAST fmt7")
))
if (nrow(c_evidence) && anyDuplicated(c_evidence$sequence_id)) {
  stop("IgBLAST C evidence sequence_id is not unique across receptor families")
}

support_index <- match(as.character(db$sequence_id), as.character(native$sequence_id))
if (anyNA(support_index)) {
  stop("MakeDb-pass contains sequence_id absent from TRUST4 AIRR")
}

# Interpret shared TRA/TRD V hits before auditing or releasing the locus.
db <- resolve_shared_trav_dv_locus(db)

# Every QC audit stage from here on uses the resolved annotation locus. This row
# connects the two bases: n_units_in is what reannotation called each locus,
# n_units_out is the subset TRUST4 had called the same, so the difference is
# what moved into that locus.
qc_audit <- rbind(
  qc_audit,
  qc_audit_gate_rows(
    library_id, "igblast_locus_agreement",
    as.character(db$locus) == as.character(native$locus[support_index]),
    support[support_index], db$locus,
    analysis = "common", stage_index = 8, kind = "diagnostic", support_conserved = FALSE)
)

ig_germ <- read_optional_airr_table(ig_germ_file, "IG germline-pass")
if (nrow(ig_germ) && anyDuplicated(ig_germ$sequence_id)) {
  stop("IG germline-pass sequence_id is not unique")
}
with_germline <- attach_d_masked_germline(db, ig_germ)
with_evidence <- attach_constant_region_evidence(with_germline$records, c_evidence)
if (with_evidence$n_call_without_evidence > 0L) {
  log_info(paste0("MakeDb c_call lacks first best-ranked C-alignment evidence ",
                  "for %d record(s); leaving c_informative_positions NA so the ",
                  "isotype-call qualification stage masks these calls"),
           with_evidence$n_call_without_evidence)
}

# Best-score ties must fit below params.igblast.num_alignments_* before the
# visible call can establish a complete candidate list. This flag accompanies
# the raw allele list; classification never silently chooses its first gene.
counts <- bind_frames_on_column_union(lapply(c(ig_fmt7_file, tr_fmt7_file), function(path) {
  parse_igblast_candidate_counts(if (nzchar(path)) readLines(path, warn = FALSE) else character())
}))
caps <- lapply(c("v", "d", "j", "c"), function(segment) script_options[[paste0("num_alignments_", segment)]])
names(caps) <- c("v", "d", "j", "c")
with_evidence$records <- attach_candidate_completeness(with_evidence$records, counts, caps)
identity_classes <- read_germline_identity_classes(reference_metadata_file)
log_info("Read %d germline allele(s) from %s, in %d identity class(es)",
         length(identity_classes), reference_metadata_file,
         length(unique(identity_classes)))

reannotated <- build_reannotated_airr_frame(with_evidence$records, library_id,
                                            support[support_index],
                                            identity_classes)
write_airr_table(reannotated, output_file)
log_info("Wrote %d all-locus MakeDb-pass records to %s", nrow(reannotated), output_file)
write_qc_audit(qc_audit, qc_audit_file)
log_info("IG germlines attached: %d/%d", with_germline$n_attached, nrow(reannotated))
log_info("First best-ranked C-alignment evidence attached: %d/%d",
         with_evidence$n_attached, nrow(reannotated))

finish_stage_log(stage_log)
