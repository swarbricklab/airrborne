#!/usr/bin/env Rscript
# Map native TRUST4 contig vote tallies to each common-QC IGH record in IMGT V
# coordinates. Check mapped bases against the observed sequence and retain
# missing values for unmapped positions. Reads the library's common-QC AIRR
# table and TRUST4 final.out, and writes the IGH consensus-vote table that
# build_consensus_vote_rows() in workflow/lib/trust4_consensus_votes.R builds,
# including the junction votes located in each contig.

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
  "..", "lib", "runtime.R"))
load_pipeline_modules(c(
  "table_io", "gene_calls", "assembly_backbone",
  "trust4_consensus_votes"))
stage_log <- start_stage_log("build_igh_consensus_votes")

suppressPackageStartupMessages({
  library(optparse)
})

script_options <- parse_script_arguments(list(
  make_option("--input_file", type = "character",
              help = "Common-QC all-locus AIRR table for one library. Required."),
  make_option("--final_out_file", type = "character",
              help = paste("TRUST4's final.out for the same library, holding",
                           "the per-position read votes. Required.")),
  make_option("--output_file", type = "character",
              help = "Path to write the IGH consensus-vote table to. Required."),
  make_option("--library_id", type = "character",
              help = "Library identifier every IGH row must carry. Required.")
))
input_file <- script_options$input_file
final_out_file <- script_options$final_out_file
output_file <- script_options$output_file
library_id <- script_options$library_id

db <- read_pipeline_table(input_file)
require_columns(db, c("sequence_id", "library_id", "locus", "sequence",
                      "sequence_alignment", "v_cigar", "v_sequence_start",
                      "v_sequence_end", "junction", "rev_comp"), "common AIRR")
db <- db[as.character(db$locus) == "IGH", , drop = FALSE]
if (nrow(db) && anyDuplicated(db$sequence_id)) stop("IGH sequence_id is not unique")
if (nrow(db) && any(as.character(db$library_id) != library_id)) {
  stop("common AIRR contains an IGH library_id other than ", library_id)
}
log_info("Read %d QC-pass IGH record(s) from %s", nrow(db), input_file)

backbones <- trust4_assembly_backbone_id(db$sequence_id)
votes <- read_backbone_votes(final_out_file, backbones)

rows <- build_consensus_vote_rows(db, votes, IMGT_V_POSITIONS)
if (nrow(rows) != nrow(db)) {
  stop("built ", nrow(rows), " vote row(s) for ", nrow(db), " IGH record(s)")
}

# Not write_airr_table(): these are vote tallies, not AIRR fields, and that
# writer would also rewrite their booleans.
write_pipeline_table(rows, output_file)
log_info("Wrote %d IGH consensus-vote row(s) to %s", nrow(rows), output_file)

if (nrow(rows)) {
  mapped <- rows$imgt_v_positions_mapped
  log_info(
    paste0("IMGT V positions with votes: median %s, min %d, max %d of %d; ",
           "records with a base rewritten by the substituted CDR3: %d, of which ",
           "the substitution reached 5' of the conserved Cys104 codon at IMGT ",
           "%d: %d"),
    stats::median(mapped), min(mapped), max(mapped), IMGT_V_POSITIONS,
    sum(rows$imgt_v_positions_consensus_mismatch > 0L, na.rm = TRUE),
    IMGT_CYS104_FIRST_POSITION,
    sum(!is.na(rows$imgt_v_first_consensus_mismatch_position) &
          rows$imgt_v_first_consensus_mismatch_position <
            IMGT_CYS104_FIRST_POSITION)
  )
  # Records IgBLAST found on the opposite strand, which Change-O then reports
  # reverse complemented. They share no position with their contig, so they carry
  # no votes at all; the count is reported next to the AIRR flag it tracks.
  opposite <- sum(rows$junction_vote_status == "opposite_contig_orientation")
  if (opposite) {
    log_info(paste0("Records emptied because their sequence is the reverse ",
                    "complement of their assembly contig: %d; records reporting ",
                    "AIRR rev_comp true: %d"),
             opposite, sum(parse_airr_boolean(db$rev_comp) == "true"))
  }
  status_counts <- table(rows$junction_vote_status)
  log_info("Junction vote status: %s",
           paste(names(status_counts), status_counts, sep = "=", collapse = ", "))
}

finish_stage_log(stage_log)
