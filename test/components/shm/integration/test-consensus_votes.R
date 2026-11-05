# Run the native vote handoff and retain its full IMGT V output width.
# Run: python test/run_tests.py --component shm --level integration

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "trust4_consensus_votes"))

work_dir <- tempfile("trust4-votes-")
dir.create(work_dir)

# --------------------------------------------------------------- fixtures ---

# Put each position's votes on its consensus base in a four-row native matrix.
vote_matrix_for <- function(consensus, consensus_votes) {
  bases <- strsplit(consensus, "", fixed = TRUE)[[1]]
  n <- length(bases)
  m <- matrix(0L, nrow = 4L, ncol = n, dimnames = list(VOTE_BASES, NULL))
  row_of <- match(bases, VOTE_BASES)
  for (j in seq_len(n)) {
    if (!is.na(row_of[j])) m[row_of[j], j] <- as.integer(consensus_votes[j])
  }
  m
}

# Write count lines exactly as TRUST4 does, with a trailing space before the newline.
count_lines_from_matrix <- function(m) {
  vapply(seq_len(nrow(m)),
         function(k) paste0(paste(m[k, ], collapse = " "), " "),
         character(1))
}

# A whole `_final.out` file. Each record is the six lines TRUST4 writes.
write_final_out <- function(path, records) {
  lines <- unlist(lapply(records, function(record) {
    c(paste0(">", record$id, " IGHV1-8*01"),
      record$consensus,
      count_lines_from_matrix(record$matrix))
  }), use.names = FALSE)
  writeLines(lines, path)
}

# One QC-pass IGH record, carrying every column build_consensus_vote_rows needs
# plus the rest of the reannotated schema so the file reads as AIRR.
airr_row <- function(sequence_id, sequence, sequence_alignment, v_cigar,
                     v_sequence_start, v_sequence_end, junction = "",
                     rev_comp = "F", library_id = "GUSS_TEST") {
  data.frame(
    sequence_id = sequence_id, library_id = library_id, locus = "IGH",
    sequence = sequence, rev_comp = rev_comp,
    v_call = "IGHV1-8*01", d_call = "IGHD3-10*01", j_call = "IGHJ5*02",
    c_call = "IGHM", c_informative_positions = 300L,
    junction = junction, junction_aa = "", cdr3 = "", cdr3_aa = "",
    sequence_alignment = sequence_alignment, germline_alignment = "",
    germline_alignment_d_mask = "",
    v_cigar = v_cigar, d_cigar = "", j_cigar = "",
    v_sequence_start = v_sequence_start, v_sequence_end = v_sequence_end,
    v_germline_start = 1L, v_germline_end = 320L,
    v_identity = 0.98, j_identity = 0.99, productive = "T", vj_in_frame = "T",
    stop_codon = "F", junction_length = nchar(junction), consensus_count = 12,
    stringsAsFactors = FALSE
  )
}

write_airr <- function(path, rows) {
  write.table(rows, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

run_script <- function(airr_path, final_out_path, out_path) {
  run_pipeline_script("build_igh_consensus_votes.R", c(
    paste0("--input_file=", airr_path),
    paste0("--final_out_file=", final_out_path),
    paste0("--output_file=", out_path),
    "--library_id=GUSS_TEST"
  ))
}

# --- The script wires the layers together ----------------------------------

# Real IMGT V positions, so the emitted vectors carry the released token count.
wide_consensus <- paste(rep("ACGT", 100L), collapse = "")
wide <- list(id = "assemble0", consensus = wide_consensus,
             matrix = vote_matrix_for(wide_consensus,
                                      rep(7L, nchar(wide_consensus))))
wide_final_out <- file.path(work_dir, "wide_final.out")
write_final_out(wide_final_out, list(wide))
wide_airr <- file.path(work_dir, "wide_airr.tsv")
write_airr(wide_airr, airr_row("assemble0_0", wide_consensus, wide_consensus,
                               paste0(nchar(wide_consensus), "="), 1L,
                               nchar(wide_consensus)))
wide_out <- file.path(work_dir, "wide_votes.tsv")
wide_run <- run_script(wide_airr, wide_final_out, wide_out)
expect_script_ok(wide_run)

test_that("the released vote table carries the full IMGT V width", {
  wide_table <- read.delim(wide_out, stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(names(wide_table), c(
    "library_id", "sequence_id", "trust4_backbone_id", "consensus_length",
    "imgt_v_positions_mapped", "imgt_v_positions_consensus_mismatch",
    "imgt_v_first_consensus_mismatch_position", "total_votes_by_imgt_position", "non_consensus_votes_by_imgt_position",
    "junction_start_in_consensus", "total_votes_by_junction_position",
    "non_consensus_votes_by_junction_position", "junction_vote_status"
  ))
  expect_identical(nrow(wide_table), 1L)
  expect_identical(wide_table$imgt_v_positions_mapped, IMGT_V_POSITIONS)
  expect_identical(length(parse_vote_vector(wide_table$total_votes_by_imgt_position,
                                            IMGT_V_POSITIONS)),
                   IMGT_V_POSITIONS)
  expect_identical(wide_table$junction_vote_status, "junction_missing")
})

# An empty IGH set still writes a header, so downstream reads do not special-case it.
empty_airr <- file.path(work_dir, "empty_airr.tsv")
empty_rows <- airr_row("assemble0_0", "ACGT", "ACGT", "4=", 1L, 4L)
empty_rows$locus <- "TRB"
write_airr(empty_airr, empty_rows)
empty_out <- file.path(work_dir, "empty_votes.tsv")
empty_run <- run_script(empty_airr, wide_final_out, empty_out)
expect_script_ok(empty_run)

test_that("an empty IGH set still writes a header", {
  empty_table <- read.delim(empty_out, stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(nrow(empty_table), 0L)
  expect_identical(ncol(empty_table), 13L)
})
