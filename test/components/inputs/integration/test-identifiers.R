# Run literal identifiers through QC, consensus votes and SHM readers.
# Run: python test/run_tests.py --component inputs --level integration

load_test_modules(c("table_io", "airr_io"))

test_that("library identifiers survive common QC, consensus votes and SHM", {
  directory <- tempfile("library-identifiers-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- function(name) file.path(directory, paste0(name, ".tsv"))
  germline <- paste0(strrep("GCT", 104), "TGTGCTTGG")
  sequence <- germline
  substr(sequence, 4, 4) <- "A"
  bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  writeLines(c(">assemble0", sequence,
    vapply(c("A", "C", "G", "T"), function(base) {
      paste(ifelse(bases == base, 7, 0), collapse = " ")
    }, character(1))), path("final_out"))
  record <- data.frame(sequence_id = "assemble0_0", library_id = "",
    locus = "IGH", sequence = sequence, rev_comp = FALSE,
    v_call = "IGHV1-8*01", d_call = "", j_call = "IGHJ5*02",
    v_gene_status = "resolved", j_gene_status = "resolved",
    v_gene = "IGHV1-8", j_gene = "IGHJ5",
    sequence_alignment = sequence, germline_alignment_d_mask = germline,
    junction = "TGTGCTTGG", junction_length = 9, productive = TRUE,
    vj_in_frame = TRUE, stop_codon = FALSE, consensus_count = 7,
    v_cigar = "312=", v_sequence_start = 1, v_sequence_end = 312)
  for (library in c("NA", "None", "0460", "nan")) {
    record$library_id <- library
    nonproductive <- record
    nonproductive$sequence_id <- "assemble1_0"
    nonproductive$productive <- FALSE
    write_airr_table(rbind(record, nonproductive), path("annotated"))
    filter_args <- c(paste0("--input_file=", path("annotated")),
      paste0("--output_file=", path("qc")), paste0("--library_id=", library),
      "--max_n_fraction=0.1", paste0("--qc_audit_file=", path("qc_audit")))
    expect_script_ok(run_pipeline_script("filter_airr.R", filter_args))
    retained <- read_pipeline_table(path("qc"))
    expect_identical(retained$library_id, library)
    expect_identical(retained$sequence_id, "assemble0_0")
    expect_identical(retained$productive, "T")
    expect_identical(retained$rev_comp, "F")
    expect_identical(retained$consensus_count, "7")
    vote_args <- c(paste0("--input_file=", path("qc")),
      paste0("--final_out_file=", path("final_out")),
      paste0("--output_file=", path("votes")), paste0("--library_id=", library))
    expect_script_ok(run_pipeline_script("build_igh_consensus_votes.R", vote_args))
    votes <- read_pipeline_table(path("votes"))
    expect_identical(votes$library_id, library)
    expect_identical(votes$imgt_v_positions_mapped, "312")
    write_pipeline_table(data.frame(sequence_id = "assemble0_0", library_id = library,
                                    qualified_c_call = ""), path("calls"))
    feature_args <- c(paste0("--input_file=", path("qc")),
      paste0("--votes_file=", path("votes")),
      paste0("--isotype_calls_file=", path("calls")),
      paste0("--output_file=", path("features")), paste0("--library_id=", library),
      "--min_position_votes=2", paste0("--qc_audit_file=", path("shm_audit")))
    expect_script_ok(run_pipeline_script("calculate_igh_features.R", feature_args))
    features <- read_pipeline_table(path("features"))
    expect_identical(features$library_id, library)
    expect_identical(features$sequence_id, "assemble0_0")
    expect_equal(as.numeric(features$shm_frequency), 1 / 312)
    expect_identical(features$informative_supported_fwr1, "78")
  }
  # Missing or mismatched identifiers must still stop all three entry points.
  for (invalid in c("", "different_library")) {
    record$library_id <- invalid
    write_airr_table(record, path("annotated"))
    write_airr_table(record, path("qc"))
    for (stage in c("filter_airr.R", "build_igh_consensus_votes.R", "calculate_igh_features.R")) {
      arguments <- switch(stage, "filter_airr.R" = filter_args,
                           "build_igh_consensus_votes.R" = vote_args,
                           "calculate_igh_features.R" = feature_args)
      run <- run_pipeline_script(stage, arguments)
      expect_false(identical(run$status, 0L))
      expect_true(any(grepl("library_id other than", run$output, fixed = TRUE)))
    }
  }
})
