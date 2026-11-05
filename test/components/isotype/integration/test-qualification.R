# Run call qualification without deleting common-QC records or their support.
# Run: python test/run_tests.py --component isotype --level integration

load_test_modules(c("table_io", "gene_calls", "reannotated_airr"))

test_that("query evidence masks calls while common-QC support remains available", {
  directory <- tempfile("isotype-boundary-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- function(name) file.path(directory, paste0(name, ".tsv"))

  # Each aligned subject has 41 known bases. The queries deliberately contain
  # 39, 40 and 41 known bases; the final record has a grouping-ineligible junction.
  ids <- paste0("assemble", 1:4, "_0")
  queries <- c(paste0(strrep("A", 39), "N-"),
               paste0(strrep("A", 40), "N"), rep(strrep("A", 41), 2))
  subjects <- rep(strrep("A", 41), 4)
  known_positions <- function(sequence) {
    vapply(strsplit(sequence, "", fixed = TRUE),
      function(bases) sum(bases %in% c("A", "C", "G", "T")), integer(1))
  }
  expect_equal(nchar(queries), rep(41L, 4))
  expect_identical(known_positions(queries), c(39L, 40L, 41L, 41L))
  expect_identical(known_positions(subjects), rep(41L, 4))
  calls <- c("IGHG1*01", "IGHG1*01", "IGHA1*01", "IGHM*01")
  fmt7 <- vapply(seq_along(ids), function(i) paste(c(
    "C", ids[i], calls[i], c("95.122", "97.561", "100", "100")[i], "41",
    c(1, 1, 0, 0)[i], c(1, 0, 0, 0)[i], c(1, 0, 0, 0)[i],
    "1", c(40, 41, 41, 41)[i], "1", "41", "1e-20", "100",
    queries[i], subjects[i], c("39NA-A", "40NA", "41", "41")[i]
  ), collapse = "\t"), character(1))
  writeLines(fmt7, path("c_alignment"))
  evidence <- read_optional_c_evidence(path("c_alignment"), "Synthetic C alignment")
  expect_identical(evidence$sequence_id, ids)
  expect_identical(evidence$c_informative_positions, c(39L, 40L, 41L, 41L))

  # Distinct constant-region evidence accompanies distinct complete sequences.
  sequences <- paste0(c(rep("AAATGTGCTTGGAAA", 3), "AAATGTNCTTGGAAA"), queries)
  records <- data.frame(
    sequence_id = ids, library_id = "0040", locus = "IGH",
    sequence = sequences, sequence_alignment = sequences,
    germline_alignment_d_mask = "AAATGTGCTTGGAAA",
    junction = c(rep("TGTGCTTGG", 3), "TGTNCTTGG"),
    junction_aa = "CAW", junction_length = 9L,
    v_call = "IGHV1-2*01", d_call = "", j_call = "IGHJ4*01",
    v_gene = "IGHV1-2", j_gene = "IGHJ4",
    v_gene_status = "resolved", j_gene_status = "resolved",
    c_call = calls, c_informative_positions = evidence$c_informative_positions,
    productive = "T", consensus_count = c(3, 5, 7, 11),
    stringsAsFactors = FALSE
  )
  write_pipeline_table(records, path("annotated"))
  expect_script_ok(run_pipeline_script("filter_airr.R", c(
    paste0("--input_file=", path("annotated")),
    paste0("--output_file=", path("common_qc")),
    paste0("--qc_audit_file=", path("common_audit")),
    "--library_id=0040", "--max_n_fraction=0.1"
  )))
  before <- readLines(path("common_qc"))
  retained <- read_pipeline_table(path("common_qc"))
  expect_identical(retained$sequence_id, ids)
  expect_equal(as.numeric(retained$consensus_count), c(3, 5, 7, 11))

  run_stage_with_options("prepare_igh_grouping.R", list(
    input_file = path("common_qc"), library_id = "0040", max_n = 0,
    grouping_file = path("grouping"), qc_audit_file = path("grouping_audit")))
  grouping <- read_pipeline_table(path("grouping"))
  expect_identical(grouping$grouping_eligible, c("TRUE", "TRUE", "TRUE", "FALSE"))

  run_stage_with_options("qualify_isotype_calls.R", list(
    input_file = path("common_qc"), library_id = "0040", min_c_positions = 40,
    qc_audit_file = path("call_audit"), output_file = path("calls")))
  qualified <- read_pipeline_table(path("calls"))
  expect_identical(qualified$sequence_id, ids)
  expect_identical(qualified$c_call, calls)
  expect_identical(qualified$qualified_c_call, c("", calls[2:4]))
  expect_identical(qualified$c_call_callable, c("FALSE", "TRUE", "TRUE", "TRUE"))
  # The masked call has no isotype. The others are named at both levels.
  expect_identical(qualified$isotype_class, c("", "IgG", "IgA", "IgM"))
  expect_identical(qualified$isotype_subclass, c("", "IgG1", "IgA1", "IgM"))

  run_stage_with_options("calculate_isotype.R", list(
    input_file = path("common_qc"), isotype_calls_file = path("calls"),
    library_id = "0040", qc_audit_file = path("summary_audit"),
    isotype_file = path("composition"), summary_file = path("summary")))
  composition <- read.delim(path("composition"), stringsAsFactors = FALSE)
  summary <- read.delim(path("summary"), stringsAsFactors = FALSE)
  expect_identical(names(composition), c("library_id", "isotype_level", "isotype",
                                          "consensus_count", "callable_consensus_count",
                                          "support_fraction"))
  classes <- composition[composition$isotype_level == "class", ]
  expect_identical(classes$isotype, c("IgM", "IgD", "IgG", "IgA", "IgE"))
  expect_equal(classes$consensus_count, c(11, 0, 5, 7, 0))
  expect_equal(classes$support_fraction, c(11, 0, 5, 7, 0) / 23)
  # Each record's sidecar class is the class whose support counts it.
  class_support <- vapply(classes$isotype, function(label) {
    sum(as.numeric(retained$consensus_count)[qualified$isotype_class == label])
  }, numeric(1), USE.NAMES = FALSE)
  expect_equal(class_support, classes$consensus_count)
  expect_equal(summary$common_qc_consensus_count, 26)
  expect_equal(summary$class_callable_consensus_count, 23)
  expect_equal(summary$class_uncallable_consensus_count, 3)
  expect_equal(summary$subclass_callable_consensus_count, 23)
  expect_equal(summary$log2_igg_iga, log2(5 / 7))
  expect_identical(summary$igg_iga_status, "available")

  audit <- read.delim(path("call_audit"), stringsAsFactors = FALSE)
  callable <- audit[audit$stage == "c_call_callable", ]
  expect_equal(callable$n_units_in, 4L)
  expect_equal(callable$n_units_out, 3L)
  expect_equal(callable$consensus_count_in, 26)
  expect_equal(callable$consensus_count_out, 23)
  expect_identical(readLines(path("common_qc")), before)
})
