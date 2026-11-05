# Check complete-reference V/J entities and candidate-list completeness.
# Run: python test/run_tests.py --component annotation --level unit

load_test_modules(c(
  "table_io", "gene_calls", "isotype_calls", "reannotated_airr", "rearrangement_keys"))

identity_fixture <- function() {
  rows <- data.frame(
    allele = c("IGHV1-2*01", "IGHV1-2*02", "IGHV1-3*01", "IGHV1-3*02", "IGHV1-2D*01"),
    gene = c("IGHV1-2", "IGHV1-2", "IGHV1-3", "IGHV1-3", "IGHV1-2D"),
    region = "V-REGION", functionality = c("F", "F", "ORF", "P", "F"),
    completeness = c("complete", "complete", "complete", "partial in 5'", "complete"),
    identity_id = c("VSEQ_A", "VSEQ_C", "VSEQ_A", "", "VSEQ_C"),
    identity_label = c("IGHV1-2,IGHV1-3", "IGHV1-2,IGHV1-2D", "IGHV1-2,IGHV1-3", "", "IGHV1-2,IGHV1-2D"),
    orphon = "F", stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".tsv")
  write_pipeline_table(rows, path)
  read_germline_identity_classes(path)
}

test_that("named alleles decide multi-gene identity before gene reduction", {
  classes <- identity_fixture()
  calls <- c("IGHV1-2*01,IGHV1-3*01", "IGHV1-2*02,IGHV1-3*01",
             "IGHV1-2*01,IGHV1-3*02", "IGHV1-2*01,IGHV1-3*99",
             "IGHV1-2*01,IGHV1-2*02")
  expect_identical(vj_tie_kind(calls, classes),
    c("assigned_group", "ambiguous", "ambiguous", "ambiguous", "resolved"))
  expect_identical(vj_canonical_token(calls, classes), c("VSEQ_A", "", "", "", "IGHV1-2"))
  expect_identical(vj_entity_label(calls[1], classes), "IGHV1-2,IGHV1-3")
})

test_that("composites remain separate from members and overlapping composites", {
  classes <- identity_fixture()
  calls <- c("IGHV1-2*01", "IGHV1-2*01,IGHV1-3*01", "IGHV1-2*02,IGHV1-2D*01", "IGHV1-2D*01")
  tokens <- vj_canonical_token(calls, classes)
  expect_identical(tokens, c("IGHV1-2", "VSEQ_A", "VSEQ_C", "IGHV1-2D"))
  data <- data.frame(v_gene = tokens, j_gene = "IGHJ1", v_call = calls)
  tool <- igh_inference_calls(data)
  expect_length(unique(tool$v_inference_call), 4L)
  expect_identical(tool$v_call, calls)
  expect_identical(exact_rearrangement_group_key(rep("ACGT", 4), tokens, rep("IGHJ1", 4)),
    paste("ACGT", tokens, "IGHJ1", sep = .GROUP_SEP))
})

test_that("incomplete lists cannot establish single-gene or composite resolution", {
  classes <- identity_fixture()
  calls <- c("IGHV1-2*01", "IGHV1-2*01,IGHV1-3*01")
  expect_identical(vj_tie_kind(calls, classes, FALSE), rep("incomplete_candidates", 2))
  expect_identical(vj_canonical_token(calls, classes, FALSE), c("", ""))
  expect_true(all(vj_calls_ambiguous(vj_tie_kind(calls, classes, FALSE), "resolved")))
})

test_that("reference functionality retains exact allele and mixed-status evidence", {
  classes <- identity_fixture()
  expect_identical(reference_call_functionality("IGHV1-3*01,IGHV1-3*02", classes),
    "IGHV1-3*01:ORF;IGHV1-3*02:P")
  expect_identical(reference_call_functionality("IGHV1-3*99", classes), "IGHV1-3*99:unknown")
  expect_true(is.na(reference_call_functionality("", classes)))
  expect_identical(reference_call_orphon(c("IGKV1/OR2-0*01", "IGHV1-2*01", "")), c(TRUE, FALSE, NA))
})

test_that("candidate completeness counts best-score ties and preserves segment caps", {
  hit <- function(segment, allele, score) paste(c(segment, "seq1", allele, "100", "8", "0", "0", "0", "1", "8", "1", "8", "0", score, "ACGTACGT", "ACGTACGT", "8"), collapse = "\t")
  counts <- parse_igblast_candidate_counts(c(hit("V", "IGHV1-2*01", 50), hit("V", "IGHV1-3*01", 50),
    hit("V", "IGHV1-4*01", 40), hit("J", "IGHJ1*01", 20)))
  records <- data.frame(sequence_id = "seq1", v_call = "IGHV1-2*01,IGHV1-3*01",
                        d_call = "", j_call = "IGHJ1*01", c_call = "")
  result <- attach_candidate_completeness(records, counts, list(v = 2, d = 3, j = 2, c = 125))
  expect_false(result$v_call_candidates_complete)
  expect_true(result$j_call_candidates_complete)
  result <- attach_candidate_completeness(records, counts, list(v = 3, d = 3, j = 2, c = 125))
  expect_true(result$v_call_candidates_complete)
})

test_that("raw capped calls remain incomplete when hit evidence is absent", {
  records <- data.frame(sequence_id = "seq1", v_call = "IGHV1-2*01,IGHV1-3*01",
                        d_call = "", j_call = "IGHJ1*01", c_call = "")
  counts <- parse_igblast_candidate_counts(character())
  result <- attach_candidate_completeness(records, counts, list(v = 2, d = 3, j = 20, c = 125))
  expect_false(result$v_call_candidates_complete)
})

test_that("an incomplete C candidate list cannot qualify an isotype call", {
  flags <- c_call_qc_flags("IGHG1*01", "IGH", 100, 20, candidates_complete = FALSE)
  expect_true(flags$enough_evidence)
  expect_false(flags$callable)
})
