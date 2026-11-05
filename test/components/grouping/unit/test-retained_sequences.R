# Check backbone choice, exact DNA merging and support transfer at all loci.
# Run: python test/run_tests.py --component grouping --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys", "igh_clones",
  "igh_grouping", "qc_audit", "rearrangement_groups"))

analysis_example <- function(locus = "IGH", sequence_id = c("assemble1_0", "assemble1_1"),
                             support = c(9, 1), junction = c("TGTGCGAAA", "TGTGCGAAT")) {
  data.frame(
    library_id = "0460", sequence_id = sequence_id, locus = locus,
    junction = junction, junction_aa = "CAK",
    v_call = paste0(locus, "V1*01"), j_call = paste0(locus, "J1*01"),
    v_gene_status = "resolved", j_gene_status = "resolved",
    v_gene = paste0(locus, "V1"), j_gene = paste0(locus, "J1"),
    consensus_count = support, v_identity = 0.99,
    sequence_alignment = paste0("ACGT", junction, "ACGT"),
    sequence = paste0("ACGT", junction, "ACGT"),
    stringsAsFactors = FALSE
  )
}

test_that("every source row remains and source support is not changed", {
  records <- analysis_example()
  before <- records
  members <- prepare_grouping_fixture(records, "0460", 0L)
  expect_identical(names(members), GROUPING_COLUMNS)
  expect_identical(members$sequence_id, records$sequence_id)
  expect_identical(members$library_id, c("0460", "0460"))
  expect_identical(members$assembly_backbone_id, c("assemble1", "assemble1"))
  expect_identical(members$grouping_eligible, c(TRUE, TRUE))
  expect_identical(members$inference_representative_sequence_id, rep("assemble1_0", 2))
  expect_identical(members$is_inference_representative, c(TRUE, FALSE))
  expect_identical(members$inference_consensus_count, c(10, 0))
  expect_identical(records, before)
  expect_false("clone_id" %in% names(members))
})

test_that("an ineligible dominant IGH variant excludes its minor eligible sibling", {
  records <- analysis_example(support = c(100, 3))
  records$v_gene_status[1] <- "ambiguous"
  records$v_call[1] <- "IGHV1*01,IGHV2*01"
  records$v_gene[1] <- ""
  members <- prepare_grouping_fixture(records, "0460", 0L)
  expect_identical(members$grouping_eligible, c(FALSE, FALSE))
  expect_identical(members$grouping_status, c("v_call_ambiguous", "ineligible_backbone_representative"))
  expect_identical(members$inference_representative_sequence_id, c(NA_character_, NA_character_))
  expect_identical(members$is_inference_representative, c(FALSE, FALSE))
  expect_identical(members$inference_consensus_count, c(0, 0))
})

test_that("an ineligible dominant light-chain variant excludes its backbone", {
  for (locus in c("IGK", "IGL")) {
    records <- analysis_example(locus, support = c(100, 3))
    records$j_gene[1] <- ""
    records$j_gene_status[1] <- "ambiguous"
    members <- prepare_grouping_fixture(records, "0460", 0L)
    expect_identical(members$grouping_status,
                     c("j_call_ambiguous", "ineligible_backbone_representative"))
    expect_identical(members$inference_representative_sequence_id, c(NA_character_, NA_character_))
    expect_identical(members$inference_consensus_count, c(0, 0))
  }
})

test_that("the configured IGH missing-base allowance changes the same gate as inference", {
  records <- analysis_example(junction = c("TGTNCGAAA", "TGTGCGAAT"))
  no_missing <- prepare_grouping_fixture(records, "0460", 0L)
  one_missing <- prepare_grouping_fixture(records, "0460", 1L)
  expect_identical(no_missing$grouping_status, c("ineligible_junction", "ineligible_backbone_representative"))
  expect_identical(no_missing$inference_consensus_count, c(0, 0))
  expect_identical(one_missing$grouping_eligible, c(TRUE, TRUE))
  expect_identical(one_missing$inference_consensus_count, c(10, 0))
  records$locus <- "TRB"
  expect_identical(prepare_grouping_fixture(records, "0460", 1L)$grouping_eligible, c(FALSE, TRUE))
})

test_that("IGH junction and ambiguous-call failures remain visible together", {
  records <- analysis_example(junction = c("TGTNCGAAA", "TGTGCGAAT"))
  records$v_gene_status[1] <- "ambiguous"
  members <- prepare_grouping_fixture(records, "0460", 0L)
  expect_identical(members$grouping_status[1], "ineligible_junction;v_call_ambiguous")
})

test_that("missing grouping entities cannot be rescued by a minor sibling", {
  records <- analysis_example()
  records$v_gene_status <- c("missing", "assigned_group")
  records$v_gene <- c("", "VSEQ_TEST")
  for (locus in c("IGH", "IGK", "IGL")) {
    records$locus <- locus
    members <- prepare_grouping_fixture(records, "0460", 0L)
    expect_false(any(members$grouping_eligible))
    expect_identical(members$grouping_status[2], "ineligible_backbone_representative")
  }
})

test_that("backbone ties use content, then sequence ID, after V identity and transfer support once", {
  # assemble1_9 carries the junction TGTGCGAAA, which sorts before TGTGCGAAT.
  records <- analysis_example(sequence_id = c("assemble1_9", "assemble1_1"), support = c(4, 4))
  members <- prepare_grouping_fixture(records, "0460", 0L)
  expect_identical(members$inference_representative_sequence_id, rep("assemble1_9", 2))
  expect_identical(members$inference_consensus_count, c(8, 0))
  records$v_identity <- c(0.9, 1)
  members <- prepare_grouping_fixture(records, "0460", 0L)
  expect_identical(members$inference_representative_sequence_id,
                   rep("assemble1_1", 2))
  expect_identical(members$inference_consensus_count, c(0, 8))
  reordered <- prepare_grouping_fixture(records[2:1, ], "0460", 0L)
  expect_identical(reordered$inference_representative_sequence_id, rep("assemble1_1", 2))
  same_content <- analysis_example(sequence_id = c("assemble1_9", "assemble1_1"),
                                   support = c(4, 4), junction = rep("TGTGCGAAA", 2))
  members <- prepare_grouping_fixture(same_content, "0460", 0L)
  expect_identical(members$inference_representative_sequence_id, rep("assemble1_1", 2))
  expect_identical(members$inference_consensus_count, c(0, 8))
})

# TRUST4 numbers the CDR3 variants of a contig differently between runs. When the
# tied variants differ in eligibility, the choice decides whether the contig
# enters grouping, so it must not follow the trailing index.
test_that("renumbering tied variants leaves the backbone choice unchanged", {
  records <- analysis_example(support = c(4, 4))
  records$v_gene_status[2] <- "missing"
  records$v_gene[2] <- ""
  renumbered <- records
  renumbered$sequence_id <- rev(records$sequence_id)
  for (locus in c("IGH", "IGK", "IGL")) {
    records$locus <- locus
    renumbered$locus <- locus
    first <- prepare_grouping_fixture(records, "0460", 0L)
    second <- prepare_grouping_fixture(renumbered, "0460", 0L)
    expect_identical(first$grouping_eligible, c(TRUE, FALSE))
    expect_identical(second$grouping_eligible, first$grouping_eligible)
    expect_identical(second$grouping_status, first$grouping_status)
    expect_identical(second$inference_consensus_count, first$inference_consensus_count)
  }
})

test_that("all four TCR loci retain each eligible variant as its own observation", {
  for (locus in c("TRA", "TRB", "TRG", "TRD")) {
    records <- analysis_example(locus)
    members <- prepare_grouping_fixture(records, "0460", 0L)
    expect_identical(members$inference_representative_sequence_id, records$sequence_id)
    expect_true(all(members$is_inference_representative))
    expect_identical(members$inference_consensus_count, records$consensus_count)
  }
})

test_that("a mixed library preserves all seven loci and never collapses across loci", {
  loci <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
  records <- do.call(rbind, lapply(seq_along(loci), function(index) {
    analysis_example(loci[index], paste0("assemble1_", index), index, "TGTGCGAAA")
  }))
  members <- prepare_grouping_fixture(records, "0460", 0L)
  expect_identical(members$locus, loci)
  expect_true(all(members$is_inference_representative))
  expect_identical(members$inference_representative_sequence_id, records$sequence_id)
  expect_equal(members$inference_consensus_count, records$consensus_count)
})

test_that("eligible light-chain observations reproduce exact group counts and support", {
  records <- analysis_example("IGK", sequence_id = c("assemble1_0", "assemble1_1", "assemble2_0"),
                               support = c(9, 1, 4), junction = rep("TGTGCGAAA", 3))
  records$sequence[3] <- paste0("T", records$sequence[3])
  records$sequence_alignment[3] <- records$sequence[3]
  members <- prepare_grouping_fixture(records, "0460", 0L)
  groups <- build_rearrangement_group_table(records, "0460", "IGK")$groups
  expect_identical(sum(members$inference_consensus_count), 14)
  expect_identical(sum(members$is_inference_representative), 2L)
  expect_equal(groups$consensus_count, sum(members$inference_consensus_count))
  expect_identical(groups$n_contigs, 2L)
})

test_that("every locus conserves exactly the support eligible for its own inference", {
  records <- rbind(analysis_example("IGH"),
                   analysis_example("IGK", c("assemble2_0", "assemble2_1")),
                   analysis_example("TRB", c("assemble3_0", "assemble3_1")))
  records$v_gene_status[1] <- "ambiguous"
  records$v_gene[1] <- ""
  records$junction[5] <- "TGTNCGAAA"
  members <- prepare_grouping_fixture(records, "0460", 0L)
  for (locus in unique(records$locus)) {
    included <- members$locus == locus & members$grouping_eligible
    expect_equal(sum(members$inference_consensus_count[members$locus == locus]),
                 sum(records$consensus_count[included]))
  }
  expect_true(all(members$inference_consensus_count[!members$is_inference_representative] == 0))
})

test_that("literal identifiers remain text through sidecar serialization", {
  records <- analysis_example(sequence_id = c("0460", "NA"))
  records$library_id <- "NA"
  members <- prepare_grouping_fixture(records, "NA", 0L)
  output <- tempfile(fileext = ".tsv")
  on.exit(unlink(output), add = TRUE)
  write_pipeline_table(members, output)
  roundtrip <- read_pipeline_table(output)
  expect_identical(roundtrip$library_id, c("NA", "NA"))
  expect_identical(roundtrip$sequence_id, c("0460", "NA"))
  expect_identical(roundtrip$inference_representative_sequence_id, c("0460", "NA"))
})

test_that("an empty or entirely ineligible library still has its declared row contract", {
  records <- analysis_example()
  empty <- prepare_grouping_fixture(records[0, ], "0460", 0L)
  expect_identical(empty, make_empty_grouping())
  expect_type(empty$grouping_eligible, "logical")
  expect_type(empty$inference_consensus_count, "double")
  records$junction <- ""
  members <- prepare_grouping_fixture(records, "0460", 0L)
  expect_identical(members$sequence_id, records$sequence_id)
  expect_false(any(members$grouping_eligible))
  expect_true(all(is.na(members$inference_representative_sequence_id)))
  expect_identical(members$inference_consensus_count, c(0, 0))
})

test_that("invalid source identity, locus or common-QC support is rejected", {
  records <- analysis_example()
  expect_error(prepare_grouping_fixture(records, 460, 0L), "string library_id")
  expect_error(prepare_grouping_fixture(records, "other", 0L), "another library")
  duplicated <- records
  duplicated$sequence_id[2] <- duplicated$sequence_id[1]
  expect_error(prepare_grouping_fixture(duplicated, "0460", 0L), "unique, nonempty")
  invalid <- records
  invalid$locus[1] <- "IG"
  expect_error(prepare_grouping_fixture(invalid, "0460", 0L), "seven receptor loci")
  records$consensus_count[1] <- 0
  expect_error(prepare_grouping_fixture(records, "0460", 0L), "non-positive")
  expect_error(prepare_grouping_fixture(analysis_example()[0, ], "0460", -1), "max_n")
})
test_that("prepared IGH rejects support on a non-representative row", {
  records <- rbind(analysis_example("IGH"),
                   analysis_example("IGK", c("assemble2_0", "assemble2_1")),
                   analysis_example("TRB", c("assemble3_0", "assemble3_1")))
  records$junction[1] <- "NNN"
  records$v_gene[3] <- ""
  igh <- prepare_igh_grouping(records[records$locus == "IGH", ], "0460", 0)
  nonigh <- lapply(c("IGK", "TRB"), function(locus) {
    build_rearrangement_group_table(records[records$locus == locus, ], "0460", locus)$mapping
  })
  actual <- do.call(rbind, c(list(igh), nonigh))
  expect_equal(actual$inference_consensus_count, c(0, 0, 0, 0, 9, 1))
  reps <- prepared_igh_representatives(records[records$locus == "IGH", ], igh, TRUE)
  expect_identical(as.numeric(reps$consensus_count), numeric())
  bad <- igh
  bad$inference_consensus_count[2] <- 2
  expect_error(validate_grouping_mapping(records[records$locus == "IGH", ], bad),
               "Invalid prepared representative links or support")
})

test_that("audit contracts reject lost support and changed checkpoint counts", {
  expect_error(qc_audit_rows("NA", "collapse", 2, 1, 10, 9, analysis = "grouping",
                            stage_index = 1, kind = "grouping", support_conserved = TRUE),
               "preserve defined support")
  expect_error(qc_audit_rows("NA", "retained", 2, 1, analysis = "common", stage_index = 1),
               "checkpoint changed")
  rows <- qc_audit_rows("NA", "retained", 2, support_in = 10, analysis = "common", stage_index = 1)
  expect_identical(names(rows), QC_AUDIT_COLUMNS)
  expect_true(is.na(rows$scope_id))
  expect_identical(rows$library_id, "NA")
  expect_error(validate_qc_audit(rbind(rows, rows)), "Duplicate")
})

test_that("branch grouping counts omit inapplicable TCR collapses", {
  records <- analysis_example("TRB")
  mapping <- build_rearrangement_group_table(records, "0460", "TRB")$mapping
  rows <- qc_audit_grouping_rows("0460", records, mapping, "TRB")
  expect_setequal(rows$stage, c("grouping_eligible", "exact_sequence_collapse", "inference_representatives"))
  expect_false("all" %in% rows$locus)
  expect_true(all(rows$n_units_out == 2))
  expect_true(all(rows$consensus_count_out == 10))
})

test_that("common criteria overlap and never add counterfactual rows", {
  flags <- data.frame(productive = c(FALSE, TRUE), low_n = c(FALSE, TRUE))
  rows <- qc_audit_qc_rows("0460", flags, c(FALSE, TRUE), c(4, 2), c("IGH", "IGH"), 1)
  expect_setequal(rows$stage, c("qc_productive", "qc_low_n", "qc_pass_total"))
  expect_true(all(rows$decision_status[rows$stage == "qc_low_n"] == "disabled"))
  at <- rows[rows$locus == "IGH", ]
  expect_equal(at$n_units_in, rep(2, 3))
  expect_equal(at$n_units_out, rep(1, 3))
  expect_equal(at$consensus_count_out, rep(2, 3))
})

# An order-dependence counterexample: two backbones retain X, and the second
# backbone also has a lower-support Y. Resolve the alternatives before merging X.
test_that("backbone-first exact merging is invariant to every source row order", {
  permutations <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
                       c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  for (locus in c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")) {
    records <- analysis_example(locus, c("assemble1_0", "assemble2_0", "assemble2_1"),
      c(100, 100, 10), c("TGTGCGAAA", "TGTGCGAAA", "TGTGCGAAT"))
    bcr <- locus %in% c("IGH", "IGK", "IGL")
    for (indices in permutations) {
      shuffled <- records[indices, ]
      mapping <- prepare_grouping_fixture(shuffled, "0460", 0)
      mapping <- mapping[match(records$sequence_id, mapping$sequence_id), ]
      expect_identical(mapping$backbone_representative_sequence_id,
        c("assemble1_0", "assemble2_0", if (bcr) "assemble2_0" else "assemble2_1"))
      expect_identical(mapping$inference_representative_sequence_id,
        c("assemble1_0", "assemble1_0", if (bcr) "assemble1_0" else "assemble2_1"))
      expect_identical(mapping$inference_consensus_count,
        if (bcr) c(210, 0, 0) else c(200, 0, 10))
      expect_identical(sum(mapping$is_inference_representative), if (bcr) 1L else 2L)
    }
  }
})

test_that("identical DNA cannot inflate a variant's support before backbone ranking", {
  records <- analysis_example("IGH", c("assemble1_0", "assemble2_0", "assemble2_1"),
    c(100, 40, 50), c("TGTGCGAAA", "TGTGCGAAA", "TGTGCGAAT"))
  mapping <- prepare_igh_grouping(records, "0460", 0)
  expect_identical(mapping$backbone_representative_sequence_id,
    c("assemble1_0", "assemble2_1", "assemble2_1"))
  expect_identical(mapping$inference_representative_sequence_id,
    mapping$backbone_representative_sequence_id)
  expect_identical(mapping$inference_consensus_count, c(100, 0, 90))
})

test_that("exact merging conserves eligible support and original backbone counts", {
  records <- analysis_example("IGK", c("assemble1_0", "assemble2_0", "assemble3_0"),
    c(6, 5, 100), rep("TGTGCGAAA", 3))
  records$v_gene_status[3] <- "ambiguous"
  records$v_gene[3] <- ""
  result <- build_rearrangement_group_table(records, "0460", "IGK")
  expect_equal(result$groups$n_contigs, 2L)
  expect_equal(result$groups$consensus_count, 11)
  expect_identical(result$mapping$inference_consensus_count, c(11, 0, 0))
  expect_true(is.na(result$mapping$backbone_representative_sequence_id[3]))
  audit <- qc_audit_grouping_rows("0460", records, result$mapping, "IGK")
  exact <- audit[audit$stage == "exact_sequence_collapse", ]
  expect_equal(as.numeric(exact[c("n_units_in", "n_units_out",
    "consensus_count_in", "consensus_count_out")]), c(2, 1, 11, 11))
  expect_identical(exact$input_unit, "backbone_representative")
  expect_identical(exact$locus_basis, "annotated")
})

test_that("identical sequences with conflicting analytical annotations are refused", {
  records <- analysis_example("TRB", c("assemble1_0", "assemble2_0"), c(6, 5),
    rep("TGTGCGAAA", 2))
  records$j_call[2] <- "TRBJ2*01"
  expect_error(prepare_grouping_fixture(records, "0460", 0), "Identical sequences disagree in j_call")
})

test_that("C-locus conflicts exclude grouping at every locus without changing raw calls", {
  for (locus in c("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")) {
    records <- analysis_example(locus, sequence_id = c("assemble1_0", "assemble2_0"))
    records$c_call <- c(if (locus == "IGH") "TRBC1*01" else "IGHM*01", "")
    original <- records
    mapping <- prepare_grouping_fixture(records, "0460", 0)
    expect_identical(mapping$grouping_eligible, c(FALSE, TRUE))
    expect_identical(mapping$grouping_status[1], "c_locus_conflict")
    expect_identical(records, original)
  }
})

test_that("eligible dominant BCR support excludes ineligible sibling support", {
  for (locus in c("IGH", "IGK", "IGL")) {
    records <- analysis_example(locus, support = c(9, 1))
    records$v_gene_status[2] <- "ambiguous"
    records$v_gene[2] <- ""
    mapping <- prepare_grouping_fixture(records, "0460", 0)
    expect_identical(mapping$grouping_eligible, c(TRUE, FALSE))
    expect_identical(mapping$inference_consensus_count, c(9, 0))
    expect_identical(mapping$inference_representative_sequence_id, c("assemble1_0", NA_character_))
  }
})

test_that("non-IGH membership gives every retained representative exactly one clonotype", {
  for (locus in c("IGK", "IGL", "TRA", "TRB", "TRD", "TRG")) {
    records <- analysis_example(locus)
    built <- build_rearrangement_group_table(records, "0460", locus)
    members <- built$membership
    expect_identical(names(members), CLONOTYPE_MEMBER_COLUMNS)
    expect_setequal(members$sequence_id, built$mapping$sequence_id[built$mapping$is_inference_representative])
    expect_true(all(members$clonotype_id %in% built$groups$clonotype_id))
    expect_equal(sum(members$member_consensus_count), sum(built$groups$consensus_count))
    expect_equal(sum(members$member_consensus_count), 10)
    expect_false(anyDuplicated(members$sequence_id) > 0)
  }
})

test_that("transferred support parsed from text stays numeric at large counts", {
  records <- analysis_example(support = c("60000", "40000"))
  mapping <- prepare_igh_grouping(records, "0460", 0)
  expect_identical(mapping$inference_consensus_count, c(100000, 0))
})

variant_backbones <- function() {
  data.frame(
    sequence_id = c("assemble1_0", "assemble1_1", "assemble2_0", "assemble2_1"),
    productive = "T",
    junction = c("TGTGCGAAA", "TGTGCGAAT", "TGTGCGCCC", "TGTGCGCCA"),
    v_call = c("IGKV1-5*01", "IGKV1-5*03", "IGKV1-5*01", "IGKV1-6*01"),
    j_call = "IGKJ1*01",
    consensus_count = c(100, 2, 20, 1),
    stringsAsFactors = FALSE)
}

describe("two backbones, one with an allele difference and one with a V-gene difference", {
  variants <- variant_backbones()

  it("maps every row of a backbone onto its best-supported row, however the calls differ",
     expect_identical(trust4_assembly_representative_index(variants),
                      c(1L, 1L, 3L, 3L)))
})

# Pooled identifiers carry a library prefix; only the trailing CDR3 index is
# removed, so one contig number cannot merge two libraries.
test_that("library-namespaced backbones stay distinct", {
  pooled <- data.frame(
    sequence_id = c("WTS_A~assemble1_0", "WTS_A~assemble1_1", "WTS_B~assemble1_0"),
    v_call = "IGHV3-23*01", j_call = "IGHJ4*02",
    junction = c("TGTGCGAAA", "TGTGCGAAT", "TGTGCGAAA"),
    consensus_count = c(9, 1, 5),
    stringsAsFactors = FALSE)
  expect_identical(trust4_assembly_representative_index(pooled), c(1L, 1L, 3L))
})

# Alternative CDR3 rows can receive different or missing V/J annotations.
# Their assembly ID still defines the backbone, and original read support
# chooses its representative. Row 1 has support 100 and row 2 has support 2.
describe("annotation disagreement between siblings", {
  disagreement <- function(v_call, j_call) {
    rows <- variant_backbones()[1:2, ]
    rows$v_call <- v_call
    rows$j_call <- j_call
    trust4_assembly_representative_index(rows)
  }

  it("collapses when one sibling has no V call at all",
     expect_identical(disagreement(c("", "IGKV1-5*01"), "IGKJ1*01"), c(1L, 1L)))
  it("collapses when both siblings report the same paralog ambiguity",
     expect_identical(disagreement("IGHV3-30*01,IGHV3-30-3*01",
                                   "IGHJ4*02,IGHJ5*02"), c(1L, 1L)))
  it("collapses when one sibling's call is narrower than the other's",
     expect_identical(disagreement("IGHV1-69*01,IGHV1-69D*01",
                                   c("IGHJ4*02,IGHJ5*02", "IGHJ4*02")),
                      c(1L, 1L)))
  it("collapses when the two J calls name entirely different genes",
     expect_identical(disagreement("IGHV1-69*01,IGHV1-69D*01",
                                   c("IGHJ4*02", "IGHJ5*02")), c(1L, 1L)))
})

# The ranking is a total order, so input order cannot change the winner.
test_that("the representative choice is independent of row order", {
  variants <- variant_backbones()
  reversed <- variants[rev(seq_len(nrow(variants))), ]
  selected <- unique(trust4_assembly_representative_index(reversed))
  expect_identical(sort(reversed$sequence_id[selected]),
                   c("assemble1_0", "assemble2_0"))
})

describe("a frame with no rows", {
  empty <- variant_backbones()[0, ]

  it("has no representatives to map",
     expect_identical(trust4_assembly_representative_index(empty), integer(0)))
})

test_that("a missing required column is an error", {
  variants <- variant_backbones()
  expect_error(trust4_assembly_representative_index(
    variants[, setdiff(names(variants), "v_call")]))
})

# The vectorised ranking must pick the same representative as the shared
# per-group ladder in rearrangement_groups.R. Exercise every criterion in the ladder,
# including a pure tie that falls through to the lexicographic sequence_id.
describe("the vectorised ranking against the shared per-group ladder", {
  ladder <- data.frame(
    sequence_id = c(
      "bbA_0", "bbA_1", "bbA_2",   # most read support wins despite the worst identity
      "bbB_0", "bbB_1",            # equal support, so an ambiguous V call loses
      "bbC_0", "bbC_1",            # equal support and both clean: v_identity decides
      "bbD_0", "bbD_1",            # equal on everything ranked: seq_id decides
      "bbE_9", "bbE_1",            # same, with the lower id written second
      "bbF_0", "bbF_1",            # equal support, so a missing junction base loses
      "bbG_0", "bbG_1",            # support outranks cleanliness when they conflict
      "bbH_0", "bbH_1",            # equal on quality: the lower junction wins
      "bbI_0", "bbI_1"             # same junction: the lower full sequence wins
    ),
    junction = c(rep("TGTGCGAAA", 11L), "TGTGCGAAN", rep("TGTGCGAAA", 3L),
                 "TGTGCGAAT", "TGTGCGAAA", rep("TGTGCGAAA", 2L)),
    v_call = c("IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01",
               "IGHV1-1*01,IGHV2-1*01", "IGHV1-1*01",
               "IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01",
               "IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*01",
               "IGHV1-1*01,IGHV2-1*01", "IGHV1-1*01", rep("IGHV1-1*01", 4L)),
    j_call = "IGHJ4*01",
    v_identity = c("0.90", "0.99", "0.95", "0.99", "0.99",
                   "0.90", "0.99", "0.95", "0.95", "0.95", "0.95",
                   "0.95", "0.95", "0.95", "0.95", rep("0.95", 4L)),
    sequence_alignment = "ACGTACGTA",
    sequence = c(rep("ACGTACGTA", 17L), "TTTTTTTTT", "AAAAAAAAA"),
    consensus_count = c(500, 10, 5, 7, 7, 4, 4, 3, 3, 3, 3, 6, 6, 400, 1, 2, 2, 2, 2),
    stringsAsFactors = FALSE)
  vectorised <- trust4_assembly_representative_index(ladder)

  it("agrees with select_rearrangement_group_representative on every backbone", {
    ladder_backbone <- sub("_[0-9]+$", "", ladder$sequence_id)
    per_group_winner <- vapply(split(seq_len(nrow(ladder)), ladder_backbone),
                               function(i) {
                                 i[select_rearrangement_group_representative(
                                   ladder[i, , drop = FALSE])$index]
                               }, integer(1))
    expect_identical(as.integer(sort(unique(vectorised))),
                     as.integer(sort(unname(per_group_winner))))
  })
  it("prefers read support, then unambiguous, then identity, then content, then sequence_id",
     expect_identical(ladder$sequence_id[unique(vectorised)],
                      c("bbA_0", "bbB_1", "bbC_1", "bbD_0", "bbE_1", "bbF_1",
                        "bbG_0", "bbH_1", "bbI_1")))
})

# ----------------------------------- deterministic representative preference ---

test_that("the best supported member represents the group even when it is the least clean", {
  # Row 2 carries an ambiguous junction base and row 3 names two V genes, but row
  # 2 has 200x the reads of the tidiest member. An ambiguous call reports what
  # IgBLAST could separate, not how well the member was observed, so it does not
  # unseat the evidence.
  group <- data.frame(
    sequence_id = c("assemble2_0", "assemble1_0", "assemble3_0"),
    junction = c("TGTGCGAAA", "TGTNCGAAA", "TGTGCGAAA"),
    v_call = c("IGKV1-5*01", "IGKV1-5*01", "IGKV1-5*01,IGKV1-9*01"),
    j_call = c("IGKJ1*01", "IGKJ1*01", "IGKJ1*01"),
    consensus_count = c(5, 999, 50),
    stringsAsFactors = FALSE)
  chosen <- select_rearrangement_group_representative(group)
  expect_identical(group$sequence_id[chosen$index], "assemble1_0")
})

test_that("cleanliness decides only once read support is tied", {
  # Same three members, now equally supported: row 1 is the only one with both an
  # unambiguous junction and a single-gene V call, so it wins the tie-break.
  group <- data.frame(
    sequence_id = c("assemble2_0", "assemble1_0", "assemble3_0"),
    junction = c("TGTGCGAAA", "TGTNCGAAA", "TGTGCGAAA"),
    v_call = c("IGKV1-5*01", "IGKV1-5*01", "IGKV1-5*01,IGKV1-9*01"),
    j_call = c("IGKJ1*01", "IGKJ1*01", "IGKJ1*01"),
    consensus_count = c(50, 50, 50),
    stringsAsFactors = FALSE)
  chosen <- select_rearrangement_group_representative(group)
  expect_identical(group$sequence_id[chosen$index], "assemble2_0")
})

test_that("members equal on every quality gate are separated by sequence_id", {
  group <- data.frame(
    sequence_id = c("assemble9_0", "assemble1_0"),
    productive = c("T", "T"), junction = c("TGTGCG", "TGTGCG"),
    v_call = c("IGKV1-5*01", "IGKV1-5*01"),
    j_call = c("IGKJ1*01", "IGKJ1*01"),
    consensus_count = c(10, 10), stringsAsFactors = FALSE)
  chosen <- select_rearrangement_group_representative(group)
  expect_identical(group$sequence_id[chosen$index], "assemble1_0")
})
