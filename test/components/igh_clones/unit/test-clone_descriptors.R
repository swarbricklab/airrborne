# Check donor clone descriptors, observed isotypes and gene agreement.
# Run: python test/run_tests.py --component igh_clones --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys", "igh_clones",
  "isotype_resolution", "pooled_records", "qc_audit", "shm_summaries", "gene_usage",
  "observed_representatives", "igh_summaries"))

MEMBERSHIP <- data.frame(
  donor_id = c("D1", "D1", "D1", "D1", "D1", "D2", "D1", "D1", "D1", "D1"),
  library_id = c("S1", "S1", "S1", "S1", "S2", "S3", "S1", "S1", "S1", "S1"),
  sequence_id = paste0("S", c(1, 1, 1, 1, 2, 3, 1, 1, 1, 1), "~seq", 1:10),
  clone_id = c("L1", "L1", "L2", "L3", "L1", "L4", "L5", "L5", "L6", "L6"),
  consensus_count = c("3", "4", "5", "6", "7", "8", "2", "2", "2", "2"),
  qualified_c_call = c("IGHG1*01", "IGHG3*01", "IGHM*01", "IGHG1,IGHA1", "IGHM*01", "",
             "IGHM*01", "IGHM*01", "IGHM*01", "IGHM*01"),
  v_call = c("IGHV1-2*01", "IGHV1-2*01,IGHV4-59*01", "IGHV3-23*01",
             "IGHV3-23*01,IGHV3-30*01", "IGHV4-59*01", "",
             "IGHV3-23*01", "IGHV3-23*01",
             "IGHV3-23*01", "IGHV3-30*01"),
  j_call = c("IGHJ4*02", "IGHJ4*02", "IGHJ4*02", "IGHJ5*02", "IGHJ4*02",
             "IGHJ6*02", "IGHJ4*02", "IGHJ6*02", "IGHJ4*02", "IGHJ4*02"),
  v_gene = c("IGHV1-2", "IGHV1-2", "IGHV3-23", "VSEQ_x", "IGHV1-2", "IGHV1-69",
             "IGHV3-23", "IGHV3-23", "VSEQ_x", "VSEQ_x"),
  j_gene = c("IGHJ4", "IGHJ4", "IGHJ4", "IGHJ5", "IGHJ4", "IGHJ6", "IGHJ4",
             "IGHJ4", "IGHJ4", "IGHJ4"),
  junction_length = c(42L, 42L, 45L, 48L, 42L, 39L, 45L, 45L, 51L, 51L),
  stringsAsFactors = FALSE
)

MEMBERSHIP$junction <- rep("TGTGCC", nrow(MEMBERSHIP))
MEMBERSHIP$sequence <- rep("ACGTTGTGCC", nrow(MEMBERSHIP))

# A composite entity reads as its member genes; a single gene is its own label.
LABELS <- c(VSEQ_x = "IGHV3-23,IGHV3-30")
FAMILIES <- build_clone_family_table(MEMBERSHIP, "scope_x", LABELS)

test_that("the released columns are declared in order", {
  expect_identical(names(FAMILIES), CLONE_FAMILY_COLUMNS)
})

describe("a clone spanning two of a donor's libraries", {
  l1 <- FAMILIES[FAMILIES$clone_id == "L1", , drop = FALSE]

  it("counts every member, whichever library it came from",
     expect_identical(l1$n_members, 3L))
  it("counts the libraries its members came from",
     expect_identical(l1$n_libraries, 2L))
  it("reports the junction length its members share",
     expect_identical(l1$junction_length, 42L))
  # Read support is a transcript-consensus_count weight summed over the members: 3+4+7.
  it("sums read support across the members",
     expect_equal(l1$consensus_count, 14))
  it("lists the classes it was seen expressing, sorted",
     expect_identical(l1$isotype_classes, "IgG,IgM"))
  it("flags a clone holding both switch states",
     expect_true(l1$has_unswitched_and_switched_members))
  it("stamps the scope and the donor", {
    expect_identical(l1$scope_id, "scope_x")
    expect_identical(l1$donor_id, "D1")
  })
})

describe("a clone with one unswitched member", {
  l2 <- FAMILIES[FAMILIES$clone_id == "L2", , drop = FALSE]

  it("reports the one class", expect_identical(l2$isotype_classes, "IgM"))
  it("is not flagged as switched",
     expect_false(l2$has_unswitched_and_switched_members))
})

# IGHG1 and IGHA1 are two classes but one switch state, so the class is
# unresolvable while the switch question still has an answer.
describe("a clone whose member names two classes at once", {
  l3 <- FAMILIES[FAMILIES$clone_id == "L3", , drop = FALSE]

  it("reports no class", expect_identical(l3$isotype_classes, ""))
  it("still answers the switching question",
     expect_false(l3$has_unswitched_and_switched_members))
})

# Absence of evidence is not evidence the clone did not switch, which is why this
# is NA rather than FALSE.
describe("a clone with no usable constant-region call", {
  l4 <- FAMILIES[FAMILIES$clone_id == "L4", , drop = FALSE]

  it("reports no class", expect_identical(l4$isotype_classes, ""))
  it("leaves the switching question undefined",
     expect_identical(l4$has_unswitched_and_switched_members, NA))
})

# Clone inference groups within exact V and J entities, so a family carries the
# one entity its members share. The fixture supplies those entities directly.
describe("the V and J entity of a clone", {
  family <- function(lineage) FAMILIES[FAMILIES$clone_id == lineage, , drop = FALSE]

  it("carries the entity every member shares",
     expect_identical(family("L1")$v_gene_key, "IGHV1-2"))
  it("names a composite by its member genes",
     expect_identical(family("L6")$v_gene_names, "IGHV3-23,IGHV3-30"))
  it("keeps the composite identifier as the key",
     expect_identical(family("L6")$v_gene_key, "VSEQ_x"))
  it("labels a single gene as itself",
     expect_identical(family("L5")$v_gene_names, "IGHV3-23"))
  it("names the J entity the same way",
     expect_identical(family("L3")$j_gene_names, "IGHJ5"))
})

test_that("a clone whose members name two V entities is rejected", {
  broken <- MEMBERSHIP
  broken$v_gene[2] <- "IGHV4-59"
  expect_error(build_clone_family_table(broken, "scope_x", LABELS),
               "Clone L1 has members with more than one V entity")
})

test_that("a clone whose members have two junction lengths is rejected", {
  broken <- MEMBERSHIP
  broken$junction_length[5] <- 45L
  expect_error(build_clone_family_table(broken, "scope_x", LABELS),
               "Clone L1 has members with more than one junction length: 42, 45")
})

test_that("families are ordered by donor and then lineage", {
  expect_identical(FAMILIES$donor_id, c("D1", "D1", "D1", "D1", "D1", "D2"))
  expect_identical(FAMILIES$clone_id, c("L1", "L2", "L3", "L5", "L6", "L4"))
})

test_that("a scope with no members still declares the released header", {
  empty <- build_clone_family_table(MEMBERSHIP[0, , drop = FALSE], "scope_x")
  expect_identical(names(empty), CLONE_FAMILY_COLUMNS)
  expect_identical(nrow(empty), 0L)
})

test_that("a scope holding one member reports one family", {
  one <- build_clone_family_table(MEMBERSHIP[3, , drop = FALSE], "scope_x")
  expect_identical(nrow(one), 1L)
  expect_identical(one$n_members, 1L)
  expect_identical(one$n_libraries, 1L)
  expect_equal(one$consensus_count, 5)
})

# Read support is the consensus_count every IGH metric shares, so an unreadable value is
# a contract violation rather than something to treat as zero.
test_that("unreadable read support is rejected", {
  broken <- MEMBERSHIP
  broken$consensus_count[2] <- "lots"
  expect_error(build_clone_family_table(broken, "scope_x"),
               "consensus_count contains 1 invalid value")
})

test_that("the donor representative uses one observed member without pooling support", {
  records <- MEMBERSHIP[MEMBERSHIP$clone_id == "L1", , drop = FALSE]
  result <- build_clone_family_table(records, "scope_x", LABELS)
  expect_identical(result$representative_sequence_id, "S2~seq5")
  expect_identical(result$representative_library_id, "S2")
  # Two identical lower-support observations do not replace the best member.
  records$consensus_count <- c(4, 4, 7)
  result <- build_clone_family_table(records, "scope_x", LABELS)
  expect_identical(result$representative_sequence_id, "S2~seq5")
})
