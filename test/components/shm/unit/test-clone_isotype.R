# Check per-library clone isotype rows and their callable denominators.
# Run: python test/run_tests.py --component shm --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys", "igh_clones",
  "isotype_resolution", "isotype_metrics", "pooled_records", "qc_audit", "shm_summaries",
  "gene_usage", "observed_representatives", "igh_summaries"))

# Read down the clones:
#
#   L1  S1 saw IGHG1 and IGHG3, and S2 saw IGHM, so each library's rows describe
#       its own members and the pooled clone table holds both states
#   L2  one member, IGHM
#   L3  one member naming IgG and IgA at once -> no class and no subclass
#   L4  one member with no call -> no rows
#   L5  two IGHM members, which share one row at each level
#   L7  an IgM beside an IgA in one library
#   L8  an IGHG1 member beside one tying IGHG1 and IGHG3 -> both are IgG, but only
#       one resolves to a subclass, so the two levels have different denominators
MEMBERSHIP <- data.frame(
  donor_id = c("D1", "D1", "D1", "D1", "D1", "D2", "D1", "D1", "D1", "D1", "D1", "D1"),
  library_id = c("S1", "S1", "S1", "S1", "S2", "S3", "S1", "S1", "S1", "S1", "S1", "S1"),
  sequence_id = paste0("S", c(1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1), "~seq", 1:12),
  clone_id = c("L1", "L1", "L2", "L3", "L1", "L4", "L5", "L5", "L7", "L7", "L8", "L8"),
  consensus_count = c("3", "4", "5", "6", "7", "8", "2", "2", "1", "3", "3", "5"),
  qualified_c_call = c("IGHG1*01", "IGHG3*01", "IGHM*01", "IGHG1,IGHA1", "IGHM*01", "",
                       "IGHM*01", "IGHM*01", "IGHM*01", "IGHA1*01", "IGHG1*01", "IGHG1,IGHG3"),
  v_gene = "IGHV3-23", j_gene = "IGHJ4", junction_length = 45L,
  stringsAsFactors = FALSE
)

MANIFEST <- data.frame(
  library_id = c("S1", "S2", "S3"), donor_id = c("D1", "D1", "D2"),
  tissue_id = c("T1", "T2", "T3"), stringsAsFactors = FALSE
)

BY_LIBRARY <- build_lineage_isotype_table(MEMBERSHIP, MANIFEST, "scope_x")

cell <- function(library_id, clone_id, level, isotype = NULL) {
  selected <- BY_LIBRARY$library_id == library_id & BY_LIBRARY$clone_id == clone_id &
    BY_LIBRARY$isotype_level == level
  if (!is.null(isotype)) selected <- selected & BY_LIBRARY$isotype == isotype
  BY_LIBRARY[selected, , drop = FALSE]
}

test_that("the released columns are declared in order", {
  expect_identical(names(BY_LIBRARY), LINEAGE_ISOTYPE_COLUMNS)
  expect_identical(LINEAGE_ISOTYPE_COLUMNS,
                   c("scope_id", "library_id", "tissue_id", "donor_id", "clone_id", "isotype_level",
                     "isotype", "n_members", "consensus_count", "callable_consensus_count",
                     "support_fraction"))
})

test_that("each library, clone, level and isotype gets one row", {
  expect_identical(nrow(BY_LIBRARY), 15L)
  expect_false(anyDuplicated(BY_LIBRARY[c("library_id", "clone_id", "isotype_level", "isotype")]) > 0)
  expect_identical(BY_LIBRARY$scope_id, rep("scope_x", nrow(BY_LIBRARY)))
})

# Scope collection merges donor shards in this order. IgM precedes IgA because
# ISOTYPE_COMPOSITION_LEVELS lists it first, although text order would reverse them.
test_that("rows follow library and clone, then class before subclass in level order", {
  expect_identical(unique(BY_LIBRARY$library_id), c("S1", "S2"))
  expect_identical(unique(BY_LIBRARY$clone_id[BY_LIBRARY$library_id == "S1"]),
                   c("L1", "L2", "L5", "L7", "L8"))
  l7 <- BY_LIBRARY[BY_LIBRARY$clone_id == "L7", ]
  expect_identical(l7$isotype_level, c("class", "class", "subclass", "subclass"))
  expect_identical(l7$isotype, c("IgM", "IgA", "IgM", "IgA1"))
  l1 <- BY_LIBRARY[BY_LIBRARY$clone_id == "L1" & BY_LIBRARY$library_id == "S1", ]
  expect_identical(l1$isotype, c("IgG", "IgG1", "IgG3"))
})

test_that("the vocabulary is the one ISOTYPE_COMPOSITION_LEVELS declares", {
  for (level in names(ISOTYPE_COMPOSITION_LEVELS)) {
    expect_true(all(BY_LIBRARY$isotype[BY_LIBRARY$isotype_level == level] %in%
                      ISOTYPE_COMPOSITION_LEVELS[[level]]))
  }
  expect_setequal(BY_LIBRARY$isotype_level, names(ISOTYPE_COMPOSITION_LEVELS))
})

test_that("the tissue is joined from the scope manifest", {
  expect_identical(unique(BY_LIBRARY$tissue_id[BY_LIBRARY$library_id == "S2"]), "T2")
})

describe("the support of one isotype within a clone", {
  it("sums the transferred support of the members resolving to it", {
    igg <- cell("S1", "L1", "class", "IgG")
    expect_identical(igg$n_members, 2L)
    expect_equal(igg$consensus_count, 7)
  })
  it("divides by the support resolving at that level", {
    igg1 <- cell("S1", "L1", "subclass", "IgG1")
    igg3 <- cell("S1", "L1", "subclass", "IgG3")
    expect_equal(igg1$callable_consensus_count, 7)
    expect_equal(igg1$support_fraction, 3 / 7)
    expect_equal(igg3$support_fraction, 4 / 7)
  })
  it("gives two members one row when they resolve alike", {
    igm <- cell("S1", "L5", "class", "IgM")
    expect_identical(igm$n_members, 2L)
    expect_equal(igm$support_fraction, 1)
  })
  it("separates minor from dominant support", {
    expect_equal(cell("S1", "L7", "class", "IgM")$support_fraction, 1 / 4)
    expect_equal(cell("S1", "L7", "class", "IgA")$support_fraction, 3 / 4)
  })
})

# IGHG1,IGHG3 is unambiguously IgG and deliberately nothing below it.
describe("a member resolving at one level and not the other", {
  it("enters the class row", {
    igg <- cell("S1", "L8", "class", "IgG")
    expect_identical(igg$n_members, 2L)
    expect_equal(igg$consensus_count, 8)
    expect_equal(igg$callable_consensus_count, 8)
  })
  it("stays out of the subclass denominator", {
    igg1 <- cell("S1", "L8", "subclass", "IgG1")
    expect_identical(nrow(cell("S1", "L8", "subclass")), 1L)
    expect_identical(igg1$n_members, 1L)
    expect_equal(igg1$consensus_count, 3)
    expect_equal(igg1$callable_consensus_count, 3)
    expect_equal(igg1$support_fraction, 1)
  })
})

# A lineage is pooled across a donor's libraries; each library's rows describe only
# its own members.
describe("a clone split across two of a donor's libraries", {
  it("reports what each library saw, not the donor's isotypes", {
    expect_identical(cell("S1", "L1", "class")$isotype, "IgG")
    expect_identical(cell("S2", "L1", "class")$isotype, "IgM")
  })
  it("counts and weights only that library's members", {
    expect_equal(cell("S2", "L1", "class")$callable_consensus_count, 7)
    expect_identical(cell("S2", "L1", "class")$n_members, 1L)
  })
  it("keeps the donor-pooled switch flag on the clone table", {
    families <- build_clone_family_table(MEMBERSHIP, "scope_x")
    expect_true(families$has_unswitched_and_switched_members[families$clone_id == "L1"])
    expect_identical(families$isotype_classes[families$clone_id == "L1"], "IgG,IgM")
  })
})

test_that("a clone with no callable member has no rows", {
  expect_false(any(BY_LIBRARY$clone_id %in% c("L3", "L4")))
})

test_that("the fractions sum to one within each library, clone and level", {
  totals <- tapply(BY_LIBRARY$support_fraction,
                   paste(BY_LIBRARY$library_id, BY_LIBRARY$clone_id, BY_LIBRARY$isotype_level), sum)
  expect_equal(as.numeric(totals), rep(1, length(totals)))
})

test_that("a scope with no members or no callable member still declares the released header", {
  empty <- build_lineage_isotype_table(MEMBERSHIP[0, , drop = FALSE], MANIFEST, "scope_x")
  expect_identical(names(empty), LINEAGE_ISOTYPE_COLUMNS)
  expect_identical(nrow(empty), 0L)
  uncalled <- build_lineage_isotype_table(MEMBERSHIP[MEMBERSHIP$clone_id %in% c("L3", "L4"), ],
                                          MANIFEST, "scope_x")
  expect_identical(names(uncalled), LINEAGE_ISOTYPE_COLUMNS)
  expect_identical(nrow(uncalled), 0L)
})
