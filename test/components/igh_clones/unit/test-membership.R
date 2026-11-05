# Check donor/scope clone labels, eligibility and membership support.
# Run: python test/run_tests.py --component igh_clones --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys", "igh_clones"))

JUNCTIONS <- data.frame(
  sequence_id = paste0("seq", 1:6),
  junction = c("ACGTACGT", "ACNTACGT", "ACNNACGT", "", "acgtacgt", NA_character_),
  # Every row here answers the junction question alone, so its calls each name
  # one germline entity.
  v_gene_status = "resolved",
  j_gene_status = "resolved",
  stringsAsFactors = FALSE)

describe("clone_record_eligible selects the junctions that can be placed", {
  it("keeps an unambiguous junction and nothing else at max_n = 0", {
    eligible <- clone_record_eligible(JUNCTIONS$junction, JUNCTIONS$v_gene_status,
                                      JUNCTIONS$j_gene_status, 0L)
    kept <- JUNCTIONS[eligible, , drop = FALSE]
    expect_identical(kept$sequence_id, "seq1")
  })
  # Clustering compares junctions of equal length, so a missing base is tolerable
  # only up to the allowance the operator set.
  it("admits one missing base at max_n = 1", {
    eligible <- clone_record_eligible(JUNCTIONS$junction, JUNCTIONS$v_gene_status,
                                      JUNCTIONS$j_gene_status, 1L)
    kept <- JUNCTIONS[eligible, , drop = FALSE]
    expect_identical(kept$sequence_id, c("seq1", "seq2"))
  })
  it("admits two at max_n = 2", {
    eligible <- clone_record_eligible(JUNCTIONS$junction, JUNCTIONS$v_gene_status,
                                      JUNCTIONS$j_gene_status, 2L)
    kept <- JUNCTIONS[eligible, , drop = FALSE]
    expect_identical(kept$sequence_id, c("seq1", "seq2", "seq3"))
  })
  it("never admits an empty, lower-case or absent junction", {
    eligible <- clone_record_eligible(JUNCTIONS$junction, JUNCTIONS$v_gene_status,
                                      JUNCTIONS$j_gene_status, 8L)
    kept <- JUNCTIONS[eligible, , drop = FALSE]
    expect_false(any(c("seq4", "seq5", "seq6") %in% kept$sequence_id))
  })
  it("returns an empty eligibility vector for an empty pool", {
    empty <- JUNCTIONS[0, ]
    expect_identical(clone_record_eligible(empty$junction, empty$v_gene_status,
                                          empty$j_gene_status, 0L), logical(0))
  })
})

# Which germline a receptor was rearranged from is the other half of the gate.
# The tie kinds arrive as columns, written once by the reannotation stage.
CALLS <- data.frame(
  sequence_id = paste0("seq", 1:5),
  junction = "ACGTACGT",
  v_gene_status = c("resolved", "assigned_group", "ambiguous",
                      "resolved", "missing"),
  j_gene_status = c("resolved", "resolved", "resolved",
                      "ambiguous", "resolved"),
  stringsAsFactors = FALSE)

describe("clone_record_eligible selects the records naming one germline", {
  eligible <- clone_record_eligible(CALLS$junction, CALLS$v_gene_status,
                                    CALLS$j_gene_status, 0L)
  kept <- CALLS[eligible, , drop = FALSE]

  # A call the reference reduces to one entity has been answered, so the record
  # takes part; a call standing for several different germline sequences has not.
  it("keeps a resolved call and a tie assigned to one entity", {
    expect_true(all(c("seq1", "seq2") %in% kept$sequence_id))
  })
  it("drops an unresolved tie in either call", {
    expect_false(any(c("seq3", "seq4") %in% kept$sequence_id))
  })
  # Absence is not ambiguity: nothing was reported, so there is no tie to reduce.
  it("excludes a record whose call names no gene", {
    expect_false("seq5" %in% kept$sequence_id)
  })
  it("still drops an ineligible junction", {
    ineligible <- CALLS
    ineligible$junction[1] <- "ACNTACGT"
    expect_false(clone_record_eligible(ineligible$junction, ineligible$v_gene_status,
                                       ineligible$j_gene_status, 0L)[1])
  })
})

# --------------------------------------------------------------- the name ---

test_that("a lineage is named by the pool its label was assigned within", {
  expect_identical(make_clone_id("cohort", "D1", 7L), "cohort~D1~7")
  expect_identical(make_clone_id("cohort", "D1", c(1L, 2L)),
                   c("cohort~D1~1", "cohort~D1~2"))
  # The same integer in another donor, or the same donor in another scope, is a
  # different lineage and gets a different identifier.
  expect_false(identical(make_clone_id("cohort", "D1", 1L),
                         make_clone_id("cohort", "D2", 1L)))
  expect_false(identical(make_clone_id("cohort", "D1", 1L),
                         make_clone_id("project_a", "D1", 1L)))
})

test_that("clone names distinguish tuples containing underscores and numeric IDs", {
  scope <- c("study_A", "study", "study-1", "001")
  donor <- c("B", "A_B", "D.01", "001")
  ids <- make_clone_id(scope, donor, 1L)
  expect_identical(ids, c("study_A~B~1", "study~A_B~1", "study-1~D.01~1", "001~001~1"))
  expect_identical(anyDuplicated(ids), 0L)
})

# --------------------------------------------------------- the propagation ---

POOL <- data.frame(
  sequence_id = c("g~assemble1_0", "g~assemble1_1", "g~assemble2_0"),
  stringsAsFactors = FALSE)
# Rows 1 and 2 are the two CDR3s of one contig; row 1 is its representative.
REPRESENTATIVE <- c(1L, 1L, 3L)
# Deliberately in a different order than the pool, because scoper may reorder.
CLUSTERED <- data.frame(
  sequence_id = c("g~assemble2_0", "g~assemble1_0"),
  clone_id = c(2L, 1L), stringsAsFactors = FALSE)

test_that("a sibling inherits its representative's lineage", {
  expect_identical(
    propagate_clone_id_to_siblings(POOL, REPRESENTATIVE, CLUSTERED),
    c(1L, 1L, 2L))
})

test_that("the lookup is by identifier, not by row position", {
  reordered <- CLUSTERED[c(2, 1), , drop = FALSE]
  expect_identical(
    propagate_clone_id_to_siblings(POOL, REPRESENTATIVE, reordered),
    propagate_clone_id_to_siblings(POOL, REPRESENTATIVE, CLUSTERED))
})

test_that("an empty pool propagates nothing", {
  expect_identical(
    length(propagate_clone_id_to_siblings(POOL[0, , drop = FALSE], integer(0),
                                          CLUSTERED)),
    0L)
})

# Either failure would silently drop receptors out of every downstream IGH
# metric, so both stop the stage.
test_that("a representative the clusterer never returned is refused", {
  expect_error(
    propagate_clone_id_to_siblings(POOL, REPRESENTATIVE,
                                   CLUSTERED[1, , drop = FALSE]),
    "omitted an inference representative")
})

test_that("a blank or absent label is refused", {
  for (bad in list("", NA_character_)) {
    broken <- CLUSTERED
    broken$clone_id <- as.character(broken$clone_id)
    broken$clone_id[2] <- bad
    expect_error(propagate_clone_id_to_siblings(POOL, REPRESENTATIVE, broken),
                 "did not receive a propagated clone_id")
  }
})

# ----------------------------------------------------------- the membership ---
#
# Two CDR3s of one contig and one further contig. Read support arrives as text,
# because every table is read as character.

CLONE_PASS <- data.frame(
  sequence_id = c("g~assemble1_0", "g~assemble1_1", "g~assemble2_0"),
  library_id = "g", locus = "IGH",
  clone_id = c(1L, 1L, 2L),
  consensus_count = c("10", "1", "4"),
  v_call = "IGHV1-1*01",
  j_call = "IGHJ4*01",
  # The germline entity each call resolves to, as the reannotation stage wrote
  # it. The record shard's match key is built from these and not from the calls.
  v_gene = "IGHV1-1",
  j_gene = "IGHJ4",
  junction = c("ACGTACGTACGT", "ACGTACGTACGA", "TTTTACGTACGT"),
  junction_aa = c("TRTR", "TRTS", "FTRT"),
  junction_length = "12",
  sequence = paste0("GGGG", c("A", "B", "C")),
  # IMGT-gapped and equal-length, as IgBLAST returns them; carried through to the
  # record shard so a lineage can be handed to a tree without realignment.
  sequence_alignment = paste0("...ACGT", c("A", "B", "C")),
  germline_alignment_d_mask = "...ACGTN",
  stringsAsFactors = FALSE)

membership_mapping <- function(rows) {
  prepare_assembly_mapping(rows, "g", rep(TRUE, nrow(rows)), rep("eligible", nrow(rows)), TRUE)
}

describe("build_clone_membership_frame makes one row per receptor", {
  members <- build_clone_membership_frame(CLONE_PASS, membership_mapping(CLONE_PASS), "cohort", "D1")

  it("collapses the alternative CDR3s of one contig to their representative",
     expect_identical(members$sequence_id,
                      c("g~assemble1_0", "g~assemble2_0")))
  # The representative absorbs its siblings' support, so summed consensus_count over a
  # set of backbones is unchanged by the collapse.
  it("carries the sibling's read support over",
     expect_equal(members$consensus_count, c(11, 4)))
  it("keeps the total read support of the pool",
     expect_equal(sum(members$consensus_count),
                  sum(as.numeric(CLONE_PASS$consensus_count))))
  it("namespaces the lineage by scope and donor",
     expect_identical(members$clone_id, c("cohort~D1~1", "cohort~D1~2")))
  it("stamps the scope and the donor on every row", {
    expect_identical(unique(members$scope_id), "cohort")
    expect_identical(unique(members$donor_id), "D1")
  })
  it("keeps the representative's contig, not a sibling's",
     expect_identical(members$sequence, c("GGGGA", "GGGGC")))
  it("declares both tables' columns", {
    expect_true(all(CLONE_MEMBER_COLUMNS %in% names(members)))
    expect_true(all(CLONE_RECORD_COLUMNS %in% names(members)))
  })
  # The exact-nucleotide key is the record shard's cross-library match key, not
  # the clonal unit: two members of one lineage have different keys.
  # Raw UTF-8 separator bytes keep the expected key unlabelled under C collation.
  it("adds the library-independent match key", {
    expect_identical(
      members$clonotype_key,
      paste(c("ACGTACGTACGT", "TTTTACGTACGT"), "IGHV1-1", "IGHJ4",
            sep = "~"))
    expect_false(members$clonotype_key[1] ==
                   members$clonotype_key[2])
  })
})

# A member whose V or J call names several germline entities has no key, which is
# the same rule the non-IGH groups apply. It is not an exclusion here: the member
# keeps its lineage and its read support, and only the cross-dataset match key is
# withheld.
test_that("a member naming no single germline entity gets no match key", {
  unresolved <- CLONE_PASS
  unresolved$v_gene <- c("", "", "IGHV1-1")
  members <- build_clone_membership_frame(unresolved, membership_mapping(unresolved), "cohort", "D1")

  expect_identical(is.na(members$clonotype_key), c(TRUE, FALSE))
  expect_equal(members$consensus_count, c(11, 4))
})

test_that("the row order is deterministic and does not depend on the input order", {
  forwards <- build_clone_membership_frame(CLONE_PASS, membership_mapping(CLONE_PASS), "cohort", "D1")
  backwards <- build_clone_membership_frame(CLONE_PASS[3:1, ], membership_mapping(CLONE_PASS[3:1, ]), "cohort", "D1")
  expect_identical(forwards$sequence_id, backwards$sequence_id)
  expect_equal(forwards$consensus_count, backwards$consensus_count)
})

# A V-call disagreement between alternative CDR3 rows does not split their
# assembly backbone. Its representative receives the eligible source support.
test_that("a V-discordant backbone still collapses to one row", {
  discordant <- CLONE_PASS
  discordant$v_call <- c("IGHV1-1*01", "IGHV2-1*01", "IGHV1-1*01")
  members <- build_clone_membership_frame(discordant, membership_mapping(discordant), "cohort", "D1")
  expect_identical(members$sequence_id, c("g~assemble1_0", "g~assemble2_0"))
  expect_equal(sort(members$consensus_count), c(4, 11))
})

test_that("unreadable read support is rejected rather than treated as zero", {
  broken <- CLONE_PASS
  broken$consensus_count[1] <- "many"
  expect_error(build_clone_membership_frame(broken, membership_mapping(broken), "cohort", "D1"),
               "consensus_count contains 1 invalid value")
})

# --------------------------------------------- IGH lineage junction eligibility ---

describe("clone_junction_eligible with no missing symbols allowed", {
  cases <- list(
    list(junction = "ACGT", expected = TRUE,
         why = "an unambiguous junction is eligible"),
    list(junction = "ACGN", expected = FALSE,
         why = "one N is one too many at max_n = 0"),
    list(junction = "", expected = FALSE,
         why = "an empty junction is not eligible"),
    list(junction = NA_character_, expected = FALSE,
         why = "a missing junction is not eligible")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      clone_junction_eligible(case$junction, max_n = 0), case$expected))
  }
})

# The accepted alphabet mirrors scoper::hierarchicalClones(IUPAC = FALSE): only
# uppercase ATCG plus the two symbols it treats as missing.
describe("clone_junction_eligible with one missing symbol allowed", {
  cases <- list(
    list(junction = "ACGN", expected = TRUE,
         why = "N is an accepted missing symbol and one is within budget"),
    list(junction = "ACG?", expected = TRUE,
         why = "? is the other accepted missing symbol"),
    list(junction = "ACNN", expected = FALSE,
         why = "two missing symbols exceed the budget"),
    list(junction = "ACGR", expected = FALSE,
         why = "an IUPAC ambiguity code is outside the accepted alphabet"),
    list(junction = "acgt", expected = FALSE,
         why = "lowercase is outside the accepted alphabet")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      clone_junction_eligible(case$junction, max_n = 1), case$expected))
  }
})
