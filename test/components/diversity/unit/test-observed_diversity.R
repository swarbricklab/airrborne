# Check diversity formulas and released group support vectors independently.
# Run: python test/run_tests.py --component diversity --level unit

load_test_modules(c("table_io", "diversity"))

EXACT <- 1e-12

test_that("the output column order is fixed", {
  row <- compute_diversity_row(
    data.frame(consensus_count = rep(1, 12)), "sample_equal", "TRB", "clonotype"
  )
  expect_identical(names(row), DIVERSITY_COLS)
})

# Twelve equally supported clones make every metric analytically transparent and
# exercise both cumulative-rank cutoffs.
describe("an evenly supported repertoire of twelve clones", {
  equal <- compute_diversity_row(
    data.frame(consensus_count = rep(1, 12)), "sample_equal", "TRB", "clonotype"
  )

  it("sums read support as consensus_count",
     expect_equal(equal$consensus_count, 12, tolerance = EXACT))
  it("counts every clone",
     expect_identical(equal$richness, 12L))
  it("has Shannon entropy log(n)",
     expect_equal(equal$shannon_entropy, log(12), tolerance = EXACT))
  it("has inverse Simpson equal to n",
     expect_equal(equal$inverse_simpson, 12, tolerance = EXACT))
  it("has zero clonality",
     expect_equal(equal$clonality, 0, tolerance = EXACT))
  it("splits support evenly at rank 1",
     expect_equal(equal$top_1_support_fraction, 1 / 12, tolerance = EXACT))
  it("splits support evenly at rank 5",
     expect_equal(equal$top_5_support_fraction, 5 / 12, tolerance = EXACT))
  it("splits support evenly at rank 10",
     expect_equal(equal$top_10_support_fraction, 10 / 12, tolerance = EXACT))
})

# Uneven positive counts independently exercise inverse Simpson and cumulative
# dominance. Every clone here is supported, because the callers that aggregate
# clones are what decide which ones count as observed.
describe("an unevenly supported repertoire", {
  uneven <- compute_diversity_row(
    data.frame(consensus_count = c(4, 2, rep(1, 6))),
    "sample_uneven", "IGH", "clone"
  )

  it("sums read support as consensus_count",
     expect_equal(uneven$consensus_count, 12, tolerance = EXACT))
  it("counts every clone it was given",
     expect_identical(uneven$richness, 8L))
  it("has inverse Simpson 72/13",
     expect_equal(uneven$inverse_simpson, 72 / 13, tolerance = EXACT))
  it("gives the top clone a third of the support",
     expect_equal(uneven$top_1_support_fraction, 1 / 3, tolerance = EXACT))
  it("gives the top five three quarters",
     expect_equal(uneven$top_5_support_fraction, 3 / 4, tolerance = EXACT))
  # Top-N sums up to N groups. These eight groups contain all support.
  it("reports all support when the repertoire has fewer than N groups",
     expect_equal(uneven$top_10_support_fraction, 1))
})

# The three indices are hand-written, so they are also pinned to the Hill-number
# parameterisation that Hill 1973 (doi:10.2307/1934352) defines and that the
# alakazam diversity vignette uses. The Hill number of order q is the number of
# equally supported clones a repertoire of that diversity would hold, and the
# forms below are arithmetically different from the ones the metric computes: a
# product instead of a sum of logarithms at q = 1, and a reciprocal power at
# q = 2. That is what makes them a check rather than a restatement. In
# particular, exp(shannon_entropy) is the q = 1 Hill number only when the
# entropy is in nats, so a switch to log2 or log10 would fail here.
#
# alakazam::alphaDiversity() computes the same profile but by resampling the
# counts to a uniform depth, so it cannot be used to pin exact values; the
# identities are asserted directly instead.
describe("the indices against the Hill-number definitions", {
  support <- c(37, 21, 13, 8, 5, 3, 2, 1, 1, 1)
  row <- compute_diversity_row(
    data.frame(consensus_count = support), "sample_hill", "TRB", "clonotype"
  )
  p <- support / sum(support)

  it("reports the q = 0 Hill number as richness",
     expect_identical(row$richness, length(support)))
  it("has exp(Shannon) equal to the q = 1 Hill number, so the entropy is in nats", {
    expect_equal(exp(row$shannon_entropy), prod(p^-p), tolerance = EXACT)
  })
  it("has inverse Simpson equal to the q = 2 Hill number", {
    expect_equal(row$inverse_simpson, sum(p^2)^(1 / (1 - 2)), tolerance = EXACT)
  })
  # Pielou evenness is the ratio of the q = 1 and q = 0 Hill numbers on a log
  # scale, and clonality is one minus it.
  it("has clonality equal to one minus Pielou evenness", {
    expect_equal(row$clonality,
                 1 - log(prod(p^-p)) / log(length(support)),
                 tolerance = EXACT)
  })
  it("orders the three Hill numbers q = 0 >= q = 1 >= q = 2", {
    expect_true(row$richness >= exp(row$shannon_entropy))
    expect_true(exp(row$shannon_entropy) >= row$inverse_simpson)
  })
})

# Rank fractions include all groups when fewer than N groups are observed.
# Check the five- and ten-group boundaries and one group below each.
describe("a repertoire at each rank cutoff", {
  fraction_at <- function(n_clones, column) {
    compute_diversity_row(
      data.frame(consensus_count = rep(1, n_clones)),
      "sample_boundary", "TRB", "clonotype"
    )[[column]]
  }

  it("sums up to five groups at five clones and at four", {
    expect_equal(fraction_at(5, "top_5_support_fraction"), 1, tolerance = EXACT)
    expect_equal(fraction_at(4, "top_5_support_fraction"), 1)
  })
  it("sums up to ten groups at ten clones and at nine", {
    expect_equal(fraction_at(10, "top_10_support_fraction"), 1, tolerance = EXACT)
    expect_equal(fraction_at(9, "top_10_support_fraction"), 1)
  })
  it("keeps the largest-group fraction defined throughout",
     expect_equal(fraction_at(1, "top_1_support_fraction"), 1, tolerance = EXACT))
})

# Clonality is documented and consumed as a proportion, so it has to stay inside
# [0, 1] for every repertoire and not merely for most of them. Shannon is summed
# independently of log(richness), so a perfectly even repertoire lands a fraction
# of an ulp either side of zero depending on the clone count, and the wide scipen
# setting writes that out in full rather than rounding it away.
describe("clonality across the range of shapes a repertoire can take", {
  clonality_of <- function(support) {
    compute_diversity_row(data.frame(consensus_count = support),
                          "sample_bounds", "TRG", "clonotype"
    )$clonality
  }

  it("is zero at equal abundances, and never below it, for every richness", {
    for (n_clones in c(2L, 3L, 5L, 7L, 12L, 100L)) {
      value <- clonality_of(rep(1, n_clones))
      expect_equal(value, 0, tolerance = EXACT)
      expect_true(value >= 0)
    }
  })

  # Five equally supported clones is the shape that put a negative value in a
  # released table, because Shannon sums to fractionally more than log(5).
  it("returns an exact zero where the arithmetic falls fractionally short",
     expect_identical(clonality_of(rep(1, 5)), 0))

  it("stays inside the unit interval however the support is distributed", {
    for (support in list(rep(1, 5), c(1000, 1, 1), c(1, 1, 1, 1, 500),
                         c(2, 2, 2, 2), c(1e6, 1))) {
      value <- clonality_of(support)
      expect_true(value >= 0 && value <= 1)
    }
  })

  it("rises towards one as a single clone takes over", {
    expect_true(clonality_of(c(1e6, 1, 1, 1)) > clonality_of(c(4, 3, 2, 1)))
    expect_true(clonality_of(c(4, 3, 2, 1)) > clonality_of(rep(1, 4)))
  })

  it("is undefined below two clones",
     expect_identical(clonality_of(5), NA_real_))
})

describe("a repertoire holding one clone", {
  single <- compute_diversity_row(
    data.frame(consensus_count = 5), "sample_single", "IGK", "clonotype"
  )

  it("has inverse Simpson one",
     expect_equal(single$inverse_simpson, 1, tolerance = EXACT))
  # Clonality is one minus Pielou evenness, and evenness needs at least two
  # groups before it means anything.
  it("leaves clonality undefined",
     expect_identical(single$clonality, NA_real_))
  it("gives the one clone all the support",
     expect_equal(single$top_1_support_fraction, 1, tolerance = EXACT))
  it("gives the wider cumulative-rank proportions all support", {
    expect_identical(
      unname(unlist(single[c("top_5_support_fraction",
                             "top_10_support_fraction")])),
      rep(1, 2)
    )
  })
})

describe("a repertoire with no eligible clones", {
  empty <- make_empty_diversity_row("sample_empty", "TRG", "clonotype")

  it("has zero consensus_count", expect_equal(empty$consensus_count, 0, tolerance = EXACT))
  it("has zero richness", expect_identical(empty$richness, 0L))
  it("leaves every other metric undefined", {
    undefined <- c("shannon_entropy", "inverse_simpson",
                   "clonality", "top_1_support_fraction",
                   "top_5_support_fraction", "top_10_support_fraction")
    expect_equal(unname(unlist(empty[undefined])), rep(NA_real_, 6))
  })
})

# Callers select the observed groups and supply positive support for each one.
# Diversity calculation rejects unsupported groups to keep richness consistent
# with the supplied population.
describe("a clone with no read support", {
  it("stops rather than being left out of the metric", {
    expect_error(
      compute_diversity_row(
        data.frame(consensus_count = c(5, 0)), "sample_one_zero", "TRD",
        "clonotype"),
      regexp = "finite positive read support")
  })

  it("stops even when every clone is unsupported", {
    expect_error(
      compute_diversity_row(
        data.frame(consensus_count = c(0, 0)), "sample_all_zero", "TRD",
        "clonotype"),
      regexp = "finite positive read support")
  })
})

# ------------------------- one library's share of a donor-pooled lineage ---
#
# IGH lineages span every sample a donor gave to a scope, so a library's
# repertoire is its own read support within lineages it may share with another
# library. These cases fix which members enter that sum and what the row that
# comes out looks like. Read support arrives as text, because every table is
# read as character.

# One lineage whose members voted a gene and one whose resolved members named
# two, so both a counted and an excluded lineage reach the aggregation.
POOLED_CALLS <- data.frame(
  clone_id = c("cohort_D1_L1", "cohort_D1_L2"),
  v_gene = c("IGHV1-69", ""),
  j_gene = c("IGHJ6", ""),
  v_vote_status = c("assigned", "resolved_members_disagree"),
  j_vote_status = c("assigned", "resolved_members_disagree"),
  stringsAsFactors = FALSE
)

POOLED_MEMBERSHIP <- data.frame(
  library_id = c("WTS_1_S1", "WTS_1_S1", "WTS_1_S1", "WTS_2_S1", "WTS_1_S1"),
  clone_id = c("cohort_D1_L2", "cohort_D1_L2", "cohort_D1_L1",
                 "cohort_D1_L1", ""),
  consensus_count = c("3", "4", "5", "99", "7"),
  junction = c("TGTGCGTGG", "TGTGCAAGG", "TGTACCTGG", "TGTAAATGG", "TGTCCCTGG"),
  stringsAsFactors = FALSE
)

# What aggregate_clone_support_per_library() returns, asserted in one place so
# a column added to it cannot pass by being added to each expectation too.
AGGREGATED_LINEAGE_COLUMNS <- c("clone_id", "consensus_count",
                                "n_members", "min_junction", "v_gene",
                                "j_gene", "v_vote_status",
                                "j_vote_status")

describe("a library that contributes to two pooled lineages", {
  lineages <- aggregate_clone_support_per_library(
    POOLED_MEMBERSHIP, "WTS_1_S1", POOLED_CALLS)

  it("reports one row per lineage the library reached",
     expect_identical(lineages$clone_id, c("cohort_D1_L1", "cohort_D1_L2")))
  # 5 for L1; 3 + 4 for L2, with the other library's 99 reads never in scope.
  it("sums only this library's member support",
     expect_equal(lineages$consensus_count, c(5, 7), tolerance = EXACT))
  it("carries the scope-wide gene vote over", {
    expect_identical(lineages$v_gene, c("IGHV1-69", ""))
    expect_identical(lineages$v_vote_status,
                     c("assigned", "resolved_members_disagree"))
  })
  # L1 came from one member, L2 from two, so a lineage's read support and the
  # number of receptors carrying it are separately reported.
  it("counts the receptors each lineage was seen as",
     expect_identical(lineages$n_members, c(1L, 2L)))
  # L2's members in this library carry TGTGCGTGG and TGTGCAAGG; the other
  # library's L1 junction does not enter.
  it("takes each lineage's smallest junction in this library",
     expect_identical(lineages$min_junction, c("TGTACCTGG", "TGTGCAAGG")))
  it("fixes the column order", {
    expect_identical(names(lineages), AGGREGATED_LINEAGE_COLUMNS)
  })
})

describe("a library holding one eligible member", {
  lineages <- aggregate_clone_support_per_library(
    POOLED_MEMBERSHIP, "WTS_2_S1", POOLED_CALLS)

  it("reports the one lineage", expect_identical(nrow(lineages), 1L))
  it("reports its support",
     expect_equal(lineages$consensus_count, 99, tolerance = EXACT))
})

test_that("several foreign-library members cannot enter the requested support", {
  foreign <- POOLED_MEMBERSHIP[c(1, 3, 4), ]
  foreign$library_id <- c("foreign_a", "foreign_b", "foreign_a")
  foreign$consensus_count <- c("100", "1000", "10000")
  members <- rbind(POOLED_MEMBERSHIP, foreign)
  for (order in list(seq_len(nrow(members)), rev(seq_len(nrow(members))))) {
    actual <- aggregate_clone_support_per_library(members[order, ], "WTS_1_S1", POOLED_CALLS)
    expect_identical(actual$clone_id, c("cohort_D1_L1", "cohort_D1_L2"))
    expect_equal(actual$consensus_count, c(5, 7), tolerance = 0)
    expect_identical(actual$n_members, c(1L, 2L))
  }
})

describe("a library with nothing to report", {
  for (case in list(
    list(what = "not in the membership table at all",
         membership = POOLED_MEMBERSHIP, library_id = "WTS_9_S1"),
    list(what = "no membership table at all",
         membership = data.frame(), library_id = "WTS_1_S1"),
    list(what = "every member unassigned to a lineage",
         membership = POOLED_MEMBERSHIP[5, , drop = FALSE],
         library_id = "WTS_1_S1"))) {
    lineages <- aggregate_clone_support_per_library(
      case$membership, case$library_id, POOLED_CALLS)

    it(sprintf("returns no rows [%s]", case$what),
       expect_identical(nrow(lineages), 0L))
    it(sprintf("still declares the aggregated columns [%s]", case$what), {
      expect_identical(names(lineages), AGGREGATED_LINEAGE_COLUMNS)
    })
    it(sprintf("hands compute_diversity_row an empty repertoire [%s]", case$what), {
      expect_equal(
        compute_diversity_row(lineages, case$library_id, "IGH", "clone"),
        make_empty_diversity_row(case$library_id, "IGH", "clone"),
        tolerance = EXACT)
    })
  }
})

# A lineage with no entry in the scope-wide summary leaves the gene and its vote
# missing rather than dropping the lineage out of the repertoire.
test_that("a lineage the scope summary does not know keeps its support", {
  lineages <- aggregate_clone_support_per_library(
    POOLED_MEMBERSHIP, "WTS_1_S1", POOLED_CALLS[1, , drop = FALSE])

  expect_equal(lineages$consensus_count, c(5, 7), tolerance = EXACT)
  expect_identical(lineages$v_gene, c("IGHV1-69", NA_character_))
  expect_identical(lineages$v_vote_status, c("assigned", NA_character_))
})

test_that("unreadable read support is rejected rather than treated as zero", {
  broken <- POOLED_MEMBERSHIP
  broken$consensus_count[1] <- "many"
  expect_error(
    aggregate_clone_support_per_library(broken, "WTS_1_S1", POOLED_CALLS),
    "consensus_count contains 1 invalid value")
})

# ------------------------- the vector a diversity row was computed from ---
#
# The point of releasing the clone-size vector is that a reader can reproduce the
# published indices from it and then compute one this pipeline does not. That
# only holds if the vector is the same object the index was taken over, so the
# assertions below recompute the released numbers from the emitted columns.

ABUNDANCE_INPUT <- data.frame(
  group_key = c("K1", "K2", "K3"),
  n_members = c(3L, 1L, 1L),
  consensus_count = c("60", "30", "10"),
  stringsAsFactors = FALSE
)

describe("build_group_abundance_rows beside compute_diversity_row", {
  rows <- build_group_abundance_rows(ABUNDANCE_INPUT, "WTS_1_S1", "TRB",
                                     "clonotype", "group_key")
  diversity <- compute_diversity_row(ABUNDANCE_INPUT, "WTS_1_S1", "TRB",
                                     "clonotype")

  it("declares the released columns in order",
     expect_identical(names(rows), CLONE_ABUNDANCE_COLS))
  it("emits one row per observed clone",
     expect_identical(nrow(rows), diversity$richness))
  it("names the unit the key identifies",
     expect_identical(unique(rows$group_unit), "clonotype"))
  it("carries the key through under one column name",
     expect_identical(rows$group_key, c("K1", "K2", "K3")))
  it("coerces the read support the metric coerced",
     expect_equal(rows$consensus_count, c(60, 30, 10), tolerance = EXACT))
  it("keeps the member count as a separate unit",
     expect_identical(rows$n_members, c(3L, 1L, 1L)))
  it("ranks the groups by support",
     expect_identical(rows$group_rank, c(1L, 2L, 3L)))
  it("takes shares over the consensus_count the diversity row reports",
     expect_equal(sum(rows$consensus_count), diversity$consensus_count,
                  tolerance = EXACT))
  it("emits shares that sum to one",
     expect_equal(sum(rows$support_fraction), 1, tolerance = EXACT))
  # The three identities a reader re-deriving the released table relies on.
  it("reproduces the published Shannon entropy",
     expect_equal(-sum(rows$support_fraction * log(rows$support_fraction)),
                  diversity$shannon_entropy, tolerance = EXACT))
  it("reproduces the published inverse Simpson",
     expect_equal(1 / sum(rows$support_fraction^2), diversity$inverse_simpson,
                  tolerance = EXACT))
  it("reproduces the published top-1 dominance",
     expect_equal(max(rows$support_fraction), diversity$top_1_support_fraction,
                  tolerance = EXACT))
})

test_that("a library with no clone at a locus still declares the columns", {
  rows <- build_group_abundance_rows(ABUNDANCE_INPUT[0, , drop = FALSE],
                                     "WTS_1_S1", "TRD", "clonotype",
                                     "group_key")
  expect_identical(nrow(rows), 0L)
  expect_identical(names(rows), CLONE_ABUNDANCE_COLS)
  expect_identical(vapply(rows[c("group_rank", "n_members", "consensus_count",
                                 "support_fraction")], class, character(1)),
                   c(group_rank = "integer", n_members = "integer",
                     consensus_count = "numeric", support_fraction = "numeric"))
})

# Tied groups share the smaller rank and follow their content, so the order is
# the same however the caller supplied them. Sorted as text, 9 would precede
# 100 and 10; the rank compares numbers.
describe("group_rank with tied support", {
  tied <- data.frame(
    group_key = c("K_c", "K_a", "K_d", "K_f", "K_b", "K_e"),
    n_members = 1L,
    consensus_count = c("9", "9", "100", "1", "9", "10"),
    stringsAsFactors = FALSE
  )
  for (rows_in in list(seq_len(6), c(6, 3, 1, 4, 2, 5))) {
    rows <- build_group_abundance_rows(tied[rows_in, ], "WTS_1_S1", "TRB",
                                       "clonotype", "group_key")
    it("shares the smaller rank among equal support",
       expect_identical(rows$group_rank, c(1L, 2L, 3L, 3L, 3L, 6L)))
    it("orders tied groups by key",
       expect_identical(rows$group_key, c("K_d", "K_e", "K_a", "K_b", "K_c", "K_f")))
  }
  # Ranks 1 to 5 hold five groups here, and the top-5 fraction sums the same
  # support whichever tied group a reader would pick.
  it("leaves the dominance fractions to the support values", {
    rows <- build_group_abundance_rows(tied, "WTS_1_S1", "TRB", "clonotype",
                                       "group_key")
    diversity <- compute_diversity_row(tied, "WTS_1_S1", "TRB", "clonotype")
    expect_equal(sum(rows$support_fraction[rows$group_rank <= 5]),
                 diversity$top_5_support_fraction, tolerance = EXACT)
  })
})

# Both arms reach the same function, so a lineage vector and a rearrangement-group
# vector cannot end up with different columns or a different denominator.
test_that("the IGH arm produces the same shape from a lineage aggregation", {
  lineages <- aggregate_clone_support_per_library(
    POOLED_MEMBERSHIP, "WTS_1_S1", POOLED_CALLS)
  rows <- build_group_abundance_rows(lineages, "WTS_1_S1", "IGH", "clone",
                                     "clone_id")

  expect_identical(names(rows), CLONE_ABUNDANCE_COLS)
  # L2 carries 7 and L1 5, so L2 ranks first.
  expect_identical(rows$group_key, c("cohort_D1_L2", "cohort_D1_L1"))
  expect_identical(rows$group_rank, c(1L, 2L))
  expect_identical(rows$n_members, c(2L, 1L))
  expect_equal(rows$support_fraction, c(7, 5) / 12, tolerance = EXACT)
})

# Equal IGH support follows the smallest member junction, which does not
# depend on SCOPer's cluster numbering, and clone_id only after that.
test_that("tied IGH clones follow their smallest junction, then clone_id", {
  lineages <- data.frame(
    clone_id = c("cohort~D1~1", "cohort~D1~2", "cohort~D1~10"),
    consensus_count = c(4, 4, 4), n_members = 1L,
    min_junction = c("TGTCCC", "TGTAAA", "TGTAAA"),
    stringsAsFactors = FALSE
  )
  rows <- build_group_abundance_rows(lineages, "WTS_1_S1", "IGH", "clone",
                                     "clone_id",
                                     order_columns = c("min_junction", "clone_id"))
  expect_identical(rows$group_rank, c(1L, 1L, 1L))
  expect_identical(rows$group_key, c("cohort~D1~10", "cohort~D1~2", "cohort~D1~1"))
})

test_that("target values equal their curve and support minima leave diagnostics", {
  groups <- data.frame(consensus_count = c(2, 1))
  result <- build_repertoire_diversity(groups, "0460", "TRB", "clonotype", 2)
  target <- result$curve[result$curve$is_target, ]
  expect_equal(result$detail$rarefied_richness, 5 / 3)
  expect_equal(result$detail$rarefied_inverse_simpson, 1.5)
  expect_equal(result$detail$rarefied_hill_shannon, 2^(2 / 3))
  expect_equal(result$detail$rarefied_inverse_simpson, target$inverse_simpson)
  expect_equal(result$detail$top_10_support_fraction, 1)
  expect_identical(result$detail$library_id, "0460")
  expect_identical(names(result$detail), DIVERSITY_DETAIL_COLUMNS)
  expect_equal(result$detail$hill_shannon, exp(result$detail$shannon_entropy))
  expect_identical(result$detail$min_consensus_count_for_diversity, 0)
  gated <- build_repertoire_diversity(groups, "0460", "TRB", "clonotype", 2, 4)
  expect_equal(gated$curve, result$curve)
  expect_equal(gated$detail$richness, 2)
  expect_equal(gated$detail$grouped_consensus_count, 3)
  expect_identical(gated$detail$min_consensus_count_for_diversity, 4)
  # The asymptote is context, not gated by the minimum. The corrected
  # concentration is (2 * 1) / (3 * 2) = 1/3, so the asymptote is 3.
  expect_equal(gated$detail$inverse_simpson_asymptote, 3)
  expect_identical(gated$detail$asymptote_status, "available")
  expect_identical(gated$detail$diversity_status, "below_minimum_support")
  expect_true(all(is.na(gated$detail[c("hill_shannon", "inverse_simpson",
    "top_1_support_fraction", "top_5_support_fraction", "top_10_support_fraction",
    "rarefied_richness", "rarefied_hill_shannon", "rarefied_inverse_simpson")])) )
})

test_that("below-target and empty repertoires retain honest observed context", {
  groups <- data.frame(consensus_count = c(2, 1))
  result <- build_repertoire_diversity(groups, "NA", "IGH", "clone", 1000)
  expect_identical(result$detail$rarefaction_status, "below_target_support")
  expect_equal(result$detail$inverse_simpson, 9 / 5)
  expect_equal(max(result$curve$support_depth), 3)
  expect_false(any(result$curve$is_target))
  empty <- build_repertoire_diversity(groups[FALSE, , drop = FALSE], "NA", "IGH", "clone", 1000)
  expect_equal(empty$detail$richness, 0)
  expect_equal(empty$detail$grouped_consensus_count, 0)
  expect_identical(empty$detail$asymptote_status, "insufficient_support")
  expect_true(is.na(empty$detail$hill_shannon))
  expect_true(is.na(empty$detail$top_5_support_fraction))
  expect_equal(empty$curve$support_depth, 0)
  expect_equal(nrow(build_repertoire_diversity(groups, "a", "TRA", "clonotype")$curve), 0)
})
