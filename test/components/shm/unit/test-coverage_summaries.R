# Check regional evidence selection and equally weighted clone SHM summaries.
# Run: python test/run_tests.py --component shm --level unit

load_test_modules(c(
  "table_io", "gene_calls", "isotype_resolution", "assembly_backbone", "qc_audit",
  "shm_summaries"))

EXACT <- 1e-12

MIN_SUBREGION_FRACTION <- 1.0
MIN_GERMLINE_OCCUPIED_V <- 270L

# These constructed regions match the IMGT order. Their occupancy below
# separates summary eligibility from the measurement-stage reference.
REGION_LEVELS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

# What each sub-region's germline occupies in this fixture. They sum to 278, which
# clears the whole-V occupancy guard, so that guard is exercised separately rather
# than shadowing every other case.
OCCUPIED <- c(fwr1 = 74L, cdr1 = 30L, fwr2 = 30L, cdr2 = 30L, fwr3 = 114L)

# Membership as the stage sees it: the feature sidecar's calls and counts already
# joined on. Read down the members --
#   seq1  covers everything            IGHG1        -> IgG and IgG1
#   seq2  covers everything but FWR1   IGHG3        -> IgG and IgG3, no whole_v
#   seq3  covers everything            IGHM         -> IgM, and no subclass row
#   seq4  covers everything            IGHG1,IGHG3  -> IgG, but no subclass: this is
#                                                     a C alignment too short to
#                                                     separate the IgG subclasses
#   seq5  covers almost nothing        IGHM         -> measured, never eligible
#   seq6  covers everything, unmutated no call      -> all_IGH only
#   seq7  covers everything            IGHA2        -> IgA and IgA2
JOINED_MEMBERSHIP <- data.frame(
  donor_id = c("D1", "D1", "D1", "D1", "D1", "D2", "D1"),
  library_id = c("S1", "S1", "S1", "S1", "S2", "S3", "S1"),
  sequence_id = paste0("S", c(1, 1, 1, 1, 2, 3, 1), "~seq", 1:7),
  clone_id = c("L1", "L1", "L2", "L3", "L1", "L4", "L5"),
  consensus_count = c("3", "4", "5", "6", "7", "8", "9"),
  qualified_c_call = c("IGHG1*01", "IGHG3*01", "IGHM*01", "IGHG1,IGHG3", "IGHM*01", "",
             "IGHA2*01"),
  shm_frequency = c(0.10, 0.30, 0.20, 0.60, 0.90, 0.00, 0.50),
  informative_v = c(278L, 214L, 278L, 278L, 100L, 278L, 278L),
  stringsAsFactors = FALSE
)
# Each member's vote support differs between sub-regions, and the five offsets are
# deliberately lopsided. Every sub-region carrying the same value would make the
# whole-V rule untestable: the mean of five equal numbers, their median, their
# smallest and their occupancy-weighted mean are all that same number, so no
# assertion could tell those rules apart. These offsets sum to +10, so the mean of
# a member's five is its base plus 2, while their median is the base itself and
# their occupancy-weighted mean is about the base plus 5.5.
VOTE_BASE <- c(4, 4, 10, 20, 6, 8, 30)
VOTE_OFFSET <- c(fwr1 = -10, cdr1 = -5, fwr2 = 0, cdr2 = 5, fwr3 = 20)

# Every sub-region is fully covered and fully supported, then the two exceptions are
# written over the top: seq2 reached only 10 of FWR1's 74 germline positions, and
# seq5 reached nothing but part of FWR3.
for (region in REGION_LEVELS) {
  JOINED_MEMBERSHIP[[paste0("shm_frequency_", region)]] <- JOINED_MEMBERSHIP$shm_frequency
  JOINED_MEMBERSHIP[[paste0("mutations_", region)]] <- c(5L, 5L, 5L, 5L, 5L, 0L, 5L)
  JOINED_MEMBERSHIP[[paste0("informative_", region)]] <- OCCUPIED[[region]]
  JOINED_MEMBERSHIP[[paste0("informative_supported_", region)]] <- OCCUPIED[[region]]
  JOINED_MEMBERSHIP[[paste0("germline_occupied_", region)]] <- OCCUPIED[[region]]
  JOINED_MEMBERSHIP[[paste0("median_majority_votes_", region)]] <-
    VOTE_BASE + VOTE_OFFSET[[region]]
}
JOINED_MEMBERSHIP$informative_supported_fwr1[2] <- 10L
JOINED_MEMBERSHIP$informative_fwr1[2] <- 10L
for (region in setdiff(REGION_LEVELS, "fwr3")) {
  JOINED_MEMBERSHIP[[paste0("informative_supported_", region)]][5] <- 0L
  JOINED_MEMBERSHIP[[paste0("informative_", region)]][5] <- 0L
}
JOINED_MEMBERSHIP$informative_supported_fwr3[5] <- 100L
JOINED_MEMBERSHIP$informative_fwr3[5] <- 100L

MANIFEST <- data.frame(
  library_id = c("S1", "S2", "S3"),
  donor_id = c("D1", "D1", "D2"),
  tissue_id = c("T1", "T2", "T3"),
  stringsAsFactors = FALSE
)

V_REGIONS <- shm_v_regions(REGION_LEVELS)

eligible_for <- function(v_region, membership = JOINED_MEMBERSHIP) {
  select_shm_eligible_members(membership, v_region, MIN_SUBREGION_FRACTION,
                              MIN_GERMLINE_OCCUPIED_V, REGION_LEVELS)
}

# --------------------------------------------------- the regions reported ---

test_that("the whole V region is reported first, then the sub-regions in order", {
  expect_identical(V_REGIONS, c("whole_v", REGION_LEVELS))
})

test_that("class rows precede subclass rows without requiring a partition", {
  expect_identical(SHM_ISOTYPE_GROUPS,
                   c("all_IGH", "IgM", "IgG", "IgA",
                     "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2"))
})

describe("the columns each region is measured from", {
  it("reads the whole V region off all five sub-regions", {
    columns <- shm_region_columns("whole_v", REGION_LEVELS)
    expect_identical(columns$rate, "shm_frequency")
    expect_identical(columns$supported,
                     paste0("informative_supported_", REGION_LEVELS))
    expect_identical(columns$germline_occupied,
                     paste0("germline_occupied_", REGION_LEVELS))
  })
  it("reads a sub-region off itself alone", {
    columns <- shm_region_columns("cdr2", REGION_LEVELS)
    expect_identical(columns$rate, "shm_frequency_cdr2")
    expect_identical(columns$supported, "informative_supported_cdr2")
    expect_identical(columns$germline_occupied, "germline_occupied_cdr2")
  })
})

# ------------------------------------------------------- who is eligible ---

describe("the members that can carry a whole-V mutation frequency", {
  eligible <- eligible_for("whole_v")

  # seq2 measured four sub-regions of five, which is what the conjunction excludes.
  it("keeps only the members that covered every sub-region",
     expect_identical(eligible$sequence_id,
                      c("S1~seq1", "S1~seq3", "S1~seq4", "S3~seq6", "S1~seq7")))
  it("carries the fields the summary reads", {
    expect_identical(names(eligible),
                     c("library_id", "sequence_id", "clone_id", "qualified_c_call",
                       "shm_frequency", "mutations", "median_trust4_consensus_base_votes",
                       "isotype_class", "isotype_subclass"))
  })
  it("sums the mutation count over the five sub-regions",
     expect_identical(eligible$mutations, c(25, 25, 25, 0, 25)))
  it("resolves an unambiguous call at both levels", {
    expect_identical(eligible$isotype_class[1], "IgG")
    expect_identical(eligible$isotype_subclass[1], "IgG1")
  })
  # IGHM is the only subclass of its class, so it has a class row and no subclass
  # row rather than two rows over one set of sequences.
  it("gives IgM no subclass row", {
    expect_identical(eligible$isotype_class[2], "IgM")
    expect_identical(eligible$isotype_subclass[2], NA_character_)
  })
  # A constant-region alignment too short to separate IGHG1 from IGHG3 ties across
  # both. That is still unambiguously IgG, and deliberately nothing below it.
  it("resolves a within-class tie to the class and no further", {
    expect_identical(eligible$isotype_class[3], "IgG")
    expect_identical(eligible$isotype_subclass[3], NA_character_)
  })
  it("leaves an absent call outside every group", {
    expect_identical(eligible$isotype_class[4], NA_character_)
    expect_identical(eligible$isotype_subclass[4], NA_character_)
  })
})

describe("the members that can carry a sub-region mutation frequency", {
  it("excludes the member that covered a tenth of the sub-region",
     expect_false("S1~seq2" %in% eligible_for("fwr1")$sequence_id))
  # The same member is summarised in FWR3, which it did cover. This is the reason
  # the sub-regions are reported separately: a five-prime truncated assembly
  # measures the three-prime end fully and the five-prime end barely.
  it("admits that member where it did cover the sub-region",
     expect_true("S1~seq2" %in% eligible_for("fwr3")$sequence_id))
  it("admits every member that covered the three-prime sub-region",
     expect_identical(eligible_for("fwr3")$sequence_id,
                      c("S1~seq1", "S1~seq2", "S1~seq3", "S1~seq4", "S3~seq6",
                        "S1~seq7")))
  it("declines a sub-region whose germline occupies nothing", {
    unoccupied <- JOINED_MEMBERSHIP
    unoccupied$germline_occupied_fwr1 <- 0L
    expect_identical(nrow(eligible_for("fwr1", unoccupied)), 0L)
  })
})

# Observing a position is not the same as observing it on reads. These two members
# measured the whole sub-region, and one of them measured it on single votes.
test_that("an observed position still has to be read-supported", {
  thin <- JOINED_MEMBERSHIP
  thin$informative_supported_cdr2 <- c(30L, 30L, 29L, 30L, 30L, 30L, 30L)
  expect_false("S1~seq3" %in% eligible_for("cdr2", thin)$sequence_id)
  expect_true("S1~seq1" %in% eligible_for("cdr2", thin)$sequence_id)
  # And it removes the same member from the whole V region, because that is the
  # conjunction of the five and CDR2 is one of them.
  expect_false("S1~seq3" %in% eligible_for("whole_v", thin)$sequence_id)
})

# Without the guard a record whose own germline alignment is truncated satisfies
# "covered all of what the germline occupies" on almost nothing.
test_that("a truncated germline cannot satisfy whole-V coverage trivially", {
  truncated <- JOINED_MEMBERSHIP
  for (region in REGION_LEVELS) {
    truncated[[paste0("germline_occupied_", region)]][1] <- 20L
    truncated[[paste0("informative_supported_", region)]][1] <- 20L
  }
  expect_false("S1~seq1" %in% eligible_for("whole_v", truncated)$sequence_id)
  # The sub-regions do not carry the guard: it is about the whole window's width.
  expect_true("S1~seq1" %in% eligible_for("cdr2", truncated)$sequence_id)
})

test_that("whole-V occupancy includes the exact configured boundary", {
  records <- JOINED_MEMBERSHIP[rep(1, 3), , drop = FALSE]
  # The other four regions occupy 164 positions. Complete FWR3 coverage at
  # 105, 106 and 107 positions gives whole-V totals of 269, 270 and 271.
  records$germline_occupied_fwr3 <- c(105L, 106L, 107L)
  records$informative_supported_fwr3 <- c(105L, 106L, 107L)
  for (region in REGION_LEVELS) {
    expect_true(all(shm_region_eligible(records, region, 1, 270, REGION_LEVELS)))
  }
  expect_identical(shm_region_eligible(records, "whole_v", 1, 270, REGION_LEVELS),
                   c(FALSE, TRUE, TRUE))
})

test_that("zero SHM coverage floors retain finite rates and occupied denominators", {
  records <- JOINED_MEMBERSHIP[rep(1, 4), , drop = FALSE]
  for (region in REGION_LEVELS) {
    records[[paste0("informative_supported_", region)]] <- 0L
    records[[paste0("germline_occupied_", region)]] <- 1L
  }
  records$germline_occupied_cdr1[2] <- 0L
  records$informative_supported_cdr1[3] <- NA_integer_
  records$shm_frequency[4] <- NA_real_
  expect_identical(shm_region_eligible(records, "whole_v", 0, 0, REGION_LEVELS),
                   c(TRUE, FALSE, FALSE, FALSE))
  expect_false(any(shm_region_eligible(records, "whole_v", 0, 270, REGION_LEVELS)))
})

test_that("a scope with no members yields nothing to summarise", {
  expect_identical(nrow(eligible_for("whole_v", data.frame())), 0L)
})

# ------------------------------------- the columns released per sequence ---
#
# The donor member shard and clone_member_shm.tsv carry the same mutation load per
# receptor that the rows above summarise. The point of these cases is that it is the
# same load under the same rule, so that a reader grouping the member rows reproduces
# the summary rather than approximating it.

describe("the mutation count of one V region", {
  it("sums the five sub-regions for the whole V gene", {
    counts <- region_mutation_counts(JOINED_MEMBERSHIP, "whole_v", REGION_LEVELS)
    # Every member carries five per sub-region except the unmutated seq6.
    expect_identical(counts, c(25, 25, 25, 25, 25, 0, 25))
  })
  it("reads a sub-region off its own column", {
    expect_identical(region_mutation_counts(JOINED_MEMBERSHIP, "cdr2",
                                            REGION_LEVELS),
                     c(5, 5, 5, 5, 5, 0, 5))
  })
  # A member that measured none of a sub-region has no count for it, and the whole
  # V gene it is part of therefore has none either -- not a smaller total.
  it("does not treat an uncounted sub-region as zero", {
    missing <- JOINED_MEMBERSHIP
    missing$mutations_cdr1[3] <- NA_integer_
    counts <- region_mutation_counts(missing, "whole_v", REGION_LEVELS)
    expect_true(is.na(counts[3]))
    expect_identical(counts[1], 25)
  })
})

describe("the donor shard's per-sequence mutation columns", {
  it("names the whole V gene, then each sub-region, then every verdict", {
    expect_identical(member_shm_columns(),
                     c("shm_frequency_whole_v", "mutations_whole_v",
                       "shm_frequency_cdr2", "shm_frequency_fwr3",
                       "mutations_cdr2", "mutations_fwr3",
                       paste0("shm_eligible_", c("whole_v", REGION_LEVELS)),
                       paste0("shm_evidence_sequence_id_", c("whole_v", REGION_LEVELS))))
  })
  it("reports the whole V gene and the two sub-regions the member shard keeps",
     expect_identical(CLONE_SHM_V_REGIONS, c("whole_v", "cdr2", "fwr3")))

  widened <- add_member_shm_columns(JOINED_MEMBERSHIP, MIN_SUBREGION_FRACTION,
                                    MIN_GERMLINE_OCCUPIED_V, REGION_LEVELS)

  it("leaves every shard column present",
     expect_true(all(member_shm_columns() %in% names(widened))))
  it("keeps the whole-V count as the sum of the sub-region counts", {
    expect_identical(widened$mutations_whole_v,
                     widened$mutations_cdr2 + widened$mutations_fwr3 +
                       as.numeric(JOINED_MEMBERSHIP$mutations_fwr1) +
                       as.numeric(JOINED_MEMBERSHIP$mutations_cdr1) +
                       as.numeric(JOINED_MEMBERSHIP$mutations_fwr2))
  })
  # The verdict is the summary table's own rule, not a second reading of it, so
  # the members it admits are exactly the ones the summary counted.
  it("gives each region the eligibility rule the summary rows use", {
    for (v_region in CLONE_SHM_V_REGIONS) {
      expect_identical(
        widened[[paste0("shm_eligible_", v_region)]],
        shm_region_eligible(JOINED_MEMBERSHIP, v_region, MIN_SUBREGION_FRACTION,
                            MIN_GERMLINE_OCCUPIED_V, REGION_LEVELS)
      )
      expect_identical(
        widened$sequence_id[widened[[paste0("shm_eligible_", v_region)]]],
        eligible_for(v_region)$sequence_id
      )
    }
  })
  # seq2 measured 10 of FWR1's 74 germline positions and all of everything else,
  # which is the case the whole-V conjunction exists to catch.
  it("can admit a member for a sub-region and refuse it the whole V gene", {
    seq2 <- widened$sequence_id == "S1~seq2"
    expect_false(widened$shm_eligible_whole_v[seq2])
    expect_true(widened$shm_eligible_fwr3[seq2])
  })
  # A rate is reported for the window a member did measure whether or not the
  # coverage rule will summarise it, so the two columns are read together.
  it("reports a rate beside a refused verdict", {
    seq2 <- widened$sequence_id == "S1~seq2"
    expect_equal(widened$shm_frequency[seq2], 0.30, tolerance = EXACT)
  })
  it("never leaves a verdict undecided",
     expect_false(any(vapply(paste0("shm_eligible_", CLONE_SHM_V_REGIONS),
                             function(column) anyNA(widened[[column]]),
                             logical(1)))))
})

# ------------------------------------------- the mean of lineage medians ---

ELIGIBLE <- eligible_for("whole_v")
MEASURED <- attr(ELIGIBLE, "measured")

# build_scope_shm_table() selects a cell's members before summarising them, so the
# helper does the same two steps: shm_group_mask() picks the isotype group and the
# library is picked by name. Both are asserted separately below, and the counts the
# summary reports are checked against them.
summarise <- function(library_id, group, v_region = "whole_v",
                      members = ELIGIBLE, measured = MEASURED) {
  in_group <- members[shm_group_mask(members, group), , drop = FALSE]
  cell <- in_group[which(in_group$library_id == library_id), , drop = FALSE]
  n_measured <- count_measured_members(measured, group, MANIFEST$library_id)
  summarise_lineage_median_shm(cell, library_id, group, v_region,
                               if (library_id %in% names(n_measured)) {
                                 n_measured[[library_id]]
                               } else {
                                 NA_integer_
                               })
}

describe("one library summarised over every eligible lineage", {
  row <- summarise("S1", "all_IGH")

  # Five of S1's members had a whole-V rate; four covered every sub-region.
  it("reports what it measured beside what it used", {
    expect_identical(row$n_members_measured, 5L)
    expect_identical(row$n_members_eligible, 4L)
  })
  it("counts the lineages they fell into", expect_identical(row$n_clones_eligible, 4L))
  # L1 = 0.10, L2 = 0.20, L3 = 0.60, L5 = 0.50.
  it("averages the per-lineage medians, not the sequences",
     expect_equal(row$mean_clone_median_shm, 0.35, tolerance = EXACT))
  it("reports the median of those medians", {
    expect_equal(row$median_clone_median_shm, 0.35, tolerance = EXACT)
  })
  it("reports the quartiles of those medians and their difference", {
    expect_equal(row$q25_clone_median_shm, 0.175, tolerance = EXACT)
    expect_equal(row$q75_clone_median_shm, 0.525, tolerance = EXACT)
    expect_equal(row$iqr_clone_median_shm, 0.35, tolerance = EXACT)
  })
  # The four summarised sequences have base votes 4, 10, 20 and 30. Each one's
  # whole-V value is the mean of its five sub-region medians, i.e. its base plus
  # 2, giving 6, 12, 22 and 32; the row reports the median of those, 17. Their
  # median-of-five would give 15 and their occupancy-weighted mean about 21.5, so
  # this value identifies which of the three rules is in force.
  it("averages a sequence's five sub-region medians, then takes the median
      over sequences",
     expect_equal(row$median_trust4_consensus_base_votes, 17, tolerance = EXACT))
  it("names the library, the group and the region it summarised", {
    expect_identical(row$library_id, "S1")
    expect_identical(row$isotype_group, "all_IGH")
    expect_identical(row$v_region, "whole_v")
  })
})

describe("the same library restricted to one isotype class", {
  igg <- summarise("S1", "IgG")
  igm <- summarise("S1", "IgM")

  # seq1 and seq4 are eligible IgG; seq2 is IgG but did not cover FWR1.
  it("counts the measured members it could not use", {
    expect_identical(igg$n_members_measured, 3L)
    expect_identical(igg$n_members_eligible, 2L)
  })
  it("averages over the lineages those sequences fell into",
     expect_equal(igg$mean_clone_median_shm, 0.35, tolerance = EXACT))
  # L1 answers for IgG here and L2 for IgM: one lineage can be summarised once in
  # every isotype its members express.
  it("lets a different lineage answer for a different isotype", {
    expect_equal(igm$mean_clone_median_shm, 0.20, tolerance = EXACT)
    expect_identical(igm$n_members_eligible, 1L)
  })
  it("leaves the median and quartiles of a single lineage equal to it", {
    expect_equal(igm$median_clone_median_shm, 0.20, tolerance = EXACT)
    expect_equal(igm$q25_clone_median_shm, 0.20, tolerance = EXACT)
    expect_equal(igm$q75_clone_median_shm, 0.20, tolerance = EXACT)
    expect_identical(igm$iqr_clone_median_shm, 0)
  })
})

describe("the subclass rows beneath a class row", {
  it("summarises the subclass its call named",
     expect_equal(summarise("S1", "IgG1")$mean_clone_median_shm, 0.10,
                  tolerance = EXACT))
  # seq4 is IgG but names no single subclass, so IgG holds two sequences and the
  # subclass rows beneath it hold one between them. They need not sum to it.
  it("does not require the subclasses to sum to their class", {
    expect_identical(summarise("S1", "IgG")$n_members_eligible, 2L)
    expect_identical(summarise("S1", "IgG1")$n_members_eligible, 1L)
    expect_identical(summarise("S1", "IgG3")$n_members_eligible, 0L)
    expect_identical(summarise("S1", "IgG4")$n_members_eligible, 0L)
  })
  # seq2 is IgG3 and is eligible in FWR3 but not whole_v, so a subclass row follows
  # the region's own eligibility exactly as a class row does.
  it("follows the region's eligibility", {
    fwr3 <- eligible_for("fwr3")
    expect_identical(summarise("S1", "IgG3", "fwr3", fwr3,
                               attr(fwr3, "measured"))$n_members_eligible, 1L)
  })
  it("summarises an IgA subclass",
     expect_equal(summarise("S1", "IgA2")$mean_clone_median_shm, 0.50,
                  tolerance = EXACT))
})

# A sub-region row reports the same statistic read off one column, so the two row
# types are asserted against each other rather than only the whole-V one being
# checked.
describe("vote support on a sub-region row", {
  it("takes the median over sequences of that sub-region's own median", {
    fwr3 <- eligible_for("fwr3")
    # S1's five FWR3-eligible sequences carry base votes 4, 4, 10, 20 and 30,
    # each raised by the FWR3 offset of 20, so their median is 30.
    expect_equal(summarise("S1", "all_IGH", "fwr3", fwr3,
                           attr(fwr3, "measured"))$median_trust4_consensus_base_votes,
                 30, tolerance = EXACT)
  })

  it("differs between sub-regions, so no sub-region stands in for another", {
    fwr1 <- eligible_for("fwr1")
    cdr2 <- eligible_for("cdr2")
    # FWR1 excludes seq2, which reached only 10 of its 74 positions, so its four
    # sequences carry bases 4, 10, 20 and 30 less the FWR1 offset of 10: median 5.
    expect_equal(summarise("S1", "all_IGH", "fwr1", fwr1,
                           attr(fwr1, "measured"))$median_trust4_consensus_base_votes,
                 5, tolerance = EXACT)
    # CDR2 keeps all five, so bases 4, 4, 10, 20 and 30 plus its offset of 5.
    expect_equal(summarise("S1", "all_IGH", "cdr2", cdr2,
                           attr(cdr2, "measured"))$median_trust4_consensus_base_votes,
                 15, tolerance = EXACT)
  })
})

describe("the zero-mutation counts", {
  # S3's only member carries no mutation anywhere.
  it("counts a lineage whose median member has no mutation",
     expect_identical(summarise("S3", "all_IGH")$n_clones_zero_median_shm, 1L))
  it("counts none where every lineage is mutated",
     expect_identical(summarise("S1", "all_IGH")$n_clones_zero_median_shm, 0L))
  # The sequence-level count answers a different question from the lineage one and
  # takes the other denominator in the row, so the two are reported separately.
  it("counts the unmutated sequences against their own denominator", {
    row <- summarise("S3", "all_IGH")
    expect_identical(row$n_members_zero_shm, row$n_members_eligible)
    expect_identical(summarise("S1", "all_IGH")$n_members_zero_shm, 0L)
  })
  # It is decided on the integer mutation count, so a rate that rounds to zero on a
  # long window is not silently counted as germline.
  it("does not read a very small rate as unmutated", {
    tiny <- JOINED_MEMBERSHIP
    tiny$shm_frequency[6] <- 1e-18
    for (region in REGION_LEVELS) tiny[[paste0("mutations_", region)]][6] <- 1L
    members <- eligible_for("whole_v", tiny)
    expect_identical(
      summarise("S3", "all_IGH", "whole_v", members,
                attr(members, "measured"))$n_clones_zero_median_shm,
      0L)
  })
})

describe("a group with nothing in it", {
  for (case in list(
    list(what = "an isotype the library never expressed",
         library_id = "S1", group = "IgA1"),
    list(what = "a library whose only member covered too little",
         library_id = "S2", group = "all_IGH"),
    list(what = "a library outside the membership entirely",
         library_id = "S9", group = "all_IGH"))) {
    row <- summarise(case$library_id, case$group)

    it(sprintf("counts nothing [%s]", case$what), {
      expect_identical(row$n_members_eligible, 0L)
      expect_identical(row$n_members_zero_shm, 0L)
      expect_identical(row$n_clones_eligible, 0L)
      expect_identical(row$n_clones_zero_median_shm, 0L)
    })
    it(sprintf("leaves every statistic undefined [%s]", case$what), {
      expect_identical(row$mean_clone_median_shm, NA_real_)
      expect_identical(row$median_clone_median_shm, NA_real_)
      expect_identical(row$q25_clone_median_shm, NA_real_)
      expect_identical(row$q75_clone_median_shm, NA_real_)
      expect_identical(row$iqr_clone_median_shm, NA_real_)
      expect_identical(row$median_trust4_consensus_base_votes, NA_real_)
    })
  }
})

# S2 measured one member and could use none of it. That is the distinction the two
# counts exist to draw: nothing was summarised, but something was there.
test_that("a library that measured something but used none of it says so", {
  row <- summarise("S2", "all_IGH")
  expect_identical(row$n_members_measured, 1L)
  expect_identical(row$n_members_eligible, 0L)
})

test_that("missing rates do not enter clone medians or their equal weighting", {
  # Every member satisfies coverage. Clone A has rates 0.1, missing and 0.3;
  # clone B has 0.6; clone C has no rate. The eligible clone medians are 0.2
  # and 0.6, whose equally weighted mean is 0.4, regardless of read support.
  records <- JOINED_MEMBERSHIP[rep(1, 5), ]
  records$sequence_id <- paste0("S1~missing_case", seq_len(5))
  records$clone_id <- c("A", "A", "A", "B", "C")
  records$consensus_count <- c(1000, 1000, 1, 1, 1000)
  rates <- c(0.1, NA_real_, 0.3, 0.6, NA_real_)
  for (column in c("shm_frequency", paste0("shm_frequency_", REGION_LEVELS))) {
    records[[column]] <- rates
  }
  for (region in c("whole_v", "cdr1")) {
    eligible <- eligible_for(region, records)
    measured <- count_measured_members(attr(eligible, "measured"), "all_IGH", "S1")
    row <- summarise_lineage_median_shm(eligible, "S1", "all_IGH", region, measured[["S1"]])
    expect_identical(row$n_members_measured, 3L)
    expect_identical(row$n_members_eligible, 3L)
    expect_identical(row$n_clones_eligible, 2L)
    expect_equal(row$mean_clone_median_shm, 0.4, tolerance = EXACT)
    expect_equal(row$median_clone_median_shm, 0.4, tolerance = EXACT)
    # R type-7 quartiles interpolate between the two clone medians.
    expect_equal(row$q25_clone_median_shm, 0.3, tolerance = EXACT)
    expect_equal(row$q75_clone_median_shm, 0.5, tolerance = EXACT)
    expect_equal(row$iqr_clone_median_shm, 0.2, tolerance = EXACT)
  }
  for (column in c("shm_frequency", paste0("shm_frequency_", REGION_LEVELS))) {
    records[[column]] <- NA_real_
  }
  for (region in c("whole_v", "cdr1")) {
    eligible <- eligible_for(region, records)
    measured <- count_measured_members(attr(eligible, "measured"), "all_IGH", "S1")
    row <- summarise_lineage_median_shm(eligible, "S1", "all_IGH", region, measured[["S1"]])
    expect_identical(row$n_members_measured, 0L)
    expect_identical(row$n_members_eligible, 0L)
    expect_identical(row$n_clones_eligible, 0L)
    expect_identical(row$mean_clone_median_shm, NA_real_)
    expect_identical(row$q25_clone_median_shm, NA_real_)
    expect_identical(row$q75_clone_median_shm, NA_real_)
  }
})

test_that("a scope with no eligible member still reports a row", {
  row <- summarise_lineage_median_shm(data.frame(), "S1", "all_IGH", "whole_v")
  expect_identical(row$library_id, "S1")
  expect_identical(row$n_members_eligible, 0L)
  expect_identical(row$mean_clone_median_shm, NA_real_)
})

# ------------------------------------------------------- the whole table ---

build_table <- function(membership, manifest = MANIFEST) {
  build_scope_shm_table(membership, manifest, "scope_x", V_REGIONS,
                        MIN_SUBREGION_FRACTION, MIN_GERMLINE_OCCUPIED_V,
                        REGION_LEVELS, SHM_ISOTYPE_GROUPS)
}

ROWS_PER_LIBRARY <- length(SHM_ISOTYPE_GROUPS) * length(V_REGIONS)

describe("the released table for a scope", {
  shm <- build_table(JOINED_MEMBERSHIP)

  it("declares the released columns in order",
     expect_identical(names(shm), SHM_COLUMNS))
  it("gives every library one row per isotype group per region",
     expect_identical(nrow(shm), nrow(MANIFEST) * ROWS_PER_LIBRARY))
  it("stamps the scope on every row",
     expect_true(all(shm$scope_id == "scope_x")))
  it("carries the manifest's donor and tissue over", {
    s3 <- shm[shm$library_id == "S3", , drop = FALSE]
    expect_true(all(s3$donor_id == "D2"))
    expect_true(all(s3$tissue_id == "T3"))
  })
  it("orders by library, then group, then region", {
    expect_identical(shm$library_id, rep(c("S1", "S2", "S3"), each = ROWS_PER_LIBRARY))
    expect_identical(shm$isotype_group[seq_len(ROWS_PER_LIBRARY)],
                     rep(SHM_ISOTYPE_GROUPS, each = 6L))
    expect_identical(shm$v_region[1:6], V_REGIONS)
  })
  # A library that contributed nothing is still reported, so a reader can tell an
  # empty repertoire from a library that was never in the scope.
  it("reports a library with no eligible sequence as zeros", {
    empty <- shm[shm$library_id == "S2" & shm$v_region == "whole_v", , drop = FALSE]
    expect_true(all(empty$n_members_eligible == 0L))
    expect_true(all(is.na(empty$mean_clone_median_shm)))
  })
  it("holds the same numbers the per-row summary produced", {
    whole_v <- shm[shm$library_id == "S1" & shm$v_region == "whole_v" &
                     shm$isotype_group == "all_IGH", , drop = FALSE]
    expect_equal(whole_v$mean_clone_median_shm, 0.35, tolerance = EXACT)
  })
  # The counts differ per region because coverage does, which is the whole point of
  # reporting the regions separately, and whole_v is the strictest because it is all
  # five at once.
  it("counts each region over its own eligible members", {
    all_igh <- shm[shm$isotype_group == "all_IGH" & shm$library_id == "S1", ,
                   drop = FALSE]
    counts <- stats::setNames(all_igh$n_members_eligible, all_igh$v_region)
    expect_identical(unname(counts[["fwr1"]]), 4L)
    expect_identical(unname(counts[["fwr3"]]), 5L)
    expect_identical(unname(counts[["whole_v"]]), 4L)
  })
  it("never summarises more than it measured",
     expect_true(all(shm$n_members_eligible <= shm$n_members_measured)))
})

test_that("a scope with no eligible member reports every library as zero", {
  shm <- build_table(data.frame())
  expect_identical(nrow(shm), nrow(MANIFEST) * ROWS_PER_LIBRARY)
  expect_true(all(shm$n_members_eligible == 0L))
  expect_true(all(is.na(shm$mean_clone_median_shm)))
})

# DNA and its measured mutation load are identical. Only the source assembly's
# vote support differs, so ranking must use coverage, then original support.
evidence_fixture <- function() {
  candidates <- JOINED_MEMBERSHIP[rep(1, 3), ]
  candidates$sequence_id <- paste0("S1~assemble", c(3, 2, 1), "_0")
  candidates$consensus_count <- c(100, 20, 20)
  for (region in REGION_LEVELS) {
    candidates[[paste0("median_majority_votes_", region)]] <- c(30, 8, 9)
  }
  members <- candidates[1, ]
  members$consensus_count <- 999  # Transferred support must never rank evidence.
  list(members = members, candidates = candidates)
}

choose_fixture_evidence <- function(fixture, fraction = 1, floor = 270) {
  choose_shm_evidence(fixture$members, fixture$candidates,
    rep(fixture$members$sequence_id, nrow(fixture$candidates)),
    fraction, floor, REGION_LEVELS)
}

test_that("regional evidence first passes its gate, then ranks coverage and source support", {
  fixture <- evidence_fixture()
  fixture$candidates$informative_supported_fwr1 <- c(10, 74, 74)
  fixture$candidates$informative_supported_fwr3 <- c(114, 100, 90)
  evidence <- choose_fixture_evidence(fixture, .8)
  source_for <- function(region) evidence$evidence_sequence_id[evidence$v_region == region]
  expect_identical(source_for("whole_v"), "S1~assemble2_0")
  expect_identical(source_for("fwr1"), "S1~assemble1_0") # stable ID breaks original support tie
  expect_identical(source_for("fwr3"), "S1~assemble3_0")
  expect_identical(source_for("cdr2"), "S1~assemble3_0") # own support breaks coverage tie
  for (indices in list(3:1, c(2, 3, 1))) {
    shuffled <- fixture
    shuffled$candidates <- shuffled$candidates[indices, ]
    expect_identical(choose_fixture_evidence(shuffled, .8), evidence)
  }
  released <- add_member_shm_columns(fixture$members, .8, 270, REGION_LEVELS, evidence)
  expect_identical(released$shm_evidence_sequence_id_fwr1, "S1~assemble1_0")
  expect_identical(released$shm_evidence_sequence_id_whole_v, "S1~assemble2_0")
  selected <- select_shm_eligible_members(fixture$members, "whole_v", .8, 270,
    REGION_LEVELS, evidence = evidence)
  expect_equal(selected$median_trust4_consensus_base_votes, 8)
})

test_that("complementary sources cannot jointly satisfy whole-V coverage", {
  fixture <- evidence_fixture()
  fixture$candidates <- fixture$candidates[1:2, ]
  fixture$candidates$informative_supported_fwr1 <- c(0, 74)
  fixture$candidates$informative_supported_fwr3 <- c(114, 0)
  evidence <- choose_fixture_evidence(fixture)
  expect_false(evidence$eligible[evidence$v_region == "whole_v"])
  expect_true(is.na(evidence$evidence_sequence_id[evidence$v_region == "whole_v"]))
  expect_true(all(evidence$eligible[evidence$v_region != "whole_v"]))
  member <- add_member_shm_columns(fixture$members, 1, 270, REGION_LEVELS, evidence)
  expect_identical(member$shm_frequency_whole_v, fixture$members$shm_frequency)
  expect_false(member$shm_eligible_whole_v)
  expect_true(is.na(member$shm_evidence_sequence_id_whole_v))
  expect_equal(nrow(select_shm_eligible_members(fixture$members, "whole_v", 1, 270,
    REGION_LEVELS, evidence = evidence)), 0L)
})

test_that("SHM evidence cannot cross libraries or change sequence-derived measurements", {
  fixture <- evidence_fixture()
  fixture$candidates$library_id[2] <- "S2"
  expect_error(choose_fixture_evidence(fixture), "crosses libraries")
  fixture <- evidence_fixture()
  fixture$candidates$mutations_fwr3[2] <- 99
  expect_error(choose_fixture_evidence(fixture), "disagrees in mutations_fwr3")
})

# The supported count travels with the source decision, so a coverage verdict can be
# read from the row rather than rebuilt from the unreleased source records.
member_shm_fixture <- function(fixture) {
  fixture$members$scope_id <- "scope_x"
  fixture$members$tissue_id <- "T1"
  fixture
}

test_that("a member's supported count comes from the source the ranking selects", {
  fixture <- member_shm_fixture(evidence_fixture())
  fixture$candidates <- fixture$candidates[1:2, ]
  # At a fraction of 0.8 the first source fails FWR1 and supports more of FWR3.
  fixture$candidates$informative_supported_fwr1 <- c(10, 74)
  fixture$candidates$informative_supported_fwr3 <- c(114, 100)
  evidence <- choose_fixture_evidence(fixture, .8)
  rows <- build_clone_member_shm_table(fixture$members, .8, 270, REGION_LEVELS, evidence = evidence)
  supported <- stats::setNames(rows$n_vote_supported_positions, rows$v_region)
  chosen <- stats::setNames(evidence$evidence_sequence_id, evidence$v_region)
  # Whole V takes the second source's 74 + 30 + 30 + 30 + 100, although the first
  # source's FWR3 is better, and FWR3 alone takes the first source's 114.
  expect_identical(chosen[["whole_v"]], "S1~assemble2_0")
  expect_equal(supported[["whole_v"]], 264)
  expect_identical(chosen[["fwr1"]], "S1~assemble2_0")
  expect_equal(supported[["fwr1"]], 74)
  expect_identical(chosen[["fwr3"]], "S1~assemble3_0")
  expect_equal(supported[["fwr3"]], 114)
  expect_true(all(rows$shm_eligible))
})

test_that("an ineligible member reports the supported count of the source ranked first", {
  fixture <- member_shm_fixture(evidence_fixture())
  fixture$candidates <- fixture$candidates[1:2, ]
  fixture$candidates$informative_supported_fwr1 <- c(0, 74)
  fixture$candidates$informative_supported_fwr3 <- c(114, 0)
  evidence <- choose_fixture_evidence(fixture)
  rows <- build_clone_member_shm_table(fixture$members, 1, 270, REGION_LEVELS, evidence = evidence)
  whole <- rows[rows$v_region == "whole_v", ]
  # Neither source passes, so the larger total ranks first: 0 + 90 + 114 against
  # 74 + 90 + 0. The votes stay empty, and the measured rate stays in this shard.
  expect_false(whole$shm_eligible)
  expect_equal(whole$n_vote_supported_positions, 204)
  expect_true(is.na(whole$median_trust4_consensus_base_votes))
  expect_equal(whole$shm_frequency, fixture$members$shm_frequency)
  expect_equal(whole$n_mutations, 25)
  # A source with no count for a subregion has no finite total to report.
  fixture$candidates$informative_supported_cdr1 <- NA_real_
  evidence <- choose_fixture_evidence(fixture)
  rows <- build_clone_member_shm_table(fixture$members, 1, 270, REGION_LEVELS, evidence = evidence)
  expect_true(is.na(rows$n_vote_supported_positions[rows$v_region == "whole_v"]))
  expect_true(is.na(rows$n_vote_supported_positions[rows$v_region == "cdr1"]))
  expect_equal(rows$n_vote_supported_positions[rows$v_region == "cdr2"], 30)
})

describe("the member SHM rows of a scope", {
  records <- JOINED_MEMBERSHIP
  records$scope_id <- "scope_x"
  records$tissue_id <- MANIFEST$tissue_id[match(records$library_id, MANIFEST$library_id)]
  # seq3 has no CDR1 rate and seq7 no whole-V rate, so each loses that region's row.
  records$shm_frequency_cdr1[3] <- NA_real_
  records$shm_frequency[7] <- NA_real_
  rows <- build_clone_member_shm_table(records, MIN_SUBREGION_FRACTION, MIN_GERMLINE_OCCUPIED_V,
                                       REGION_LEVELS)

  it("declares the released columns in order",
     expect_identical(names(rows), CLONE_MEMBER_SHM_COLUMNS))
  it("holds one row per member and region with a finite rate", {
    expect_identical(nrow(rows), 7L * 6L - 2L)
    expect_false(any(rows$sequence_id == "S1~seq3" & rows$v_region == "cdr1"))
    expect_false(any(rows$sequence_id == "S1~seq7" & rows$v_region == "whole_v"))
    expect_false(anyDuplicated(rows[c("sequence_id", "v_region")]) > 0)
  })
  it("counts n_members_measured of the library table", {
    shm <- build_scope_shm_table(records, MANIFEST, "scope_x", V_REGIONS, MIN_SUBREGION_FRACTION,
                                 MIN_GERMLINE_OCCUPIED_V, REGION_LEVELS)
    all_igh <- shm[shm$isotype_group == "all_IGH", ]
    for (index in seq_len(nrow(all_igh))) {
      expect_identical(sum(rows$library_id == all_igh$library_id[index] &
                             rows$v_region == all_igh$v_region[index]),
                       all_igh$n_members_measured[index])
    }
  })
  # seq5 measured 100 positions of FWR3 alone, so its whole-V denominator is 100 while
  # its germline occupies all 278.
  it("reads the whole-V denominator and sums germline occupancy", {
    seq5 <- rows[rows$sequence_id == "S2~seq5" & rows$v_region == "whole_v", ]
    expect_equal(seq5$n_informative_positions, 100)
    expect_equal(seq5$n_germline_occupied_positions, 278)
    expect_equal(seq5$n_vote_supported_positions, 100)
    expect_false(seq5$shm_eligible)
  })
  it("reads a subregion's counts off its own columns", {
    seq2 <- rows[rows$sequence_id == "S1~seq2" & rows$v_region == "fwr1", ]
    expect_equal(seq2$n_informative_positions, 10)
    expect_equal(seq2$n_germline_occupied_positions, 74)
    expect_equal(seq2$n_vote_supported_positions, 10)
    expect_equal(seq2$n_mutations, 5)
  })
  it("uses the eligibility the summary rows use", {
    for (region in V_REGIONS) {
      expect_setequal(rows$sequence_id[rows$v_region == region & rows$shm_eligible],
                      eligible_for(region, records)$sequence_id)
    }
  })
  it("leaves a header alone when there are no members or no measured rate", {
    empty <- build_clone_member_shm_table(records[0, ], MIN_SUBREGION_FRACTION,
                                          MIN_GERMLINE_OCCUPIED_V, REGION_LEVELS)
    expect_identical(names(empty), CLONE_MEMBER_SHM_COLUMNS)
    expect_identical(nrow(empty), 0L)
    unmeasured <- records
    for (column in c("shm_frequency", paste0("shm_frequency_", REGION_LEVELS))) {
      unmeasured[[column]] <- NA_real_
    }
    none <- build_clone_member_shm_table(unmeasured, MIN_SUBREGION_FRACTION,
                                         MIN_GERMLINE_OCCUPIED_V, REGION_LEVELS)
    expect_identical(names(none), CLONE_MEMBER_SHM_COLUMNS)
    expect_identical(nrow(none), 0L)
  })
})

# The complete audit uses the same evidence owner as the sparse measurements.
test_that("the coverage audit retains unmeasured and ineligible member regions", {
  regions <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
  records <- data.frame(scope_id = "scope", library_id = "S1", tissue_id = "T1",
    donor_id = "D1", clone_id = "c1", sequence_id = c("S1~a", "S1~b", "S1~c"),
    germline_alignment_d_mask = c(NA, "ACGT", "ACGT"),
    shm_frequency = c(NA, 0, 0.25), stringsAsFactors = FALSE)
  for (region in regions) records[[paste0("shm_frequency_", region)]] <- c(NA, 0, 0.25)
  evidence <- do.call(rbind, lapply(c("whole_v", regions), function(region) {
    data.frame(sequence_id = records$sequence_id, v_region = region,
      eligible = c(FALSE, FALSE, TRUE), evidence_sequence_id = c(NA, NA, "S1~c"))
  }))
  audit <- build_member_shm_coverage(records, evidence, regions)
  expect_equal(nrow(audit), 18L)
  expect_equal(table(audit$measurement_status), table(rep(c("unmeasured", "measured_ineligible", "eligible"), 6)))
  expect_true(all(is.na(audit$evidence_sequence_id[audit$measurement_status != "eligible"])))
  expect_identical(unique(audit$germline_reconstruction_status[audit$sequence_id == "S1~a"]), "unavailable_other")
})
