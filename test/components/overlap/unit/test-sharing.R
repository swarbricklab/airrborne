# Check pairwise sharing, endpoint denominators and directional comparisons.
# Run: python test/run_tests.py --component overlap --level unit

load_test_modules(c("table_io", "diversity", "gene_calls", "observed_representatives", "overlap_metrics"))

describe("make_empty_overlap_tables against the declared schemas", {
  tables <- make_empty_overlap_tables()

  it("gives the pair table exactly the declared pair columns, in order",
     expect_identical(names(tables$donor_overlap_pairs), OVERLAP_SCHEMAS$donor_overlap_pairs))
  it("gives the donor table exactly the declared donor columns, in order",
     expect_identical(names(tables$donor_overlap_summary), OVERLAP_SCHEMAS$donor_overlap_summary))
  it("gives the shared-clone table exactly the declared columns, in order",
     expect_identical(names(tables$shared_groups),
                      OVERLAP_SCHEMAS$shared_groups))
  it("returns the three tables and nothing else",
     expect_identical(names(tables),
                      c("donor_overlap_pairs", "donor_overlap_summary", "shared_groups")))
  it("returns the pair table with no rows",
     expect_identical(nrow(tables$donor_overlap_pairs), 0L))
  it("returns the donor table with no rows",
     expect_identical(nrow(tables$donor_overlap_summary), 0L))
  it("returns the shared-clone table with no rows",
     expect_identical(nrow(tables$shared_groups), 0L))
  it("types the clone counts as integers",
     expect_true(all(vapply(tables$donor_overlap_pairs[c("n_groups_1", "n_groups_2", "n_groups_shared")],
                            is.integer, logical(1)))))
  it("types the similarity indices as numeric",
     expect_true(all(vapply(tables$donor_overlap_pairs[c("jaccard", "morisita_horn")],
                            is.numeric, logical(1)))))
  it("types the identifiers as character, not factor",
     expect_true(all(vapply(tables$donor_overlap_pairs[c("donor_id", "locus", "library_id_1")],
                            is.character, logical(1)))))
  # The pair table leaves out every value that follows from its own columns or
  # from diversity.tsv.
  it("releases no derivable pair measurement", {
    derivable <- c("overlap_status", "simpson_concentration_difference",
                   "simpson_difference_status", "shared_higher_in_2_fraction",
                   "shared_higher_in_1_fraction", "shared_equal_fraction",
                   "shared_difference_balance", "shared_difference_status",
                   "support_fraction_only_in_2", "support_fraction_only_in_1")
    expect_length(intersect(derivable, OVERLAP_SCHEMAS$donor_overlap_pairs), 0L)
  })

})

# Every expected value below is an exact identity -- 1/3, 0.5, 1/6 -- so the
# comparisons should hold to close to machine precision rather than to
# testthat's default of about 1.5e-8.
EXACT <- 1e-12

describe("jaccard_index", {
  it("is shared over union on two presence masks", {
    expect_equal(jaccard_index(c(TRUE, TRUE, FALSE, FALSE),
                               c(TRUE, FALSE, TRUE, FALSE)),
                 1 / 3, tolerance = EXACT)
  })
  it("is one for two libraries detecting the same clones", {
    expect_equal(jaccard_index(c(TRUE, TRUE), c(TRUE, TRUE)), 1,
                 tolerance = EXACT)
  })
  it("is zero for two libraries sharing nothing", {
    expect_equal(jaccard_index(c(TRUE, FALSE), c(FALSE, TRUE)), 0,
                 tolerance = EXACT)
  })
  it("counts a single shared clone", {
    expect_equal(jaccard_index(TRUE, TRUE), 1, tolerance = EXACT)
  })
  it("is undefined when neither library detected anything", {
    expect_identical(jaccard_index(c(FALSE, FALSE), c(FALSE, FALSE)), NA_real_)
  })
  it("is undefined when there are no keys at all", {
    expect_identical(jaccard_index(logical(0), logical(0)), NA_real_)
  })
  it("is symmetric in its arguments", {
    expect_identical(jaccard_index(c(TRUE, TRUE, FALSE), c(TRUE, FALSE, TRUE)),
                     jaccard_index(c(TRUE, FALSE, TRUE), c(TRUE, TRUE, FALSE)))
  })
})

describe("morisita_horn_index", {
  it("is one for two libraries with identical proportions", {
    p <- c(0.5, 0.3, 0.2)
    expect_equal(morisita_horn_index(p, p), 1, tolerance = EXACT)
  })
  it("is zero for two libraries sharing no clone", {
    expect_equal(morisita_horn_index(c(1, 0), c(0, 1)), 0, tolerance = EXACT)
  })
  it("is one for two single-clone libraries holding the same clone", {
    expect_equal(morisita_horn_index(1, 1), 1, tolerance = EXACT)
  })
  # 2*(0.5*0.25 + 0.5*0.75) / (0.5 + 0.625) = 1 / 1.125.
  it("weighs each clone by its share of read support", {
    expect_equal(morisita_horn_index(c(0.5, 0.5), c(0.25, 0.75)), 1 / 1.125,
                 tolerance = EXACT)
  })
  it("is undefined when both sides are all zero", {
    expect_identical(morisita_horn_index(c(0, 0), c(0, 0)), NA_real_)
  })
  it("is undefined when there are no keys at all", {
    expect_identical(morisita_horn_index(numeric(0), numeric(0)), NA_real_)
  })
})

describe("build_sharing_units", {
  it("sums read support per key and sorts by key", {
    expect_equal(
      build_sharing_units(c("b", "a", "b"), c(1, 2, 3)),
      data.frame(group_key = c("a", "b"), consensus_count = c(2, 4),
                 stringsAsFactors = FALSE),
      tolerance = EXACT
    )
  })
  it("keeps a library with one clone", {
    expect_equal(
      build_sharing_units("a", 7),
      data.frame(group_key = "a", consensus_count = 7, stringsAsFactors = FALSE),
      tolerance = EXACT
    )
  })
  it("drops a record with no key", {
    expect_identical(build_sharing_units(c(NA, "", "a"), c(9, 9, 1))$group_key, "a")
  })
  it("drops a record with no positive support", {
    expect_identical(build_sharing_units(c("a", "b", "c"), c(0, -1, 2))$group_key, "c")
  })
  it("returns no rows when nothing survives", {
    expect_identical(nrow(build_sharing_units(c("a", "b"), c(0, 0))), 0L)
  })
  it("returns no rows for a library with no records", {
    expect_identical(nrow(build_sharing_units(character(0), numeric(0))), 0L)
  })
  it("always returns the two sharing-unit columns", {
    expect_identical(names(build_sharing_units(character(0), numeric(0))),
                     c("group_key", "consensus_count"))
  })
  it("refuses keys and abundances of different lengths", {
    expect_error(build_sharing_units(c("a", "b", "c"), c(1, 2)))
  })
})

# The oracle's toy repertoires, rebuilt as sharing-unit frames.
oracle_reference <- function() {
  utils::read.delim(file.path(MODULE_ROOT, "test", "data", "expected", "overlap_reference.tsv"),
                    comment.char = "#", stringsAsFactors = FALSE,
                    check.names = FALSE)
}
oracle_units <- function(keys, counts) {
  build_sharing_units(strsplit(keys, ",", fixed = TRUE)[[1]],
                      as.numeric(strsplit(counts, ",", fixed = TRUE)[[1]]))
}

# Recorded immunarch results independently check the published sharing indices.
describe("the two indices against immunarch::repOverlap", {
  reference <- oracle_reference()

  it("reads five recorded pairs", expect_identical(nrow(reference), 5L))

  for (index in seq_len(nrow(reference))) {
    row <- reference[index, ]
    first <- oracle_units(row$keys_first, row$counts_first)
    second <- oracle_units(row$keys_second, row$counts_second)
    values <- calculate_repertoire_pair_overlap(first, second)

    it(paste("reproduces repOverlap jaccard for", row$pair), {
      expect_equal(values$jaccard, row$jaccard, tolerance = EXACT)
    })
    it(paste("reproduces repOverlap morisita for", row$pair), {
      expect_equal(values$morisita_horn, row$morisita, tolerance = EXACT)
    })
  }
})

# A partially overlapping pair, small enough to check by hand: the union is
# K1..K4, so two of four keys are shared.
describe("calculate_repertoire_pair_overlap", {
  first <- build_sharing_units(c("K1", "K2", "K3"), c(10, 5, 1))
  second <- build_sharing_units(c("K2", "K3", "K4"), c(20, 4, 6))
  values <- calculate_repertoire_pair_overlap(first, second)

  it("reports similarity and endpoint support on the same aligned union",
     expect_identical(names(values), c("pair_status", "min_consensus_count_per_library",
       "n_groups_shared", "shared_consensus_count_1", "shared_consensus_count_2",
       "shared_support_fraction_1", "shared_support_fraction_2",
       "jaccard", "morisita_horn", "bray_curtis_proportion",
       "n_shared_higher_in_1", "n_shared_higher_in_2", "n_shared_equal")))
  it("counts the shared keys", expect_identical(values$n_groups_shared, 2L))
  it("takes jaccard over the union of both key sets",
     expect_equal(values$jaccard, 0.5, tolerance = EXACT))
  it("is symmetric in its arguments", {
    reversed <- calculate_repertoire_pair_overlap(second, first)
    symmetric <- c("n_groups_shared", "jaccard", "morisita_horn")
    expect_equal(unlist(reversed[symmetric]), unlist(values[symmetric]), tolerance = EXACT)
    expect_equal(reversed$shared_support_fraction_1, values$shared_support_fraction_2)
    expect_equal(reversed$shared_support_fraction_2, values$shared_support_fraction_1)
  })
  it("scores a pair sharing nothing as zero on both indices", {
    apart <- calculate_repertoire_pair_overlap(
      first, build_sharing_units(c("K8", "K9"), c(1, 1)))
    expect_equal(c(apart$jaccard, apart$morisita_horn), c(0, 0),
                 tolerance = EXACT)
  })
  it("scores a library against itself as one on both indices", {
    same <- calculate_repertoire_pair_overlap(first, first)
    expect_equal(c(same$jaccard, same$morisita_horn), c(1, 1),
                 tolerance = EXACT)
  })
})

describe("build_pairwise_overlap_rows", {
  tissue_by_sample <- c(S1 = "T0", S2 = "T1", S3 = "T2")
  repertoires <- list(
    S1 = build_sharing_units(c("K1", "K2", "K3"), c(10, 5, 1)),
    S2 = build_sharing_units(c("K2", "K3", "K4"), c(20, 4, 6)),
    S3 = build_sharing_units(c("K5", "K6"), c(1, 1))
  )
  rows <- build_pairwise_overlap_rows(repertoires, "D1", "TRB", "clonotype", tissue_by_sample)

  it("writes the released columns in the declared order",
     expect_identical(names(rows), OVERLAP_SCHEMAS$donor_overlap_pairs))
  it("enumerates every unordered pair once",
     expect_identical(nrow(rows), 3L))
  it("orders the pairs as combn does", {
    expect_identical(paste(rows$library_id_1, rows$library_id_2),
                     c("S1 S2", "S1 S3", "S2 S3"))
  })
  it("carries each library's tissue across",
     expect_identical(rows$tissue_id_1, c("T0", "T0", "T1")))
  it("counts each side's clones",
     expect_identical(rows$n_groups_1, c(3L, 3L, 3L)))
  it("names the unit the counts hold",
     expect_identical(rows$group_unit, rep("clonotype", 3)))
  it("sums each side's read support",
     expect_equal(rows$grouped_consensus_count_1, c(16, 16, 30), tolerance = EXACT))
  it("returns the released empty table for one contributing library", {
    expect_identical(
      build_pairwise_overlap_rows(repertoires["S1"], "D1", "TRB", "clonotype", tissue_by_sample),
      make_empty_overlap_tables()$donor_overlap_pairs
    )
  })
  it("returns the released empty table when nothing contributed", {
    expect_identical(
      build_pairwise_overlap_rows(list(), "D1", "TRB", "clonotype", tissue_by_sample),
      make_empty_overlap_tables()$donor_overlap_pairs
    )
  })
})

describe("summarise_donor_locus_occupancy", {
  tissue_by_sample <- c(S1 = "T0", S2 = "T1", S3 = "T2")
  repertoires <- list(
    S1 = build_sharing_units(c("K1", "K2", "K3"), c(10, 5, 1)),
    S2 = build_sharing_units(c("K2", "K3", "K4"), c(20, 4, 6)),
    S3 = build_sharing_units(c("K5", "K6"), c(1, 1))
  )
  rows <- build_pairwise_overlap_rows(repertoires, "D1", "TRB", "clonotype", tissue_by_sample)
  summary_row <- summarise_donor_locus_occupancy(repertoires, rows, "D1", "TRB", "clonotype", 3L)

  it("writes the released columns in the declared order",
     expect_identical(names(summary_row), OVERLAP_SCHEMAS$donor_overlap_summary))
  it("is one row", expect_identical(nrow(summary_row), 1L))
  it("counts the union of every contributing library's clones",
     expect_identical(summary_row$n_groups_union, 6L))
  it("counts the clones seen in two or more libraries",
     expect_identical(summary_row$n_groups_shared_2plus, 2L))
  it("counts the clones seen in all of them",
     expect_identical(summary_row$n_groups_shared_all, 0L))
  it("names the unit the counts hold",
     expect_identical(summary_row$group_unit, "clonotype"))
  it("counts every unordered pair and those available", {
    expect_identical(summary_row$n_pairs_total, 3L)
    expect_identical(summary_row$n_pairs_available, 3L)
  })
  it("leaves below-minimum pairs out of the available count", {
    shallow <- build_pairwise_overlap_rows(repertoires, "D1", "TRB", "clonotype", tissue_by_sample,
                                           minimum_consensus_count = 17)
    row <- summarise_donor_locus_occupancy(repertoires, shallow, "D1", "TRB", "clonotype", 3L)
    expect_identical(row$n_pairs_total, 3L)
    expect_identical(row$n_pairs_available, 0L)
  })
  it("averages jaccard over the pairwise rows",
     expect_equal(summary_row$mean_jaccard, 1 / 6,
                  tolerance = EXACT))

  # The means are the coverage-gated pair of columns, so a dropped-out library
  # must suppress them even though the surviving pairs are perfectly comparable.
  it("withholds both means when a library contributed nothing", {
    incomplete <- summarise_donor_locus_occupancy(repertoires, rows, "D1", "TRB", "clonotype", 4L)
    expect_identical(
      unname(unlist(incomplete[c("mean_jaccard",
                                 "mean_morisita_horn")])),
      c(NA_real_, NA_real_)
    )
  })
  it("still reports the donor's and the nonempty library counts when it does", {
    incomplete <- summarise_donor_locus_occupancy(repertoires, rows, "D1", "TRB", "clonotype", 4L)
    expect_identical(incomplete$n_libraries_in_donor, 4L)
    expect_identical(incomplete$n_libraries_nonempty, 3L)
    expect_identical(incomplete$n_pairs_total, 6L)
  })
  it("counts a clone as shared by all only when all expected libraries hold it", {
    everywhere <- list(
      S1 = build_sharing_units("K1", 1),
      S2 = build_sharing_units("K1", 2)
    )
    expect_identical(
      summarise_donor_locus_occupancy(
        everywhere, build_pairwise_overlap_rows(everywhere, "D1", "TRB", "clonotype",
                                                tissue_by_sample),
        "D1", "TRB", "clonotype", 2L)$n_groups_shared_all,
      1L
    )
  })
  it("withholds both means for a donor with no comparable pair", {
    lone <- repertoires["S1"]
    alone <- summarise_donor_locus_occupancy(
      lone, make_empty_overlap_tables()$donor_overlap_pairs, "D1", "TRB", "clonotype", 1L)
    expect_identical(
      unname(unlist(alone[c("mean_jaccard",
                            "mean_morisita_horn")])),
      c(NA_real_, NA_real_)
    )
  })
  it("reports an empty locus as zero clones", {
    none <- summarise_donor_locus_occupancy(
      list(), make_empty_overlap_tables()$donor_overlap_pairs, "D1", "IGH", "clone", 2L)
    expect_identical(none$n_groups_union, 0L)
    expect_identical(none$n_libraries_nonempty, 0L)
    expect_identical(none$group_unit, "clone")
    expect_identical(none$n_pairs_available, 0L)
  })
})

# --------------------------------------------- which clones are persistent ---
#
# The listing and the two summaries are produced by one calculation, so what is
# asserted here is that they cannot disagree: the clones listed are exactly the
# clones counted, and each library's support for one of them is the support the
# indices were computed over.

SHARED_REPERTOIRES <- list(
  S1 = build_sharing_units(c("K1", "K2", "K3"), c(60, 30, 10)),
  S2 = build_sharing_units(c("K1", "K4"), c(20, 80))
)
SHARED_TISSUES <- c(S1 = "T1", S2 = "T2", S3 = "T3")

SHARED_RECORDS <- data.frame(
  group_key = c("K1", "K1", "K2", "K3", "K4"),
  consensus_count = c(60, 20, 30, 10, 80),
  record_id = c("r1", "r2", "r3", "r4", "r5"),
  representative_junction = c("ACGTAAA", "ACGTAAG", "TTTTAAA", "GGGGAAA", "CCCCAAA"),
  representative_junction_aa = c("TK", "TK", "FK", "GK", "PK"),
  v_gene_key = "IGHV1-69",
  v_gene_names = "IGHV1-69",
  j_gene_key = "IGHJ4",
  j_gene_names = "IGHJ4",
  representative_sequence_id = paste0("S", c(1, 2, 1, 1, 2), "~r", 1:5),
  sequence = paste0("CONTIG", 1:5),
  stringsAsFactors = FALSE
)

describe("summarise_shared_clone_descriptors", {
  descriptors <- summarise_shared_clone_descriptors(SHARED_RECORDS)

  it("describes each clone once",
     expect_identical(descriptors$group_key, c("K1", "K2", "K3", "K4")))
  # K1 was seen at 60 reads in S1 and 20 in S2; the better-supported observation
  # is the one the row shows.
  it("takes the best-supported observation of a clone", {
    expect_identical(descriptors$representative_junction[1], "ACGTAAA")
    expect_identical(descriptors$representative_sequence_id[1], "S1~r1")
  })
  it("carries the descriptive columns the schema declares",
     expect_identical(names(descriptors),
                      c("group_key", SHARED_CLONE_DESCRIPTOR_COLUMNS)))
})

test_that("a tie between two observations is broken by content, not record ID", {
  tied <- SHARED_RECORDS[1:2, ]
  tied$consensus_count <- c(20, 20)
  forwards <- summarise_shared_clone_descriptors(tied)
  backwards <- summarise_shared_clone_descriptors(tied[2:1, ])
  expect_identical(forwards$representative_junction, backwards$representative_junction)
  expect_identical(forwards$representative_junction, "ACGTAAA")
  relabelled <- tied
  relabelled$record_id <- rev(tied$record_id)
  expect_identical(summarise_shared_clone_descriptors(relabelled)$representative_junction, "ACGTAAA")
  # With one junction and one pair of entities, the full sequence decides.
  same_junction <- tied
  same_junction$representative_junction <- "ACGTAAA"
  same_junction$sequence <- c("TTTT", "AAAA")
  expect_identical(summarise_shared_clone_descriptors(same_junction)$representative_sequence_id, "S2~r2")
  expect_false("sequence" %in% names(summarise_shared_clone_descriptors(same_junction)))
})

describe("build_donor_locus_shared_clone_rows", {
  descriptors <- summarise_shared_clone_descriptors(SHARED_RECORDS)
  rows <- build_donor_locus_shared_clone_rows(
    SHARED_REPERTOIRES, descriptors, "D1", "TRB", "clonotype",
    c("S1", "S2"), SHARED_TISSUES)

  it("declares the released columns in order",
     expect_identical(names(rows), OVERLAP_SCHEMAS$shared_groups))
  # Only K1 is in both libraries, and it gets a row for each of them.
  it("lists only the clones detected more than once",
     expect_identical(unique(rows$group_key), "K1"))
  it("gives the clone a row for every library of the donor",
     expect_identical(rows$library_id, c("S1", "S2")))
  it("reports the support each library saw it at",
     expect_equal(rows$consensus_count, c(60, 20), tolerance = EXACT))
  it("takes each share over its own library's repertoire",
     expect_equal(rows$support_fraction, c(60 / 100, 20 / 100),
                  tolerance = EXACT))
  it("repeats the descriptor so the row stands alone",
     expect_identical(unique(rows$representative_junction), "ACGTAAA"))
  it("says how many libraries detected it and how many the donor has", {
    expect_identical(unique(rows$n_libraries_detected), 2L)
    expect_identical(unique(rows$n_libraries_in_donor), 2L)
  })
  it("names the unit the key identifies",
     expect_identical(unique(rows$group_unit), "clonotype"))
  it("attaches each library's tissue",
     expect_identical(rows$tissue_id, c("T1", "T2")))
})

# The listing and the summary count the same thing, so this is the identity a
# reader is entitled to rely on when they filter the released table themselves.
test_that("the clones listed are exactly the clones the summary counts", {
  descriptors <- summarise_shared_clone_descriptors(SHARED_RECORDS)
  rows <- build_donor_locus_shared_clone_rows(
    SHARED_REPERTOIRES, descriptors, "D1", "TRB", "clonotype",
    c("S1", "S2"), SHARED_TISSUES)
  summary_row <- summarise_donor_locus_occupancy(
    SHARED_REPERTOIRES,
    build_pairwise_overlap_rows(SHARED_REPERTOIRES, "D1", "TRB", "clonotype", SHARED_TISSUES),
    "D1", "TRB", "clonotype", 2L)

  expect_identical(length(unique(rows$group_key)),
                   summary_row$n_groups_shared_2plus)
})

# A library that contributed nothing at a locus is not the same as one where the
# clone was looked for and not found, so its share is missing rather than zero.
test_that("a library that contributed nothing gets zero support and no share", {
  descriptors <- summarise_shared_clone_descriptors(SHARED_RECORDS)
  rows <- build_donor_locus_shared_clone_rows(
    SHARED_REPERTOIRES, descriptors, "D1", "TRB", "clonotype",
    c("S1", "S2", "S3"), SHARED_TISSUES)

  empty_library <- rows[rows$library_id == "S3", ]
  expect_identical(nrow(empty_library), 1L)
  expect_equal(empty_library$consensus_count, 0, tolerance = EXACT)
  expect_identical(empty_library$support_fraction, NA_real_)
  # Occupancy is over nonempty libraries; the donor count is over all of them,
  # which is what makes an empty library visible.
  expect_identical(empty_library$n_libraries_detected, 2L)
  expect_identical(empty_library$n_libraries_in_donor, 3L)
})

test_that("an undetected shared clone has zero share in a contributing library", {
  repertoires <- c(SHARED_REPERTOIRES, list(S3 = build_sharing_units("K9", 7)))
  rows <- build_donor_locus_shared_clone_rows(
    repertoires, summarise_shared_clone_descriptors(SHARED_RECORDS), "D1", "TRB",
    "clonotype", c("S1", "S2", "S3"), SHARED_TISSUES)
  undetected <- rows[rows$library_id == "S3", ]
  expect_identical(undetected$group_key, "K1")
  expect_equal(undetected$consensus_count, 0)
  expect_equal(undetected$support_fraction, 0)
  expect_identical(undetected$n_libraries_detected, 2L)
  expect_identical(undetected$n_libraries_in_donor, 3L)
})

# A reader scanning the table meets the most supported shared clones first, and
# every library's view of a clone stays together.
test_that("shared clones run by total support, then key, then library", {
  repertoires <- list(
    S2 = build_sharing_units(c("K1", "K2", "K3"), c(5, 1, 6)),
    S1 = build_sharing_units(c("K1", "K2", "K3"), c(5, 30, 4))
  )
  rows <- build_donor_locus_shared_clone_rows(
    repertoires, summarise_shared_clone_descriptors(SHARED_RECORDS), "D1", "TRB",
    "clonotype", c("S2", "S1"), SHARED_TISSUES)
  # K2 carries 31 support units in all, K1 and K3 10 each.
  expect_identical(rows$group_key, c("K2", "K2", "K1", "K1", "K3", "K3"))
  expect_identical(rows$library_id, rep(c("S1", "S2"), 3))
  expect_identical(rownames(rows), as.character(seq_len(6)))
})

test_that("a donor sharing no clone at a locus still declares the columns", {
  disjoint <- list(
    S1 = build_sharing_units("K1", 5),
    S2 = build_sharing_units("K9", 5)
  )
  rows <- build_donor_locus_shared_clone_rows(
    disjoint, summarise_shared_clone_descriptors(SHARED_RECORDS), "D1", "TRB",
    "clonotype", c("S1", "S2"), SHARED_TISSUES)

  expect_identical(nrow(rows), 0L)
  expect_identical(names(rows), OVERLAP_SCHEMAS$shared_groups)
})

# Constructed longitudinal cases distinguish redistribution from the number of
# redetected groups changing direction. Expectations are analytic. The
# concentration difference, the shared-group fractions and balance and the
# support only at one endpoint follow from these values and diversity.tsv, and
# the report derives them; scientific_values.test.js checks that derivation.
describe("longitudinal repertoire measurements", {
  units <- function(counts) build_sharing_units(names(counts), unname(counts))
  compare_pair <- function(a, b, ...) {
    calculate_repertoire_pair_overlap(units(a), units(b), ...)
  }

  it("recognises identical proportions at different support depths", {
    row <- compare_pair(c(a = 1, b = 2, c = 3), c(a = 10, b = 20, c = 30))
    expect_equal(row$bray_curtis_proportion, 0)
    expect_equal(row$n_shared_higher_in_2, 0)
    expect_equal(row$n_shared_higher_in_1, 0)
    expect_equal(row$n_shared_equal, 3)
  })

  it("distinguishes replacement from increased dominance", {
    row <- compare_pair(c(a = 8, b = 2), c(a = 2, b = 8))
    expect_equal(row$bray_curtis_proportion, .6)
    expect_equal(row$n_shared_higher_in_2, 1)
    expect_equal(row$n_shared_higher_in_1, 1)
  })

  it("allows most groups to decrease while dominance increases", {
    before <- setNames(rep(10, 10), letters[1:10])
    after <- setNames(c(110, rep(10, 9)), letters[1:10])
    row <- compare_pair(before, after)
    expect_equal(row$bray_curtis_proportion, .45)
    expect_equal(row$n_shared_higher_in_2, 1)
    expect_equal(row$n_shared_higher_in_1, 9)
  })

  it("preserves union denominators and reverses all directional counts", {
    a <- c(a = 8, b = 2)
    b <- c(a = 4, c = 6)
    row <- compare_pair(a, b)
    reverse <- compare_pair(b, a)
    expect_equal(row$bray_curtis_proportion, .6)
    expect_equal(row$n_shared_higher_in_1, 1)
    # 1 - shared_support_fraction_k is endpoint k's support only there.
    expect_equal(1 - row$shared_support_fraction_2, .6)
    expect_equal(1 - row$shared_support_fraction_1, .2)
    expect_equal(reverse$bray_curtis_proportion, row$bray_curtis_proportion)
    expect_equal(reverse$n_shared_higher_in_2, row$n_shared_higher_in_1)
    expect_equal(reverse$n_shared_higher_in_1, row$n_shared_higher_in_2)
    expect_equal(reverse$n_shared_equal, row$n_shared_equal)
    expect_equal(reverse$shared_support_fraction_2, row$shared_support_fraction_1)
  })

  it("separates disjoint nonempty repertoires from empty endpoints", {
    row <- compare_pair(c(a = 1), c(b = 1))
    expect_equal(row$bray_curtis_proportion, 1)
    expect_equal(row$n_groups_shared, 0L)
    expect_equal(row$n_shared_higher_in_2, 0)
    expect_equal(row$n_shared_higher_in_1, 0)
    expect_equal(row$n_shared_equal, 0)
    expect_equal(row$pair_status, "available")
    for (inputs in list(list(numeric(), c(a = 1)), list(c(a = 1), numeric()),
                        list(numeric(), numeric()))) {
      empty <- compare_pair(inputs[[1]], inputs[[2]])
      for (field in c("bray_curtis_proportion", "jaccard", "morisita_horn",
                      "n_shared_higher_in_2", "n_shared_higher_in_1", "n_shared_equal")) {
        expect_true(is.na(empty[[field]]), info = field)
      }
    }
    expect_equal(compare_pair(numeric(), c(a = 1))$pair_status, "empty_first")
    expect_equal(compare_pair(c(a = 1), numeric())$pair_status, "empty_second")
    expect_equal(compare_pair(numeric(), numeric())$pair_status, "empty_both")
  })

  it("flags pair support without masking scores", {
    a <- c(a = 8, b = 2)
    b <- c(a = 20, b = 80)
    base <- compare_pair(a, b)
    boundary <- compare_pair(a, b, minimum_consensus_count = 10)
    low <- compare_pair(a, b, minimum_consensus_count = 11)
    expect_equal(boundary$pair_status, "available")
    expect_equal(low$pair_status, "below_minimum_support")
    expect_equal(low$min_consensus_count_per_library, 11)
    expect_equal(low$bray_curtis_proportion, base$bray_curtis_proportion)
    expect_equal(low$n_shared_higher_in_2, base$n_shared_higher_in_2)
    expect_equal(compare_pair(b, a, minimum_consensus_count = 11)$pair_status,
                 "below_minimum_support")
    for (minimum in list(-1, 1.5, NA_real_, Inf)) {
      expect_error(compare_pair(a, b, minimum_consensus_count = minimum), "non-negative integer")
    }
  })

  it("counts small nonzero changes without an effect-size tolerance", {
    row <- compare_pair(c(a = 1000000, b = 1000000), c(a = 1000001, b = 1000000))
    expect_equal(row$n_shared_higher_in_2, 1)
    expect_equal(row$n_shared_higher_in_1, 1)
    expect_equal(row$n_shared_equal, 0)
    expect_gt(row$bray_curtis_proportion, 0)
  })
})

test_that("sharing keeps endpoint denominators and unavailable pairs", {
  first <- build_sharing_units(c("a", "b"), c(90, 10))
  second <- build_sharing_units(c("a", "c"), c(5, 5))
  pair <- calculate_repertoire_pair_overlap(first, second)
  expect_equal(pair$shared_consensus_count_1, 90)
  expect_equal(pair$shared_consensus_count_2, 5)
  expect_equal(pair$shared_support_fraction_1, .9)
  expect_equal(pair$shared_support_fraction_2, .5)
  expect_equal(pair$morisita_horn, .9 / 1.32)
  # Scaling both counts preserves the observed proportions. A shallow draw
  # that misses the minor shared support can still change the overlap.
  scaled <- build_sharing_units(c("a", "c"), c(50, 50))
  expect_equal(calculate_repertoire_pair_overlap(first, scaled)$morisita_horn,
               pair$morisita_horn)
  missed <- build_sharing_units("c", 1)
  expect_equal(calculate_repertoire_pair_overlap(first, missed)$morisita_horn, 0)
  rows <- build_pairwise_overlap_rows(list(a = first, b = second), "0460", "IGH", "clone",
                                      c(a = "t1", b = "t2", c = "t3"), c("a", "b", "c"))
  expect_equal(nrow(rows), 3)
  expect_identical(rows$pair_status, c("available", "empty_second", "empty_second"))
  expect_true(all(is.na(rows$morisita_horn[2:3])))
  expect_true(all(is.na(rows$shared_support_fraction_2[2:3])))
  expect_equal(rows$shared_support_fraction_1[2:3], c(0, 0))
})
