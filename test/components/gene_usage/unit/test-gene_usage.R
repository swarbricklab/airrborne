# Check gene usage populations and clone gene-vote outcomes.
# Run: python test/run_tests.py --component gene_usage --level unit

load_test_modules(c("table_io", "gene_calls", "gene_usage"))

EXACT <- 1e-12

# A composite entity's reference facts. IGHV3-30 and IGHV3-30-3 share one
# allele, so both genes and the composite form one equivalence group. Single
# genes outside the lookup name their own group.
REFERENCE <- list(
  labels = c(VSEQ_example = "IGHV3-30,IGHV3-30-3"),
  equivalence_groups = c(`IGHV3-30` = "IGHV3-30,IGHV3-30-3", `IGHV3-30-3` = "IGHV3-30,IGHV3-30-3",
                         VSEQ_example = "IGHV3-30,IGHV3-30-3"),
  alleles = c(VSEQ_example = "IGHV3-30*04,IGHV3-30-3*03")
)

# Seven units with V genes and one with none. IGHV3-30 and the composite share a
# group whose two entities together outnumber IGHV1-69, so group order is
# exercised as well as the order inside a group and a tie between groups.
UNITS <- data.frame(
  v_gene = c("IGHV3-23", "IGHV1-69", "IGHV3-23", "VSEQ_example",
             "IGHV1-69", "IGHV3-30", "VSEQ_example", ""),
  support = c(5, 40, 3, 2, 1, 7, 2, 9),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------- count_gene_usage ---

describe("a typical set of units", {
  usage <- count_gene_usage(UNITS$v_gene, UNITS$support, REFERENCE)

  it("declares the released columns", {
    expect_identical(names(usage), GENE_USAGE_COUNT_COLUMNS)
  })
  it("orders equivalence groups by their unit count and keeps a group's rows together", {
    expect_identical(usage$gene_equivalence_group,
                     c("IGHV3-30,IGHV3-30-3", "IGHV3-30,IGHV3-30-3", "IGHV1-69", "IGHV3-23"))
    expect_identical(usage$gene_key, c("VSEQ_example", "IGHV3-30", "IGHV1-69", "IGHV3-23"))
  })
  it("counts every unit carrying each entity", {
    expect_identical(usage$n_groups, c(2L, 1L, 2L, 2L))
  })
  # Seven of the eight units had a V gene; the eighth is in neither denominator.
  it("takes the group share over the units that had a gene", {
    expect_equal(usage$group_fraction, c(2, 1, 2, 2) / 7, tolerance = EXACT)
  })
  it("sums support and takes its share over the units that had a gene", {
    expect_equal(usage$consensus_count, c(4, 7, 41, 8), tolerance = EXACT)
    expect_equal(usage$support_fraction, c(4, 7, 41, 8) / 60, tolerance = EXACT)
  })
  it("names a composite's genes and identical alleles, and a gene by itself", {
    expect_identical(usage$gene_names[1:2], c("IGHV3-30,IGHV3-30-3", "IGHV3-30"))
    expect_identical(usage$identical_alleles[1:2], c("IGHV3-30*04,IGHV3-30-3*03", ""))
  })
  it("sums each fraction to one over the equivalence groups", {
    expect_equal(sum(tapply(usage$group_fraction, usage$gene_equivalence_group, sum)), 1, tolerance = EXACT)
    expect_equal(sum(tapply(usage$support_fraction, usage$gene_equivalence_group, sum)), 1, tolerance = EXACT)
  })
})

describe("one unit", {
  usage <- count_gene_usage("IGHV3-23", 12, REFERENCE)

  it("counts it once", expect_identical(usage$n_groups, 1L))
  it("gives it the whole of both shares", {
    expect_equal(usage$group_fraction, 1, tolerance = EXACT)
    expect_equal(usage$support_fraction, 1, tolerance = EXACT)
  })
  it("makes an unlinked gene its own group", expect_identical(usage$gene_equivalence_group, "IGHV3-23"))
})

describe("nothing to count", {
  for (case in list(list(what = "no units at all", sets = character(0), support = numeric(0)),
                    list(what = "no unit with a gene", sets = c("", ""), support = c(1, 2)))) {
    usage <- count_gene_usage(case$sets, case$support, REFERENCE)

    it(sprintf("returns no rows [%s]", case$what),
       expect_identical(nrow(usage), 0L))
    it(sprintf("still declares the released columns and their types [%s]", case$what), {
      expect_identical(
        vapply(usage, class, character(1)),
        c(gene_key = "character", gene_names = "character", gene_equivalence_group = "character",
          identical_alleles = "character", n_groups = "integer", group_fraction = "numeric",
          consensus_count = "numeric", support_fraction = "numeric"))
    })
  }
})

test_that("a composite entity missing from the reference stops the count", {
  expect_error(count_gene_usage("VSEQ_unknown", 1, REFERENCE), "lacks the composite entities")
})

# ---------------------------------------------------- count_vj_pair_usage ---

# Seven units with both entities. IGHV3-23 pairs with two J genes, and the
# composite entity counts under its own key beside the gene IGHV3-30 it lists.
PAIR_UNITS <- data.frame(
  v_gene = c("IGHV3-23", "IGHV1-69", "IGHV3-23", "VSEQ_example",
             "IGHV1-69", "IGHV3-30", "VSEQ_example"),
  j_gene = c("IGHJ4", "IGHJ4", "IGHJ6", "IGHJ4", "IGHJ4", "IGHJ4", "IGHJ4"),
  support = c(5, 40, 3, 2, 1, 7, 2),
  stringsAsFactors = FALSE
)

describe("a typical set of pairs", {
  usage <- count_vj_pair_usage(PAIR_UNITS$v_gene, PAIR_UNITS$j_gene, PAIR_UNITS$support, REFERENCE)

  it("declares the released columns", {
    expect_identical(names(usage), VJ_PAIR_USAGE_COUNT_COLUMNS)
  })
  it("orders pairs by their unit count, then by V key and J key", {
    expect_identical(usage$v_gene_key,
                     c("IGHV1-69", "VSEQ_example", "IGHV3-23", "IGHV3-23", "IGHV3-30"))
    expect_identical(usage$j_gene_key, c("IGHJ4", "IGHJ4", "IGHJ4", "IGHJ6", "IGHJ4"))
  })
  it("counts every unit carrying each pair", {
    expect_identical(usage$n_groups, c(2L, 2L, 1L, 1L, 1L))
  })
  it("takes the group share over the locus's units", {
    expect_equal(usage$group_fraction, c(2, 2, 1, 1, 1) / 7, tolerance = EXACT)
  })
  it("sums support and takes its share over the locus's support", {
    expect_equal(usage$consensus_count, c(41, 4, 5, 3, 7), tolerance = EXACT)
    expect_equal(usage$support_fraction, c(41, 4, 5, 3, 7) / 60, tolerance = EXACT)
  })
  it("counts a composite entity once under its own key and names its genes", {
    composite <- usage[usage$v_gene_key == "VSEQ_example", ]
    expect_identical(composite$v_gene_names, "IGHV3-30,IGHV3-30-3")
    expect_identical(composite$n_groups, 2L)
    expect_identical(usage$n_groups[usage$v_gene_key == "IGHV3-30"], 1L)
  })
  it("names a gene by itself", {
    expect_identical(usage$j_gene_names, usage$j_gene_key)
  })
  for (segment in c("v", "j")) {
    it(sprintf("gives gene usage when summed over the other segment [%s]", segment), {
      single <- count_gene_usage(PAIR_UNITS[[paste0(segment, "_gene")]], PAIR_UNITS$support, REFERENCE)
      keys <- usage[[paste0(segment, "_gene_key")]]
      for (column in c("n_groups", "group_fraction", "consensus_count", "support_fraction")) {
        expect_equal(as.numeric(tapply(usage[[column]], keys, sum)[single$gene_key]),
                     as.numeric(single[[column]]), tolerance = EXACT, info = column)
      }
    })
  }
})

describe("a unit lacking one segment", {
  usage <- count_vj_pair_usage(c("IGHV1-69", "", "IGHV3-23"), c("IGHJ4", "IGHJ4", ""),
                               c(6, 50, 50), REFERENCE)

  it("has no pair and leaves both denominators", {
    expect_identical(usage$n_groups, 1L)
    expect_equal(usage$group_fraction, 1, tolerance = EXACT)
    expect_equal(usage$support_fraction, 1, tolerance = EXACT)
  })
})

describe("no pair to count", {
  for (case in list(list(what = "no units at all", v = character(0), j = character(0), support = numeric(0)),
                    list(what = "no unit with both", v = c("IGHV1-69", ""), j = c("", "IGHJ4"),
                         support = c(1, 2)))) {
    usage <- count_vj_pair_usage(case$v, case$j, case$support, REFERENCE)

    it(sprintf("returns no rows [%s]", case$what), expect_identical(nrow(usage), 0L))
    it(sprintf("still declares the released columns and their types [%s]", case$what), {
      expect_identical(
        vapply(usage, class, character(1)),
        c(v_gene_key = "character", v_gene_names = "character", j_gene_key = "character",
          j_gene_names = "character", n_groups = "integer", group_fraction = "numeric",
          consensus_count = "numeric", support_fraction = "numeric"))
    })
  }
})

# ------------------------------------------------ read_vj_entity_labels ---

write_reference_metadata <- function(rows) {
  path <- tempfile(fileext = ".tsv")
  utils::write.table(rows, path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}

describe("reference entity labels", {
  rows <- data.frame(
    allele = c("IGLJ2*01", "IGLJ3*01", "IGLJ3*02", "IGHA1*01"),
    identity_id = c("JSEQ_a", "JSEQ_a", "JSEQ_b", ""),
    identity_label = c("IGLJ2,IGLJ3", "IGLJ2,IGLJ3", "IGLJ3", ""),
    stringsAsFactors = FALSE)

  it("names each identifier once and skips alleles without one", {
    expect_identical(read_vj_entity_labels(write_reference_metadata(rows)),
                     c(JSEQ_a = "IGLJ2,IGLJ3", JSEQ_b = "IGLJ3"))
  })
  it("rejects an identifier with two labels", {
    rows$identity_label[2] <- "IGLJ3"
    expect_error(read_vj_entity_labels(write_reference_metadata(rows)),
                 "more than one identity_label")
  })
})

describe("reference gene equivalence groups", {
  # IGLJ2 and IGLJ9 share no allele but are joined through IGLJ3. The constant
  # region row is outside V and J.
  rows <- data.frame(
    allele = c("IGLJ2*01", "IGLJ3*01", "IGLJ3*02", "IGLJ3*03", "IGLJ9*01", "IGLJ7*01", "IGHA1*01"),
    gene = c("IGLJ2", "IGLJ3", "IGLJ3", "IGLJ3", "IGLJ9", "IGLJ7", "IGHA1"),
    region = c(rep("J-REGION", 6), "C-REGION"),
    identity_id = c("JSEQ_a", "JSEQ_a", "JSEQ_b", "JSEQ_d", "JSEQ_d", "JSEQ_c", ""),
    identity_label = c("IGLJ2,IGLJ3", "IGLJ2,IGLJ3", "IGLJ3", "IGLJ3,IGLJ9", "IGLJ3,IGLJ9", "IGLJ7", ""),
    stringsAsFactors = FALSE)
  reference <- read_vj_entity_reference(write_reference_metadata(rows))

  it("joins genes linked directly or through a third gene", {
    expect_identical(unname(reference$equivalence_groups[c("IGLJ2", "IGLJ3", "IGLJ9", "JSEQ_a", "JSEQ_d")]),
                     rep("IGLJ2,IGLJ3,IGLJ9", 5))
  })
  it("makes a gene sharing no sequence its own group", {
    expect_identical(unname(reference$equivalence_groups["IGLJ7"]), "IGLJ7")
    expect_false("IGHA1" %in% names(reference$equivalence_groups))
  })
  it("lists the identical alleles of composite entities only", {
    expect_identical(reference$alleles,
                     c(JSEQ_a = "IGLJ2*01,IGLJ3*01", JSEQ_d = "IGLJ3*03,IGLJ9*01"))
  })
})

test_that("an unassigned entity keeps an empty label", {
  expect_identical(vj_entity_display_label(c("", "VSEQ_example"), REFERENCE$labels),
                   c("", "IGHV3-30,IGHV3-30-3"))
})

# --------------------------------------- summarise_clone_vj_call_sets ---

# The canonical gene columns are written out rather than derived, because they are
# what the reannotation stage puts on each row and the vote reads them as given.
# Here no two of these genes are one germline, so each is its own entity and the
# allele suffix and the duplicate-gene D alias are all that collapse.
MEMBERSHIP <- data.frame(
  library_id = c("WTS_1_S1", "WTS_1_S1", "WTS_2_S1", "WTS_2_S1", "WTS_2_S1"),
  clone_id = c("cohort_D1_L2", "cohort_D1_L2", "cohort_D1_L1",
                 "cohort_D1_L2", ""),
  v_call = c("IGHV3-23*01", "IGHV3-23D*01", "IGHV1-69", "IGHV4-34", "IGHV1-8"),
  j_call = c("IGHJ4*02", "IGHJ5", "IGHJ6", "IGHJ4", "IGHJ4"),
  v_gene = c("IGHV3-23", "IGHV3-23", "IGHV1-69", "IGHV4-34",
                       "IGHV1-8"),
  j_gene = c("IGHJ4", "IGHJ5", "IGHJ6", "IGHJ4", "IGHJ4"),
  stringsAsFactors = FALSE
)

describe("a scope's lineage genes", {
  calls <- summarise_clone_vj_call_sets(MEMBERSHIP)

  it("gives one row per lineage",
     expect_identical(calls$clone_id, c("", "cohort_D1_L1", "cohort_D1_L2")))
  # The vote is taken over every in-scope member, whichever library it came
  # from, so both of WTS_2_S1's and WTS_1_S1's contributions to L2 are counted.
  # Two of L2's members are IGHV3-23 under two spellings and one is IGHV4-34, so
  # its resolved members name two genes.
  it("excludes a lineage whose resolved members name different V genes", {
    expect_identical(calls$v_vote_status[calls$clone_id == "cohort_D1_L2"],
                     "resolved_members_disagree")
    expect_identical(calls$v_gene[calls$clone_id == "cohort_D1_L2"],
                     "")
  })
  it("excludes it on J for the same reason", {
    expect_identical(calls$j_vote_status[calls$clone_id == "cohort_D1_L2"],
                     "resolved_members_disagree")
  })
  it("assigns a lineage whose one member names one gene", {
    expect_identical(calls$v_vote_status[calls$clone_id == "cohort_D1_L1"],
                     "assigned")
    expect_identical(calls$v_gene[calls$clone_id == "cohort_D1_L1"],
                     "IGHV1-69")
  })
})

# An allele and an IMGT duplicate-gene alias are one gene, so members spelling it
# either way agree rather than disagreeing.
test_that("an allele and a duplicate-gene alias vote for one gene", {
  spellings <- MEMBERSHIP[1:2, ]
  calls <- summarise_clone_vj_call_sets(spellings)
  expect_identical(calls$v_vote_status, "assigned")
  expect_identical(calls$v_gene, "IGHV3-23")
})

describe("a lineage with one member", {
  calls <- summarise_clone_vj_call_sets(MEMBERSHIP[3, , drop = FALSE])

  it("reports that member's genes", {
    expect_identical(calls, data.frame(clone_id = "cohort_D1_L1",
                                       v_gene = "IGHV1-69",
                                       j_gene = "IGHJ6",
                                       v_vote_status = "assigned",
                                       j_vote_status = "assigned",
                                       stringsAsFactors = FALSE))
  })
})

# The four statuses, one lineage each, on the smallest frame that separates them.
# A member whose call spans two germlines has an empty canonical gene and
# abstains, which is what makes a lineage holding nothing else
# "members_ambiguous" and leaves one holding a decided member assigned to that
# member's germline.
describe("the vote status of each kind of lineage", {
  cases <- list(
    list(lineage = "L_assigned",
         v_call = c("IGHV3-23*01", "IGHV3-23*04"),
         v_token = c("IGHV3-23", "IGHV3-23"),
         status = "assigned", gene = "IGHV3-23",
         why = "members resolving to one germline assign it"),
    list(lineage = "L_none",
         v_call = c("", NA_character_),
         v_token = c("", ""),
         status = "no_callable_member", gene = "",
         why = "no member with a callable gene assigns nothing"),
    list(lineage = "L_ambiguous",
         v_call = c("IGHV3-23*01,IGHV1-8*01", "IGHV4-34*01,IGHV1-8*01"),
         v_token = c("", ""),
         status = "members_ambiguous", gene = "",
         why = "callable members that resolve to no germline assign nothing"),
    list(lineage = "L_disagree",
         v_call = c("IGHV3-23*01", "IGHV1-8*01"),
         v_token = c("IGHV3-23", "IGHV1-8"),
         status = "resolved_members_disagree", gene = "",
         why = "members resolving to different germlines assign nothing"),
    list(lineage = "L_mixed",
         v_call = c("IGHV3-23*01", "IGHV4-34*01,IGHV1-8*01"),
         v_token = c("IGHV3-23", ""),
         status = "assigned", gene = "IGHV3-23",
         why = "an undecided member abstains and the decided one carries it"),
    # IGHV3-30-22 and IGHV3-30-62 spell the same germline, so the reference gives
    # both members one token whether the call named one of them or both. Voting on
    # the name instead would read this lineage as holding two germlines.
    list(lineage = "L_one_germline_two_names",
         v_call = c("IGHV3-30-22*01,IGHV3-30-62*01", "IGHV3-30-22*01"),
         v_token = c("IGHV3-30-22,IGHV3-30-62", "IGHV3-30-22,IGHV3-30-62"),
         status = "assigned", gene = "IGHV3-30-22,IGHV3-30-62",
         why = "two IMGT names for one germline are one entity, not a disagreement")
  )
  for (case in cases) {
    calls <- summarise_clone_vj_call_sets(data.frame(
      clone_id = case$lineage, v_call = case$v_call, j_call = "IGHJ4*02",
      v_gene = case$v_token, j_gene = "IGHJ4",
      stringsAsFactors = FALSE))
    it(case$why, {
      expect_identical(calls$v_vote_status, case$status)
      expect_identical(calls$v_gene, case$gene)
    })
  }
})

# Every lineage lands in exactly one status, which is what lets the released
# n_groups column sum to the library's lineage count across its rows.
test_that("every lineage gets exactly one status", {
  calls <- summarise_clone_vj_call_sets(MEMBERSHIP)
  for (gene_type in c("v", "j")) {
    status <- calls[[paste0(gene_type, "_vote_status")]]
    expect_true(all(status %in% CLONE_VOTE_STATUSES))
    expect_identical(sum(table(factor(status, levels = CLONE_VOTE_STATUSES))),
                     length(unique(MEMBERSHIP$clone_id)))
  }
})

# A gene is named when and only when the vote was decided, so a reader can filter
# on either column and get the same rows.
test_that("a gene is named for exactly the assigned lineages", {
  calls <- summarise_clone_vj_call_sets(MEMBERSHIP)
  expect_identical(nzchar(calls$v_gene),
                   calls$v_vote_status == "assigned")
  expect_identical(nzchar(calls$j_gene),
                   calls$j_vote_status == "assigned")
})

describe("a scope with no membership at all", {
  calls <- summarise_clone_vj_call_sets(data.frame())

  it("returns no rows", expect_identical(nrow(calls), 0L))
  it("still declares the five columns", {
    expect_identical(names(calls),
                     c("clone_id", "v_gene", "j_gene",
                       "v_vote_status", "j_vote_status"))
  })
})

# A member that was never assigned to a lineage cannot be looked up by one, so
# it gets no row. The per-library aggregation drops those members anyway.
test_that("an unassigned member contributes no lineage", {
  unassigned <- MEMBERSHIP
  unassigned$clone_id <- NA_character_
  expect_identical(nrow(summarise_clone_vj_call_sets(unassigned)), 0L)
})
