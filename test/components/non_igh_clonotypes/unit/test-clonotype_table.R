# Check non-IGH nucleotide junction and V/J groups with conserved support.
# Run: python test/run_tests.py --component non_igh_clonotypes --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys",
  "rearrangement_groups"))

# --------------------------------------------- the rearrangement-group key ---

describe("exact_rearrangement_group_key", {
  # The key takes the germline entity each call resolved to, not the call, so
  # these are the canonical tokens the reannotation stage wrote.
  key_for <- function(junction) {
    exact_rearrangement_group_key(junction, "IGHV1-69", "IGHJ4")
  }
  # The three fields are joined by ~.
  group_separator <- "~"
  cases <- list(
    list(junction = "tgtgcgaaa",
         expected = paste("TGTGCGAAA", "IGHV1-69", "IGHJ4",
                          sep = group_separator),
         why = "uppercases the junction and joins it to the canonical V and J genes"),
    list(junction = "TGTNCGAAA", expected = NA_character_,
         why = "rejects a junction carrying an ambiguous base"),
    list(junction = "   ", expected = NA_character_,
         why = "rejects a blank junction"),
    list(junction = NA_character_, expected = NA_character_,
         why = "rejects a missing junction"),
    list(junction = "TGT-GCGAAA", expected = NA_character_,
         why = "rejects a junction carrying an alignment gap")
  )
  for (case in cases) {
    it(case$why, expect_identical(key_for(case$junction), case$expected))
  }
  # An empty token is what a call naming several germline entities produces, so
  # the drop is the existing eligibility test rather than a separate rule.
  it("rejects a record whose V names no single entity", {
    expect_identical(
      exact_rearrangement_group_key("TGTGCGAAA", "", "IGKJ1"), NA_character_)
  })
  it("rejects a record whose J names no single entity", {
    expect_identical(
      exact_rearrangement_group_key("TGTGCGAAA", "IGKV1-5", ""), NA_character_)
  })
})

# Nucleotide clonotype semantics, on pairs of records that differ in one way.
describe("what does and does not split two records into separate groups", {
  it("splits two distinct nucleotide junctions", {
    # The key compares nucleotide junctions directly.
    keys <- exact_rearrangement_group_key(
      c("TGTGCGAAA", "TGTGCGAAC"),
      c("IGKV1-5", "IGKV1-5"), c("IGKJ1", "IGKJ1"))
    expect_identical(length(unique(keys)), 2L)
  })
  it("splits two distinct V entities", {
    keys <- exact_rearrangement_group_key(
      c("TGTGCGAAA", "TGTGCGAAA"),
      c("IGKV1-5", "IGKV1-6"), c("IGKJ1", "IGKJ1"))
    expect_identical(length(unique(keys)), 2L)
  })
  # A complete-reference composite keeps its own analytical boundary.
  it("keeps a single gene separate from a composite", {
    calls <- c("TRBV6-2*01", "TRBV6-2*01,TRBV6-3*01")
    identity_classes <- c("TRBV6-2*01" = "VSEQ_TEST", "TRBV6-3*01" = "VSEQ_TEST")
    keys <- exact_rearrangement_group_key(
      c("TGTGCGAAA", "TGTGCGAAA"),
      vj_canonical_token(calls, identity_classes), c("TRBJ1-1", "TRBJ1-1"))
    expect_identical(length(unique(keys)), 2L)
  })
})

describe("make_rearrangement_group_id", {
  # The numeric suffix is the rank among SORTED unique keys, so "z_key" gets the
  # higher number even though it appears first.
  ids <- make_rearrangement_group_id("WTS_1_S1", "IGK",
                                     c("z_key", "a_key", "z_key", NA))

  it("namespaces the label with the library and locus and ranks sorted keys",
     expect_identical(ids[1:3], c("WTS_1_S1_IGK_ct00002", "WTS_1_S1_IGK_ct00001",
                                  "WTS_1_S1_IGK_ct00002")))
  it("passes a missing key through as missing",
     expect_identical(ids[[4]], NA_character_))
})

# The id assigned to a given key must not depend on the order the rows arrived
# in, or two libraries would label the same repertoire differently.
describe("make_rearrangement_group_id is invariant to row order", {
  forward <- make_rearrangement_group_id("g", "IGL", c("a_key", "z_key"))
  backward <- make_rearrangement_group_id("g", "IGL", c("z_key", "a_key"))

  it("gives the first key of one ordering the label the other gives it second",
     expect_identical(forward[[1]], backward[[2]]))
  it("gives the second key of one ordering the label the other gives it first",
     expect_identical(forward[[2]], backward[[1]]))
})

describe("make_rearrangement_group_id is collision-free", {
  ids <- make_rearrangement_group_id("g", "IGL", c("k1", "k2", "k1"))

  it("gives one key one label", expect_identical(ids[[1]], ids[[3]]))
  it("gives distinct keys distinct labels", expect_false(ids[[1]] == ids[[2]]))
})


EXACT <- 1e-12

# The four columns the reannotation stage writes beside each call, added here the
# way that stage adds them so the fixtures cannot drift from the classifier the
# pipeline runs. The empty lookup makes every gene its own germline entity, which
# is true of all but 12 V and J genes in the pinned human reference, which form
# six two-gene entities.
with_germline_entities <- function(df, identity_classes = character(0)) {
  df$v_gene_status <- vj_tie_kind(df$v_call, identity_classes)
  df$j_gene_status <- vj_tie_kind(df$j_call, identity_classes)
  df$v_gene <- vj_canonical_token(df$v_call, identity_classes)
  df$j_gene <- vj_canonical_token(df$j_call, identity_classes)
  df
}

# Seven records at one locus. The first two share a junction and a V and J
# entity, so they collapse into one group; the third and fourth differ in the
# junction and in the V gene respectively; the last three are ineligible, one for
# an ambiguous base in the junction, one for having no V call at all, and one for
# a V call naming two germline entities.
LOCUS_RECORDS <- with_germline_entities(data.frame(
  sequence_id = c("assemble0_0", "assemble1_0", "assemble2_0", "assemble3_0",
                  "assemble4_0", "assemble5_0", "assemble6_0"),
  productive = "T",
  junction = c("TGTGCGAAA", "tgtgcgaaa", "TGTGCGAAC", "TGTGCGAAC", "TGTNCGAAA",
               "TGTGCGAAT", "TGTGCGAAA"),
  junction_aa = c("CAK", "CAK", "CAN", "CAN", "CAX", "CAY", "CAK"),
  v_call = c("IGKV1-5*01", "IGKV1-5*03", "IGKV1-5*01", "IGKV1-6*01",
             "IGKV1-5*01", "", "IGKV1-5*01,IGKV1-12*01"),
  j_call = c("IGKJ1*01", "IGKJ1*02", "IGKJ1*01", "IGKJ1*01", "IGKJ1*01",
             "IGKJ1*01", "IGKJ1*01"),
  consensus_count = c("10", "3", "7", "1", "5", "2", "4"),
  sequence = c("AAACCCTTT", "AAACCCTTG", "TTTGGG", "TTTGGA", "CCCCCC",
               "GGGGGG", "AAACCCTTT"),
  stringsAsFactors = FALSE
))

describe("a locus with four eligible records in three groups", {
  built <- build_rearrangement_group_table(LOCUS_RECORDS, "WTS_1_S1", "IGK")
  groups <- built$groups

  it("writes the released columns in order", {
    expect_identical(names(groups), REARRANGEMENT_GROUP_COLUMNS)
  })
  it("reports the records it was given", expect_identical(built$n_records, 7L))
  it("reports the records with no usable key", expect_identical(built$n_dropped, 3L))
  # A record whose call names several germline entities is excluded for a
  # different reason than one with no call or an unreadable junction, and the two
  # are reported apart.
  it("reports separately the records naming several germline entities",
     expect_identical(built$n_dropped_ambiguous_vj, 1L))
  it("reports the records that reached the collapse",
     expect_identical(built$n_eligible, 4L))
  it("collapses them to one row per distinct rearrangement",
     expect_identical(nrow(groups), 3L))
  it("labels the groups in sorted-key order", {
    expect_identical(groups$clonotype_id,
                     paste0("WTS_1_S1_IGK_ct", c("00001", "00002", "00003")))
  })
  it("stamps every row with the library and locus", {
    expect_identical(unique(groups$library_id), "WTS_1_S1")
    expect_identical(unique(groups$locus), "IGK")
  })
  # 10 + 3 for the merged group, then 7 and 1.
  it("sums its members' read support",
     expect_equal(groups$consensus_count, c(13, 7, 1), tolerance = EXACT))
  # The merged group was assembled twice; the other two once each.
  it("counts the contigs its members came from",
     expect_identical(groups$n_contigs, c(2L, 1L, 1L)))
  # The merged group's two records differ in full DNA, so both remain members.
  it("counts its distinct retained member sequences", {
    expect_identical(groups$n_members, c(2L, 1L, 1L))
    expect_identical(as.vector(table(built$membership$clonotype_id)[groups$clonotype_id]),
                     groups$n_members)
  })
  # Each member carries its transferred support, which sums to its group's.
  it("releases each member's support under its member name", {
    expect_identical(names(built$membership), CLONOTYPE_MEMBER_COLUMNS)
    expect_equal(as.vector(tapply(built$membership$member_consensus_count,
                                  built$membership$clonotype_id, sum)[groups$clonotype_id]),
                 groups$consensus_count, tolerance = EXACT)
  })
  it("takes each group's share of the locus's eligible support",
     expect_equal(groups$support_fraction, c(13, 7, 1) / 21, tolerance = EXACT))
  # Two alleles of one gene are one entity, so the second record joins the first;
  # the displayed set is the entity the group is keyed on.
  it("displays the group's V entity",
     expect_identical(groups$v_gene_key,
                      c("IGKV1-5", "IGKV1-5", "IGKV1-6")))
  it("displays the group's J entity",
     expect_identical(groups$j_gene_key, rep("IGKJ1", 3)))
  # Descriptive fields come from one representative. Record 1 outranks record 2
  # on read support, so the uppercase junction and record 1's ID are what the row
  # shows.
  it("takes the junction from the representative",
     expect_identical(groups$junction[1], "TGTGCGAAA"))
  it("gives the junction's nucleotide length",
     expect_identical(groups$junction_length, nchar(groups$junction)))
  # build_clonotypes.R qualifies these library-local IDs with the library.
  it("names the representative's record, one of the group's members", {
    expect_identical(groups$representative_sequence_id,
                     c("assemble0_0", "assemble2_0", "assemble3_0"))
    expect_identical(built$membership$clonotype_id[match(
      groups$representative_sequence_id, built$membership$sequence_id)],
      groups$clonotype_id)
  })
  it("takes the amino-acid junction from the representative",
     expect_identical(groups$junction_aa, c("CAK", "CAN", "CAN")))
  # The key is content, so two libraries reporting the same rearrangement agree
  # on it; the label above is local to this library.
  # Raw UTF-8 separator bytes keep the expected key unlabelled under C collation.
  it("keys each group on its own content", {
    expect_identical(groups$clonotype_key,
                     paste(c("TGTGCGAAA", "TGTGCGAAC", "TGTGCGAAC"),
                           c("IGKV1-5", "IGKV1-5", "IGKV1-6"), "IGKJ1",
                           sep = "~"))
  })
})

describe("a locus with one eligible record", {
  built <- build_rearrangement_group_table(LOCUS_RECORDS[1, , drop = FALSE],
                                           "WTS_1_S1", "IGK")

  it("makes it one group", expect_identical(nrow(built$groups), 1L))
  it("gives it the whole support share",
     expect_equal(built$groups$support_fraction, 1, tolerance = EXACT))
  it("drops nothing", expect_identical(built$n_dropped, 0L))
  it("counts no germline-entity exclusion",
     expect_identical(built$n_dropped_ambiguous_vj, 0L))
})

# Two members alike on productivity, ambiguity and read support fall through to
# the full-sequence tie-break, so the row is reproducible whichever order the
# records arrived in and however TRUST4 numbered them.
test_that("a tie between two members is broken deterministically", {
  tied <- with_germline_entities(data.frame(
    sequence_id = c("assemble9_0", "assemble1_0"),
    productive = c("T", "T"), junction = c("TGTGCGAAA", "TGTGCGAAA"),
    junction_aa = c("CAK", "CAK"),
    v_call = c("IGKV1-5*01", "IGKV1-5*01"), j_call = c("IGKJ1*01", "IGKJ1*01"),
    consensus_count = c("4", "4"),
    sequence = c("AAAAAA", "GGGGGG"), stringsAsFactors = FALSE
  ))
  built <- build_rearrangement_group_table(tied, "WTS_1_S1", "IGK")

  expect_identical(nrow(built$groups), 1L)
  expect_identical(built$groups$representative_sequence_id, "assemble9_0")
  expect_equal(built$groups$consensus_count, 8, tolerance = EXACT)
})

describe("a locus with no group to report", {
  for (case in list(
    list(what = "no records at all", records = LOCUS_RECORDS[0, , drop = FALSE],
         n_records = 0L, n_dropped = 0L),
    list(what = "no record with a usable key",
         records = LOCUS_RECORDS[5:7, , drop = FALSE],
         n_records = 3L, n_dropped = 3L))) {
    built <- build_rearrangement_group_table(case$records, "WTS_1_S1", "IGK")

    it(sprintf("returns no rows [%s]", case$what),
       expect_identical(nrow(built$groups), 0L))
    it(sprintf("still declares the released columns [%s]", case$what), {
      expect_identical(names(built$groups), REARRANGEMENT_GROUP_COLUMNS)
    })
    it(sprintf("reports what it saw [%s]", case$what), {
      expect_identical(built$n_records, case$n_records)
      expect_identical(built$n_dropped, case$n_dropped)
      expect_identical(built$n_eligible, 0L)
    })
  }
  it("still separates the germline-entity exclusions", {
    built <- build_rearrangement_group_table(LOCUS_RECORDS[5:7, , drop = FALSE],
                                             "WTS_1_S1", "IGK")
    expect_identical(built$n_dropped_ambiguous_vj, 1L)
  })
  it("declares the counting columns as numbers", {
    empty <- make_empty_rearrangement_group_table()
    expect_identical(
      vapply(empty[c("junction_length", "n_contigs", "n_members", "consensus_count", "support_fraction")],
             class, character(1)),
      c(junction_length = "integer", n_contigs = "integer", n_members = "integer", consensus_count = "numeric",
        support_fraction = "numeric"))
  })
})

# TRUST4 reports alternative CDR3s on one contig as separate AIRR rows.
# n_contigs counts contributing contigs, so those rows count once.
test_that("alternative CDR3 rows count once per contributing assembly", {
  siblings <- with_germline_entities(data.frame(
    sequence_id = c("assemble7_0", "assemble7_1", "assemble8_0"),
    productive = c("T", "T", "T"),
    junction = rep("TGTGCGAAA", 3), junction_aa = rep("CAK", 3),
    v_call = rep("IGKV1-5*01", 3), j_call = rep("IGKJ1*01", 3),
    consensus_count = c("6", "2", "4"),
    sequence = rep("AAACCCTTT", 3), stringsAsFactors = FALSE
  ))
  built <- build_rearrangement_group_table(siblings, "WTS_1_S1", "IGK")

  expect_identical(nrow(built$groups), 1L)
  expect_identical(built$groups$n_contigs, 2L)
  # The two assemblies retain identical full DNA, which is one member.
  expect_identical(built$groups$n_members, 1L)
  expect_equal(built$groups$consensus_count, 12, tolerance = EXACT)
})

# Two spellings of one contig whose junctions differ. At the light chains they
# share one assembly representative, chosen by original read support;
# at the TCR loci each reported CDR3 is its own rearrangement. The two records
# differ only in the locus they are handed to and the gene names that go with it.
sibling_records <- function(v_gene, j_gene) {
  with_germline_entities(data.frame(
    sequence_id = c("assemble7_0", "assemble7_1"),
    productive = c("T", "T"),
    junction = c("TGTGCGAAA", "TGTGCGAAC"),
    junction_aa = c("CAK", "CAN"),
    v_call = paste0(v_gene, "*01"), j_call = paste0(j_gene, "*01"),
    consensus_count = c("6", "2"),
    sequence = c("AAATGTGCGAAATTT", "AAATGTGCGAACTTT"), stringsAsFactors = FALSE
  ))
}

describe("two spellings of one light-chain contig", {
  built <- build_rearrangement_group_table(
    sibling_records("IGKV1-5", "IGKJ1"), "WTS_1_S1", "IGK")

  it("become one group", expect_identical(nrow(built$groups), 1L))
  it("count as one contig", expect_identical(built$groups$n_contigs, 1L))
  it("keep every read between them",
     expect_equal(built$groups$consensus_count, 8, tolerance = EXACT))
  it("are represented by the better-supported spelling",
     expect_identical(built$groups$junction, "TGTGCGAAA"))
  # The collapse is not an exclusion: both rows were eligible and are counted.
  it("are both counted as having reached the collapse",
     expect_identical(built$n_eligible, 2L))
})

describe("the same two spellings at a T-cell locus", {
  built <- build_rearrangement_group_table(
    sibling_records("TRBV6-1", "TRBJ1-1"), "WTS_1_S1", "TRB")

  it("stay two groups", expect_identical(nrow(built$groups), 2L))
  it("keep their own read support",
     expect_equal(sort(built$groups$consensus_count), c(2, 6),
                  tolerance = EXACT))
  it("still sum to the same read support as the light-chain case",
     expect_equal(sum(built$groups$consensus_count), 8, tolerance = EXACT))
})

# The amino-acid description is optional; full DNA is required for grouping.
test_that("a record set without the optional descriptive columns still builds", {
  spare <- LOCUS_RECORDS[, setdiff(names(LOCUS_RECORDS),
                                   "junction_aa")]
  built <- build_rearrangement_group_table(spare, "WTS_1_S1", "IGK")

  expect_identical(names(built$groups), REARRANGEMENT_GROUP_COLUMNS)
  expect_identical(built$groups$junction_aa, rep(NA_character_, 3))
  expect_true(all(!is.na(built$groups$representative_sequence_id)))
  expect_error(build_rearrangement_group_table(spare[, setdiff(names(spare), "sequence")],
    "WTS_1_S1", "IGK"), "sequence")
})

# Read support at or below zero means the common QC gate upstream let something
# through that every consensus_count metric would then divide by, so it is refused
# rather than dropped.
test_that("non-positive read support is refused", {
  unsupported <- LOCUS_RECORDS
  unsupported$consensus_count[1] <- "0"
  expect_error(
    build_rearrangement_group_table(unsupported, "WTS_1_S1", "IGK"),
    "[IGK] common QC input contains non-positive read support", fixed = TRUE)
})

test_that("unreadable read support is refused", {
  unreadable <- LOCUS_RECORDS
  unreadable$consensus_count[2] <- "lots"
  expect_error(
    build_rearrangement_group_table(unreadable, "WTS_1_S1", "IGK"),
    "consensus_count contains 1 invalid value")
})

# Every downstream stage reads the germline entity rather than re-deriving it, so
# a table arriving without those columns is a broken contract and not something
# to work around.
test_that("records with no germline-entity columns are refused", {
  bare <- LOCUS_RECORDS[, setdiff(names(LOCUS_RECORDS), "v_gene")]
  expect_error(build_rearrangement_group_table(bare, "WTS_1_S1", "IGK"),
               "v_gene")
})
