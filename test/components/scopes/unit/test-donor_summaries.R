# Check donor enrichment, manifest identity and scientific summary populations.
# Run: python test/run_tests.py --component scopes --level unit

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys", "igh_clones",
  "isotype_resolution", "isotype_metrics", "pooled_records", "qc_audit", "shm_summaries",
  "gene_usage", "observed_representatives", "igh_summaries"))

DONOR_TEST_REGIONS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

# Two lineages across two libraries, plus a library with no lineage assignment.
# The second member has partial FWR1 coverage; the fourth has no measured SHM.
donor_summary_inputs <- function() {
  libraries <- c("0460_a", "0460_a", "0460_b", "0460_b")
  sequences <- paste0("s", 1:4)
  records <- data.frame(
    scope_id = "project", donor_id = "0460", library_id = libraries,
    sequence_id = make_pooled_sequence_id(libraries, sequences),
    clone_id = c("project~0460~1", "project~0460~1", "project~0460~1", "project~0460~2"),
    clonotype_key = paste0("key", 1:4),
    v_call = "IGHV1-2*01", j_call = "IGHJ4*01",
    v_gene = "IGHV1-2", j_gene = "IGHJ4",
    junction = "ATG", junction_aa = "M", consensus_count = c(4, 8, 12, 16),
    sequence = "ATG", junction_length = 3,
    sequence_alignment = "ATG", germline_alignment_d_mask = "ATG",
    stringsAsFactors = FALSE
  )
  features <- data.frame(
    library_id = libraries, sequence_id = sequences,
    qualified_c_call = c("IGHM", "IGHG1", "IGHM", "IGHA1"),
    shm_frequency = c(0, 0.1, 0.2, NA), informative_v = c(300, 300, 300, NA)
  )
  for (region in DONOR_TEST_REGIONS) {
    features[[paste0("shm_frequency_", region)]] <- c(0, 0.1, 0.2, NA)
    features[[paste0("mutations_", region)]] <- c(0, 6, 12, NA)
    features[[paste0("informative_", region)]] <- c(60, 60, 60, NA)
    features[[paste0("informative_supported_", region)]] <- c(60, 60, 60, NA)
    features[[paste0("germline_occupied_", region)]] <- c(60, 60, 60, NA)
    features[[paste0("median_majority_votes_", region)]] <- c(4, 5, 6, NA)
  }
  features$informative_supported_fwr1[2] <- 30
  manifest <- data.frame(
    library_id = c("0460_b", "0460_empty", "0460_a", "0500_a"),
    donor_id = c("0460", "0460", "0460", "0500"),
    tissue_id = c("0002", "0003", "0001", "0004")
  )
  list(records = records, membership = records[, CLONE_MEMBER_COLUMNS],
       features = features[4:1, ], manifest = manifest)
}

test_that("one donor enrichment preserves assignment support and library identities", {
  inputs <- donor_summary_inputs()
  joined <- enrich_donor_igh_records(inputs$records, inputs$membership, inputs$features,
                                     inputs$manifest, "project", "0460", DONOR_TEST_REGIONS)
  expect_identical(joined$manifest$library_id, c("0460_a", "0460_b", "0460_empty"))
  expect_identical(joined$records$sequence_id, inputs$membership$sequence_id)
  expect_identical(joined$records$consensus_count, c(4, 8, 12, 16))
  expect_identical(joined$records$tissue_id, c("0001", "0001", "0002", "0002"))
  expect_identical(joined$records$qualified_c_call, c("IGHM", "IGHG1", "IGHM", "IGHA1"))
  expect_equal(joined$records$shm_frequency, c(0, 0.1, 0.2, NA))
})

test_that("the six donor outputs use the same members and existing SHM weighting", {
  inputs <- donor_summary_inputs()
  joined <- enrich_donor_igh_records(inputs$records, inputs$membership, inputs$features,
                                     inputs$manifest, "project", "0460", DONOR_TEST_REGIONS)
  tables <- build_donor_igh_tables(joined$records, joined$manifest, "project", 1, 270,
                                    DONOR_TEST_REGIONS)
  expect_identical(names(tables), c("clone_shm", "shm", "clones", "clone_isotype", "members",
                                    "clone_member_shm", "clone_member_shm_coverage"))
  expect_identical(names(tables$shm), SHM_COLUMNS)
  expect_identical(names(tables$clones), CLONE_FAMILY_COLUMNS)
  expect_identical(names(tables$clone_isotype), LINEAGE_ISOTYPE_COLUMNS)
  expect_identical(names(tables$members), donor_member_export_columns())
  expect_identical(names(tables$clone_member_shm), CLONE_MEMBER_SHM_COLUMNS)
  expect_identical(tables$members$locus, rep("IGH", 4))
  whole <- tables$shm[tables$shm$isotype_group == "all_IGH" &
                       tables$shm$v_region == "whole_v", ]
  expect_identical(whole$library_id, c("0460_a", "0460_b", "0460_empty"))
  expect_equal(whole$mean_clone_median_shm, c(0, 0.2, NA))
  expect_equal(whole$n_members_eligible, c(1, 1, 0))
  expect_equal(whole$n_members_measured, c(2, 1, 0))
  expect_equal(whole$n_clones_eligible, c(1, 1, 0))
  expect_equal(tables$clones$n_members, c(3, 1))
  expect_equal(tables$clones$n_libraries, c(2, 1))
  expect_equal(tables$clones$junction_length, c(3, 3))
  expect_equal(tables$clones$consensus_count, c(24, 16))
  # Clone 1 holds IgM and IgG1 in the first library and IgM in the second, and clone 2
  # holds IgA1 in the second. Each isotype has a row at both levels.
  expect_equal(nrow(tables$clone_isotype), 8)
  expect_identical(tables$clone_isotype$isotype,
                   c("IgM", "IgG", "IgM", "IgG1", "IgM", "IgM", "IgA", "IgA1"))
  expect_equal(tables$clone_isotype$support_fraction, c(1 / 3, 2 / 3, 1 / 3, 2 / 3, 1, 1, 1, 1))
  # The SHM shards list each library's rows by isotype group, then V region, in
  # reporting order, and clone SHM rows by clone within those.
  regions <- shm_v_regions(DONOR_TEST_REGIONS)
  expect_identical(tables$shm$isotype_group[seq_len(length(SHM_ISOTYPE_GROUPS) * 6)],
                   rep(SHM_ISOTYPE_GROUPS, each = 6))
  expect_identical(tables$shm$v_region[1:6], regions)
  library_b <- tables$clone_shm[tables$clone_shm$library_id == "0460_b", ]
  expect_identical(unique(tables$clone_shm$library_id), c("0460_a", "0460_b"))
  expect_identical(unique(library_b$isotype_group), c("all_IGH", "IgM", "IgA", "IgA1"))
  expect_identical(library_b$v_region[1:4], rep(regions[1:2], each = 2))
  expect_identical(library_b$clone_id[1:2], c("project~0460~1", "project~0460~2"))
  expect_identical(tables$members$shm_eligible_whole_v, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(tables$members$shm_frequency_whole_v, c(0, 0.1, 0.2, NA))
  expect_equal(sum(tables$members$consensus_count), sum(tables$clones$consensus_count))
  # The member shard keeps the evidence IDs and names each entity.
  expect_identical(tables$members$v_gene_names, rep("IGHV1-2", 4))
  expect_true("shm_evidence_sequence_id_fwr1" %in% names(tables$members))
  expect_false("informative_supported_fwr1" %in% names(tables$members))
  # A member SHM row exists where the region has a rate, so its rows count
  # n_members_measured, and the shard follows library, clone, sequence and region.
  member_shm <- tables$clone_member_shm
  for (region in regions) {
    measured <- tables$shm[tables$shm$isotype_group == "all_IGH" & tables$shm$v_region == region, ]
    expect_equal(as.integer(table(factor(member_shm$library_id[member_shm$v_region == region],
                                         levels = measured$library_id))),
                 measured$n_members_measured)
  }
  expect_false(any(member_shm$sequence_id == inputs$records$sequence_id[4]))
  expect_identical(order(member_shm$library_id, member_shm$clone_id, member_shm$sequence_id,
                         match(member_shm$v_region, regions)), seq_len(nrow(member_shm)))
  # The second member fails FWR1 coverage with 30 of its 60 positions supported.
  # The shard keeps its measured count and rate, and its position counts explain why.
  partial <- member_shm[member_shm$sequence_id == inputs$records$sequence_id[2] &
                          member_shm$v_region == "whole_v", ]
  expect_false(partial$shm_eligible)
  expect_equal(partial$n_informative_positions, 300)
  expect_equal(partial$n_germline_occupied_positions, 300)
  expect_equal(partial$n_vote_supported_positions, 270)
  expect_equal(partial$n_mutations, 30)
  expect_equal(partial$shm_frequency, 0.1)
  expect_true(is.na(partial$median_trust4_consensus_base_votes))
  expect_true(member_shm$shm_eligible[member_shm$sequence_id == inputs$records$sequence_id[2] &
                                        member_shm$v_region == "cdr2"])
  expect_true("informative_supported_fwr1" %in% names(joined$records))
  relaxed <- build_donor_igh_tables(joined$records, joined$manifest, "project", 0.5, 270,
                                     DONOR_TEST_REGIONS)
  first <- relaxed$shm$library_id == "0460_a" & relaxed$shm$isotype_group == "all_IGH" &
    relaxed$shm$v_region == "whole_v"
  expect_equal(relaxed$shm$mean_clone_median_shm[first], 0.05)
  expect_identical(relaxed$members$shm_eligible_whole_v, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("an empty donor has SHM rows for every library and header-only member tables", {
  inputs <- donor_summary_inputs()
  joined <- enrich_donor_igh_records(inputs$records[FALSE, ], inputs$membership[FALSE, ],
                                     inputs$features[FALSE, ], inputs$manifest,
                                     "project", "0460", DONOR_TEST_REGIONS)
  tables <- build_donor_igh_tables(joined$records, joined$manifest, "project", 1, 270,
                                    DONOR_TEST_REGIONS)
  expect_equal(nrow(tables$shm), 3 * 6 * length(SHM_ISOTYPE_GROUPS))
  expect_true(all(tables$shm$n_members_eligible == 0))
  expect_true(all(is.na(tables$shm$mean_clone_median_shm)))
  expect_equal(nrow(tables$clones), 0)
  expect_equal(nrow(tables$clone_isotype), 0)
  expect_equal(nrow(tables$members), 0)
  expect_identical(names(tables$members), donor_member_export_columns())
  expect_equal(nrow(tables$clone_member_shm), 0)
  expect_identical(names(tables$clone_member_shm), CLONE_MEMBER_SHM_COLUMNS)
})

test_that("donor enrichment rejects mismatched assignments or out-of-scope features", {
  inputs <- donor_summary_inputs()
  changed <- inputs$records
  changed$consensus_count[1] <- changed$consensus_count[1] + 1
  expect_error(enrich_donor_igh_records(changed, inputs$membership, inputs$features,
                                        inputs$manifest, "project", "0460", DONOR_TEST_REGIONS),
                 "disagree in consensus_count")
  expect_error(enrich_donor_igh_records(inputs$records[-1, ], inputs$membership, inputs$features,
                                        inputs$manifest, "project", "0460", DONOR_TEST_REGIONS),
                 "different sequence_id keys")
  expect_error(enrich_donor_igh_records(inputs$records, inputs$membership, inputs$features[-1, ],
                                        inputs$manifest, "project", "0460", DONOR_TEST_REGIONS),
                 "missing from IGH features")
  changed <- inputs$features
  changed$library_id[1] <- "0500_a"
  expect_error(enrich_donor_igh_records(inputs$records, inputs$membership, changed,
                                        inputs$manifest, "project", "0460", DONOR_TEST_REGIONS),
                 "outside the donor")
  expect_error(enrich_donor_igh_records(inputs$records, inputs$membership, inputs$features,
                                        inputs$manifest, "other", "0460", DONOR_TEST_REGIONS),
                 "another scope or donor")
})

test_that("SHM audit counts use the region masks used by the numerical summaries", {
  inputs <- donor_summary_inputs()
  joined <- enrich_donor_igh_records(inputs$records, inputs$membership, inputs$features,
                                     inputs$manifest, "project", "0460", DONOR_TEST_REGIONS)
  tables <- build_donor_igh_tables(joined$records, joined$manifest, "project", 1, 270,
                                    DONOR_TEST_REGIONS, audit_counts = TRUE)
  audit <- attr(tables, "qc_audit")
  rows <- audit[audit$stage == "shm_eligible_members" & audit$isotype_group == "all_IGH" &
                  audit$v_region == "whole_v", ]
  expect_equal(rows$n_units_in, c(2, 2, 0))
  expect_equal(rows$n_units_out, c(1, 1, 0))
  expect_equal(rows$consensus_count_in, c(12, 28, 0))
  expect_equal(rows$consensus_count_out, c(4, 12, 0))
  expect_equal(rows$n_units_out, tables$shm$n_members_eligible[tables$shm$isotype_group == "all_IGH" &
                                                        tables$shm$v_region == "whole_v"])
})

# Manifest joins preserve scientific rows and require a unique crosswalk.

MANIFEST <- data.frame(
  scope_id = rep("luminal", 3),
  library_id = c("S2", "S1", "S3"),
  donor_id = c("D2", "D1", "D1"),
  tissue_id = c("T2", "T1", "T3"),
  stringsAsFactors = FALSE
)

RECORDS <- data.frame(
  scope_id = rep("luminal", 3),
  donor_id = c("D1", "D1", "D2"),
  library_id = c("S1", "S1", "S2"),
  sequence_id = c("S1~seq1", "S1~seq2", "S2~seq3"),
  stringsAsFactors = FALSE
)

describe("a scope whose libraries are all in the manifest", {
  joined <- add_scope_manifest_identifiers(RECORDS, MANIFEST)

  it("cannot change the row count",
     expect_identical(nrow(joined), nrow(RECORDS)))
  it("looks the tissue up rather than assuming manifest order",
     expect_identical(joined$tissue_id, c("T1", "T1", "T2")))
  it("leaves the record columns as they were",
     expect_identical(joined$sequence_id, RECORDS$sequence_id))
})

# A table whose shards carry neither identifier asks for both in one lookup, so
# the donor and the tissue cannot come from two differently-ordered joins.
test_that("both manifest identifiers can be attached at once", {
  bare <- RECORDS[, c("scope_id", "library_id"), drop = FALSE]
  joined <- add_scope_manifest_identifiers(bare, MANIFEST,
                                           columns = c("donor_id", "tissue_id"))
  expect_identical(joined$donor_id, c("D1", "D1", "D2"))
  expect_identical(joined$tissue_id, c("T1", "T1", "T2"))
})

test_that("a manifest lacking a requested column stops the scope", {
  expect_error(
    add_scope_manifest_identifiers(RECORDS, MANIFEST[, c("library_id", "tissue_id")],
                                   columns = c("donor_id", "tissue_id")),
    "donor_id"
  )
})

test_that("no records leaves an empty tissue column", {
  joined <- add_scope_manifest_identifiers(RECORDS[0, , drop = FALSE], MANIFEST)
  expect_identical(nrow(joined), 0L)
  expect_true("tissue_id" %in% names(joined))
})

test_that("a library outside the scope manifest stops the scope", {
  expect_error(
    add_scope_manifest_identifiers(
      rbind(RECORDS, data.frame(scope_id = "luminal", donor_id = "D9",
                                library_id = "S9", sequence_id = "S9~seq1",
                                stringsAsFactors = FALSE)),
      MANIFEST
    ),
    "absent from the scope manifest: S9"
  )
})

test_that("a manifest listing a library twice stops the scope", {
  expect_error(
    add_scope_manifest_identifiers(RECORDS, rbind(MANIFEST, MANIFEST[2, , drop = FALSE])),
    "scope manifest has a duplicate library_id"
  )
})
