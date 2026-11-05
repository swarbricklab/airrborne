# Run regional SHM summaries with source evidence and complete library rows.
# Run: python test/run_tests.py --component shm --level integration

load_test_modules(c(
  "table_io", "gene_calls", "assembly_backbone", "rearrangement_keys", "igh_clones",
  "igh_grouping", "isotype_resolution", "pooled_records", "qc_audit", "shm_summaries",
  "observed_representatives", "igh_summaries", "shm_regions"))

work_dir <- tempfile("scope-shm-")
dir.create(work_dir)

MIN_SUBREGION_FRACTION <- 1.0
MIN_GERMLINE_OCCUPIED_V <- 270L

feature_file <- file.path(work_dir, "features.tsv")
membership_file <- file.path(work_dir, "membership.tsv")
records_file <- file.path(work_dir, "records.tsv")
manifest_file <- file.path(work_dir, "manifest.tsv")
shm_file <- file.path(work_dir, "shm.tsv")
clones_file <- file.path(work_dir, "clones.tsv")
lineage_isotype_file <- file.path(work_dir, "clone_isotype.tsv")
members_file <- file.path(work_dir, "members.tsv")
member_shm_file <- file.path(work_dir, "clone_member_shm.tsv")

features <- data.frame(
  library_id = c(rep("sample_a", 6), "sample_b"),
  sequence_id = paste0("seq", seq_len(7)),
  qualified_c_call = c("IGHG1*01", "IGHG1*01", "IGHA1*01", "IGHG2*01",
             "IGHM*01", "IGHD*01", "IGHG1*01"),
  v_call_status = "resolved",
  v_germline_allele = "IGHV1-1*01",
  shm_frequency = c(0.1, 0.3, 0.4, 0.6, 0.0, 0.8, NA),
  informative_v = c(rep(300L, 6), NA_integer_),
  stringsAsFactors = FALSE
)
# Each record's vote support differs between sub-regions, by offsets that sum to
# +10 rather than cancelling. Giving every sub-region the same number would leave
# the whole-V pooling rule untestable, since the mean, median and occupancy-weighted
# mean of five equal values are the same value.
VOTE_OFFSET <- c(fwr1 = -10, cdr1 = -5, fwr2 = 0, cdr2 = 5, fwr3 = 20)

# Every sub-region carries the record's own rate; what varies is coverage. The sixth
# record reaches only a fraction of FWR1, so it must drop out of FWR1 and out of
# whole_v while remaining in the other four sub-regions. The fifth is unmutated.
for (region in V_REGION_LEVELS) {
  features[[paste0("shm_frequency_", region)]] <- features$shm_frequency
  features[[paste0("mutations_", region)]] <-
    c(6L, 18L, 24L, 36L, 0L, 48L, NA_integer_)
  features[[paste0("informative_", region)]] <- c(rep(60L, 6), NA_integer_)
  features[[paste0("informative_supported_", region)]] <-
    c(rep(60L, 6), NA_integer_)
  features[[paste0("germline_occupied_", region)]] <- 60L
  features[[paste0("median_majority_votes_", region)]] <-
    c(rep(10 + VOTE_OFFSET[[region]], 6), NA_real_)
}
features$informative_fwr1 <- c(60L, 60L, 60L, 60L, 60L, 5L, NA_integer_)
features$informative_supported_fwr1 <- c(60L, 60L, 60L, 60L, 60L, 5L, NA_integer_)

membership <- data.frame(
  scope_id = "scope",
  donor_id = "donor",
  library_id = features$library_id,
  sequence_id = paste(features$library_id, features$sequence_id, sep = "~"),
  clone_id = c("lineage_1", "lineage_1", "lineage_1", "lineage_2",
                 "lineage_2", "lineage_4", "lineage_5"),
  consensus_count = 1,
  v_call = "IGHV1-1*01",
  j_call = "IGHJ4*01",
  # Each call names one gene, and no other gene spells the same germline, so the
  # entity each resolves to is that gene.
  v_gene = "IGHV1-1",
  j_gene = "IGHJ4",
  junction = "TGTGCGTGG",
  junction_aa = "CAW",
  stringsAsFactors = FALSE
)
manifest <- data.frame(
  library_id = c("sample_a", "sample_b"),
  donor_id = "donor",
  tissue_id = c("tissue_a", "tissue_b"),
  stringsAsFactors = FALSE
)

# The donor producer receives the wider projection of these same members.
# Feature values above remain the numerical fixture; the record columns supply
# the sequence evidence carried through to the member export.
records <- membership
records$clonotype_key <- exact_rearrangement_group_key(
  records$junction, records$v_gene, records$j_gene)
records$sequence <- paste0(vapply(seq_len(nrow(records)), function(i)
  paste(rep("A", i), collapse = ""), character(1)), records$junction)
records$junction_length <- nchar(records$junction)
records$sequence_alignment <- records$sequence
records$germline_alignment_d_mask <- records$sequence
records <- records[, CLONE_RECORD_COLUMNS, drop = FALSE]

# Add a second backbone with the same full sequence as seq1. The higher-support
# source remains the sequence representative, but its FWR1 votes do not pass.
# The duplicate supplies whole-V/FWR1 evidence; all numerical SHM expectations
# below remain those of the seven unique sequences. A discarded different-DNA
# sibling has full coverage and a different rate; it must never supply evidence.
source_records <- records
source_records$locus <- "IGH"
source_records$v_gene_status <- source_records$j_gene_status <- "resolved"
source_records$consensus_count[1] <- 5
source_duplicate <- source_records[1, ]
source_duplicate$sequence_id <- "sample_a~seq8"
source_duplicate$consensus_count <- 2
source_sibling <- source_records[1, ]
source_sibling$sequence_id <- "sample_a~seq1_1"
source_sibling$sequence <- paste0("T", source_sibling$sequence)
source_sibling$sequence_alignment <- source_sibling$sequence
source_sibling$consensus_count <- 1
source_records <- rbind(source_records, source_duplicate, source_sibling)
feature_duplicate <- features[1, ]
feature_duplicate$sequence_id <- "seq8"
features$informative_supported_fwr1[1] <- 0L
feature_sibling <- feature_duplicate
feature_sibling$sequence_id <- "seq1_1"
feature_sibling$shm_frequency <- .9
features <- rbind(features, feature_duplicate, feature_sibling)
membership$consensus_count[1] <- records$consensus_count[1] <- 8
mapping <- do.call(rbind, lapply(unique(source_records$library_id), function(library) {
  prepare_igh_grouping(source_records[source_records$library_id == library, ], library, 0)
}))
mapping <- mapping[match(source_records$sequence_id, mapping$sequence_id), ]
pooled_file <- file.path(work_dir, "pooled.tsv")
grouping_file <- file.path(work_dir, "grouping.tsv")
write_pipeline_table(source_records, pooled_file)
write_pipeline_table(mapping, grouping_file)

write.table(features, feature_file, sep = "\t", row.names = FALSE, quote = FALSE,
            na = "")
write.table(membership, membership_file, sep = "\t", row.names = FALSE,
            quote = FALSE)
write.table(records, records_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(manifest, manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)
# No composite entities in this fixture, so the reference names none.
reference_metadata_file <- file.path(work_dir, "reference_metadata.tsv")
writeLines("allele\tidentity_id\tidentity_label", reference_metadata_file)

run <- run_pipeline_script("calculate_donor_igh_summaries.R", c(
  paste0("--qc_audit_file=", tempfile(fileext = ".tsv")),
  paste0("--out_shm=", shm_file),
  paste0("--out_clone_shm=", file.path(work_dir, "clone_shm.tsv")),
  paste0("--out_clones=", clones_file),
  paste0("--out_lineage_isotype=", lineage_isotype_file),
  paste0("--out_members=", members_file),
  paste0("--out_clone_member_shm=", member_shm_file),
  paste0("--out_clone_member_shm_coverage=", file.path(work_dir, "coverage.tsv")),
  "--scope_id=scope",
  "--donor_id=donor",
  paste0("--manifest_file=", manifest_file),
  paste0("--feature_files=", feature_file),
  paste0("--membership_file=", membership_file),
  paste0("--records_file=", records_file),
  paste0("--pooled_file=", pooled_file),
  paste0("--grouping_file=", grouping_file),
  paste0("--min_subregion_fraction=", MIN_SUBREGION_FRACTION),
  paste0("--min_germline_occupied_v=", MIN_GERMLINE_OCCUPIED_V),
  paste0("--reference_metadata_file=", reference_metadata_file)
))
expect_script_ok(run)

shm <- read.delim(shm_file, stringsAsFactors = FALSE, check.names = FALSE)
clone_shm <- read.delim(file.path(work_dir, "clone_shm.tsv"), stringsAsFactors = FALSE,
                        check.names = FALSE)
V_REGIONS <- shm_v_regions(V_REGION_LEVELS)
# Summary rows follow isotype and region order. Clone detail puts library-local
# abundance first and retains that order during scope collection.
test_that("SHM summaries retain strata and clone detail follows local abundance", {
  expected <- expand.grid(v_region = V_REGIONS, isotype_group = SHM_ISOTYPE_GROUPS,
                          library_id = manifest$library_id, stringsAsFactors = FALSE)
  for (column in c("library_id", "isotype_group", "v_region")) {
    expect_identical(shm[[column]], expected[[column]])
  }
  expect_identical(order(clone_shm$library_id, clone_shm$sort_group_support,
                         clone_shm$sort_group_junction, clone_shm$clone_id,
                         match(clone_shm$isotype_group, SHM_ISOTYPE_GROUPS),
                         match(clone_shm$v_region, V_REGIONS)),
                   seq_len(nrow(clone_shm)))
  expect_setequal(unique(clone_shm$isotype_group[clone_shm$library_id == "sample_a"]),
                   c("all_IGH", "IgM", "IgG", "IgA", "IgG1", "IgG2", "IgA1"))
})
test_that("every library gets one row per isotype group per V region", {
  expect_identical(unique(shm$isotype_group), SHM_ISOTYPE_GROUPS)
  expect_identical(unique(shm$v_region), V_REGIONS)
  expect_equal(nrow(shm),
               nrow(manifest) * length(SHM_ISOTYPE_GROUPS) * length(V_REGIONS))
  expect_true(!anyDuplicated(shm[c("library_id", "isotype_group", "v_region")]))
  expect_identical(names(shm), SHM_COLUMNS)
})

row_for <- function(sample, group, v_region = "whole_v") {
  shm[shm$library_id == sample & shm$isotype_group == group &
        shm$v_region == v_region, , drop = FALSE]
}
all_a <- row_for("sample_a", "all_IGH")
igg_a <- row_for("sample_a", "IgG")
iga_a <- row_for("sample_a", "IgA")
igm_a <- row_for("sample_a", "IgM")
# lineage_1 = median(0.1, 0.3, 0.4) = 0.3; lineage_2 = median(0.6, 0.0) = 0.3.
# The sixth record covers a twelfth of FWR1, so it is absent from whole_v.
test_that("each isotype group averages its own lineage medians", {
  expect_equal(all_a$n_members_eligible, 5L)
  expect_equal(all_a$n_clones_eligible, 2L)
  expect_true(isTRUE(all.equal(all_a$mean_clone_median_shm, 0.3)))
  expect_equal(igg_a$n_members_eligible, 3L)
  expect_equal(igg_a$n_clones_eligible, 2L)
  expect_true(isTRUE(all.equal(igg_a$mean_clone_median_shm,
                               mean(c(0.2, 0.6)))))
  expect_equal(iga_a$n_members_eligible, 1L)
  expect_equal(iga_a$n_clones_eligible, 1L)
  expect_true(isTRUE(all.equal(iga_a$mean_clone_median_shm, 0.4)))
  expect_equal(igm_a$n_members_eligible, 1L)
  expect_equal(igm_a$n_clones_eligible, 1L)
  expect_true(isTRUE(all.equal(igm_a$mean_clone_median_shm, 0.0)))
})

# The subclass rows hold only the sequences whose call named one subclass, so they
# need not sum to the class row above them.
test_that("a subclass row summarises the subclass its call named", {
  igg1 <- row_for("sample_a", "IgG1")
  igg2 <- row_for("sample_a", "IgG2")
  iga1 <- row_for("sample_a", "IgA1")
  expect_equal(igg1$n_members_eligible, 2L)
  expect_true(isTRUE(all.equal(igg1$mean_clone_median_shm, 0.2)))
  expect_equal(igg2$n_members_eligible, 1L)
  expect_true(isTRUE(all.equal(igg2$mean_clone_median_shm, 0.6)))
  expect_equal(iga1$n_members_eligible, 1L)
  # IgG4 and IgA2 were never expressed, so their rows exist and hold nothing.
  expect_equal(row_for("sample_a", "IgG4")$n_members_eligible, 0L)
  expect_true(is.na(row_for("sample_a", "IgA2")$mean_clone_median_shm))
})

# The sub-region rows are not copies of the whole-V row: each answers over the
# records that measured that sub-region, which is the point of reporting them, and
# whole_v is the strictest because it is all five at once.
test_that("a sub-region is summarised over the records that covered it", {
  fwr1 <- row_for("sample_a", "all_IGH", "fwr1")
  fwr3 <- row_for("sample_a", "all_IGH", "fwr3")
  expect_equal(fwr1$n_members_eligible, 5L)
  expect_equal(fwr1$n_clones_eligible, 2L)
  expect_true(isTRUE(all.equal(fwr1$mean_clone_median_shm, 0.3)))
  expect_equal(fwr3$n_members_eligible, 6L)
  expect_equal(fwr3$n_clones_eligible, 3L)
  expect_true(isTRUE(all.equal(fwr3$mean_clone_median_shm,
                               mean(c(0.3, 0.3, 0.8)))))
  # Every record had a rate, whichever region admitted it, so the measured count is
  # the same across regions while the summarised count is not.
  expect_equal(fwr1$n_members_measured, 6L)
  expect_equal(fwr3$n_members_measured, 6L)
})

test_that("the released spread and support statistics survive the file round trip", {
  expect_true(isTRUE(all.equal(all_a$q25_clone_median_shm, 0.3)))
  expect_true(isTRUE(all.equal(all_a$q75_clone_median_shm, 0.3)))
  expect_true(isTRUE(all.equal(all_a$iqr_clone_median_shm, 0)))
  # Every record carries 10 votes plus its sub-region's offset, so a record's
  # whole-V value is the mean of those five, 12. A median or a smallest-of-five
  # rule would give 10, and an occupancy-weighted mean over equal-width
  # sub-regions would also give 12, which the unit suite separates on unequal
  # widths.
  expect_true(isTRUE(all.equal(all_a$median_trust4_consensus_base_votes, 12)))
  # A sub-region row reads its own column, so FWR3 sits 20 above the base.
  expect_true(isTRUE(all.equal(
    row_for("sample_a", "all_IGH", "fwr3")$median_trust4_consensus_base_votes, 30)))
  expect_true(isTRUE(all.equal(
    row_for("sample_a", "all_IGH", "fwr1")$median_trust4_consensus_base_votes, 0)))
})

# The unmutated member is the only one in its group, so its lineage's median member
# carries no mutation and the lineage counts as unmutated.
test_that("an unmutated lineage is counted as one", {
  expect_equal(igm_a$n_clones_zero_median_shm, 1L)
  expect_equal(igg_a$n_clones_zero_median_shm, 0L)
  # Mixed lineage_2 holds one mutated and one unmutated member, so its median member
  # is mutated and all_IGH counts no unmutated lineage.
  expect_equal(all_a$n_clones_zero_median_shm, 0L)
})

test_that("a library with no callable record still reports its zeros", {
  empty_rows <- shm[shm$library_id == "sample_b", , drop = FALSE]
  expect_true(all(empty_rows$n_members_eligible == 0L))
  expect_true(all(empty_rows$n_clones_eligible == 0L))
  expect_true(all(is.na(empty_rows$mean_clone_median_shm)))
})

test_that("a lineage holding both unswitched and switched members is flagged", {
  families <- read.delim(clones_file, stringsAsFactors = FALSE,
                         check.names = FALSE)
  expect_true("has_unswitched_and_switched_members" %in% names(families))
  expect_true(families$has_unswitched_and_switched_members[
    families$clone_id == "lineage_2"
  ] == TRUE)
  # Each of the four clones holds members of one library, and every member of
  # this fixture carries the nine-nucleotide junction.
  expect_identical(families$n_libraries, c(1L, 1L, 1L, 1L))
  expect_identical(families$junction_length, rep(9L, 4))
})

# ------------------------- the same lineages as one library saw them ---
#
# The third table the stage writes. lineage_1 holds the transferred support of seq1,
# so its IgG and IgA rows separate dominant from minor support, and lineage_2 is an
# IgM beside an IgG2 in one library.

clone_isotype <- read.delim(lineage_isotype_file, stringsAsFactors = FALSE,
                              check.names = FALSE)

test_that("the lineage shard keeps its payload and collection sort keys", {
  expect_identical(names(clone_isotype), c(LINEAGE_ISOTYPE_COLUMNS, "sort_group_support", "sort_group_junction"))
})

test_that("the rows land in the order scope collection merges", {
  expect_identical(nrow(clone_isotype), 12L)
  expect_identical(clone_isotype$clone_id,
                   rep(c("lineage_1", "lineage_2", "lineage_4", "lineage_5"), c(4, 4, 2, 2)))
  expect_identical(clone_isotype$library_id, c(rep("sample_a", 10), rep("sample_b", 2)))
  expect_identical(clone_isotype$isotype_level[1:4], c("class", "class", "subclass", "subclass"))
  expect_identical(clone_isotype$isotype[1:8],
                   c("IgG", "IgA", "IgG1", "IgA1", "IgM", "IgG", "IgM", "IgG2"))
})

test_that("each isotype row carries its members, support and callable denominator", {
  isotype_row <- function(clone_id, level, isotype) {
    clone_isotype[clone_isotype$clone_id == clone_id & clone_isotype$isotype_level == level &
                    clone_isotype$isotype == isotype, , drop = FALSE]
  }
  igg <- isotype_row("lineage_1", "class", "IgG")
  expect_equal(igg$n_members, 2L)
  expect_equal(igg$consensus_count, 9)
  expect_equal(igg$callable_consensus_count, 10)
  expect_true(isTRUE(all.equal(igg$support_fraction, 0.9)))
  expect_true(isTRUE(all.equal(isotype_row("lineage_1", "class", "IgA")$support_fraction, 0.1)))
  expect_true(isTRUE(all.equal(isotype_row("lineage_2", "subclass", "IgG2")$support_fraction, 0.5)))
  expect_equal(isotype_row("lineage_4", "subclass", "IgD")$n_members, 1L)
  expect_identical(clone_isotype$tissue_id[clone_isotype$library_id == "sample_b"], rep("tissue_b", 2))
})

test_that("post-merge SHM counts the sequence once and exports regional source IDs", {
  released <- read_pipeline_table(members_file)
  member <- released[released$sequence_id == "sample_a~seq1", ]
  expect_equal(nrow(released), 7L)
  expect_equal(as.numeric(member$consensus_count), 8)
  expect_identical(member$shm_evidence_sequence_id_whole_v, "sample_a~seq8")
  expect_identical(member$shm_evidence_sequence_id_fwr1, "sample_a~seq8")
  expect_identical(member$shm_evidence_sequence_id_fwr3, "sample_a~seq1")
  expect_identical(as.numeric(member$shm_frequency_whole_v), .1)
})

# seq1's FWR1 votes do not pass, so whole V and FWR1 read seq8, the identical source.
# seq6 reaches 5 of FWR1's 60 positions and seq7 has no rate at all.
test_that("the member SHM shard lands in key order with the selected source's counts", {
  rows <- read_pipeline_table(member_shm_file)
  expect_identical(names(rows), c(CLONE_MEMBER_SHM_COLUMNS, "sort_group_support", "sort_group_junction"))
  expect_equal(nrow(rows), 6L * 6L)
  expect_false(any(rows$sequence_id == "sample_b~seq7"))
  expect_identical(order(rows$library_id, as.numeric(rows$sort_group_support),
                         rows$sort_group_junction, rows$clone_id, rows$sequence_id,
                         match(rows$v_region, V_REGIONS)), seq_len(nrow(rows)))
  expect_equal(sum(rows$v_region == "fwr1"), row_for("sample_a", "all_IGH", "fwr1")$n_members_measured)
  cell <- function(sequence_id, region) {
    rows[rows$sequence_id == sequence_id & rows$v_region == region, , drop = FALSE]
  }
  expect_equal(as.numeric(cell("sample_a~seq1", "whole_v")$n_vote_supported_positions), 300)
  expect_equal(as.numeric(cell("sample_a~seq1", "fwr1")$n_vote_supported_positions), 60)
  expect_identical(as.character(cell("sample_a~seq1", "whole_v")$shm_eligible), "TRUE")
  ineligible <- cell("sample_a~seq6", "whole_v")
  expect_identical(as.character(ineligible$shm_eligible), "FALSE")
  expect_equal(as.numeric(ineligible$n_informative_positions), 300)
  expect_equal(as.numeric(ineligible$n_germline_occupied_positions), 300)
  expect_equal(as.numeric(ineligible$n_vote_supported_positions), 245)
  expect_equal(as.numeric(ineligible$n_mutations), 240)
  expect_equal(as.numeric(ineligible$shm_frequency), 0.8)
  expect_true(is.na(ineligible$median_trust4_consensus_base_votes) || ineligible$median_trust4_consensus_base_votes == "")
})
