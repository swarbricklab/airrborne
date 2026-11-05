# Check clone mutation states, missing evidence and library-specific medians.
# Run: python test/run_tests.py --component shm --level unit

load_test_modules(c(
  "table_io", "gene_calls", "isotype_resolution", "assembly_backbone", "qc_audit",
  "shm_summaries"))

CLONE_SHM_TEST_REGIONS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

# The first clone has counts [0, 0, 8]. Other clones cover every observed state,
# partial coverage, ambiguous class/subclass evidence and a second library.
clone_shm_records <- function() {
  mutations <- c(0, 0, 8, 8, 8, 0, 0, NA, 8, 0, 8)
  measured <- is.finite(mutations)
  records <- data.frame(
    scope_id = "project", donor_id = "0460",
    tissue_id = c(rep("0001", 10), "0002"),
    library_id = c(rep("0460_a", 10), "0460_b"),
    sequence_id = paste0("s", seq_along(mutations)),
    clone_id = c(rep("mixed_clone", 3), rep("positive_clone", 2),
                   rep("zero_clone", 2), "missing_clone", "partial_clone",
                   "uncalled_clone", "mixed_clone"),
    qualified_c_call = c("IGHM", "IGHM", "IGHG1", "IGHG1", "IGHG1", "IGHA1", "IGHA1",
               "IGHG1", "IGHG1,IGHG3", "", "IGHM"),
    consensus_count = c(rep(1, 7), 1000000, 10, 5, 2),
    shm_frequency = mutations / 300,
    informative_v = ifelse(measured, 300, NA_real_)
  )
  for (region in CLONE_SHM_TEST_REGIONS) {
    counts <- if (region == "cdr2") mutations else ifelse(measured, 0, NA_real_)
    records[[paste0("mutations_", region)]] <- counts
    records[[paste0("shm_frequency_", region)]] <- counts / 60
    records[[paste0("informative_", region)]] <- ifelse(measured, 60, NA_real_)
    records[[paste0("informative_supported_", region)]] <- ifelse(measured, 60, NA_real_)
    records[[paste0("germline_occupied_", region)]] <- ifelse(measured, 60, NA_real_)
    records[[paste0("median_majority_votes_", region)]] <- ifelse(measured, 10, NA_real_)
  }
  records$informative_supported_fwr1[9] <- 30
  records
}

test_that("a median-zero mixture stays distinct from all-zero members", {
  result <- build_clone_shm_table(clone_shm_records(), 1, 270, CLONE_SHM_TEST_REGIONS)
  expect_identical(names(result), CLONE_SHM_COLUMNS)
  expect_identical(names(result)[6:7], c("isotype_group", "v_region"))
  # Rows follow library, then isotype group and V region in reporting order,
  # then clone.
  expect_identical(order(result$library_id, match(result$isotype_group, SHM_ISOTYPE_GROUPS),
                         match(result$v_region, shm_v_regions(CLONE_SHM_TEST_REGIONS)),
                         result$clone_id),
                   seq_len(nrow(result)))
  expect_identical(result$isotype_group[1:2], c("all_IGH", "all_IGH"))
  expect_identical(result$v_region[1:2], c("whole_v", "whole_v"))
  selected <- result$library_id == "0460_a" & result$clone_id == "mixed_clone" &
    result$v_region == "whole_v"
  mixed <- result[selected & result$isotype_group == "all_IGH", ]
  expect_equal(mixed$consensus_count, 3)
  expect_equal(mixed$eligible_consensus_count, 3)
  expect_equal(mixed$n_members, 3)
  expect_equal(mixed$n_members_eligible, 3)
  expect_equal(mixed$n_members_zero_shm, 2)
  expect_equal(mixed$n_members_positive_shm, 1)
  expect_equal(mixed$median_shm_frequency, 0)
  expect_equal(mixed$mutation_state, "mixed")
  igm <- result[selected & result$isotype_group == "IgM", ]
  igg <- result[selected & result$isotype_group == "IgG", ]
  expect_equal(igm$mutation_state, "all_zero")
  expect_equal(igm$consensus_count, 2)
  expect_equal(igg$mutation_state, "all_positive")
  expect_equal(igg$median_shm_frequency, 8 / 300)
  other_library <- result[result$library_id == "0460_b" &
                            result$clone_id == "mixed_clone" &
                            result$v_region == "whole_v" &
                            result$isotype_group == "all_IGH", ]
  expect_equal(other_library$mutation_state, "all_positive")
  expect_identical(other_library$tissue_id, "0002")
  expect_identical(other_library$donor_id, "0460")
})

# The donor shard keeps these rows; scope collection leaves them out of the
# released clone_shm.tsv, and group_abundance.tsv keeps the clone's support.
test_that("dominant clones without evaluable SHM retain support and member counts", {
  result <- build_clone_shm_table(clone_shm_records(), 1, 270, CLONE_SHM_TEST_REGIONS)
  missing <- result[result$clone_id == "missing_clone" & result$v_region == "whole_v" &
                      result$isotype_group == "all_IGH", ]
  expect_equal(missing$consensus_count, 1000000)
  expect_equal(missing$eligible_consensus_count, 0)
  expect_equal(missing$n_members, 1)
  expect_equal(missing$n_members_eligible, 0)
  expect_equal(missing$n_members_zero_shm, 0)
  expect_equal(missing$n_members_positive_shm, 0)
  expect_true(is.na(missing$median_shm_frequency))
  expect_equal(missing$mutation_state, "no_evaluable_members")
  expect_equal(missing$consensus_count, max(result$consensus_count))
  expect_true(all(result$n_members_zero_shm + result$n_members_positive_shm ==
                    result$n_members_eligible))
})

test_that("coverage and class ambiguity follow the existing SHM populations", {
  records <- clone_shm_records()
  result <- build_clone_shm_table(records, 1, 270, CLONE_SHM_TEST_REGIONS)
  partial <- result[result$clone_id == "partial_clone", ]
  expect_setequal(partial$isotype_group, c("all_IGH", "IgG"))
  whole <- partial[partial$v_region == "whole_v" & partial$isotype_group == "IgG", ]
  cdr2 <- partial[partial$v_region == "cdr2" & partial$isotype_group == "IgG", ]
  expect_equal(whole$consensus_count, 10)
  expect_equal(whole$eligible_consensus_count, 0)
  expect_equal(whole$n_members, 1)
  expect_equal(whole$mutation_state, "no_evaluable_members")
  expect_equal(cdr2$eligible_consensus_count, 10)
  expect_equal(cdr2$mutation_state, "all_positive")
  expect_equal(cdr2$median_shm_frequency, 8 / 60)
  uncalled <- result[result$clone_id == "uncalled_clone", ]
  expect_setequal(uncalled$isotype_group, "all_IGH")
  expect_true(all(uncalled$mutation_state == "all_zero"))
  relaxed <- build_clone_shm_table(records, 0.5, 270, CLONE_SHM_TEST_REGIONS)
  selected <- relaxed$clone_id == "partial_clone" & relaxed$v_region == "whole_v" &
    relaxed$isotype_group == "all_IGH"
  expect_equal(relaxed$mutation_state[selected], "all_positive")
})

# Two members of one clone differ in eligibility, so eligible support is a subset
# of the clone's support in the region that separates them.
test_that("eligible support counts only the members eligible in the region", {
  records <- clone_shm_records()
  records$clone_id[9] <- "positive_clone"
  records$qualified_c_call[9] <- "IGHG1"
  result <- build_clone_shm_table(records, 1, 270, CLONE_SHM_TEST_REGIONS)
  row_for <- function(region) {
    result[result$clone_id == "positive_clone" & result$isotype_group == "all_IGH" &
             result$v_region == region, ]
  }
  expect_equal(row_for("whole_v")$consensus_count, 12)
  expect_equal(row_for("whole_v")$eligible_consensus_count, 2)
  expect_equal(row_for("fwr1")$eligible_consensus_count, 2)
  expect_equal(row_for("cdr2")$eligible_consensus_count, 12)
})

test_that("clone medians reproduce the existing equally weighted library summaries", {
  records <- clone_shm_records()
  detail <- build_clone_shm_table(records, 1, 270, CLONE_SHM_TEST_REGIONS)
  manifest <- unique(records[, c("library_id", "donor_id", "tissue_id")])
  summary <- build_scope_shm_table(records, manifest, "project",
                                   shm_v_regions(CLONE_SHM_TEST_REGIONS), 1, 270,
                                   CLONE_SHM_TEST_REGIONS)
  for (index in seq_len(nrow(summary))) {
    row <- summary[index, ]
    clones <- detail[detail$library_id == row$library_id &
                       detail$isotype_group == row$isotype_group &
                       detail$v_region == row$v_region &
                       detail$n_members_eligible > 0, ]
    expect_equal(nrow(clones), row$n_clones_eligible)
    expect_equal(sum(clones$n_members_eligible), row$n_members_eligible)
    expect_equal(sum(clones$median_shm_frequency == 0), row$n_clones_zero_median_shm)
    expected <- if (nrow(clones)) mean(clones$median_shm_frequency) else NA_real_
    expect_equal(expected, row$mean_clone_median_shm, tolerance = 1e-12)
  }
  whole <- detail[detail$library_id == "0460_a" & detail$v_region == "whole_v" &
                    detail$isotype_group == "all_IGH" & detail$n_members_eligible > 0, ]
  expect_equal(mean(whole$median_shm_frequency), 1 / 150, tolerance = 1e-12)
})

test_that("empty and absent isotype subsets yield a complete header without invented clones", {
  records <- clone_shm_records()
  empty <- build_clone_shm_table(records[FALSE, ], 1, 270, CLONE_SHM_TEST_REGIONS)
  expect_identical(names(empty), CLONE_SHM_COLUMNS)
  expect_equal(nrow(empty), 0)
  absent <- build_clone_shm_table(records, 1, 270, CLONE_SHM_TEST_REGIONS,
                                   isotype_groups = "IgA2")
  expect_identical(names(absent), CLONE_SHM_COLUMNS)
  expect_equal(nrow(absent), 0)
})
