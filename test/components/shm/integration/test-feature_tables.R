# Run IGH feature tables and compare regional mutation rates with pinned SHazaM.
# Run: python test/run_tests.py --component shm --level integration

load_test_modules(c(
  "table_io", "gene_calls", "shm_regions", "assembly_backbone",
  "trust4_consensus_votes", "igh_features", "isotype_calls", "qc_audit"))

work_dir <- tempfile("igh-features-")
dir.create(work_dir)

fixture <- shm_feature_fixture()
run_features <- function(rows, label, votes = fixture$vote_rows(rows$sequence_id)) {
  input_file <- file.path(work_dir, paste0(label, "_airr.tsv"))
  votes_file <- file.path(work_dir, paste0(label, "_votes.tsv"))
  calls_file <- file.path(work_dir, paste0(label, "_isotype_calls.tsv"))
  output_file <- file.path(work_dir, paste0(label, "_features.tsv"))
  calls <- build_qualified_isotype_calls(rows, "GUSS_TEST", fixture$MIN_C_POSITIONS)
  for (pair in list(list(rows, input_file), list(votes, votes_file),
                    list(calls, calls_file))) {
    write.table(pair[[1]], pair[[2]], sep = "\t", row.names = FALSE,
                quote = FALSE, na = "")
  }
  run <- run_pipeline_script("calculate_igh_features.R", c(
    paste0("--qc_audit_file=", file.path(work_dir, paste0(label, "_audit.tsv"))),
    paste0("--input_file=", input_file),
    paste0("--votes_file=", votes_file),
    paste0("--isotype_calls_file=", calls_file),
    paste0("--min_position_votes=", fixture$MIN_POSITION_VOTES),
    paste0("--output_file=", output_file),
    "--library_id=GUSS_TEST"
  ))
  expect_script_ok(run)
  read.delim(output_file, stringsAsFactors = FALSE, check.names = FALSE)
}

# --- A full-length record: the sub-regions partition the released rate -------
# One mutation each in FWR1, CDR1 and FWR3, none in FWR2 or CDR2.

full_sequence <- fixture$mutate_positions(fixture$GERMLINE, c(10L, 100L, 250L))
full <- run_features(fixture$airr_row("assemble0_0", full_sequence, fixture$GERMLINE), "full")

test_that("SHM audit separates measured zero mutation from zero informative positions", {
  partial <- fixture$GERMLINE
  substr(partial, 79, 114) <- strrep(".", 36)
  rows <- rbind(fixture$airr_row("unmutated", fixture$GERMLINE, fixture$GERMLINE),
                fixture$airr_row("unobserved", strrep("N", 312), fixture$GERMLINE),
                fixture$airr_row("missing_germline", fixture$GERMLINE, ""),
                fixture$airr_row("partial", partial, fixture$GERMLINE))
  result <- run_features(rows, "measured_zero")
  result <- result[match(rows$sequence_id, result$sequence_id), ]
  expect_equal(result$shm_frequency[1], 0)
  expect_equal(result$informative_v[2], 0)
  expect_true(is.na(result$shm_frequency[2]))
  expect_true(is.na(result$informative_v[3]))
  expect_equal(result$informative_cdr1[4], 0)
  expect_equal(result$informative_v[4], 276)
  audit <- read_qc_audit(file.path(work_dir, "measured_zero_audit.tsv"))
  measured <- audit[audit$stage == "shm_measured_records", ]
  expect_true(all(measured$n_units_in == 4))
  expect_equal(measured$n_units_out, ifelse(measured$v_region == "cdr1", 1, 2))
  expect_true(all(measured$consensus_count_in == 48))
  expect_equal(measured$consensus_count_out, ifelse(measured$v_region == "cdr1", 12, 24))
})

test_that("a full-length record releases the fixed column set", {
  expect_identical(names(full), fixture$EXPECTED_COLUMNS)
  expect_identical(ncol(full), 37L)
  expect_identical(nrow(full), 1L)
  expect_identical(full$qualified_c_call, "IGHM*01")
})

# Every position of this germline is occupied, and the record covers all of them.
test_that("a full-length record covers every occupied germline position", {
  expect_identical(fixture$region_values(full, "germline_occupied_"),
                   unname(fixture$REGION_WIDTH[V_REGION_LEVELS]))
  expect_identical(fixture$region_values(full, "informative_"),
                   unname(fixture$REGION_WIDTH[V_REGION_LEVELS]))
  expect_identical(full$informative_v, 312L)
  expect_identical(full$informative_v, sum(fixture$region_values(full, "informative_")))
})

test_that("the sub-regions partition the released rate", {
  expect_true(isTRUE(all.equal(full$shm_frequency_fwr1, 1 / 78)))
  expect_true(isTRUE(all.equal(full$shm_frequency_cdr1, 1 / 36)))
  expect_true(isTRUE(all.equal(full$shm_frequency_fwr2, 0)))
  expect_true(isTRUE(all.equal(full$shm_frequency_cdr2, 0)))
  expect_true(isTRUE(all.equal(full$shm_frequency_fwr3, 1 / 117)))
  expect_true(isTRUE(all.equal(full$shm_frequency, 3 / 312)))
})

# The numerator is released so a reader can re-aggregate in integers rather than
# multiplying a rate back out, and the whole-V numerator is its exact sum.
test_that("each rate ships the mutation count it was divided from", {
  expect_identical(fixture$region_values(full, "mutations_"), c(1L, 1L, 0L, 0L, 1L))
  expect_identical(sum(fixture$region_values(full, "mutations_")), 3L)
})

test_that("the germline the rate was measured against is named", {
  expect_identical(full$v_call_status, "resolved")
  expect_identical(full$v_germline_allele, "IGHV1-8*01")
})

# The contract: the released whole-V number does not move.
test_that("the released whole-V number does not move", {
  expect_true(isTRUE(all.equal(
    full$shm_frequency, fixture$shazam_whole_v_mu_freq(full_sequence, fixture$GERMLINE),
    tolerance = fixture$SERIALISATION_TOLERANCE
  )))
})

# --- Read support is taken over the positions the rate was taken over --------

test_that("each sub-region reports the read support behind its own positions", {
  expect_identical(as.numeric(full[paste0("median_majority_votes_",
                                          V_REGION_LEVELS)][1, ]),
                   as.numeric(fixture$REGION_MAJORITY_VOTES[V_REGION_LEVELS]))
})

# A position with no mapped coordinate, and one with a total but no base to agree
# with, both leave the median with one fewer value rather than pulling it towards
# zero. Counting them as zero is the mistake this column exists to make visible.
test_that("an absent vote is not a zero vote", {
  expect_identical(fixture$region_values(full, "informative_supported_"),
                   unname(fixture$VOTED_WIDTH[V_REGION_LEVELS]))
  expect_identical(full$informative_supported_cdr2, 20L)
  expect_identical(full$informative_cdr2, 30L)
  # The median is unmoved: the ten absent CDR2 positions are dropped, not scored.
  # A whole number reads back from the text table as an integer, so compare values.
  expect_equal(as.numeric(full$median_majority_votes_cdr2),
               as.numeric(fixture$REGION_MAJORITY_VOTES[["cdr2"]]))
})

# The support count is the one the released eligibility rule reads, so it has to
# separate a position resting on a single agreeing read from one resting on
# several. The rate and the median are deliberately left alone: this column
# changes which sequences a summary admits, not what their mutation load is.
thin <- run_features(
  fixture$airr_row("assemble6_0", full_sequence, fixture$GERMLINE), "thin",
  votes = fixture$vote_rows("assemble6_0", fixture$THIN_VOTES$total, fixture$THIN_VOTES$non_consensus)
)

test_that("a position resting on one vote is measured but not read-supported", {
  # FWR1 still contributed every one of its positions to the rate ...
  expect_identical(thin$informative_fwr1, full$informative_fwr1)
  expect_identical(thin$shm_frequency_fwr1, full$shm_frequency_fwr1)
  # ... and none of them is supported, because the majority there is exactly one.
  expect_identical(thin$informative_supported_fwr1, 0L)
  expect_equal(as.numeric(thin$median_majority_votes_fwr1), 1)
  # Every other sub-region is untouched, so the threshold is per position and not
  # a property of the record.
  expect_identical(fixture$region_values(thin, "informative_supported_")[-1],
                   unname(fixture$VOTED_WIDTH[V_REGION_LEVELS])[-1])
})

# The three counts nest by construction, which is what lets one comparison against
# germline_occupied assert both halves of the coverage rule at once.
test_that("supported is a subset of informative is a subset of occupied", {
  for (record in list(full, thin)) {
    expect_true(all(fixture$region_values(record, "informative_supported_") <=
                      fixture$region_values(record, "informative_")))
    expect_true(all(fixture$region_values(record, "informative_") <=
                      fixture$region_values(record, "germline_occupied_")))
  }
})

# A record the vote stage emptied because it was assembled in the opposite
# orientation to its contig still has a mutation rate; what it has no answer for
# is how well voted those positions were.
emptied <- run_features(
  fixture$airr_row("assemble5_0", full_sequence, fixture$GERMLINE), "emptied",
  votes = fixture$vote_rows("assemble5_0", fixture$EMPTY_VOTES, fixture$EMPTY_VOTES)
)
test_that("a record with no votes at all keeps its rate and loses its support", {
  expect_true(isTRUE(all.equal(emptied$shm_frequency, 3 / 312)))
  expect_identical(fixture$region_values(emptied, "informative_supported_"), rep(0L, 5L))
  expect_true(all(is.na(as.numeric(
    emptied[paste0("median_majority_votes_", V_REGION_LEVELS)][1, ]))))
})

# --- A five-prime truncated record reports what it actually covered ----------
# The germline still occupies every position; the record does not.

truncated_sequence <- fixture$blank_positions(
  fixture$mutate_positions(fixture$GERMLINE, 250L), 1L, 200L
)
truncated <- run_features(
  fixture$airr_row("assemble1_0", truncated_sequence, fixture$GERMLINE), "truncated"
)

test_that("a truncated record reports what it actually covered", {
  expect_identical(fixture$region_values(truncated, "informative_"),
                   c(0L, 0L, 0L, 0L, 112L))
  expect_identical(truncated$informative_v, 112L)
  # Not covered is not the same as not there: the germline occupancy is unchanged.
  expect_identical(fixture$region_values(truncated, "germline_occupied_"),
                   unname(fixture$REGION_WIDTH[V_REGION_LEVELS]))
  # A sub-region with no contributing position has no rate, rather than a zero.
  expect_true(is.na(truncated$shm_frequency_fwr1))
  expect_true(isTRUE(all.equal(truncated$shm_frequency_fwr3, 1 / 112)))
})

# Feature measurement reports the rate and its measured denominator.
# SHM summaries apply the regional coverage requirements separately.
test_that("a thinly covered record still receives its whole-V rate", {
  expect_true(isTRUE(all.equal(truncated$shm_frequency, 1 / 112)))
  expect_true(isTRUE(all.equal(truncated$shm_frequency,
                               fixture$shazam_whole_v_mu_freq(truncated_sequence, fixture$GERMLINE),
                               tolerance = fixture$SERIALISATION_TOLERANCE)))
})

# Votes exist at every FWR1 position, but the record measured none of them, so
# there is nothing for a support median to describe.
test_that("read support is not reported for a sub-region that was not measured", {
  expect_identical(fixture$region_values(truncated, "informative_supported_"),
                   c(0L, 0L, 0L, 0L, 112L))
  expect_true(is.na(truncated$median_majority_votes_fwr1))
  expect_equal(as.numeric(truncated$median_majority_votes_fwr3),
               as.numeric(fixture$REGION_MAJORITY_VOTES[["fwr3"]]))
})

# --- IMGT spacers are not coverage a record is missing ----------------------
# A germline whose real CDR1 is shorter than the 36 columns IMGT reserves leaves
# spacers there. They must reduce the sub-region's denominator, not appear as
# positions the record failed to cover.

SPACER_FROM <- 79L
SPACER_TO <- 90L
spacer_germline <- fixture$blank_positions(fixture$GERMLINE, SPACER_FROM, SPACER_TO)
spacer_sequence <- fixture$blank_positions(
  fixture$mutate_positions(fixture$GERMLINE, 100L), SPACER_FROM, SPACER_TO
)
spacer <- run_features(
  fixture$airr_row("assemble2_0", spacer_sequence, spacer_germline), "spacer"
)

n_spacers <- SPACER_TO - SPACER_FROM + 1L
test_that("IMGT spacers are not coverage a record is missing", {
  expect_identical(spacer$germline_occupied_cdr1, fixture$REGION_WIDTH[["cdr1"]] - n_spacers)
  expect_identical(spacer$informative_cdr1, spacer$germline_occupied_cdr1)
  # Only CDR1 was shortened.
  expect_identical(spacer$germline_occupied_fwr1, fixture$REGION_WIDTH[["fwr1"]])
  expect_identical(spacer$informative_v, 312L - n_spacers)
  # The rate divides by the positions that exist, not by the reserved width.
  expect_true(isTRUE(all.equal(spacer$shm_frequency_cdr1, 1 / spacer$germline_occupied_cdr1)))
  expect_true(isTRUE(all.equal(spacer$shm_frequency,
                               fixture$shazam_whole_v_mu_freq(spacer_sequence, spacer_germline),
                               tolerance = fixture$SERIALISATION_TOLERANCE)))
})

# A spacer carries a vote count, because the position exists on the contig, but it
# is not a mutation opportunity, so it must not enter the support median either.
test_that("a spacer is not a voted position", {
  expect_identical(spacer$informative_supported_cdr1, spacer$germline_occupied_cdr1)
})

# --- A library with no IGH records still writes the full header -------------

non_igh <- fixture$airr_row("assemble3_0", full_sequence, fixture$GERMLINE)
non_igh$locus <- "TRB"
empty <- run_features(non_igh, "empty")
test_that("a library with no IGH records still writes the full header", {
  expect_identical(names(empty), fixture$EXPECTED_COLUMNS)
  expect_identical(nrow(empty), 0L)
})
