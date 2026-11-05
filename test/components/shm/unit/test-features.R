# Check regional SHM mutations, informative positions and read-supported evidence.
# Run: python test/run_tests.py --component shm --level unit

load_test_modules(c(
  "table_io", "gene_calls", "shm_regions", "assembly_backbone",
  "trust4_consensus_votes", "igh_features"))

fixture <- shm_feature_fixture()
full_sequence <- fixture$mutate_positions(fixture$GERMLINE, c(10L, 100L, 250L))
truncated_sequence <- fixture$blank_positions(
  fixture$mutate_positions(fixture$GERMLINE, 250L), 1L, 200L)
full <- build_igh_v_region_feature_table(
  fixture$airr_row("assemble0_0", full_sequence, fixture$GERMLINE),
  fixture$vote_rows("assemble0_0"), fixture$MIN_POSITION_VOTES)

# Two records the module scores in one call: one covering the whole window and
# one truncated to FWR3, so every per-region diagnostic has two different values
# to summarise. The second also carries a tied V call.
DIRECT_RECORDS <- rbind(
  fixture$airr_row("assemble0_0", full_sequence, fixture$GERMLINE),
  fixture$airr_row("assemble1_0", truncated_sequence, fixture$GERMLINE,
           v_call = "IGHV1-8*01,IGHV1-8*03")
)
DIRECT_VOTES <- fixture$vote_rows(DIRECT_RECORDS$sequence_id)
direct <- build_igh_v_region_feature_table(DIRECT_RECORDS, DIRECT_VOTES,
                                           fixture$MIN_POSITION_VOTES)

test_that("the module builds the released columns in the released order", {
  expect_identical(names(direct), fixture$EXPECTED_COLUMNS)
  expect_identical(nrow(direct), 2L)
})

# The contract the sub-region columns exist to keep, asserted without going
# through a text file: the whole-V rate is still SHazaM's own number.
test_that("the whole-V rate is SHazaM's, before any serialisation", {
  expect_identical(direct$shm_frequency[1],
                   fixture$shazam_whole_v_mu_freq(full_sequence, fixture$GERMLINE))
})

test_that("every record with a comparison to make receives a rate", {
  expect_identical(direct$informative_v, c(312L, 112L))
  expect_false(anyNA(direct$shm_frequency))
})

# Two alleles of one gene are a resolved gene call and an unmade germline choice,
# which is the distinction this column exists to carry.
test_that("a tied V call is reported as the ambiguity it is", {
  expect_identical(direct$v_call_status,
                   c("resolved", "ambiguous_multi_allele"))
  expect_identical(direct$v_germline_allele, rep("IGHV1-8*01", 2L))
})

# A record with no germline to compare against is not scored at all, which is
# distinct from being scored and found to cover little.
test_that("a record with no alignment pair is left uncounted", {
  no_germline <- fixture$airr_row("assemble4_0", full_sequence, fixture$GERMLINE)
  no_germline$germline_alignment_d_mask <- ""
  unscored <- build_igh_v_region_feature_table(no_germline,
                                               fixture$vote_rows("assemble4_0"), fixture$MIN_POSITION_VOTES)
  expect_true(is.na(unscored$informative_v))
  expect_true(is.na(unscored$shm_frequency))
  expect_true(is.na(unscored$informative_fwr3))
  # Nothing measured is not the same as measuring nothing, so the voted count is
  # absent here where an emptied vote vector would have made it zero.
  expect_true(is.na(unscored$informative_supported_fwr3))
  expect_true(is.na(unscored$median_majority_votes_fwr3))
})

# Both tables are built from the same filtered AIRR file, so a record without a
# vote row means they describe different assemblies. That pairing runs clean and
# produces plausible-looking numbers, so it has to stop here.
test_that("a record with no vote row stops the library", {
  expect_error(
    build_igh_v_region_feature_table(DIRECT_RECORDS,
                                     DIRECT_VOTES[1, , drop = FALSE], fixture$MIN_POSITION_VOTES),
    "1 IGH record(s) have no consensus-vote row", fixed = TRUE
  )
  expect_error(
    build_igh_v_region_feature_table(DIRECT_RECORDS,
                                     rbind(DIRECT_VOTES, DIRECT_VOTES[1, ]), fixture$MIN_POSITION_VOTES),
    "sequence_id is not unique in the IGH consensus-vote table", fixed = TRUE
  )
})

describe("the per-region callability diagnostics", {
  callability <- summarise_igh_feature_callability(direct)

  it("reports one row per sub-region, 5' to 3'",
     expect_identical(callability$region, V_REGION_LEVELS))
  # FWR1: 78 positions in one record and none in the other; FWR3: 117 and 112.
  it("takes the median of what was actually covered", {
    expect_equal(callability$median_informative,
                 c(39, 18, 25.5, 15, 114.5))
  })
  # Occupancy is a property of the germline, so both records report the same one
  # and the median is that value rather than something between two.
  it("reports germline occupancy beside it", {
    expect_equal(callability$median_germline_occupied,
                 unname(as.numeric(fixture$REGION_WIDTH[V_REGION_LEVELS])))
  })
  # Both records see the same votes, so the median of the medians is that value;
  # the truncated record contributes nothing where it measured nothing.
  it("reports the read support behind the measured positions", {
    expect_equal(callability$median_trust4_consensus_base_votes,
                 as.numeric(fixture$REGION_MAJORITY_VOTES[V_REGION_LEVELS]))
  })
  it("counts the records given a rate in each sub-region",
     expect_identical(callability$n_rate_called, c(1L, 1L, 1L, 1L, 2L)))
})

test_that("a library with no records still reports every sub-region", {
  callability <- summarise_igh_feature_callability(
    build_igh_v_region_feature_table(DIRECT_RECORDS[0, , drop = FALSE],
                                     DIRECT_VOTES[0, , drop = FALSE], fixture$MIN_POSITION_VOTES)
  )
  expect_identical(callability$region, V_REGION_LEVELS)
  expect_true(all(is.na(callability$median_informative)))
  expect_identical(callability$n_rate_called, rep(0L, 5L))
})

# --- the shared per-position definitions --------------------------------------
#
# The identity that matters is that the mask counts what SHazaM counts.

test_that("the informative mask counts what the informative columns count", {
  mask <- imgt_v_informative_positions(full_sequence, fixture$GERMLINE, fixture$N_POSITIONS)
  expect_identical(length(mask), fixture$N_POSITIONS)
  expect_identical(sum(mask), full$informative_v)
  for (region in V_REGION_LEVELS) {
    expect_identical(sum(mask & fixture$REGION_OF_POSITION == region),
                     full[[paste0("informative_", region)]])
  }
})

test_that("a truncated record's mask is empty where it observed nothing", {
  mask <- imgt_v_informative_positions(truncated_sequence, fixture$GERMLINE, fixture$N_POSITIONS)
  expect_identical(sum(mask), 112L)
  expect_false(any(mask[1:200]))
})

test_that("the configured vote floor changes support without changing measured SHM", {
  calls <- fixture$airr_row("assemble6_0", full_sequence, fixture$GERMLINE)
  votes <- fixture$vote_rows("assemble6_0", fixture$THIN_VOTES$total, fixture$THIN_VOTES$non_consensus)
  relaxed <- build_igh_v_region_feature_table(calls, votes, min_position_votes = 1L)
  strict <- build_igh_v_region_feature_table(calls, votes, min_position_votes = 2L)
  expect_equal(relaxed$informative_supported_fwr1, 78L)
  expect_equal(strict$informative_supported_fwr1, 0L)
  expect_identical(relaxed$shm_frequency, strict$shm_frequency)
  expect_identical(relaxed$median_majority_votes_fwr1, strict$median_majority_votes_fwr1)
  missing <- build_igh_v_region_feature_table(calls,
    fixture$vote_rows("assemble6_0", fixture$EMPTY_VOTES, fixture$EMPTY_VOTES), min_position_votes = 0L)
  expect_identical(fixture$region_values(missing, "informative_supported_"), rep(0L, 5L))
})

test_that("the vote floor includes exactly two agreeing reads", {
  majority <- position_majority_votes(c(20L, 20L, 20L), c(19L, 18L, 17L))
  support <- summarise_vote_support_by_region(
    majority, rep(TRUE, 3), 2L, rep("fwr1", 3))
  expect_identical(majority, c(1L, 2L, 3L))
  expect_equal(support$n_supported[["fwr1"]], 2L)
  expect_equal(support$median[["fwr1"]], 2)
})

describe("the per-region vote summary", {
  mask <- imgt_v_informative_positions(full_sequence, fixture$GERMLINE, fixture$N_POSITIONS)
  majority <- position_majority_votes(
    parse_vote_vector(fixture$VOTES$total, fixture$N_POSITIONS),
    parse_vote_vector(fixture$VOTES$non_consensus, fixture$N_POSITIONS)
  )
  support <- summarise_vote_support_by_region(majority, mask, fixture$MIN_POSITION_VOTES,
                                                fixture$REGION_OF_POSITION)

  it("medians each region over its own positions",
     expect_equal(unname(support$median[V_REGION_LEVELS]),
                  as.numeric(fixture$REGION_MAJORITY_VOTES[V_REGION_LEVELS])))
  it("says how many positions each median stands on",
     expect_identical(unname(support$n[V_REGION_LEVELS]),
                      unname(fixture$VOTED_WIDTH[V_REGION_LEVELS])))
  it("leaves a region with no counted position without a median", {
    nothing <- summarise_vote_support_by_region(
      majority, rep(FALSE, fixture$N_POSITIONS), fixture$MIN_POSITION_VOTES, fixture$REGION_OF_POSITION)
    expect_true(all(is.na(nothing$median)))
    expect_identical(unname(nothing$n), rep(0L, 5L))
  })
  it("refuses vectors that do not describe the same window",
     expect_error(summarise_vote_support_by_region(majority[-1], mask,
                                                   fixture$MIN_POSITION_VOTES, fixture$REGION_OF_POSITION),
                  "IMGT V window", fixed = TRUE))
})

# SHazaM scores each mismatch alone in its germline codon, and an N is not an
# observed base. The released mutations_ and informative_ columns follow both.
test_that("a mismatch that alone would make a stop, and an N, count no mutation", {
  # Germline codon 2 is TGG. The record reads TAC, which is not a stop, but G5A
  # alone would make TAG, so only G6C, making TGC, is scored.
  germline <- fixture$GERMLINE
  substr(germline, 4, 6) <- "TGG"
  sequence <- germline
  substr(sequence, 4, 6) <- "TAC"
  substr(sequence, 130, 130) <- "N"
  features <- build_igh_v_region_feature_table(
    fixture$airr_row("assemble7_0", sequence, germline), fixture$vote_rows("assemble7_0"),
    fixture$MIN_POSITION_VOTES)
  expect_identical(fixture$region_values(features, "mutations_"), c(1L, 0L, 0L, 0L, 0L))
  expect_identical(fixture$region_values(features, "informative_"),
                   c(78L, 36L, 50L, 30L, 117L))
  expect_equal(features$shm_frequency, 1 / 311)
})
