# Check library-local feature joins, string identifiers and missing evidence.
# Run: python test/run_tests.py --component scopes --level unit

load_test_modules(c("table_io", "pooled_records"))

test_that("the key namespaces a library-local identifier", {
  expect_identical(make_pooled_sequence_id("WTS_1_S1", "assemble3_0"),
                   "WTS_1_S1~assemble3_0")
})

test_that("the key is vectorised over both arguments", {
  expect_identical(
    make_pooled_sequence_id(c("WTS_1_S1", "WTS_2_S1"), c("a0_0", "a1_0")),
    c("WTS_1_S1~a0_0", "WTS_2_S1~a1_0")
  )
})

# ------------------------------------------- joining the feature sidecar on ---
#
# Feature rows are keyed per library and arrive in a different order from the
# membership they are joined onto, so the join has to be a lookup rather than a
# row-for-row assumption. Every table is read as character, so every rate and
# count arrives as text and leaves as a number.

# These constructed regions make the columns carried by the join explicit.
REGION_LEVELS <- c("fwr1", "cdr1", "fwr2", "cdr2", "fwr3")

MEMBERSHIP <- data.frame(
  donor_id = c("D1", "D1", "D2"),
  library_id = c("S1", "S1", "S2"),
  sequence_id = c("S1~seq1", "S1~seq2", "S2~seq3"),
  clone_id = c("L1", "L1", "L2"),
  consensus_count = c("3", "4", "5"),
  stringsAsFactors = FALSE
)

FEATURES <- data.frame(
  library_id = c("S2", "S1", "S1", "S1"),
  sequence_id = c("seq3", "seq2", "seq1", "seq9"),
  qualified_c_call = c("IGHM*01", "IGHG1*01", "IGHA1*01", "IGHG2*01"),
  shm_frequency = c("0.05", "0.30", "", "0.90"),
  informative_v = c("280", "270", "12", "300"),
  stringsAsFactors = FALSE
)
# Three counts nest per sub-region -- occupied, of those informative, of those
# read-supported -- because the SHM eligibility rule is a comparison between the
# first and the last. The third row is the sequence that measured almost nothing,
# and the fourth is the feature row no member claims.
for (region in REGION_LEVELS) {
  FEATURES[[paste0("shm_frequency_", region)]] <- c("0.05", "0.30", "", "0.90")
  FEATURES[[paste0("mutations_", region)]] <- c("2", "9", "0", "45")
  FEATURES[[paste0("informative_", region)]] <- c("40", "30", "2", "50")
  FEATURES[[paste0("informative_supported_", region)]] <- c("35", "30", "0", "50")
  FEATURES[[paste0("germline_occupied_", region)]] <- c("60", "60", "60", "60")
  FEATURES[[paste0("median_majority_votes_", region)]] <- c("7", "12", "1", "22")
}

attach_features <- function(membership, features) {
  attach_igh_features_to_membership(membership, features, REGION_LEVELS)
}

describe("a scope whose members all have a feature row", {
  joined <- attach_features(MEMBERSHIP, FEATURES)

  it("keeps one row per member",
     expect_identical(nrow(joined), nrow(MEMBERSHIP)))
  it("attaches each member's own constant-region call",
     expect_identical(joined$qualified_c_call, c("IGHA1*01", "IGHG1*01", "IGHM*01")))
  it("reads shm_frequency back as a number",
     expect_equal(joined$shm_frequency, c(NA, 0.30, 0.05)))
  # Every column a region's eligibility rule reads has to arrive, or the rule
  # would silently admit nothing.
  it("carries each region's rate and both its denominators", {
    expect_equal(joined$shm_frequency_cdr2, c(NA, 0.30, 0.05))
    expect_identical(joined$informative_cdr2, c(2, 30, 40))
    expect_identical(joined$germline_occupied_cdr2, c(60, 60, 60))
    expect_identical(joined$informative_v, c(12, 270, 280))
  })
  # A feature row for a record that was never lineage-callable is expected, and
  # is not an error: the sidecar covers every common-QC IGH record.
  it("ignores a feature row no member claims",
     expect_false("seq9" %in% joined$sequence_id))
  it("leaves the member columns as they were",
     expect_identical(joined$clone_id, MEMBERSHIP$clone_id))
})

# A caller projects the join's output onto a released schema, so a scope that
# contributed nothing needs the columns to exist and be empty rather than to be
# absent -- otherwise the projection fails instead of writing the header.
describe("a scope with no members", {
  joined <- attach_features(MEMBERSHIP[0, , drop = FALSE], FEATURES)

  it("has no rows", expect_identical(nrow(joined), 0L))
  it("still declares every joined column", {
    expect_identical(joined$qualified_c_call, character(0))
    expect_identical(joined$shm_frequency, numeric(0))
    expect_identical(joined$mutations_cdr2, numeric(0))
    expect_identical(joined$germline_occupied_fwr3, numeric(0))
  })
  it("keeps the member columns it was given",
     expect_true(all(names(MEMBERSHIP) %in% names(joined))))
})

test_that("no members and no features still declares the joined columns", {
  joined <- attach_features(MEMBERSHIP[0, , drop = FALSE], data.frame())
  expect_identical(nrow(joined), 0L)
  expect_identical(joined$qualified_c_call, character(0))
  expect_identical(joined$shm_frequency, numeric(0))
})

# ------------------------------------------------------------ error paths ---

test_that("a member with no feature row stops the scope", {
  expect_error(
    attach_features(MEMBERSHIP, FEATURES[1:2, , drop = FALSE]),
    "contains 1 sequence\\(s\\) missing from IGH features"
  )
})

test_that("no feature rows at all stops a scope that has members", {
  expect_error(
    attach_features(MEMBERSHIP, data.frame()),
    "contains 3 sequence\\(s\\) missing from IGH features"
  )
})

# Two libraries can each hold a sequence_id, so the uniqueness that matters is of
# the namespaced key, not of the identifier the library assigned.
test_that("a repeated pooled key stops the scope", {
  duplicated_features <- rbind(FEATURES, FEATURES[3, , drop = FALSE])
  expect_error(
    attach_features(MEMBERSHIP, duplicated_features),
    "sequence_id is not unique across IGH feature tables"
  )
})

test_that("a feature table missing a joined column is named in the message", {
  expect_error(
    attach_features(MEMBERSHIP, FEATURES[, 1:3, drop = FALSE]),
    "IGH feature table is missing required columns: shm_frequency"
  )
})


test_that("malformed joined numeric evidence names its column and pooled record", {
  for (column in c("shm_frequency", "informative_v", "mutations_cdr2")) {
    features <- FEATURES
    features[[column]][3] <- "invalid"
    expect_error(attach_features(MEMBERSHIP, features),
                 paste0("column ", column, ".*record S1~seq1.*invalid"))
  }
})

test_that("numeric evidence preserves allowed missing values and unused rows", {
  features <- FEATURES
  features$shm_frequency <- c(NA_character_, " ", "", "unused text")
  joined <- attach_features(MEMBERSHIP, features)
  expect_identical(joined$shm_frequency, rep(NA_real_, 3))
  features$shm_frequency <- c("NaN", "Inf", "0", "unused text")
  expect_equal(attach_features(MEMBERSHIP, features)$shm_frequency, c(0, Inf, NaN))
})
