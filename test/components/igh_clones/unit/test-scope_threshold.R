# Check configured, fitted and fallback clone thresholds within each scope.
# Run: python test/run_tests.py --component igh_clones --level unit

load_test_modules(c("table_io", "scope_clone_threshold"))

describe("parse_configured_threshold reads the clonal_threshold setting", {
  it("reads 'auto' as no configured number",
     expect_identical(parse_configured_threshold("auto"), NA_real_))
  it("reads a number", {
    expect_equal(parse_configured_threshold("0.08"), 0.08)
    expect_equal(parse_configured_threshold("0.15"), 0.15)
  })
  it("accepts both ends of the distance range", {
    expect_equal(parse_configured_threshold("0"), 0)
    expect_equal(parse_configured_threshold("1"), 1)
  })
  for (bad in list("1.5", "-0.1", "AUTO", "", "0.1 or so", NA_character_)) {
    it(sprintf("rejects %s", if (is.na(bad)) "a missing value" else sQuote(bad)), {
      expect_error(parse_configured_threshold(bad),
                   "must be 'auto' or a number between 0 and 1")
    })
  }
})

# --------------------------------------------------- which donors may speak ---

describe("donor_fit_is_usable selects the donor fits that inform the median", {
  usable <- donor_fit_is_usable(
    threshold = c(0.12, 0.20, NA_real_, 0.11, 1.4, -0.1, 0.0, 1.0),
    fitted    = c(TRUE, FALSE, TRUE,     NA,  TRUE, TRUE, TRUE, TRUE))

  it("keeps a fitted distance",
     expect_identical(usable[1], TRUE))
  it("drops a donor whose own stage said it could not be fitted",
     expect_identical(usable[2], FALSE))
  it("drops a fitted donor with no value",
     expect_identical(usable[3], FALSE))
  it("drops an unreadable flag rather than treating it as a pass",
     expect_identical(usable[4], FALSE))
  it("drops a value that is not a normalised distance",
     expect_identical(usable[5:6], c(FALSE, FALSE)))
  it("keeps both ends of the distance range",
     expect_identical(usable[7:8], c(TRUE, TRUE)))
})

# ------------------------------------------------------------- the decision ---

describe("choose_scope_clone_threshold with a configured cutoff", {
  choice <- choose_scope_clone_threshold(c(0.12, 0.20), configured = 0.08)

  it("uses the configured number", expect_equal(choice$threshold, 0.08))
  it("says so", expect_identical(choice$source, "configured"))
  # The donor fits are still counted, because they are still generated and
  # retained for review even when they are not operative.
  it("still counts the donors that were fitted",
     expect_identical(choice$n_donors_fitted, 2L))
})

describe("choose_scope_clone_threshold in automatic mode", {
  it("takes the one fit when there is one", {
    choice <- choose_scope_clone_threshold(0.17, configured = NA_real_)
    expect_equal(choice$threshold, 0.17)
    expect_identical(choice$source, "automatic_fit")
    expect_identical(choice$n_donors_fitted, 1L)
  })
  it("takes the median of an odd number of fits", {
    choice <- choose_scope_clone_threshold(c(0.30, 0.10, 0.20),
                                             configured = NA_real_)
    expect_equal(choice$threshold, 0.20)
    expect_identical(choice$n_donors_fitted, 3L)
  })
  # The median rather than the mean, so one donor with a pathological
  # distribution cannot pull the cutoff for everyone else: the mean here is 0.30.
  it("takes the median of an even number of fits, and ignores an outlier", {
    choice <- choose_scope_clone_threshold(c(0.10, 0.14, 0.16, 0.80),
                                             configured = NA_real_)
    expect_equal(choice$threshold, 0.15)
    expect_identical(choice$source, "automatic_fit")
  })
})

describe("choose_scope_clone_threshold with nothing to learn from", {
  choice <- choose_scope_clone_threshold(numeric(0), configured = NA_real_)

  it("uses the stated fallback",
     expect_equal(choice$threshold, AUTOMATIC_CLONE_THRESHOLD_FALLBACK))
  it("that fallback is the human bulk BCR reference distance",
     expect_equal(AUTOMATIC_CLONE_THRESHOLD_FALLBACK, 0.1))
  it("says the value was assumed rather than measured",
     expect_identical(choice$source, "automatic_fallback"))
  it("reports no fitted donor",
     expect_identical(choice$n_donors_fitted, 0L))
  it("can be overridden by the caller", {
    expect_equal(
      choose_scope_clone_threshold(numeric(0), NA_real_, fallback = 0.25)$threshold,
      0.25)
  })
})

# ---------------------------------------------------------- the input check ---

DONORS <- data.frame(
  scope_id = "cohort", donor_id = c("D1", "D2"), n_libraries = 2L,
  n_members = c(650L, 720L), n_distances = c(600L, 700L),
  threshold = c(0.11, 0.13), fit_status = "fitted", stringsAsFactors = FALSE)
BINS <- data.frame(
  scope_id = "cohort", donor_id = rep(c("D1", "D2"), each = 2L),
  distance_bin_start = c(0, 0.01, 0, 0.01), stringsAsFactors = FALSE)

test_that("consistent donor fits and diagnostics pass", {
  expect_silent(validate_donor_fit_consistency(DONORS, BINS, "cohort"))
})

test_that("an empty scope passes, because a scope may have no donor at all", {
  expect_silent(validate_donor_fit_consistency(DONORS[0, ], BINS[0, ], "cohort"))
})

test_that("a fit from another scope is refused", {
  strayed <- DONORS
  strayed$scope_id[2] <- "other_project"
  expect_error(validate_donor_fit_consistency(strayed, BINS, "cohort"),
               "donor threshold table contains a different scope_id")
})

test_that("diagnostics from another scope are refused", {
  strayed <- BINS
  strayed$scope_id[1] <- "other_project"
  expect_error(validate_donor_fit_consistency(DONORS, strayed, "cohort"),
               "distance diagnostics contain a different scope_id")
})

# One donor counted twice would move the median without being visible in the
# released row.
test_that("the same donor twice is refused", {
  doubled <- rbind(DONORS, DONORS[1, ])
  expect_error(validate_donor_fit_consistency(doubled, BINS, "cohort"),
               "duplicate donor_id")
})

test_that("a donor in one table and not the other is refused", {
  expect_error(
    validate_donor_fit_consistency(DONORS[1, ], BINS, "cohort"),
    "cover different donors")
  expect_error(
    validate_donor_fit_consistency(DONORS, BINS[BINS$donor_id == "D1", ], "cohort"),
    "cover different donors")
})

# ------------------------------------------------------- reading it back in ---

write_scope_row <- function(...) {
  path <- tempfile(fileext = ".tsv")
  write.table(data.frame(..., stringsAsFactors = FALSE), path, sep = "\t",
              row.names = FALSE, quote = FALSE)
  path
}

test_that("a scope's frozen distance is read back", {
  path <- write_scope_row(scope_id = "cohort", clone_threshold = 0.1234,
                          threshold_source = "automatic_fit",
                          n_donors = 3L, n_donors_fitted = 3L)
  expect_equal(read_scope_clone_threshold(path, "cohort"), 0.1234)
})

# Assignment runs once per donor from a path built out of wildcards, so being
# handed another scope's table has to fail rather than cluster at the wrong
# distance.
test_that("another scope's threshold table is refused", {
  path <- write_scope_row(scope_id = "other_project", clone_threshold = 0.1)
  expect_error(read_scope_clone_threshold(path, "cohort"),
               "does not match expected scope_id cohort")
})

test_that("a threshold table that is not one row is refused", {
  path <- write_scope_row(scope_id = c("cohort", "cohort"),
                          clone_threshold = c(0.1, 0.2))
  expect_error(read_scope_clone_threshold(path, "cohort"),
               "must contain exactly one row")
})

test_that("a cutoff that is not a distance is refused", {
  for (bad in list("1.5", "-0.2", "unset", "")) {
    path <- write_scope_row(scope_id = "cohort", clone_threshold = bad)
    expect_error(read_scope_clone_threshold(path, "cohort"),
                 "must be between 0 and 1")
  }
})

test_that("a table missing the cutoff column names itself in the message", {
  path <- write_scope_row(scope_id = "cohort", threshold = 0.1)
  expect_error(read_scope_clone_threshold(path, "cohort"),
               "Scope threshold table is missing required columns: clone_threshold")
})
