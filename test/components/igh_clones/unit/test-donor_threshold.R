# Check pinned donor fitting, seed 1337 and retained calibration evidence.
# Run: python test/run_tests.py --component igh_clones --level unit

load_test_modules(c(
  "table_io", "gene_calls", "donor_distance_threshold", "scope_clone_threshold"))

test_that("the pipeline seed is fixed at 1337", {
  expect_identical(PIPELINE_SEED, 1337L)
})

# ----------------------------------------------------------- the evidence gate ---

describe("distance_evidence_status says why a mixture cannot be identified", {
  it("needs the minimum number of distances", {
    expect_identical(distance_evidence_status(seq(0, 1, length.out = 500), 500L),
                     NA_character_)
    expect_identical(distance_evidence_status(seq(0, 1, length.out = 499), 500L),
                     "too_few_distances")
    expect_identical(distance_evidence_status(seq(0, 1, length.out = 501), 500L),
                     NA_character_)
  })
  # A donor whose distances are all the same value has no mixture to find,
  # however many of them there are.
  it("needs at least two distinct values", {
    expect_identical(distance_evidence_status(rep(0.1, 1000L), 500L),
                     "constant_distances")
    expect_identical(distance_evidence_status(c(rep(0.1, 999L), 0.2), 500L),
                     NA_character_)
  })
  it("reports a short constant donor as too few distances",
     expect_identical(distance_evidence_status(rep(0.1, 10L), 500L),
                      "too_few_distances"))
  it("refuses an empty donor",
     expect_identical(distance_evidence_status(numeric(0), 2L), "too_few_distances"))
})

describe("is_threshold_in_range decides whether a cutoff is a distance", {
  it("accepts a normalised distance", expect_true(is_threshold_in_range(0.13)))
  it("accepts both ends", {
    expect_true(is_threshold_in_range(0))
    expect_true(is_threshold_in_range(1))
  })
  it("refuses anything outside the range", {
    expect_false(is_threshold_in_range(1.0001))
    expect_false(is_threshold_in_range(-1e-9))
  })
  it("refuses a missing or infinite value", {
    expect_false(is_threshold_in_range(NA_real_))
    expect_false(is_threshold_in_range(Inf))
  })
})

describe("fit_threshold_status names the outcome of an attempted fit", {
  it("calls a missing cutoff a failed fit",
     expect_identical(fit_threshold_status(NA_real_), "fit_failed"))
  it("keeps an out-of-range cutoff apart from a failure", {
    expect_identical(fit_threshold_status(1.2), "out_of_range")
    expect_identical(fit_threshold_status(-0.1), "out_of_range")
  })
  it("accepts a normalised distance",
     expect_identical(fit_threshold_status(0.13), "fitted"))
})

# --------------------------------------------------------------- the two tables ---

test_that("a donor with no records still reports a row", {
  row <- build_donor_threshold_row("cohort", "D1")
  expect_identical(names(row), DONOR_THRESHOLD_COLUMNS)
  expect_identical(c(row$n_libraries, row$n_members, row$n_distances), c(0L, 0L, 0L))
  expect_identical(row$threshold, NA_real_)
  expect_identical(row$fit_status, "too_few_distances")
  expect_identical(row$n_distances_fitted, NA_integer_)
  expect_true(all(is.na(row[names(DONOR_FIT_PARAMETERS)])))
})

test_that("a fitted donor reports its cutoff and populations", {
  row <- build_donor_threshold_row("cohort", "D1", 2, 650, 600, 0.13, "fitted")
  expect_identical(c(row$n_libraries, row$n_members, row$n_distances), c(2L, 650L, 600L))
  expect_equal(row$threshold, 0.13)
  expect_identical(row$fit_status, "fitted")
  expect_error(build_donor_threshold_row("cohort", "D1", fit_status = "fittable"),
               "fit_status must be one of")
})

describe("the review surface of a donor with no distances", {
  bins <- build_distance_histogram_diagnostics("cohort", "D1")

  it("still has every bin", expect_identical(nrow(bins), 100L))
  it("declares the released columns",
     expect_identical(names(bins), DONOR_DIAGNOSTIC_COLUMNS))
  it("uses the shared 0.01-wide breaks", {
    expect_equal(bins$distance_bin_start, head(DISTANCE_HISTOGRAM_BREAKS, -1L))
    expect_equal(bins$distance_bin_end, tail(DISTANCE_HISTOGRAM_BREAKS, -1L))
    expect_equal(bins$distance_bin_midpoint, bins$distance_bin_start + 0.005)
  })
  it("counts nothing", expect_true(all(bins$n_distances_in_bin == 0L)))
  it("reports zero density rather than an undefined one",
     expect_true(all(bins$observed_density == 0)))
  it("has no fitted components", {
    expect_true(all(is.na(bins$near_component_density)))
    expect_true(all(is.na(bins$far_component_density)))
  })
})

describe("the review surface of a donor with distances but no fit", {
  # 0.01 is a break, and 1 is the last one: bins are left-closed, and the highest
  # value belongs to the last bin rather than falling off the end.
  distances <- c(0, 0.005, 0.01, 0.5, 1)
  bins <- build_distance_histogram_diagnostics("cohort", "D1", distances)

  it("counts every distance", expect_identical(sum(bins$n_distances_in_bin), 5L))
  # Bins are left-closed, so 0.01 opens the second bin rather than closing the
  # first, and 0 and 0.005 are the first bin's two members.
  it("puts a value on a break into the bin it opens",
     expect_identical(bins$n_distances_in_bin[2], 1L))
  it("puts the two below it in the first bin",
     expect_identical(bins$n_distances_in_bin[1], 2L))
  it("keeps the highest value in the last bin",
     expect_identical(bins$n_distances_in_bin[100], 1L))
  it("has a density that integrates to one",
     expect_equal(sum(bins$observed_density * diff(DISTANCE_HISTOGRAM_BREAKS)), 1))
  it("still has no fitted components",
     expect_true(all(is.na(bins$near_component_density))))
})

# ------------------------------------------------------------------- the fit ---
#
FIT_FIXTURE <- separated_distance_fixture()
FIT_FIXTURE$v_gene <- sub("\\*.*$", "", FIT_FIXTURE$v_call)
FIT_FIXTURE$j_gene <- "IGHJ4"

# The sparse case first, because it returns before the override runs and so
# leaves the shazam namespace alone.
test_that("a donor below the evidence floor is measured but not fitted", {
  sparse <- fit_donor_distance_threshold(FIT_FIXTURE[1:6, ], 500L, 0.9, 0.995, 10000)
  expect_identical(length(sparse$distances), 6L)
  expect_identical(sparse$fit_status, "too_few_distances")
  expect_identical(sparse$threshold, NA_real_)
  expect_null(sparse$fit)
  expect_null(sparse$error)
  # No fit was attempted, so no distance was used by one.
  expect_identical(sparse$n_distances_fitted, NA_integer_)
})

# The fitting checks restore the package namespace even when a test aborts.
check_seeded_fitting <- function() {
  ORIGINAL_FITTERS <- lapply(
    stats::setNames(nm = c("gmmFit", "rocSpace")),
    function(name) get(name, envir = asNamespace("shazam"), inherits = FALSE))
  on.exit({
    for (name in names(ORIGINAL_FITTERS)) {
      assignInNamespace(name, ORIGINAL_FITTERS[[name]], ns = "shazam")
    }
  }, add = TRUE)
  fitter_body_text <- function(name) {
    paste(deparse(body(get(name, envir = asNamespace("shazam"),
                           inherits = FALSE))), collapse = "\n")
  }
  count_seed_resets <- function(name) {
    positions <- gregexpr("set.seed(NULL)", fitter_body_text(name),
                          fixed = TRUE)[[1]]
    if (positions[1] == -1L) 0L else length(positions)
  }

  test_that("the installed SHazaM resets the seed exactly once in each fitter", {
    expect_identical(count_seed_resets("gmmFit"), 1L)
    expect_identical(count_seed_resets("rocSpace"), 1L)
  })

  FIT <- fit_donor_distance_threshold(FIT_FIXTURE, 500L, 0.9, 0.995, 10000)

  test_that("fitting replaced the internal reset with the pipeline seed", {
    expect_identical(count_seed_resets("gmmFit"), 0L)
    expect_identical(count_seed_resets("rocSpace"), 0L)
    expect_true(grepl(paste0("set.seed(", PIPELINE_SEED),
                      fitter_body_text("gmmFit"), fixed = TRUE))
  })

  test_that("the seed rewrite preserves every other package expression, including NULL", {
    for (name in names(ORIGINAL_FITTERS)) {
      original <- paste(deparse(body(ORIGINAL_FITTERS[[name]])), collapse = "\n")
      expected <- parse(text = gsub("set.seed(NULL)", "set.seed(1337L)",
                                    original, fixed = TRUE), keep.source = FALSE)[[1]]
      patched <- get(name, envir = asNamespace("shazam"), inherits = FALSE)
      expect_identical(body(patched), expected)
    }
  })

  # The assertion is what stops a SHazaM upgrade quietly restoring the
  # non-determinism: with the reset already replaced there is nothing left to find.
  test_that("an override that finds no reset stops the stage", {
    expect_error(override_shazam_gmm_seed(1L),
                 "Expected one internal SHazaM seed reset in gmmFit; found 0")
  })

  test_that("a separated mixture is fitted", {
    expect_identical(length(FIT$distances), 600L)
    expect_identical(FIT$fit_status, "fitted")
    expect_true(is_threshold_in_range(FIT$threshold))
    expect_null(FIT$error)
    expect_true(isS4(FIT$fit) && is(FIT$fit, "GmmThreshold"))
    # The ceiling above the fixture's distances leaves all of them in the fit.
    expect_identical(FIT$n_distances_fitted, length(FIT$fit@x))
    expect_identical(FIT$n_distances_fitted, 600L)
  })

  describe("the row of a fitted donor", {
    row <- build_donor_threshold_row("fit_scope", "fit_donor", 1L, 600L, 600L, FIT$threshold,
                                     FIT$fit_status, FIT$n_distances_fitted, FIT$fit)

    it("carries the fit's columns in their written order",
       expect_identical(names(row), DONOR_THRESHOLD_COLUMNS))
    it("reads each mixture parameter from its slot", {
      expect_equal(c(row$near_weight, row$near_shape, row$near_scale),
                   c(FIT$fit@a1, FIT$fit@b1, FIT$fit@c1))
      expect_equal(c(row$far_weight, row$far_mean, row$far_sd),
                   c(FIT$fit@a2, FIT$fit@b2, FIT$fit@c2))
      expect_equal(row$near_weight + row$far_weight, 1)
    })
    it("reads the fit statistics from their slots", {
      expect_equal(row$absolute_log_likelihood, FIT$fit@loglk)
      expect_equal(row$fit_sensitivity, FIT$fit@sensitivity)
    })
    it("reproduces both component densities of the histogram", {
      bins <- build_distance_histogram_diagnostics("fit_scope", "fit_donor",
                                                   FIT$distances, FIT$fit)
      expect_equal(bins$near_component_density,
                   row$near_weight * dgamma(bins$distance_bin_midpoint, shape = row$near_shape,
                                            scale = row$near_scale))
      expect_equal(bins$far_component_density,
                   row$far_weight * dnorm(bins$distance_bin_midpoint, mean = row$far_mean,
                                          sd = row$far_sd))
    })
  })

  # A cutoff inside [0, 1] is not the same thing as a usable cutoff. This fixture's
  # two modes sit at 1 to 4 mismatches and 14 to 18 mismatches over a 60 nt
  # junction, so a cutoff that separates them has to land in the empty band
  # between: above every near distance and below every far one. Without this, the
  # arguments to findThreshold() and distToNearest() could all change and move the
  # cutoff for every scope with the suite still green.
  test_that("the fitted cutoff falls in the gap between the two modes", {
    near_distances <- 4 / 60
    far_distances <- 14 / 60
    expect_gt(FIT$threshold, near_distances)
    expect_lt(FIT$threshold, far_distances)
  })

  # The override is one-shot within a process, which is why the reset has to be put
  # back before the fit can be repeated here. The stage itself is one Rscript
  # invocation per donor, so it overrides once and exits.
  test_that("the same input gives the same cutoff", {
    for (name in names(ORIGINAL_FITTERS)) {
      assignInNamespace(name, ORIGINAL_FITTERS[[name]], ns = "shazam")
    }
    again <- fit_donor_distance_threshold(FIT_FIXTURE, 500L, 0.9, 0.995, 10000)
    expect_identical(again$threshold, FIT$threshold)
  })

  describe("the review surface of a fitted donor", {
    bins <- build_distance_histogram_diagnostics("fit_scope", "fit_donor",
                                                 FIT$distances, FIT$fit)

    it("counts every distance",
       expect_identical(sum(bins$n_distances_in_bin), 600L))
    it("evaluates both fitted components over the same grid", {
      expect_true(any(is.finite(bins$near_component_density)))
      expect_true(any(is.finite(bins$far_component_density)))
    })
    # The near component is the same-lineage mass and the far component the
    # different-lineage mass, so the near one has to peak at the shorter distance.
    it("puts the near component's mass below the far one's", {
      expect_true(which.max(bins$near_component_density) <
                    which.max(bins$far_component_density))
    })
    it("drops the components when no fit is passed", {
      expect_true(all(is.na(
        build_distance_histogram_diagnostics("fit_scope", "fit_donor",
                                             FIT$distances)$near_component_density)))
    })
  })

  for (name in names(ORIGINAL_FITTERS)) {
    assignInNamespace(name, ORIGINAL_FITTERS[[name]], ns = "shazam")
  }

  test_that("a failed initialization returns no fit and the automatic policy falls back", {
    expect_output(
      sparse <- fit_donor_distance_threshold(FIT_FIXTURE[1:4, ], 2L,
                                             0.9, 0.995, NULL),
      "No fit found", fixed = TRUE)
    expect_identical(length(sparse$distances), 4L)
    # findThreshold found no fit and raised no error, which is still a failed fit.
    expect_identical(sparse$fit_status, "fit_failed")
    expect_null(sparse$fit)
    expect_null(sparse$error)
    # A fit was attempted on every distance, and it reports no parameters.
    expect_identical(sparse$n_distances_fitted, 4L)
    expect_true(all(is.na(build_donor_threshold_row(
      "cohort", "D1", 1L, 4L, 4L, sparse$threshold, sparse$fit_status,
      sparse$n_distances_fitted, sparse$fit)[names(DONOR_FIT_PARAMETERS)])))
    choice <- choose_scope_clone_threshold(numeric(0), NA_real_)
    expect_identical(choice$threshold, 0.1)
    expect_identical(choice$source, "automatic_fallback")
    for (name in names(ORIGINAL_FITTERS)) {
      assignInNamespace(name, ORIGINAL_FITTERS[[name]], ns = "shazam")
    }
  })

  test_that("configured fitting controls reach the fitted distribution and cutoff", {
    changed <- fit_donor_distance_threshold(FIT_FIXTURE, 500L, 0.8, 0.99, 300L)
    expect_identical(changed$fit_status, "fitted")
    expect_identical(length(changed$distances), 600L)
    expect_identical(length(changed$fit@x), 300L)
    expect_identical(changed$n_distances_fitted, 300L)
    expect_equal(changed$fit@specificity, 0.99, tolerance = 1e-6)
    for (name in names(ORIGINAL_FITTERS)) {
      assignInNamespace(name, ORIGINAL_FITTERS[[name]], ns = "shazam")
    }
  })

  # findThreshold() subsamples by position, and record order follows TRUST4
  # numbering, which changes between runs.
  test_that("a subsampled fit does not depend on record order", {
    forwards <- fit_donor_distance_threshold(FIT_FIXTURE, 500L, 0.9, 0.995, 300L)
    for (name in names(ORIGINAL_FITTERS)) {
      assignInNamespace(name, ORIGINAL_FITTERS[[name]], ns = "shazam")
    }
    reversed <- FIT_FIXTURE[rev(seq_len(nrow(FIT_FIXTURE))), , drop = FALSE]
    backwards <- fit_donor_distance_threshold(reversed, 500L, 0.9, 0.995, 300L)
    expect_identical(backwards$threshold, forwards$threshold)
    for (name in names(ORIGINAL_FITTERS)) {
      assignInNamespace(name, ORIGINAL_FITTERS[[name]], ns = "shazam")
    }
  })

  test_that("the shazam namespace is left as it was found", {
    expect_identical(count_seed_resets("gmmFit"), 1L)
    expect_identical(count_seed_resets("rocSpace"), 1L)
  })
}
check_seeded_fitting()

# --------------------------------------------------------- reading them back ---

test_that("a donor fit round-trips through its table", {
  path <- tempfile(fileext = ".tsv")
  write_pipeline_table(
    build_donor_threshold_row("cohort", "D1", 2L, 650L, 600L, 0.13, "fitted"), path)
  row <- read_donor_threshold_row(path)
  expect_identical(nrow(row), 1L)
  expect_identical(as.character(row$donor_id), "D1")
  expect_equal(as.numeric(row$threshold), 0.13)
  expect_identical(row$fit_status, "fitted")
})

test_that("a review surface round-trips through its table", {
  path <- tempfile(fileext = ".tsv")
  write_pipeline_table(
    build_distance_histogram_diagnostics("cohort", "D1", c(0.1, 0.2)), path)
  bins <- read_donor_distance_diagnostics(path)
  expect_identical(nrow(bins), 100L)
  expect_identical(names(bins), DONOR_DIAGNOSTIC_COLUMNS)
})

# A short histogram means the writer and the reader no longer agree on what one
# bin is, and the pooled scope table would otherwise mix two grids.
test_that("a truncated review surface is refused", {
  path <- tempfile(fileext = ".tsv")
  bins <- build_distance_histogram_diagnostics("cohort", "D1")
  write_pipeline_table(bins[1:99, , drop = FALSE], path)
  expect_error(read_donor_distance_diagnostics(path),
               "Invalid donor distance diagnostics table")
})

test_that("a donor fit table that is not one row is refused", {
  path <- tempfile(fileext = ".tsv")
  write_pipeline_table(rbind(
    build_donor_threshold_row("cohort", "D1"),
    build_donor_threshold_row("cohort", "D2")), path)
  expect_error(read_donor_threshold_row(path), "Invalid donor threshold table")
})
