# Check analytical rarefaction against exact draws and large-support boundaries.
# Run: python test/run_tests.py --component diversity --level unit

load_test_modules(c("table_io", "diversity"))

test_that("the three orders follow exact draws from two groups", {
  result <- interpolate_hill_diversity(c(2, 1), 2)
  expect_equal(result$richness, 5 / 3, tolerance = 1e-12)
  expect_equal(result$hill_shannon, 2^(2 / 3), tolerance = 1e-12)
  expect_equal(result$inverse_simpson, 3 / 2, tolerance = 1e-12)
  # Averaging exp(H) or inverse concentration over the three draws gives 5/3.
  expect_gt(abs(result$hill_shannon - 5 / 3), 0.01)
  expect_gt(abs(result$inverse_simpson - 5 / 3), 0.01)
})

test_that("interpolation agrees with enumeration over uneven repertoires", {
  for (counts in list(c(3, 2, 1), c(2, 2, 1, 1), c(4, 1, 1))) {
    groups <- rep(seq_along(counts), counts)
    for (depth in seq_along(groups)) {
      draws <- combn(seq_along(groups), depth, simplify = FALSE)
      observations <- vapply(draws, function(draw) {
        drawn <- tabulate(groups[draw], nbins = length(counts))
        p <- drawn[drawn > 0] / depth
        c(length(p), -sum(p * log(p)), sum(p^2))
      }, numeric(3))
      expected <- rowMeans(observations)
      result <- interpolate_hill_diversity(counts, depth)
      expect_equal(result$richness, expected[1], tolerance = 1e-12)
      expect_equal(result$hill_shannon, exp(expected[2]), tolerance = 1e-12)
      expect_equal(result$inverse_simpson, 1 / expected[3], tolerance = 1e-12)
    }
  }
})

test_that("full-depth values reproduce observed Hill diversity", {
  counts <- c(37, 21, 13, 8, 5, 3, 2, 1, 1, 1)
  p <- counts / sum(counts)
  endpoint <- interpolate_hill_diversity(counts, sum(counts))
  expect_equal(endpoint$richness, length(counts))
  expect_equal(endpoint$hill_shannon, prod(p^-p), tolerance = 1e-12)
  expect_equal(endpoint$inverse_simpson, 1 / sum(p^2), tolerance = 1e-12)
})

test_that("target equality is available independently of the diversity minimum", {
  # Target 5 separates totals 4, 5 and 6. At total 5, the only draw is the
  # complete population: richness 2, q1 = exp(H) and q2 = 25 / 13.
  target <- 5L
  below <- summarise_rarefaction(c(2, 2), target)
  equal <- summarise_rarefaction(c(3, 2), target)
  above <- summarise_rarefaction(c(4, 2), target)
  expect_identical(below$rarefaction_status, "below_target_support")
  expect_identical(equal$rarefaction_status, "available")
  expect_identical(above$rarefaction_status, "available")
  expect_equal(equal$rarefied_richness, 2, tolerance = 1e-12)
  expect_equal(equal$rarefied_hill_shannon,
               exp(-sum(c(3 / 5, 2 / 5) * log(c(3 / 5, 2 / 5)))), tolerance = 1e-12)
  expect_equal(equal$rarefied_inverse_simpson, 25 / 13, tolerance = 1e-12)
  allowed <- summarise_rarefaction(c(3, 2), target, minimum_consensus_count = 5)
  gated <- summarise_rarefaction(c(3, 2), target, minimum_consensus_count = 6)
  expect_identical(allowed$rarefaction_status, "available")
  expect_identical(gated$rarefaction_status, "below_minimum_support")
  expect_equal(gated$consensus_count, 5)
  expect_equal(gated$richness, 2)
  expect_true(is.na(gated$rarefied_inverse_simpson))
  curve <- build_rarefaction_curve(c(3, 2), target)
  expect_equal(curve$inverse_simpson[curve$is_target], 25 / 13, tolerance = 1e-12)
})

test_that("empty and singleton populations retain their defined quantities", {
  empty <- interpolate_hill_diversity(numeric(), 0)
  expect_equal(empty$richness, 0)
  expect_true(is.na(empty$hill_shannon))
  expect_true(is.na(empty$inverse_simpson))
  one_group <- interpolate_hill_diversity(12, c(1, 5, 12))
  expect_equal(one_group$richness, rep(1, 3))
  expect_equal(one_group$hill_shannon, rep(1, 3))
  expect_equal(one_group$inverse_simpson, rep(1, 3))
  singleton <- interpolate_hill_diversity(1, 1)
  expect_equal(unname(unlist(singleton)), rep(1, 4))
  all_singletons <- interpolate_hill_diversity(rep(1, 12), c(1, 5, 12))
  expect_equal(all_singletons$richness, c(1, 5, 12))
  expect_equal(all_singletons$hill_shannon, c(1, 5, 12), tolerance = 1e-12)
  expect_equal(all_singletons$inverse_simpson, c(1, 5, 12), tolerance = 1e-12)
})

test_that("q2 adjustment depends on concentration as well as support depth", {
  balanced <- interpolate_hill_diversity(c(5, 5), 5)
  dominated <- interpolate_hill_diversity(c(9, 1), 5)
  expect_equal(balanced$inverse_simpson, 1.8, tolerance = 1e-12)
  expect_equal(dominated$inverse_simpson, 25 / 21, tolerance = 1e-12)
})

test_that("curves include targets and endpoints on a bounded reproducible grid", {
  counts <- c(4000, 2000, 1000)
  curve <- build_rarefaction_curve(counts, c(400, 1000, 9000))
  expect_lte(nrow(curve), RAREFACTION_GRID_POINTS + 3L)
  expect_equal(curve$support_depth[curve$is_target], c(400, 1000))
  expect_equal(curve$support_depth[curve$is_endpoint], sum(counts))
  expect_true(all(diff(curve$support_depth) > 0))
  expect_identical(curve, build_rarefaction_curve(rev(counts), c(400, 1000, 9000)))
  expect_equal(rarefaction_depth_grid(c(2, 1), 2), c(1, 2, 3))
  empty <- build_rarefaction_curve(numeric(), 400)
  expect_equal(empty$support_depth, 0)
  expect_true(empty$is_endpoint)
  expect_false(empty$is_target)
})

test_that("minimum support gates summaries while diagnostic curves remain", {
  counts <- c(3, 2)
  below_minimum <- summarise_rarefaction(counts, 4, minimum_consensus_count = 6)
  expect_equal(below_minimum$consensus_count, 5)
  expect_equal(below_minimum$richness, 2L)
  expect_equal(below_minimum$diversity_status, "below_minimum_support")
  expect_equal(below_minimum$rarefaction_status, "below_minimum_support")
  expect_true(is.na(below_minimum$rarefied_inverse_simpson))
  curve <- build_rarefaction_curve(counts, 4)
  expect_true(all(is.finite(curve$inverse_simpson)))
  eligible <- summarise_rarefaction(counts, 4, minimum_consensus_count = 5)
  expect_equal(eligible$rarefaction_status, "available")
  expect_equal(eligible$rarefied_inverse_simpson,
               curve$inverse_simpson[curve$is_target], tolerance = 1e-12)
  below_target <- summarise_rarefaction(counts, 6)
  expect_equal(below_target$diversity_status, "available")
  expect_equal(below_target$rarefaction_status, "below_target_support")
  expect_true(is.na(below_target$rarefied_richness))
  disabled <- summarise_rarefaction(counts)
  expect_equal(disabled$rarefaction_status, "not_configured")
  expect_true(is.na(disabled$rarefaction_target))
  empty <- summarise_rarefaction(numeric(), 400)
  expect_equal(empty$diversity_status, "empty_repertoire")
  expect_equal(empty$rarefaction_status, "empty_repertoire")
})

test_that("the q2 asymptote retains its finite and undefined cases", {
  counts <- c(2, rep(1, 398))
  # Corrected concentration is 2 / (400 * 399), whose reciprocal is 79800.
  result <- q2_asymptote(counts)
  expect_identical(names(result), c("inverse_simpson_asymptote", "asymptote_status"))
  expect_equal(result$inverse_simpson_asymptote, 79800, tolerance = 1e-12)
  expect_equal(result$asymptote_status, "available")
  for (counts in list(numeric(), 1)) {
    result <- q2_asymptote(counts)
    expect_equal(result$asymptote_status, "insufficient_support")
    expect_true(is.na(result$inverse_simpson_asymptote))
  }
  all_singletons <- q2_asymptote(rep(1, 5))
  expect_equal(all_singletons$asymptote_status, "no_repeated_support")
  expect_true(is.na(all_singletons$inverse_simpson_asymptote))
  constant <- q2_asymptote(5)
  expect_equal(constant$inverse_simpson_asymptote, 1)
})

# The dictionary describes the curve through its asymptote A: inverse Simpson at
# depth m is m / (1 + (m - 1) / A), and it reaches 95 percent of A at
# 19 (A - 1), because 0.95 / (1 - 0.95) = 19.
test_that("the asymptote gives the q2 curve in closed form", {
  counts <- c(40, 25, 10, 5, 3, 1, 1)
  asymptote <- q2_asymptote(counts)$inverse_simpson_asymptote
  depths <- c(2, 7, 30, 84)
  expect_equal(interpolate_hill_diversity(counts, depths)$inverse_simpson,
               depths / (1 + (depths - 1) / asymptote), tolerance = 1e-12)
  depth_95 <- 19 * (asymptote - 1)
  expect_equal(depth_95 / (1 + (depth_95 - 1) / asymptote), 0.95 * asymptote,
               tolerance = 1e-12)
})

test_that("large integer support remains finite without binomial coefficients", {
  # choose(14597, 1000) overflows. An independently expressed q0 probability
  # product and the closed q2 expression remain available as numerical checks.
  counts <- c(14000, 500, 97)
  depth <- 1000
  result <- interpolate_hill_diversity(counts, depth)
  missing <- vapply(counts, function(count) {
    if (sum(counts) - count < depth) return(0)
    exp(sum(log1p(-count / (sum(counts) - seq.int(0, depth - 1)))))
  }, numeric(1))
  expect_true(all(is.finite(unlist(result))))
  expect_equal(result$richness, sum(1 - missing), tolerance = 1e-12)
  expect_gte(result$richness, result$hill_shannon)
  expect_gte(result$hill_shannon, result$inverse_simpson)
  repeated <- interpolate_hill_diversity(rep(100, 10000), c(400, 1000, 1000000))
  expect_true(all(is.finite(unlist(repeated))))
  expect_equal(repeated$richness[3], 10000)
  # The third depth is the full observed support, where both values are exactly
  # 10000. Each one sums a term per group over ten thousand groups, so it carries
  # more rounding than its neighbours here. The relative errors measured on one
  # host are 1.27e-12 for the Shannon value and 1.3e-13 for the q2 value, and
  # they vary between platforms. A 1e-9 bound still holds nine significant
  # figures, far finer than any difference a diversity estimate can express.
  expect_equal(repeated$hill_shannon[3], 10000, tolerance = 1e-9)
  expect_equal(repeated$inverse_simpson[3], 10000, tolerance = 1e-9)
})

test_that("chunk boundaries preserve the full probability sum", {
  # Each group can contribute more than one chunk of possible counts. At half
  # depth, a symmetric two-group distribution has independently computable H.
  counts <- c(70000, 70000)
  result <- interpolate_hill_diversity(counts, 70000)
  drawn <- seq_len(70000)
  p <- drawn / 70000
  expected_entropy <- -2 * sum(stats::dhyper(drawn, 70000, 70000, 70000) *
                                p * log(p))
  expect_equal(result$hill_shannon, exp(expected_entropy), tolerance = 1e-12)
  expect_lte(result$hill_shannon, 2)
})

test_that("invalid support and extrapolation inputs fail explicitly", {
  for (counts in list(c(1, 0), c(1, -1), c(1, 0.5), c(1, NA), c(1, Inf))) {
    expect_error(interpolate_hill_diversity(counts, 1), "consensus_count")
  }
  expect_error(interpolate_hill_diversity(2^53, 1), "consecutive-integer")
  for (depth in c(-1, 1.5, 4, NA, Inf)) {
    expect_error(interpolate_hill_diversity(c(2, 1), depth), "depths")
  }
  expect_error(summarise_rarefaction(c(2, 1), 0), "target")
  expect_error(summarise_rarefaction(c(2, 1), 2, -1), "minimum")
  expect_error(rarefaction_depth_grid(c(2, 1), 0), "targets")
  expect_error(rarefaction_depth_grid(c(2, 1), grid_points = 1), "grid_points")
})
