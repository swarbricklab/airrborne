# Check class/subclass support fractions, callable populations and ratios.
# Run: python test/run_tests.py --component isotype --level unit

load_test_modules(c("table_io", "gene_calls", "isotype_resolution", "isotype_metrics"))

test_that("class ties recover support without assigning an unknown subclass", {
  records <- data.frame(
    qualified_c_call = c("IGHM*01", "IGHD*01", "IGHG1*01", "IGHG1*01,IGHG3*01",
               "IGHA1*01", "IGHA1*01,IGHA2*01", "IGHE*01",
               "IGHG1*01,IGHA1*01", "IGHGP", ""),
    consensus_count = c(10, 5, 8, 12, 6, 4, 5, 7, 3, 2)
  )
  result <- build_isotype_metrics(records)
  composition <- result$composition
  summary <- result$summary
  expect_identical(names(composition), ISOTYPE_METRIC_COLUMNS)
  expect_identical(names(summary), ISOTYPE_SUMMARY_COLUMNS)
  expect_equal(nrow(composition), 14L)
  expect_equal(nrow(summary), 1L)
  classes <- composition[composition$isotype_level == "class", ]
  subclasses <- composition[composition$isotype_level == "subclass", ]
  expect_identical(classes$isotype, c("IgM", "IgD", "IgG", "IgA", "IgE"))
  expect_equal(classes$consensus_count, c(10, 5, 20, 10, 5))
  expect_equal(classes$support_fraction, c(10, 5, 20, 10, 5) / 50)
  expect_equal(subclasses$consensus_count, c(10, 5, 8, 0, 0, 0, 6, 0, 5))
  expect_equal(subclasses$support_fraction, subclasses$consensus_count / 34)
  expect_equal(summary$common_qc_consensus_count, 62)
  expect_equal(summary$class_callable_consensus_count, 50)
  expect_equal(summary$class_uncallable_consensus_count, 12)
  expect_equal(summary$subclass_callable_consensus_count, 34)
  expect_equal(summary$class_only_consensus_count, 16)
  expect_equal(summary$igm_consensus_count, 10)
  expect_equal(summary$igg_consensus_count, 20)
  expect_equal(summary$iga_consensus_count, 10)
  expect_equal(summary$log2_igg_iga, 1)
  expect_equal(summary$igg_iga_status, "available")
  for (table in list(classes, subclasses)) {
    expect_equal(sum(table$consensus_count), table$callable_consensus_count[1])
    expect_equal(sum(table$support_fraction), 1)
  }
  # The summary keeps the common-QC total that the uncallable support completes.
  expect_equal(classes$callable_consensus_count[1] + summary$class_uncallable_consensus_count,
               summary$common_qc_consensus_count)
  expect_equal(summary$subclass_callable_consensus_count +
                 summary$class_only_consensus_count,
               summary$class_callable_consensus_count)
})

test_that("classes without subdivisions occur in both documented resolutions", {
  records <- data.frame(qualified_c_call = c("IGHM", "IGHD", "IGHE"),
                         consensus_count = c(2, 3, 5))
  result <- build_isotype_metrics(records)
  for (isotype in c("IgM", "IgD", "IgE")) {
    rows <- result$composition[result$composition$isotype == isotype, ]
    expect_identical(rows$isotype_level, c("class", "subclass"))
    expect_equal(rows$consensus_count[1], rows$consensus_count[2])
    expect_equal(rows$support_fraction[1], rows$support_fraction[2])
  }
  expect_equal(result$summary$class_only_consensus_count, 0)
  expect_equal(result$summary$igg_iga_status, "both_zero")
})

test_that("cross-class and unfamiliar calls stay unresolved at both resolutions", {
  records <- data.frame(
    qualified_c_call = c("IGHG1,IGHA1", "IGHG1,IGHGP", "unrecognised", "", NA_character_),
    consensus_count = c(3, 5, 7, 11, 13)
  )
  result <- build_isotype_metrics(records)
  expect_equal(result$composition$consensus_count, rep(0, 14))
  expect_true(all(is.na(result$composition$support_fraction)))
  expect_equal(result$composition$callable_consensus_count, rep(0, 14))
  expect_equal(result$summary$common_qc_consensus_count, 39)
  expect_equal(result$summary$class_callable_consensus_count, 0)
  expect_equal(result$summary$class_only_consensus_count, 0)
  expect_equal(result$summary$igg_iga_status, "no_class_callable_evidence")
  expect_true(is.na(result$summary$log2_igg_iga))
})

test_that("every unavailable IgG:IgA reason retains its support inputs", {
  cases <- list(
    list(qualified_c_call = c("IGHG1", "IGHA1"), counts = c(8, 2),
         status = "available", ratio = 2),
    list(qualified_c_call = "IGHA1", counts = 3, status = "zero_igg", ratio = NA_real_),
    list(qualified_c_call = "IGHG1", counts = 3, status = "zero_iga", ratio = NA_real_),
    list(qualified_c_call = "IGHM", counts = 3, status = "both_zero", ratio = NA_real_),
    list(qualified_c_call = "", counts = 3, status = "no_class_callable_evidence", ratio = NA_real_)
  )
  for (case in cases) {
    result <- build_isotype_metrics(data.frame(qualified_c_call = case$qualified_c_call,
                                                consensus_count = case$counts))
    expect_equal(result$summary$igg_iga_status, case$status)
    expect_equal(result$summary$log2_igg_iga, case$ratio)
    expect_equal(result$summary$common_qc_consensus_count, sum(case$counts))
    classes <- result$composition[result$composition$isotype_level == "class", ]
    expect_equal(result$summary$igg_consensus_count,
                 classes$consensus_count[classes$isotype == "IgG"])
    expect_equal(result$summary$iga_consensus_count,
                 classes$consensus_count[classes$isotype == "IgA"])
  }
})

test_that("empty and zero-support inputs retain zero counts and missing fractions", {
  for (records in list(
    data.frame(qualified_c_call = character(), consensus_count = numeric()),
    data.frame(qualified_c_call = "IGHM", consensus_count = 0)
  )) {
    result <- build_isotype_metrics(records)
    expect_equal(nrow(result$composition), 14L)
    expect_equal(result$composition$consensus_count, rep(0, 14))
    expect_true(all(is.na(result$composition$support_fraction)))
    counts <- c("common_qc_consensus_count", "class_callable_consensus_count",
                "class_uncallable_consensus_count", "subclass_callable_consensus_count",
                "class_only_consensus_count", "igm_consensus_count",
                "igg_consensus_count", "iga_consensus_count")
    expect_true(all(unlist(result$summary[counts]) == 0))
    expect_equal(result$summary$igg_iga_status, "no_class_callable_evidence")
    expect_true(is.na(result$summary$log2_igg_iga))
  }
})

test_that("qualified gene and allele conventions follow the existing resolver", {
  records <- data.frame(qualified_c_call = c("IGHG4A*01", "IGHG1*01,IGHG1*02", "IGHA1,IGHA2"),
                         consensus_count = c(2, 3, 5))
  result <- build_isotype_metrics(records)
  subclasses <- result$composition[result$composition$isotype_level == "subclass", ]
  expect_equal(subclasses$consensus_count[subclasses$isotype == "IgG4"], 2)
  expect_equal(subclasses$consensus_count[subclasses$isotype == "IgG1"], 3)
  expect_equal(result$summary$igg_consensus_count, 5)
  expect_equal(result$summary$iga_consensus_count, 5)
  expect_equal(result$summary$class_only_consensus_count, 5)
  expect_equal(result$summary$log2_igg_iga, 0)
})

test_that("the finite log ratio does not overflow at unequal finite support", {
  result <- build_isotype_metrics(data.frame(qualified_c_call = c("IGHG1", "IGHA1"),
                                               consensus_count = c(1e300, 1e-100)))
  expect_equal(result$summary$igg_iga_status, "available")
  expect_true(is.finite(result$summary$log2_igg_iga))
  expect_equal(result$summary$log2_igg_iga, 400 * log2(10), tolerance = 1e-12)
})

test_that("invalid evidence inputs fail without excluding their support silently", {
  expect_error(build_isotype_metrics(data.frame(qualified_c_call = "IGHM")), "consensus_count")
  for (support in c(NA, Inf, -1)) {
    expect_error(build_isotype_metrics(data.frame(qualified_c_call = "IGHM", consensus_count = support)),
                 "consensus_count")
  }
})

test_that("callable classes and subclasses retain their own support populations", {
  records <- data.frame(qualified_c_call = c("IGHM*01", "IGHG1,IGHG3", "IGHA1"),
                        consensus_count = c(3, 5, 2))
  metrics <- build_isotype_metrics(records)
  expect_identical(!is.na(metrics$resolved$subclass), c(TRUE, FALSE, TRUE))
  expect_identical(metrics$resolved$subclass, c("IgM", NA_character_, "IgA1"))
  expect_equal(metrics$summary$subclass_callable_consensus_count, 5, tolerance = 1e-12)
  expect_equal(metrics$summary$subclass_callable_consensus_count /
                 metrics$summary$common_qc_consensus_count, 0.5, tolerance = 1e-12)
  records$qualified_c_call <- c("IGHM", "IGHG4A", "IGHA1*02")
  expect_identical(build_isotype_metrics(records)$resolved$subclass, c("IgM", "IgG4", "IgA1"))
  records$qualified_c_call <- c("", "IGHV3-23", "IGHM,IGHA1")
  expect_true(all(is.na(build_isotype_metrics(records)$resolved$subclass)))
})

# Subclass fractions use callable subclass support (24), while total IGH
# support (100) includes the uncalled record.
test_that("composition sums repeated subclasses and retains canonical empty categories", {
  records <- data.frame(qualified_c_call = c("IGHG1", "IGHM", "IGHG1", "IGHA1", "IGHG3", ""),
                        consensus_count = c(10, 4, 5, 4, 1, 76))
  result <- build_isotype_metrics(records)$composition
  result <- result[result$isotype_level == "subclass", ]
  expect_identical(result$isotype, CANONICAL_IGH_ISOTYPES)
  expect_equal(nrow(result), 9L)
  expect_equal(result$consensus_count[result$isotype == "IgG1"], 15, tolerance = 1e-12)
  expect_equal(result$consensus_count[result$isotype == "IgG2"], 0, tolerance = 1e-12)
  expect_equal(result$support_fraction[result$isotype == "IgG1"], 15 / 24, tolerance = 1e-12)
  expect_equal(sum(result$support_fraction), 1, tolerance = 1e-12)
  expect_equal(result$callable_consensus_count, rep(24, 9), tolerance = 1e-12)
})

test_that("single, complete and fractional compositions retain their support", {
  single <- build_isotype_metrics(data.frame(qualified_c_call = "IGHE", consensus_count = 3))
  expect_equal(single$composition$support_fraction[single$composition$isotype == "IgE"], c(1, 1), tolerance = 1e-12)
  calls <- c("IGHM", "IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGHE")
  flat <- build_isotype_metrics(data.frame(qualified_c_call = calls, consensus_count = 1))
  rows <- flat$composition[flat$composition$isotype_level == "subclass", ]
  expect_identical(rows$isotype, CANONICAL_IGH_ISOTYPES)
  expect_equal(rows$consensus_count, rep(1, 9), tolerance = 1e-12)
  fractional <- build_isotype_metrics(data.frame(qualified_c_call = c("IGHM", "IGHA1"), consensus_count = c(4.2, 10.5)))
  expect_equal(fractional$summary$subclass_callable_consensus_count, 14.7, tolerance = 1e-12)
  expect_equal(fractional$summary$common_qc_consensus_count, 14.7, tolerance = 1e-12)
})
