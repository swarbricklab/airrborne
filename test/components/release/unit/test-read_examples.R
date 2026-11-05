# Protect the base-R import and join examples independently of Python.
source(file.path(Sys.getenv("AIRRBORNE_ROOT"), "examples/read_release.R"))

test_that("string imports keep literal NA and zeroes and metadata cannot multiply rows", {
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path))
  writeLines(c("library_id\ttissue_id\tinput_read_pairs", "0460\tNA\t", "NA\t001\t0"), path)
  headline <- read_rows(path)
  expect_identical(headline$library_id, c("0460", "NA"))
  expect_identical(headline$tissue_id, c("NA", "001"))
  expect_identical(number(headline$input_read_pairs), c(NA_real_, 0))
  metadata <- data.frame(tissue_id = "NA", group = "observed")
  joined <- add_metadata(headline, metadata)
  expect_identical(joined$group, c("observed", NA_character_))
  expect_error(add_metadata(headline, rbind(metadata, metadata)), "Duplicate")
  expect_error(add_metadata(headline, data.frame(tissue_id = "NA", library_id = "wrong")), "conflicts")
  writeLines(c('tissue_id,group', 'NA,"with, comma"'), path)
  expect_identical(read_rows(path, sep = ",")$group, "with, comma")
})

test_that("member representatives retain transferred support exactly once", {
  members <- data.frame(scope_id = "cohort", library_id = "0460",
                        sequence_id = "0460~assemble0_0", member_consensus_count = "12")
  airr <- data.frame(library_id = "0460", sequence_id = c("0460~assemble0_0", "0460~assemble0_1"),
                     is_member_representative = c("T", "F"), member_consensus_count = c("12", "0"),
                     consensus_count = c("10", "2"), sequence = c("ACTG", "ACTA"))
  joined <- join_member_airr(members, airr)
  expect_equal(nrow(joined), 1)
  expect_identical(joined$member_consensus_count, "12")
  expect_identical(joined$source_consensus_count, "10")
})
