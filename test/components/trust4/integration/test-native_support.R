# Run TRUST4-reported support and yield with malformed-count and empty-library boundaries.
# Run: python test/run_tests.py --component trust4 --level integration

load_test_modules(c("table_io"))

test_that("TRUST4-reported support preserves paired depth and rejects malformed support", {
  directory <- tempfile("native-support-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- list(airr = file.path(directory, "native.tsv"),
                depth = file.path(directory, "depth.tsv"))
  output <- list(library_yield = file.path(directory, "yield.tsv"))
  params <- list(sample_label = "0460")
  records <- data.frame(locus = c("IGH", "TRB"), consensus_count = c("5", "3"))
  write_pipeline_table(records, input$airr)
  depth <- function(library_id, pairs) data.frame(library_id = library_id, n_fastq_file_pairs = 2,
                                                  input_read_pairs = pairs, total_bases = 30000000000)
  write_pipeline_table(depth("0460", 100), input$depth)
  execute <- function(ok = TRUE) run_stage_with_options("summarise_library_yield.R",
    list(airr_file = input$airr, depth_file = input$depth,
         library_id = params$sample_label, output_file = output$library_yield),
    expect_success = ok)
  execute()
  row <- read_pipeline_table(output$library_yield)
  expect_identical(row$library_id, "0460")
  expect_identical(names(row)[1:6], c("library_id", "n_fastq_file_pairs", "input_read_pairs",
                                      "input_bases", "yield_status", "trust4_consensus_count"))
  expect_identical(c(row$n_fastq_file_pairs, row$input_bases), c("2", "30000000000"))
  expect_equal(as.numeric(row$trust4_consensus_count), 8)
  expect_equal(as.numeric(row$input_read_pairs), 100)
  expect_equal(as.numeric(row$igh_yield_per_million_read_pairs), 50000)
  expect_equal(as.numeric(row$trb_yield_per_million_read_pairs), 30000)
  expect_identical(row$yield_status, "available")
  valid_output <- readLines(output$library_yield)
  for (token in c("", "NA", "Inf", "malformed", "-1")) {
    invalid <- records
    invalid$consensus_count[2] <- token
    write_pipeline_table(invalid, input$airr)
    run <- execute(FALSE)
    expect_gt(run$status, 0L)
    expect_match(paste(run$output, collapse = "\n"), "consensus_count.*invalid")
    expect_identical(readLines(output$library_yield), valid_output)
  }
  # Header-only receptor input is measured zero at positive input depth.
  write_pipeline_table(records[FALSE, ], input$airr)
  execute()
  row <- read_pipeline_table(output$library_yield)
  expect_equal(as.numeric(row$trust4_consensus_count), 0)
  expect_equal(as.numeric(row$igh_yield_per_million_read_pairs), 0)
  # Zero input depth leaves yield unavailable and retains a literal NA ID.
  params$sample_label <- "NA"
  write_pipeline_table(depth("NA", 0), input$depth)
  execute()
  row <- read_pipeline_table(output$library_yield)
  expect_identical(row$library_id, "NA")
  expect_identical(row$yield_status, "no_input_read_pairs")
  expect_identical(row$igh_yield_per_million_read_pairs, "")
})
