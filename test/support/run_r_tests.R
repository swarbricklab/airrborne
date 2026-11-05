# Run selected R test files in separate environments inside the analysis image.
# Each file names its production modules. A failure or skipped required case
# makes the component command fail.
root <- Sys.getenv("AIRRBORNE_ROOT")
stopifnot(nzchar(root), requireNamespace("testthat", quietly = TRUE))
paths <- commandArgs(trailingOnly = TRUE)
stopifnot(length(paths) > 0L)
failed <- FALSE
for (path in paths) {
  cat("\n", sub(paste0(root, "/test/components/"), "", path, fixed = TRUE), "\n", sep = "")
  environment <- new.env(parent = globalenv())
  source(file.path(root, "test/support/helpers.R"), local = environment)
  result <- testthat::test_file(path, env = environment, reporter = "summary")
  cases <- as.data.frame(result)
  if (nrow(cases) == 0L || any(cases$failed > 0L | cases$error | cases$skipped)) {
    failed <- TRUE
  }
}
if (failed) stop("Required R checks failed or were skipped.")
