# Check input table cells, read support and named argument values.
# Run: python test/run_tests.py --component inputs --level unit

load_test_modules(c("table_io"))

describe("ensure_output_dirs", {
  it("creates the parent directory of every output in a nested list", {
    root <- withr::local_tempdir()
    ensure_output_dirs(list(
      metrics = file.path(root, "metrics", "overlap.tsv"),
      nested  = list(file.path(root, "summary", "overlap_summary.tsv"))
    ))
    expect_true(all(dir.exists(file.path(root, c("metrics", "summary")))))
  })

  it("creates the directory but not the file itself", {
    root <- withr::local_tempdir()
    target <- file.path(root, "metrics", "overlap.tsv")
    ensure_output_dirs(list(target))
    expect_false(file.exists(target))
  })

  it("is quiet when the directory already exists", {
    root <- withr::local_tempdir()
    target <- file.path(root, "metrics", "overlap.tsv")
    ensure_output_dirs(list(target))
    expect_silent(ensure_output_dirs(list(target)))
  })

  it("accepts an empty set of outputs",
     expect_silent(ensure_output_dirs(list())))

  it("returns nothing",
     expect_null(ensure_output_dirs(list())))
})

describe("require_columns", {
  it("passes a frame that carries every required column",
     expect_silent(require_columns(
       data.frame(donor_id = "D1", locus = "TRB", jaccard = NA_real_),
       c("donor_id", "locus"), "the metrics table")))

  it("passes when nothing is required",
     expect_silent(require_columns(data.frame(a = 1), character(0), "anything")))

  it("names every missing column so the caller can see all of them at once", {
    expect_error(
      require_columns(data.frame(donor_id = "D1"),
                       c("donor_id", "locus", "jaccard"), "the metrics table"),
      regexp = "locus, jaccard"
    )
  })

  it("names the context so the caller knows which table is wrong", {
    expect_error(
      require_columns(data.frame(donor_id = "D1"), c("locus"),
                       "the donor overlap summary"),
      regexp = "the donor overlap summary"
    )
  })

  it("reports a frame with no columns at all as missing everything", {
    expect_error(
      require_columns(data.frame(), c("donor_id", "locus"), "an empty table"),
      regexp = "donor_id, locus"
    )
  })
})

# Aggregating stages bind one shard per donor or per library, and a shard written
# against an older schema is missing a column rather than malformed. The binder
# fills that gap with NA and says nothing, so what it does with a partial shard is
# worth pinning: callers that cannot tolerate the gap check each shard as they
# read it instead of relying on the bind to notice.

describe("bind_frames_on_column_union", {
  it("stacks frames that share their columns", {
    bound <- bind_frames_on_column_union(
      list(data.frame(library_id = "G1", richness = 3L),
           data.frame(library_id = "G2", richness = 5L)))
    expect_identical(bound$library_id, c("G1", "G2"))
    expect_identical(bound$richness, c(3L, 5L))
  })

  it("fills a column only one frame carries with NA rather than failing", {
    bound <- bind_frames_on_column_union(
      list(data.frame(library_id = "G1", richness = 3L),
           data.frame(library_id = "G2")))
    expect_identical(bound$richness, c(3L, NA_integer_))
  })

  it("takes the union of the columns in the order the frames present them", {
    bound <- bind_frames_on_column_union(
      list(data.frame(library_id = "G1", richness = 3L),
           data.frame(library_id = "G2", shannon = 1.1)))
    expect_identical(names(bound), c("library_id", "richness", "shannon"))
  })

  it("projects onto the requested columns and drops the rest", {
    bound <- bind_frames_on_column_union(
      list(data.frame(library_id = "G1", richness = 3L, extra = "x")),
      columns = c("library_id", "richness"))
    expect_identical(names(bound), c("library_id", "richness"))
  })

  it("NA-fills a requested column that no frame carries", {
    bound <- bind_frames_on_column_union(
      list(data.frame(library_id = "G1")),
      columns = c("library_id", "shannon"))
    expect_true(is.na(bound$shannon))
  })

  it("ignores frames with no rows", {
    bound <- bind_frames_on_column_union(
      list(data.frame(library_id = character(0), richness = integer(0)),
           data.frame(library_id = "G1", richness = 3L)))
    expect_identical(nrow(bound), 1L)
  })

  it("returns the requested header when every frame is empty", {
    bound <- bind_frames_on_column_union(
      list(data.frame(library_id = character(0))),
      columns = c("library_id", "richness"))
    expect_identical(nrow(bound), 0L)
    expect_identical(names(bound), c("library_id", "richness"))
  })

  it("returns a frame with nothing in it when there is nothing to bind",
     expect_identical(nrow(bind_frames_on_column_union(list())), 0L))
})

# The shared writer fixes null spelling, numeric notation and schema projection.
# Every stage uses these rules, including empty outputs.

describe("write_pipeline_table serialises the same way for every stage", {
  it("writes a null as the empty string, not as the literal NA", {
    path <- withr::local_tempfile()
    write_pipeline_table(
      data.frame(library_id = "G1", threshold = NA_real_, note = NA_character_,
                 stringsAsFactors = FALSE),
      path)
    expect_identical(readLines(path), c("library_id\tthreshold\tnote", "G1\t\t"))
  })

  it("writes a large round count in full rather than as 1e+05", {
    path <- withr::local_tempfile()
    write_pipeline_table(
      data.frame(consensus_count = c(100000, 1e7), stringsAsFactors = FALSE), path)
    expect_identical(readLines(path), c("consensus_count", "100000", "10000000"))
  })

  it("leaves the caller's own scipen alone", {
    path <- withr::local_tempfile()
    before <- getOption("scipen")
    write_pipeline_table(data.frame(consensus_count = 100000), path)
    expect_identical(getOption("scipen"), before)
  })

  it("projects and orders the frame to the released schema", {
    path <- withr::local_tempfile()
    write_pipeline_table(
      data.frame(b = 2, a = 1, dropped = 3), path, columns = c("a", "b"))
    expect_identical(readLines(path), c("a\tb", "1\t2"))
  })

  it("says which column is missing when the schema cannot be met", {
    path <- withr::local_tempfile()
    expect_error(
      write_pipeline_table(data.frame(a = 1), path, columns = c("a", "b")),
      regexp = "b")
  })

  it("writes unquoted rows and returns the output path invisibly", {
    path <- withr::local_tempfile()
    frame <- data.frame(library_id = c("0460", "NA"), note = c("has a space", ""))
    expect_invisible(result <- write_pipeline_table(frame, path))
    expect_identical(result, path)
    expect_identical(readLines(path), c("library_id\tnote", "0460\thas a space", "NA\t"))
    expect_identical(read_pipeline_table(path), frame)
  })

  it("writes a header-only table when no rows are available", {
    path <- withr::local_tempfile()
    write_pipeline_table(data.frame(library_id = character(), count = numeric()), path)
    expect_identical(readLines(path), "library_id\tcount")
  })
})

# Sorting a character vector uses the locale's collation, which can differ
# between runtimes, so without a declared collation the order of released rows
# and the numeric suffix of every clonotype_id would depend on which runtime ran
# the stage. set_pipeline_collation() answers it once; the test helper calls it,
# so this suite sorts the way the runtimes do.

describe("the declared collation", {
  it("is in force in this suite",
     expect_identical(Sys.getlocale("LC_COLLATE"), "C"))

  it("orders by code point, so case and punctuation are not folded away", {
    expect_identical(sort(c("b", "A", "a", "B")), c("A", "B", "a", "b"))
    expect_identical(sort(c("WTS_1", "WTSa", "WTS1")),
                     c("WTS1", "WTS_1", "WTSa"))
  })
})

test_that("path lists trim whitespace and empty components while keeping order", {
  expect_identical(split_path_list(" second.tsv, ,first.tsv,, "), c("second.tsv", "first.tsv"))
  expect_identical(split_path_list("one.tsv"), "one.tsv")
  expect_identical(split_path_list("one.tsv,one.tsv"), rep("one.tsv", 2))
  for (value in c("", " ", ",", " , , ")) {
    expect_error(split_path_list(value), "Expected at least one path")
  }
})

test_that("named values keep order, accept an empty setting and reject malformed pairs", {
  parsed <- parse_named_values("TRB=100, IGH=1000")
  expect_identical(names(parsed), c("TRB", "IGH"))
  expect_equal(unname(parsed), c(100, 1000))
  expect_equal(parse_named_values("IGH=0.5"), c(IGH = 0.5))
  empty <- parse_named_values("")
  expect_identical(names(empty), character(0))
  expect_length(empty, 0L)
  expect_identical(parse_named_values(" , "), empty)
  expect_error(parse_named_values("IGH"), "NAME=value")
  expect_error(parse_named_values("IGH=1=2"), "NAME=value")
  expect_error(parse_named_values("=1"), "NAME=value")
  expect_error(parse_named_values("IGH=many"), "NAME=number")
  expect_error(parse_named_values("IGH=1,IGH=2"), "Repeated name")
})

describe("parse_airr_boolean reads a text AIRR boolean", {
  cases <- list(
    list(token = "T", expected = "true",
         why = "reads the AIRR shorthand T as true"),
    list(token = "F", expected = "false",
         why = "reads the AIRR shorthand F as false"),
    list(token = "TRUE", expected = "true",
         why = "reads a spelled-out TRUE"),
    list(token = "FALSE", expected = "false",
         why = "reads a spelled-out FALSE"),
    list(token = "1", expected = "true",
         why = "reads 1 as true"),
    list(token = "0", expected = "false",
         why = "reads 0 as false"),
    list(token = "", expected = "missing",
         why = "reads a blank cell as missing"),
    list(token = "NA", expected = "invalid",
         why = "reports the literal string NA as invalid, since AIRR forbids it as a null"),
    list(token = NA_character_, expected = "missing",
         why = "reads a true NA as missing"),
    list(token = "yes", expected = "invalid",
         why = "reports an unrecognised token as invalid, never as false")
  )
  for (case in cases) {
    it(case$why, expect_identical(parse_airr_boolean(case$token), case$expected))
  }
})

describe("parse_airr_boolean reads a logical AIRR boolean", {
  cases <- list(
    list(value = TRUE, expected = "true", why = "maps TRUE"),
    list(value = FALSE, expected = "false", why = "maps FALSE"),
    list(value = NA, expected = "missing", why = "maps a logical NA to missing")
  )
  for (case in cases) {
    it(case$why, expect_identical(parse_airr_boolean(case$value), case$expected))
  }
})

describe("as_read_support converts a consensus_count field to a number", {
  cases <- list(
    list(value = "0", expected = 0,
         why = "accepts zero, which is an observation rather than an error"),
    list(value = "2", expected = 2,
         why = "parses a count held as text"),
    list(value = 3, expected = 3,
         why = "passes an already numeric count straight through")
  )
  for (case in cases) {
    it(case$why, expect_identical(as_read_support(case$value), case$expected))
  }
})

# An IMGT-aligned string uses gaps for two different things. The leading dots
# below are numbering spacers the germline shares, so they must stay aligned;
# the dot and dash opposite germline G and G are positions the read never
# covered, and those are not callable mutation opportunities.
