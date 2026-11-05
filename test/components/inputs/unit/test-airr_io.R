# Check AIRR boolean/null spelling and literal identifier round trips.
# Run: python test/run_tests.py --component inputs --level unit

load_test_modules(c("table_io", "airr_io"))

make_airr_frame <- function() {
  data.frame(
    sequence_id  = c("s1", "s2"),
    productive   = c(TRUE, FALSE),
    rev_comp     = c(NA, TRUE),
    note         = c("has a space", NA),
    stringsAsFactors = FALSE
  )
}

describe("normalise_airr_frame", {
  normalised <- normalise_airr_frame(make_airr_frame())

  it("keeps every column, in the order it was given",
     expect_identical(names(normalised), names(make_airr_frame())))
  it("keeps every row",
     expect_identical(nrow(normalised), 2L))
  it("returns every column as character",
     expect_true(all(vapply(normalised, is.character, logical(1)))))
  it("spells a logical column as the AIRR T and F tokens",
     expect_identical(normalised$productive, c("T", "F")))
  it("writes a missing boolean as the empty AIRR null",
     expect_identical(normalised$rev_comp, c("", "T")))
  it("writes a missing text cell as the empty AIRR null",
     expect_identical(normalised$note, c("has a space", "")))
  it("leaves a present text cell untouched",
     expect_identical(normalised$sequence_id, c("s1", "s2")))
})

describe("normalise_airr_frame on the cells it has to coerce", {
  it("spells an all-missing logical column as all nulls", {
    frame <- data.frame(flag = c(NA, NA), stringsAsFactors = FALSE)
    expect_identical(normalise_airr_frame(frame)$flag, c("", ""))
  })
  it("preserves literal strings that resemble missing-value tokens", {
    frame <- data.frame(text = c("NA", " NA ", "None", "nan", "0460"),
                        stringsAsFactors = FALSE)
    expect_identical(normalise_airr_frame(frame)$text, frame$text)
    path <- withr::local_tempfile()
    write_airr_table(frame, path)
    expect_identical(read_pipeline_table(path)$text, frame$text)
  })
  # Character columns preserve serialized AIRR tokens and identifier strings.
  # Boolean-like text must not be coerced by the writer.
  it("leaves a character column holding boolean-like words unchanged", {
    frame <- data.frame(text = c("TRUE", "maybe"), stringsAsFactors = FALSE)
    expect_identical(normalise_airr_frame(frame)$text, c("TRUE", "maybe"))
  })
})

# A numeric column is turned into text inside normalise_airr_frame(), so it is
# that function, and not just the writer, that has to prefer fixed notation. Two
# million is the smallest round count where the two spellings differ:
# as.character(2e6) is "2e+06" unless the writer's scipen window is open.
describe("normalise_airr_frame on a numeric column", {
  counts <- function() {
    data.frame(sequence_id = c("s1", "s2"),
               consensus_count = c(2e6, 3),
               stringsAsFactors = FALSE)
  }

  it("spells a large count in full rather than in scientific notation",
     expect_identical(normalise_airr_frame(counts())$consensus_count,
                      c("2000000", "3")))

  it("leaves the caller's scipen option as it found it", {
    before <- getOption("scipen")
    normalise_airr_frame(counts())
    expect_identical(getOption("scipen"), before)
  })

  it("puts the full count in the bytes on disk", {
    path <- withr::local_tempfile()
    write_airr_table(counts(), path)
    expect_identical(readLines(path),
                     c("sequence_id\tconsensus_count",
                       "s1\t2000000", "s2\t3"))
  })
})

describe("normalise_airr_frame on a frame with no rows", {
  empty <- normalise_airr_frame(
    data.frame(sequence_id = character(0), productive = logical(0),
               stringsAsFactors = FALSE)
  )

  it("keeps the column names",
     expect_identical(names(empty), c("sequence_id", "productive")))
  it("still has no rows",
     expect_identical(nrow(empty), 0L))
})

describe("write_airr_table", {
  it("writes AIRR spellings into the bytes on disk", {
    path <- withr::local_tempfile()
    write_airr_table(make_airr_frame(), path)
    expect_identical(
      readLines(path),
      c("sequence_id\tproductive\trev_comp\tnote",
        "s1\tT\t\thas a space",
        "s2\tF\tT\t")
    )
  })


  it("round trips through read.delim with the AIRR tokens intact", {
    path <- withr::local_tempfile()
    write_airr_table(make_airr_frame(), path)
    expect_identical(
      read.delim(path, colClasses = "character", na.strings = character(0)),
      data.frame(
        sequence_id  = c("s1", "s2"),
        productive   = c("T", "F"),
        rev_comp     = c("", "T"),
        note         = c("has a space", ""),
        stringsAsFactors = FALSE
      )
    )
  })

  it("writes a frame with no rows as a header line and nothing else", {
    path <- withr::local_tempfile()
    write_airr_table(
      data.frame(sequence_id = character(0), productive = logical(0),
                 stringsAsFactors = FALSE),
      path
    )
    expect_identical(readLines(path), "sequence_id\tproductive")
  })

  it("reads a header-only file back as a frame with no rows", {
    path <- withr::local_tempfile()
    write_airr_table(
      data.frame(sequence_id = character(0), productive = logical(0),
                 stringsAsFactors = FALSE),
      path
    )
    round_trip <- read.delim(path, colClasses = "character",
                             na.strings = character(0))
    expect_identical(dim(round_trip), c(0L, 2L))
  })

  it("returns the path it wrote, invisibly", {
    path <- withr::local_tempfile()
    expect_invisible(write_airr_table(make_airr_frame(), path))
  })

  it("returns the path it was given", {
    path <- withr::local_tempfile()
    expect_identical(write_airr_table(make_airr_frame(), path), path)
  })
})

PRESENT <- c("ACGT", "0", "na_token_inside_a_word", "N/A", "NAN")
ABSENT <- c(NA_character_, "", "   ", "\t")
NA_TOKEN <- c("NA", "na", "Na", " NA ")

describe("is_airr_null reads an AIRR null", {
  it("treats a value as present", expect_false(any(is_airr_null(PRESENT))))
  it("treats NA and blank as absent", expect_true(all(is_airr_null(ABSENT))))
  it("treats the literal NA as present, because it is an ordinary string",
     expect_false(any(is_airr_null(NA_TOKEN))))
  it("handles an empty vector",
     expect_identical(is_airr_null(character(0)), logical(0)))
})

describe("no writer produces a token is_airr_null() would miss", {
  # Nulls arrive as actual missing values or empty strings read from AIRR tables.
  frame <- data.frame(
    kept = c("ACGT", "ACGT"),
    real_na = c(NA_character_, "present"),
    passed_through = c("", "present"),
    stringsAsFactors = FALSE
  )

  round_trip <- function(writer) {
    path <- withr::local_tempfile(fileext = ".tsv")
    writer(frame, path)
    read_pipeline_table(path)
  }

  it("holds for the metric writer", {
    back <- round_trip(write_pipeline_table)
    expect_true(all(is_airr_null(back$real_na[1]), is_airr_null(back$passed_through[1])))
    expect_false(any(is_airr_null(back$real_na[2]), is_airr_null(back$passed_through[2])))
  })

  it("holds for the AIRR writer", {
    back <- round_trip(write_airr_table)
    expect_true(all(is_airr_null(back$real_na[1]), is_airr_null(back$passed_through[1])))
    expect_false(any(is_airr_null(back$real_na[2]), is_airr_null(back$passed_through[2])))
  })

  it("means the two writers spell a null identically", {
    metric <- withr::local_tempfile(fileext = ".tsv")
    airr <- withr::local_tempfile(fileext = ".tsv")
    write_pipeline_table(frame, metric)
    write_airr_table(frame, airr)
    expect_identical(readLines(metric), readLines(airr))
  })
})
