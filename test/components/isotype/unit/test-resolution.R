# Check canonical isotype class, subclass and switch-state resolution.
# Run: python test/run_tests.py --component isotype --level unit

load_test_modules(c("table_io", "gene_calls", "isotype_resolution"))

# ------------------------------------------------------------------ isotype ---

describe("resolve_isotype at subclass level", {
  cases <- list(
    list(call = "IGHG1*01,IGHG1*03", expected = "IgG1",
         why = "two alleles of one subclass resolve to it"),
    list(call = "IGHG1,IGHG3", expected = NA_character_,
         why = "two different subclasses are ambiguous at subclass level"),
    list(call = "IGHM", expected = "IgM",
         why = "IGHM is its own subclass")
  )
  for (case in cases) {
    it(case$why, expect_identical(resolve_isotype(case$call, "subclass"),
                                  case$expected))
  }
})

describe("resolve_isotype at class level", {
  cases <- list(
    list(call = "IGHG1,IGHG3", expected = "IgG",
         why = "two IgG subclasses agree on the class"),
    list(call = "IGHM,IGHD", expected = NA_character_,
         why = "IGHM and IGHD are separate classes, so the call is ambiguous"),
    list(call = "IGHM", expected = "IgM",
         why = "IGHM is kept as its own class rather than merged"),
    list(call = "IGHG1,IGHA1", expected = NA_character_,
         why = "two different classes are ambiguous")
  )
  for (case in cases) {
    it(case$why, expect_identical(resolve_isotype(case$call, "class"),
                                  case$expected))
  }
})

describe("resolve_isotype at switch level", {
  cases <- list(
    list(call = "IGHM,IGHD", expected = "unswitched",
         why = "IgM and IgD are both unswitched"),
    list(call = "IGHG1,IGHA1", expected = "switched",
         why = "IgG and IgA are both class-switched"),
    list(call = "IGHM,IGHG1", expected = NA_character_,
         why = "a call spanning both sides of the switch is ambiguous")
  )
  for (case in cases) {
    it(case$why, expect_identical(resolve_isotype(case$call, "switch"),
                                  case$expected))
  }
})

# Dropping an unrecognised token instead of rejecting the whole call would turn
# annotation noise into false isotype evidence.
describe("resolve_isotype rejects a call containing anything unrecognised", {
  cases <- list(
    list(call = "IGHG1,FOO", why = "one unrecognised token alongside a real one"),
    list(call = "FOO", why = "an entirely unrecognised call")
  )
  for (case in cases) {
    it(case$why, expect_identical(resolve_isotype(case$call, "switch"),
                                  NA_character_))
  }
})

describe("resolve_isotype accepts only canonical subclasses", {
  cases <- list(
    list(call = "IGHG", why = "a bare class with no subclass number"),
    list(call = "IGHG5", why = "an IgG subclass number that does not exist"),
    list(call = "IGHA9", why = "an IgA subclass number that does not exist")
  )
  for (case in cases) {
    it(case$why, expect_identical(resolve_isotype(case$call, "class"),
                                  NA_character_))
  }
})

# The allotype rule is deliberately narrow: only the one suffix known to be an
# allotype is normalised, so an unfamiliar suffix stays uncallable.
describe("resolve_isotype normalises only the known IGHG4A allotype", {
  cases <- list(
    list(call = "IGHG4A*01", expected = "IgG4",
         why = "IGHG4A is an allotype of IGHG4"),
    list(call = "IGHG1P*01", expected = NA_character_,
         why = "IGHG1P is not a known allotype and stays uncallable")
  )
  for (case in cases) {
    it(case$why, expect_identical(resolve_isotype(case$call, "subclass"),
                                  case$expected))
  }
})

test_that("CANONICAL_IGH_ISOTYPES is the shared output category set", {
  expect_identical(CANONICAL_IGH_ISOTYPES,
                   c("IgM", "IgD", "IgG1", "IgG2", "IgG3", "IgG4",
                     "IgA1", "IgA2", "IgE"))
})

test_that("ANALYSIS_IGH_ISOTYPE_GROUPS fixes the released group order", {
  expect_identical(ANALYSIS_IGH_ISOTYPE_GROUPS,
                   c("all_IGH", "IgM", "IgG", "IgA"))
})

# IgD and IgE stay inside all-IGH analyses but get no stratified endpoint of
# their own, so they resolve to no group here.
describe("resolve_analysis_igh_isotype", {
  cases <- list(
    list(call = "IGHM*01", expected = "IgM", why = "an IgM call becomes IgM"),
    list(call = "IGHG1*01,IGHG1*03", expected = "IgG",
         why = "two alleles of one IgG subclass become IgG"),
    list(call = "IGHA2", expected = "IgA", why = "an IgA subclass becomes IgA"),
    list(call = "IGHD", expected = NA_character_,
         why = "IgD has no stratified group"),
    list(call = "IGHE", expected = NA_character_,
         why = "IgE has no stratified group"),
    list(call = "IGHM,IGHG1", expected = NA_character_,
         why = "an ambiguous call belongs to no group")
  )
  for (case in cases) {
    it(case$why, expect_identical(resolve_analysis_igh_isotype(case$call),
                                  case$expected))
  }
})
