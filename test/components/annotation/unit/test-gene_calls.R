# Check locus-compatible gene calls, complete-reference entities and shared TRA/TRD V genes.
# Run: python test/run_tests.py --component annotation --level unit

load_test_modules(c("table_io", "gene_calls", "reannotated_airr", "record_qc"))

# -------------------------------------- segment-aware gene-call compatibility ---

describe("gene_call_locus_compatible on V calls", {
  cases <- list(
    list(call = "IGHV1-2*01", locus = "IGH", expected = TRUE,
         why = "accepts a plain heavy-chain V gene"),
    list(call = "TRAV38-2/DV8*01", locus = "TRA", expected = TRUE,
         why = "accepts a V gene IMGT shares between TRA and TRD, read as TRA"),
    list(call = "TRAV38-2/DV8*01", locus = "TRD", expected = TRUE,
         why = "accepts the same shared V gene read as TRD"),
    list(call = "TRDV1*01", locus = "TRD", expected = TRUE,
         why = "accepts a V gene belonging to TRD alone"),
    list(call = "TRGVA*01", locus = "TRG", expected = TRUE,
         why = "accepts TRGVA, an IMGT V gene named without a number"),
    list(call = "IGHJ4*02", locus = "IGH", expected = FALSE,
         why = "rejects a J gene offered as a V call at the right locus"),
    list(call = "IGKV1-5*01,IGLV1-44*01", locus = "IGK", expected = FALSE,
         why = "rejects a call spanning two loci, because every gene must fit"),
    list(call = "", locus = "IGH", expected = FALSE,
         why = "rejects an absent V call"),
    list(call = "IGHD3-10*01", locus = "IGH", expected = FALSE,
         why = "rejects a D segment offered as a V call")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      gene_call_locus_compatible(case$call, case$locus, "V"), case$expected))
  }
})

describe("gene_call_locus_compatible on J calls", {
  cases <- list(
    list(call = "IGHJ4*02", locus = "IGH", expected = TRUE,
         why = "accepts a plain heavy-chain J gene"),
    list(call = "TRAJ1*01", locus = "TRA", expected = TRUE,
         why = "accepts a TRA J gene"),
    list(call = "TRDJ1*01", locus = "TRD", expected = TRUE,
         why = "accepts a TRD J gene"),
    list(call = "TRBJ2-7*02", locus = "TRB", expected = TRUE,
         why = "accepts a TRB J gene"),
    list(call = "IGHV1-2*01", locus = "IGH", expected = FALSE,
         why = "rejects a V gene offered as a J call"),
    list(call = "TRDJ1*01", locus = "TRA", expected = FALSE,
         why = "rejects a TRD J gene on TRA, since J genes are not shared the way the V genes are"),
    list(call = "IGHJ4*01,IGHJ4*02", locus = "IGH", expected = TRUE,
         why = "accepts two alleles of one J gene"),
    list(call = "", locus = "IGH", expected = FALSE,
         why = "rejects an absent J call")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      gene_call_locus_compatible(case$call, case$locus, "J"), case$expected))
  }
})

describe("gene_call_locus_compatible on D calls", {
  cases <- list(
    list(call = "IGHD3-10*01", locus = "IGH", expected = TRUE,
         why = "accepts a heavy-chain D segment"),
    list(call = "IGHD1/OR15-1a*01", locus = "IGH", expected = TRUE,
         why = "accepts an orphan heavy-chain D segment"),
    list(call = "TRBD1*01", locus = "TRB", expected = TRUE,
         why = "accepts a TRB D segment"),
    list(call = "TRDD2*01", locus = "TRD", expected = TRUE,
         why = "accepts a TRD D segment"),
    # IGHD is both the delta constant gene and the prefix of every heavy-chain D
    # segment, which is why a true D segment has to carry a digit after the D.
    list(call = "IGHD*01", locus = "IGH", expected = FALSE,
         why = "rejects IGHD the constant gene, which is not a D segment"),
    list(call = "TRDD2*01", locus = "TRB", expected = FALSE,
         why = "rejects a D segment from another locus"),
    list(call = "IGHD3-10*01", locus = "IGK", expected = FALSE,
         why = "rejects a D segment on a locus that joins V straight to J"),
    list(call = "IGHD2-2*01,IGHD2-8*01", locus = "IGH", expected = TRUE,
         why = "accepts two D genes when both belong to the locus"),
    list(call = "", locus = "IGH", expected = FALSE,
         why = "reports an absent D call as incompatible, leaving absence for the caller to handle")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      gene_call_locus_compatible(case$call, case$locus, "D"), case$expected))
  }
})

describe("gene_call_locus_compatible on C calls", {
  cases <- list(
    list(call = "IGHG1*01", locus = "IGH", expected = TRUE,
         why = "accepts a heavy-chain constant gene named by its isotype letter"),
    list(call = "IGKC*01", locus = "IGK", expected = TRUE,
         why = "accepts the kappa constant gene"),
    list(call = "IGLC3*03", locus = "IGL", expected = TRUE,
         why = "accepts a lambda constant gene"),
    list(call = "IGLJ-C/OR18*01", locus = "IGL", expected = TRUE,
         why = "accepts the IGLJ-C/OR18 orphan the pinned IMGT C database contains"),
    list(call = "TRAC*01", locus = "TRA", expected = TRUE,
         why = "accepts the TRA constant gene"),
    list(call = "TRDC*01", locus = "TRD", expected = TRUE,
         why = "accepts the TRD constant gene"),
    list(call = "TRBC1*03", locus = "TRB", expected = TRUE,
         why = "accepts a numbered TRB constant gene"),
    list(call = "TRGC1*01", locus = "TRG", expected = TRUE,
         why = "accepts a numbered TRG constant gene"),
    list(call = "IGHG1*01", locus = "IGK", expected = FALSE,
         why = "rejects a heavy-chain constant gene on a light-chain record"),
    list(call = "IGHJ4*02", locus = "IGH", expected = FALSE,
         why = "rejects a J gene offered as a C call"),
    list(call = "IGHD3-10*01", locus = "IGH", expected = FALSE,
         why = "rejects a D segment offered as a C call")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      gene_call_locus_compatible(case$call, case$locus, "C"), case$expected))
  }
})

# A gene needs a complete IMGT-style name, not merely the right prefix,
# otherwise annotation noise would be read as a locus-consistent call.
describe("gene_call_locus_compatible rejects prefix-only and malformed calls", {
  cases <- list(
    list(call = "IGHV", locus = "IGH", segment = "V",
         why = "a V prefix with no gene number"),
    list(call = "IGHVgarbage", locus = "IGH", segment = "V",
         why = "a V prefix followed by nonsense"),
    list(call = "TRAVx/DVy", locus = "TRA", segment = "V",
         why = "a shared TRA/TRD V name with no gene numbers"),
    list(call = "IGHJ", locus = "IGH", segment = "J",
         why = "a J prefix with no gene number"),
    list(call = "IGHJgarbage", locus = "IGH", segment = "J",
         why = "a J prefix followed by nonsense"),
    list(call = "IGKCgarbage", locus = "IGK", segment = "C",
         why = "a kappa constant name with trailing nonsense"),
    list(call = "TRACgarbage", locus = "TRA", segment = "C",
         why = "a TRA constant name with trailing nonsense")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      gene_call_locus_compatible(case$call, case$locus, case$segment), FALSE))
  }
})

# ------------------------------------------------------- canonical gene sets ---

describe("parse_vj_call_set drops allele suffixes and preserves gene locations", {
  cases <- list(
    list(call = "IGHV3-23*01,IGHV3-23*04,IGHV3-23D*01", expected = "IGHV3-23,IGHV3-23D",
         why = "drops allele suffixes and preserves IGHV3-23D"),
    list(call = "IGKV1-39*01,IGKV1D-39*01", expected = "IGKV1-39,IGKV1D-39",
         why = "preserves the distinct IGKV1D-39 location"),
    list(call = "IGHD3-10*01", expected = "IGHD3-10",
         why = "leaves a true heavy-chain D segment's D in place"),
    list(call = "TRBD1*01", expected = "TRBD1",
         why = "leaves a true TRB D segment's D in place"),
    list(call = "TRAV38-2/DV8*01", expected = "TRAV38-2/DV8",
         why = "leaves the shared TRA/TRD V gene name intact"),
    list(call = "V1*01|V2*01", expected = "V1,V2",
         why = "accepts the pipe delimiter and returns a sorted set"),
    list(call = "", expected = "",
         why = "returns an empty string for a blank call"),
    list(call = NA_character_, expected = "",
         why = "returns an empty string for a missing call")
  )
  for (case in cases) {
    it(case$why, expect_identical(parse_vj_call_set(case$call), case$expected))
  }
})

# The intersection answers the question a tally of the members' genes cannot:
# whether one gene is consistent with every member. Members {A,B} and {A,B} name
# the same two genes between them as members {A} and {B}, but only the first pair
# has a gene in common.
describe("intersect_vj_call_sets over the calls of one group", {
  it("keeps a gene every member names", {
    expect_identical(
      intersect_vj_call_sets(c("IGHV3-23*01", "IGHV1-8*01,IGHV3-23*04")),
      "IGHV3-23"
    )
  })
  it("returns nothing where the members name no gene in common", {
    expect_identical(
      intersect_vj_call_sets(c("IGHV1-2*01", "IGHV4-59*01")), "")
  })
  it("separates two groups the gene count cannot tell apart", {
    same <- c("IGHV1-2*01,IGHV4-59*01", "IGHV1-2*01,IGHV4-59*01")
    split <- c("IGHV1-2*01", "IGHV4-59*01")
    expect_identical(count_vj_genes(same), count_vj_genes(split))
    expect_identical(intersect_vj_call_sets(same), "IGHV1-2,IGHV4-59")
    expect_identical(intersect_vj_call_sets(split), "")
  })
  it("returns several shared genes sorted", {
    expect_identical(
      intersect_vj_call_sets(c("IGHV4-59*01,IGHV1-2*01,IGHV3-23*01",
                               "IGHV3-23*02,IGHV1-2*03")),
      "IGHV1-2,IGHV3-23"
    )
  })
  it("preserves duplicated-location names when intersecting calls", {
    expect_identical(
      intersect_vj_call_sets(c("IGKV1-39*01", "IGKV1D-39*02")), "")
  })
  it("returns the member's own genes when there is only one member", {
    expect_identical(
      intersect_vj_call_sets("IGHV3-23*01,IGHV3-30*01"), "IGHV3-23,IGHV3-30")
  })
  # A member that names no gene cannot confirm one, so it empties the
  # intersection. The gene count ignores it instead, and that asymmetry is
  # deliberate.
  it("is emptied by a member with no call, where the gene count ignores one", {
    expect_identical(intersect_vj_call_sets(c("IGHV3-23*01", "")), "")
    expect_identical(intersect_vj_call_sets(c("IGHV3-23*01", NA_character_)), "")
    expect_identical(count_vj_genes(c("IGHV3-23*01", "")), 1L)
  })
  it("returns one string, not one per input call", {
    expect_length(intersect_vj_call_sets(c("IGHV3-23*01", "IGHV3-23*02")), 1L)
  })
  it("returns an empty string for no members at all", {
    expect_identical(intersect_vj_call_sets(character(0)), "")
  })
  it("never names a gene no member was called with", {
    for (calls in list(c("IGHV3-23*01", "IGHV1-8*01,IGHV3-23*04"),
                       c("IGHV1-2*01", "IGHV4-59*01"),
                       c("IGKV1-39*01", "IGKV1D-39*02"),
                       c("IGHV3-23*01", ""))) {
      shared <- .split_calls(intersect_vj_call_sets(calls))
      named <- unlist(lapply(calls, .gene_set), use.names = FALSE)
      expect_true(all(shared %in% named))
    }
  })
})

# Two genes IgBLAST could not separate are usually paralogs within one family, so
# the family count is what distinguishes an unresolvable duplicate from members
# called to genuinely different V genes.
describe("count_vj_gene_families", {
  it("counts one family for several genes that share it", {
    expect_identical(
      count_vj_gene_families(c("IGHV3-23*01", "IGHV3-30*01,IGHV3-33*01")), 1L)
  })
  it("counts each distinct family once", {
    expect_identical(
      count_vj_gene_families(c("IGHV1-2*01", "IGHV4-59*01", "IGHV1-69*01")), 2L)
  })
  it("counts a family named inside an ambiguous call", {
    expect_identical(count_vj_gene_families("IGHV1-2*01,IGHV4-59*01"), 2L)
  })
  it("treats a duplicate-D alias as its base gene's family", {
    expect_identical(count_vj_gene_families(c("IGHV3-23*01", "IGHV3-23D*01")), 1L)
  })
  it("counts a family for a gene name carrying no number suffix", {
    expect_identical(count_vj_gene_families("IGHV6*01"), 1L)
  })
  # The gene-level count is what separates two genes of one family from one gene,
  # which is exactly the distinction the family count cannot make.
  it("collapses to one family what the gene count reports as two", {
    calls <- c("IGHV3-23*01", "IGHV3-30*01")
    expect_identical(count_vj_gene_families(calls), 1L)
    expect_identical(count_vj_genes(calls), 2L)
  })
  # Heavy-chain J names carry no family level, so each J gene counts as its own.
  # That is why the released column reports V families and not J ones.
  it("counts each J gene separately, there being no J family to collapse to", {
    expect_identical(count_vj_gene_families(c("IGHJ4*02", "IGHJ6*02")), 2L)
    expect_identical(count_vj_gene_families(c("IGHJ4*01", "IGHJ4*02")), 1L)
  })
  it("counts nothing where no call names a gene", {
    expect_identical(count_vj_gene_families(c("", NA_character_)), 0L)
    expect_identical(count_vj_gene_families(character(0)), 0L)
  })
})

describe("count_vj_genes", {
  it("counts each distinct gene once across the members' calls", {
    expect_identical(count_vj_genes(c("IGHV1-2*01", "IGHV4-59*01")), 2L)
    expect_identical(count_vj_genes(c("IGHV3-23*01", "IGHV3-23*04")), 1L)
  })
  it("counts the genes named inside one ambiguous call", {
    expect_identical(count_vj_genes("IGHV1-2*01,IGHV4-59*01"), 2L)
  })
  it("counts a duplicated gene as a separate location", {
    expect_identical(count_vj_genes(c("IGHV3-23*01", "IGHV3-23D*01")), 2L)
  })
  it("counts nothing where no call names a gene", {
    expect_identical(count_vj_genes(c("", NA_character_)), 0L)
    expect_identical(count_vj_genes(character(0)), 0L)
  })
})

describe("canonical_gene_name", {
  it("strips the allele suffix and preserves the D marker", {
    expect_identical(
      canonical_gene_name(c("IGHV3-23*01", "IGHV3-23D*01", "IGKV1D-39*02",
                            "IGHD3-10*01", "TRAV38-2/DV8*01", "IGHJ4")),
      c("IGHV3-23", "IGHV3-23D", "IGKV1D-39", "IGHD3-10", "TRAV38-2/DV8", "IGHJ4"))
  })
  it("returns nothing for no tokens",
     expect_identical(canonical_gene_name(character(0)), character(0)))
})

# ------------------------------------------------- germline entity of a call ---

# Two identity classes taken from the human reference: TRBV6-2 with TRBV6-3, and
# IGLJ2 with IGLJ2A. The named IGLJ3 allele also shares that sequence. An
# unlisted allele does not affect the identity of the alleles named in a call.
# TRBV7-9 represents a gene the lookup does not name.
IDENTITY_CLASSES <- c(
  "TRBV6-2*01" = "VSEQ_TEST",
  "TRBV6-3*01" = "VSEQ_TEST",
  "IGLJ2*01"   = "JSEQ_TEST",
  "IGLJ2A*01"  = "JSEQ_TEST",
  "IGLJ3*01"   = "JSEQ_TEST"
)

describe("vj_tie_kind over the four statuses", {
  cases <- list(
    list(call = "TRBV6-2*01", expected = "resolved",
         why = "one gene is resolved"),
    list(call = "TRBV6-2*01,TRBV6-2*03", expected = "resolved",
         why = "two alleles of one gene are resolved without consulting the reference"),
    list(call = "TRBV6-2*01,TRBV6-3*01", expected = "assigned_group",
         why = "two genes the reference spells identically are one entity"),
    list(call = "TRBV6-2*01,TRBV7-9*01", expected = "ambiguous",
         why = "two genes standing for different sequence leave a real question"),
    list(call = "IGLJ2*01,IGLJ2A*01,IGLJ3*01", expected = "assigned_group",
         why = "an unlisted extra allele does not change identity of the named alleles"),
    list(call = "TRBV7-9*01,TRBV7-8*01", expected = "ambiguous",
         why = "two genes the reference does not name are two entities"),
    list(call = "", expected = "missing", why = "a blank call names nothing"),
    list(call = NA_character_, expected = "missing", why = "an NA call names nothing")
  )
  for (case in cases) {
    it(case$why, expect_identical(vj_tie_kind(case$call, IDENTITY_CLASSES),
                                  case$expected))
  }
  it("is vectorised over a whole column", {
    expect_identical(
      vj_tie_kind(c("TRBV6-2*01", "TRBV6-2*01,TRBV6-3*01",
                    "TRBV6-2*01,TRBV7-9*01", ""), IDENTITY_CLASSES),
      c("resolved", "assigned_group", "ambiguous", "missing"))
  })
})

describe("vj_canonical_token over the four statuses", {
  cases <- list(
    list(call = "TRBV6-2*01,TRBV6-3*01", expected = "VSEQ_TEST",
         why = "a resolved tie names the entity its genes share"),
    list(call = "TRBV7-9*01", expected = "TRBV7-9",
         why = "a gene the reference does not name is its own entity"),
    list(call = "TRBV6-2*01,TRBV7-9*01", expected = "",
         why = "a call spanning two entities names neither"),
    list(call = "", expected = "", why = "a blank call names nothing"),
    list(call = NA_character_, expected = "", why = "an NA call names nothing")
  )
  for (case in cases) {
    it(case$why, expect_identical(
      vj_canonical_token(case$call, IDENTITY_CLASSES), case$expected))
  }
  it("keeps single-gene calls separate from a complete-reference composite", {
    expect_identical(vj_canonical_token("TRBV6-2*01", IDENTITY_CLASSES), "TRBV6-2")
    expect_identical(vj_canonical_token("TRBV6-3*07", IDENTITY_CLASSES), "TRBV6-3")
    expect_identical(vj_canonical_token("TRBV6-2*01,TRBV6-3*01", IDENTITY_CLASSES), "VSEQ_TEST")
  })
  it("names an entity for exactly the calls vj_tie_kind can reduce", {
    calls <- c("TRBV6-2*01", "TRBV6-2*01,TRBV6-3*01", "TRBV6-2*01,TRBV7-9*01", "")
    expect_identical(
      nzchar(vj_canonical_token(calls, IDENTITY_CLASSES)),
      vj_tie_kind(calls, IDENTITY_CLASSES) %in% c("resolved", "assigned_group"))
  })
})

# Both arms read the tie kinds off a record the same way, so the reading itself
# is stated once and tested once.
describe("vj_calls_ambiguous over a record's two tie kinds", {
  KINDS <- c("resolved", "assigned_group", "ambiguous", "missing")

  it("reports an irreducible tie whichever call carries it", {
    expect_identical(
      vj_calls_ambiguous(KINDS, rep("resolved", 4L)),
      c(FALSE, FALSE, TRUE, FALSE))
    expect_identical(
      vj_calls_ambiguous(rep("resolved", 4L), KINDS),
      c(FALSE, FALSE, TRUE, FALSE))
  })
  # A call the reference reduces to one entity is answered, and a call naming no
  # gene at all asks nothing; neither is the ambiguity this reports.
  it("reports neither a reducible tie nor a missing call", {
    expect_identical(
      vj_calls_ambiguous(c("assigned_group", "missing"),
                             c("missing", "assigned_group")),
      c(FALSE, FALSE))
  })
  it("reports nothing for no records at all", {
    expect_identical(vj_calls_ambiguous(character(0), character(0)),
                     logical(0))
  })
})

# A multi-valued call is a set of equally scoring hits rather than a ranking, so
# "first" means an arbitrary but stable representative of the tie.
describe("first_vj_gene", {
  cases <- list(
    list(call = "IGLV1-44*01,IGLV1-36*01", expected = "IGLV1-44",
         why = "keeps the gene the call lists first"),
    list(call = "IGLV1-36*01,IGLV1-44*01", expected = "IGLV1-36",
         why = "keeps the gene listed first even when that is not the lower gene number"),
    list(call = "", expected = "",
         why = "returns an empty string for a blank call")
  )
  for (case in cases) {
    it(case$why, expect_identical(first_vj_gene(case$call), case$expected))
  }
})

# At gene level, a multi-allele call naming one gene is resolved. This helper
# reports multiple gene names without consulting sequence identity classes.
describe("vj_gene_status", {
  cases <- list(
    list(call = "IGHV3-23*01", expected = "resolved",
         why = "a single-allele call is resolved"),
    list(call = "IGHV3-23*01,IGHV3-23*04", expected = "resolved",
         why = "two alleles of one gene are resolved"),
    list(call = "IGHV1-69*01,IGHV1-69D*01", expected = "ambiguous_multi_gene",
         why = "a gene and its duplicated location remain two genes"),
    list(call = "", expected = "missing",
         why = "a blank call is missing"),
    list(call = NA_character_, expected = "missing",
         why = "an NA call is missing")
  )
  for (case in cases) {
    it(case$why, expect_identical(vj_gene_status(case$call), case$expected))
  }
})

# Reconstructing a germline picks one reference sequence, so it needs the allele
# the call names and not just the gene. These two are the pair the germline stage
# reads: which allele was used, and how much of a choice that was.
describe("first_vj_allele", {
  cases <- list(
    list(call = "IGHV3-23*01,IGHV3-23*04", expected = "IGHV3-23*01",
         why = "keeps the allele suffix that first_vj_gene strips"),
    list(call = "IGLV1-36*01,IGLV1-44*01", expected = "IGLV1-36*01",
         why = "takes the allele the call lists first, as CreateGermlines does"),
    list(call = "IGHV1-69", expected = "IGHV1-69",
         why = "returns a call with no allele suffix unchanged"),
    list(call = "", expected = "",
         why = "returns an empty string for a blank call")
  )
  for (case in cases) {
    it(case$why, expect_identical(first_vj_allele(case$call), case$expected))
  }
})

# Two alleles of one gene are resolved at gene level and unresolved here, which is
# exactly the distinction between grouping records and building one germline.
describe("vj_call_ambiguity_status", {
  cases <- list(
    list(call = "IGHV3-23*01", expected = "resolved",
         why = "a single-allele call is resolved"),
    list(call = "IGHV3-23*01,IGHV3-23*01", expected = "resolved",
         why = "a repeated allele is one allele"),
    list(call = "IGHV3-23*01,IGHV3-23*04", expected = "ambiguous_multi_allele",
         why = "two alleles of one gene leave the germline sequence unchosen"),
    list(call = "IGHV1-69*01,IGHV1-69D*01", expected = "ambiguous_multi_gene",
         why = "a gene and its duplicated location remain different genes"),
    list(call = "IGLV1-44*01,IGLV1-36*01", expected = "ambiguous_multi_gene",
         why = "two genes are the coarser ambiguity"),
    list(call = "", expected = "missing",
         why = "a blank call is missing"),
    list(call = NA_character_, expected = "missing",
         why = "an NA call is missing")
  )
  for (case in cases) {
    it(case$why,
       expect_identical(vj_call_ambiguity_status(case$call), case$expected))
  }
  it("is vectorised over a whole column", {
    expect_identical(
      vj_call_ambiguity_status(c("IGHV3-23*01", "IGHV3-23*01,IGHV3-23*04",
                                 "IGLV1-44*01,IGLV1-36*01", "")),
      c("resolved", "ambiguous_multi_allele", "ambiguous_multi_gene", "missing")
    )
  })
})


test_that("shared V calls use concordant J evidence without requiring D or C", {
  records <- data.frame(
    locus = c("TRA", "TRA", "TRD", "TRA"),
    v_call = c("TRAV14/DV4*01", "TRAV38-2/DV8*01",
               "TRAV14/DV4*01", "TRAV14/DV4*01,TRAV14/DV4*02"),
    j_call = c("TRDJ1*01", "TRDJ1*01,TRDJ2*01", "TRAJ1*01", "TRDJ1*01"),
    d_call = c("TRDD2*01", "", NA, ""),
    c_call = c("TRDC*01", NA, "TRAC*01", ""),
    sequence_alignment = "TGTGCCTTT", productive = "T",
    junction_length = 9L, consensus_count = c(3, 5, 7, 11))
  resolved <- resolve_shared_trav_dv_locus(records)
  expect_identical(resolved$locus, c("TRD", "TRD", "TRA", "TRD"))
  expect_identical(resolved[names(records) != "locus"], records[names(records) != "locus"])
  expect_true(all(airr_qc_pass(airr_qc_flags(resolved))))
  expect_identical(resolve_shared_trav_dv_locus(resolved), resolved)
  expect_identical(resolve_shared_trav_dv_locus(records[FALSE, ]), records[FALSE, ])
})

test_that("shared V resolution does not rescue conflicting or incomplete calls", {
  records <- data.frame(
    locus = "TRA", v_call = "TRAV14/DV4*01", j_call = "TRDJ1*01",
    d_call = "", c_call = "", sequence_alignment = "TGTGCCTTT",
    productive = "T", junction_length = 9L, consensus_count = 3)
  conflicts <- list(
    c(v_call = "TRAV14/DV4*01,bogus"),
    c(j_call = "TRDJ1*01,TRAJ1*01"),
    c(j_call = ""), c(j_call = "TRBJ1-1*01"),
    c(d_call = "TRBD1*01"), c(d_call = "TRDD2*01,TRBD1*01"),
    c(c_call = "TRAC*01"), c(c_call = "TRDC*01,TRAC*01"))
  for (conflict in conflicts) {
    input <- records
    input[names(conflict)] <- as.list(conflict)
    resolved <- resolve_shared_trav_dv_locus(input)
    expect_identical(resolved, input)
    expect_false(airr_qc_pass(airr_qc_flags(resolved)))
  }
  # C qualification remains separate for records whose locus already agrees.
  records$locus <- "TRD"
  records$c_call <- ""
  expect_true(airr_qc_pass(airr_qc_flags(resolve_shared_trav_dv_locus(records))))
})


test_that("the complete TRA/TRD V pool resolves from J with consistent D and C", {
  records <- data.frame(locus = c("TRD", "TRA", "TRA"),
    v_call = c("TRDV1*01", "TRAV1-1*01", "TRAV1-1*01,TRAV14/DV4*01"),
    j_call = c("TRAJ1*01", "TRDJ1*01", "TRDJ1*01"),
    d_call = c("", "TRDD2*01", ""), c_call = c("TRAC*01", "TRDC*01", ""),
    sequence_alignment = "TGTGCCTTT", productive = "T", junction_length = 9L, consensus_count = 3)
  resolved <- resolve_shared_trav_dv_locus(records)
  expect_identical(resolved$locus, c("TRA", "TRD", "TRD"))
  expect_true(all(airr_qc_pass(airr_qc_flags(resolved))))
})
