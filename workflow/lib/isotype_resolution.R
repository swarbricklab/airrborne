# Name antibody classes and subclasses from IGH constant-region gene calls.
# The calls have already passed the evidence checks in isotype_calls.R.
# resolve_isotype() and the gene-to-isotype tables below decide how tied calls
# resolve at each level. For example, IGHG1 and IGHG3 together identify IgG but
# not one subclass. isotype_metrics.R calculates their read-support fractions.
# Unresolved calls remain missing. These helpers use base R.

.IGH_CONSTANT_GENE_ISOTYPE <- c(
  IGHM = "IgM", IGHD = "IgD", IGHG1 = "IgG1", IGHG2 = "IgG2",
  IGHG3 = "IgG3", IGHG4 = "IgG4", IGHA1 = "IgA1", IGHA2 = "IgA2", IGHE = "IgE"
)
CANONICAL_IGH_ISOTYPES <- unname(.IGH_CONSTANT_GENE_ISOTYPE)

# Map allele-stripped calls to subclasses. IGHG4A, a duplicated copy of IGHG4
# that the IMGT reference names as a separate gene, counts as IgG4. No other
# gene name is mapped to a subclass it does not name, because guessing could
# create false isotype evidence.
.subclass_set <- function(cell) {
  genes <- sub("\\*.*$", "", .split_calls(cell))
  genes[genes == "IGHG4A"] <- "IGHG4"
  unname(.IGH_CONSTANT_GENE_ISOTYPE[unique(genes[genes != ""])])
}

# Unknown tokens, including IGHGP, invalidate the whole call at every level.
.SUBCLASS_TO_CLASS <- c(
  IgM = "IgM", IgD = "IgD", IgE = "IgE",
  IgG1 = "IgG", IgG2 = "IgG", IgG3 = "IgG", IgG4 = "IgG",
  IgA1 = "IgA", IgA2 = "IgA"
)

# The major class of each subclass token, or NA where the token is not recognised.
.iso_class <- function(subclasses) unname(.SUBCLASS_TO_CLASS[subclasses])

# IgM and IgD are classified as unswitched; IgG, IgA and IgE are switched.
# This classification does not establish a cell's differentiation state.
.CLASS_TO_SWITCH <- c(IgM = "unswitched", IgD = "unswitched",
                      IgG = "switched", IgA = "switched", IgE = "switched")

# The switch status of each major class, or NA where the class is not recognised.
.iso_switch <- function(classes) unname(.CLASS_TO_SWITCH[classes])

#' Resolve each qualified C call to one subclass, class or switch status.
#'
#' Tied subclasses may resolve at a broader level (IGHG1,IGHG3 resolves to IgG).
#' Missing, unfamiliar or still-ambiguous calls return NA. Return one character
#' value per input call in input order, including character(0) for empty input.
resolve_isotype <- function(x, level = c("subclass", "class", "switch")) {
  level <- match.arg(level)
  # Resolve distinct calls once because many records share the same call.
  calls <- unique(as.character(x))
  resolved <- vapply(calls, function(cell) {
    subclasses <- .subclass_set(cell)
    if (!length(subclasses)) return(NA_character_)
    # Reject the whole call if any token is not a recognised IGH isotype.
    # Dropping invalid tokens could turn annotation noise into false evidence.
    if (any(is.na(.iso_class(subclasses)))) return(NA_character_)
    resolved_values <- switch(level,
      subclass = subclasses,
      class    = .iso_class(subclasses),
      switch   = .iso_switch(.iso_class(subclasses)))
    resolved_values <- unique(resolved_values)
    if (length(resolved_values) == 1L) resolved_values else NA_character_
  }, character(1), USE.NAMES = FALSE)
  resolved[match(as.character(x), calls)]
}

# SHM and BASELINe report these groups in this order. all_IGH needs no C call.
ANALYSIS_IGH_ISOTYPE_GROUPS <- c("all_IGH", "IgM", "IgG", "IgA")

# SHM adds these subclasses beneath their classes; membership overlaps. An
# unresolved IgG subclass remains in IgG alone. IgM needs no duplicate subclass
# row, while IgD and IgE remain within all_IGH. BASELINe keeps class groups because
# subclass splitting leaves fewer clones contributing to each selection estimate.
ANALYSIS_IGH_SUBCLASS_GROUPS <- c("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2")

# Each stratified group uses the call resolution specified here.
.ANALYSIS_GROUP_LEVEL <- c(
  IgM = "class", IgG = "class", IgA = "class",
  IgG1 = "subclass", IgG2 = "subclass", IgG3 = "subclass", IgG4 = "subclass",
  IgA1 = "subclass", IgA2 = "subclass"
)

#' Resolve evidence-qualified C calls to the released SHM/BASELINe isotypes.
#'
#' level selects the three reported classes or six subclasses. Other isotypes
#' return NA; IgD and IgE remain in all_IGH only. Calls tied within one class
#' resolve to that class but stay missing at subclass level.
resolve_analysis_igh_isotype <- function(x, level = c("class", "subclass")) {
  level <- match.arg(level)
  resolved <- resolve_isotype(x, level)
  allowed <- names(.ANALYSIS_GROUP_LEVEL)[.ANALYSIS_GROUP_LEVEL == level]
  resolved[!resolved %in% allowed] <- NA_character_
  resolved
}

#' Whether a released isotype group is decided on the major class or the subclass.
#'
#' `group` is one of the stratified groups, so not `all_IGH`, which is decided on
#' no call at all. Returns "class" or "subclass", vectorised.
analysis_igh_isotype_level <- function(group) {
  level <- unname(.ANALYSIS_GROUP_LEVEL[group])
  if (anyNA(level)) {
    stop("no resolution level for isotype group(s): ",
         paste(unique(group[is.na(level)]), collapse = ", "))
  }
  level
}
