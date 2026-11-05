"""Generate the methods-as-run paragraphs for one scope release.

The paragraphs come from the scope's provenance.json, effective_config.yaml,
manifest and operative clone threshold, so every setting, version and identity
they state is one the release records. The release README and the offline
report print the same paragraphs. The text uses the passive voice of a methods
section. A sentence that states no setting describes fixed pipeline behaviour,
so a change to that behaviour needs a matching edit here. FIXED_METHODS records
the same kind of fixed choice as data, for effective_config.yaml.
"""

import csv
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml

# Every library receives a diversity row at each of these loci, so the methods
# state a diversity minimum for each.
DIVERSITY_LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")

# Scientific choices the code hard-codes rather than reads from configuration.
# build_scope_provenance.py records them under fixed_methods in
# effective_config.yaml. Each entry names one R function and its arguments as the
# call writes them. test_provenance.py compares every entry with the source file
# that makes the call, so a change there needs a matching edit here.
FIXED_METHODS = {
    "clone_distance": {"function": "shazam::distToNearest", "model": "ham", "normalize": "len",
                       "symmetry": "avg", "first": False, "VJthenLen": True},
    "clone_threshold_fit": {"function": "shazam::findThreshold", "method": "gmm",
                            "model": "gamma-norm", "cutoff": "user"},
    "clone_assignment": {"function": "scoper::hierarchicalClones", "method": "nt",
                         "linkage": "single", "normalize": "len"},
    "shm_counting": {"function": "shazam::calcObservedMutations", "regionDefinition": "IMGT_V_BY_REGIONS",
                     "mutationDefinition": None, "ambiguousMode": "eitherOr"},
    "selection": {
        "clone_consensus": {"function": "shazam::collapseClones", "method": "thresholdedFreq",
                            "includeAmbiguous": False, "breakTiesStochastic": False},
        "baseline": {"function": "shazam::calcBaseline", "testStatistic": "focused",
                     "targetingModel": "HH_S5F"},
    },
}


def read_provenance(path: Path, scope_id: str, library_ids: set[str]) -> dict:
    """Require the current scope schema and exact assembly-history coverage."""
    record = json.loads(path.read_text())
    if record["schema_version"] != "3.2":
        raise ValueError("Unsupported scope provenance schema")
    if record["scope_id"] != scope_id:
        raise ValueError("Release provenance belongs to another scope")
    history = record["assembly_producers"]
    known = set(history["libraries"])
    unknown = history["unknown"]["libraries"]
    if known & set(unknown) or len(set(unknown)) != len(unknown) or known | set(unknown) != library_ids:
        raise ValueError("Provenance must retain exactly the scope's assembly histories")
    if set(history["libraries"].values()) != set(history["producers"]):
        raise ValueError("Known assembly histories do not match their producer identities")
    return record


def read_methods_inputs(config_path: Path, manifest_path: Path,
                        threshold_path: Path) -> tuple[dict, list[dict], dict]:
    """Read the required inputs with identifiers preserved as strings."""
    config = yaml.safe_load(config_path.read_text())
    with manifest_path.open(newline="") as handle:
        libraries = list(csv.DictReader(handle, delimiter="\t"))
    with threshold_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError("Methods require exactly one operative scope threshold")
    return config, libraries, rows[0]


def _libraries(count: int) -> str:
    """Name a library count with the matching noun, as in "1 library"."""
    return f"{count} library" if count == 1 else f"{count} libraries"


def methods_paragraphs(record: Mapping[str, Any], config: Mapping[str, Any],
                        libraries: Sequence[Mapping[str, str]],
                        threshold: Mapping[str, Any]) -> list[str]:
    """Return the methods paragraphs for one scope, in reading order.

    `record` is the scope provenance, `config` the effective configuration,
    `libraries` the scope manifest rows and `threshold` the one operative
    threshold row. Stops when an input belongs to another scope or the manifest
    differs from the libraries the provenance records assembly histories for.
    """
    scope = record["scope_id"]
    if (config["scope_id"] != scope or threshold["scope_id"] != scope or not libraries
            or any(row["scope_id"] != scope for row in libraries)):
        raise ValueError("Methods inputs belong to another scope")
    ids = [row["library_id"] for row in libraries]
    history = record["assembly_producers"]
    if (len(set(ids)) != len(ids)
            or set(ids) != set(history["libraries"]) | set(history["unknown"]["libraries"])):
        raise ValueError("Methods manifest differs from the assembly histories")
    params = config["params"]
    source = record["source"]
    methods = [
        "Receptor sequences were recovered from bulk RNA sequencing. Read support is not a cell count. "
        "Heavy and light chains, and the two T-cell receptor chains, are unpaired.",
        f"The release records source revision {source['revision']} "
        + (f"with local changes identified by source state {source['state_sha256']}, "
           if source["dirty"] else "with no recorded local changes, ")
        + f"provenance schema {record['schema_version']} and seed {config['seed']}.",
    ]
    measured = [f"{entry['name']} {entry['version']}" for entry in record["software"]["analysis"]["tools"]
                if entry["status"] == "measured"]
    if measured:
        methods.append("The analysis image contains " + ", ".join(measured) + ".")
    unknown_tools = [entry["name"] for entry in record["software"]["analysis"]["tools"]
                     if entry["status"] != "measured"]
    if unknown_tools:
        methods.append("The versions of " + ", ".join(unknown_tools)
                       + " in the analysis image could not be queried.")
    qc, lineage, igblast = params["airr_qc"], params["lineage"], params["igblast"]
    for producer_id, producer in sorted(history["producers"].items()):
        count = sum(value == producer_id for value in history["libraries"].values())
        versions = ", ".join(f"{entry['name']} {entry['version']}" for entry in producer["software"]
                             if entry["status"] == "measured") or "TRUST4 with an unavailable version"
        settings = ", ".join(f"{key}={value}" for key, value in sorted(producer["settings"].items()))
        methods.append(f"Assembly of {_libraries(count)} was recorded with {versions}, "
                       f"reference {producer['reference']['name']} (SHA256 {producer['reference']['sha256']}) "
                       + (f"and settings {settings}." if settings else "and no recorded settings."))
    methods.append("Records whose TRUST4 CDR3 was imputed rather than observed were removed before "
                   "annotation. Source identifiers and original support were preserved in records "
                   "retained by annotation and common QC. "
                   f"IgBLAST reported at most {igblast['num_alignments_v']} V, "
                   f"{igblast['num_alignments_j']} J, {igblast['num_alignments_d']} D and "
                   f"{igblast['num_alignments_c']} C alignments per query.")
    methods.append("MakeDb V and J CIGARs and supported D CIGARs were normalised to AIRR "
                   "query/reference semantics. Insertions consume query positions and deletions "
                   "consume reference positions. Standard V endpoints were derived on the selected "
                   "ungapped reference. The released standard alignments were checked against "
                   "query and reference intervals, selected alleles and reference bases before "
                   "publication. Unsupported D CIGARs and germline endpoints were left empty, "
                   "with an availability status and retained native hit evidence. Native Change-O "
                   "reconstruction coordinates were preserved separately. Queries were stored in "
                   "V orientation with standard rev_comp set to F, and the historical MakeDb "
                   "flag was retained in changeo_rev_comp. These semantic checks "
                   "were separate from AIRR structural validation and metadata completeness.")
    methods.append("Calls from the shared TRA/TRD V pool were resolved as TRA or TRD from the J call when "
                   "every V candidate belonged to that pool and every present D or C call agreed. "
                   "Missing D and C calls were allowed. This interpretation preceded the locus "
                   "audit and common QC. Raw gene calls were retained.")
    methods.append("The common record-QC gate required that the V and J calls, and any D call, agree "
                   "with the assigned locus, that the IgBLAST/MakeDb productivity annotation be true, "
                   "that its junction length be divisible by three, that its read support be "
                   "positive, and that the aligned N "
                   f"fraction be at most {qc['max_alignment_n_fraction']}. A constant-region call was "
                   f"read as an isotype only when it was present, had a complete candidate list, named only IGH constant genes and had "
                   f"at least {qc['min_informative_c_alignment_positions']} informative alignment "
                   "positions. A call that failed was masked and its record retained.")
    methods.append("Productivity and stop-codon calls describe the annotated V(D)J sequence. "
                   "Unaligned terminal bases are outside the stop assessment. Passing common QC "
                   "does not establish functionality of the complete assembled transcript.")
    methods.append("The dominant common-QC variant per contig and locus at IGH, IGK and IGL was selected "
                   "by original record support. Ties were broken by an A/C/G/T-only junction with "
                   "single-gene V and J calls, then V identity, usable alignment length, junction, "
                   "full sequence, V, J and C calls, and sequence ID. Grouping eligibility was tested after this choice, without substituting "
                   "a minor eligible sibling for an ineligible dominant variant. TCR variants remained distinct at "
                   "this step. Identical full nucleotide sequences among retained representatives "
                   "were then merged within each library and locus, requiring agreeing analytical "
                   "annotations. Each eligible source record contributed support once. The source "
                   "records and both representative links were retained. IGH fitting and assignment "
                   "used these unique-sequence members. Non-IGH clonotypes used exact "
                   "nucleotide junctions with resolved V and J germline entities.")
    methods.append("A V or J call naming alleles of one gene was assigned that gene. "
                   "Multi-gene V/J calls were resolved only when every named allele shared one "
                   "complete, unambiguous reference-region sequence in the pinned snapshot. "
                   "Partial reference entries and potentially truncated candidate lists could not "
                   "establish equivalence. D-marker names were retained unless sequence identity "
                   "established a composite entity. Raw calls and reference metadata were retained. "
                   "Unresolved V/J ambiguity and present C-locus conflicts prevented grouping.")
    methods.append(f"The IGH threshold setting was {lineage['clonal_threshold']}. "
                   f"A donor fit required at least {lineage['min_fit_distances']} distances. "
                   f"The permitted junction N or ? count was {lineage['max_n']}.")
    methods.append(f"The scope used threshold {threshold['clone_threshold']} ({threshold['threshold_source']}).")
    methods.append("Automatic threshold selection uses the median of usable donor fits. "
                   "A usable fit has fit_status fitted: its distances are not all equal, and its "
                   "fitted cutoff is finite and between 0 and 1. "
                   "If no donor fit is usable, the automatic fallback is 0.1. "
                   "A configured numeric threshold takes precedence over the donor fits.")
    subsample = lineage["fit_subsample"] if lineage["fit_subsample"] is not None else "all available distances"
    methods.append("IGH nearest-neighbour distances (SHazaM distToNearest) and clone assignment "
                   "(SCOPer hierarchicalClones) partitioned records by exact V and J entities and "
                   "junction length, used nucleotide Hamming distance normalised by junction length, "
                   "and clustered with single linkage within each donor and scope. Composite entities "
                   "joined only themselves. The threshold constrained links, not all pairwise "
                   "distances within a clone. Donor threshold calibration fitted a "
                   "gamma-normal mixture with SHazaM findThreshold, with edge "
                   f"{lineage['fit_edge']}, specificity {lineage['fit_specificity']} and subsample "
                   f"{subsample}, on distances sorted so that the subsample did not depend on "
                   "record order. These donor fits ran whether or not they supplied the operative "
                   "threshold. SHazaM's internal seed resets in gmmFit and rocSpace were replaced "
                   f"with seed {config['seed']} before fitting. The fitting algorithm was otherwise "
                   "unchanged.")
    targets = "; ".join(f"{locus} = {value}" for locus, value in sorted(params["rarefaction"]["targets"].items()))
    minima = "; ".join(
        f"{locus} = {params['diversity']['min_consensus_count'].get(locus, 0)}"
        for locus in DIVERSITY_LOCI)
    rarefaction = (
        "Rarefaction was calculated within the observed support vector without replacement. "
        "The target-depth quantities were expected richness (q = 0), exp(E[H]) for Shannon "
        "entropy H (q = 1) and the reciprocal of expected Simpson concentration (q = 2). "
        "The q = 2 quantity was not the mean inverse Simpson across draws. "
        f"Target support depths were {targets}. "
        "The calculation is analytical and uses no simulated draws. "
        if targets else "No rarefaction target depths or curves were configured. ")
    methods.append(rarefaction + f"Minimum support for diversity estimates was {minima}.")
    shm = params["shm"]
    methods.append("SHM frequency was the replacement and silent mutation count from SHazaM "
                   "calcObservedMutations in each IMGT V subregion, divided by the positions with an "
                   "observed base in both the record and its D-masked germline. Whole V summed both "
                   "counts over FWR1, CDR1, FWR2, CDR2 and FWR3, excluding CDR3 and FWR4. "
                   f"A supported IMGT V position carried at least {shm['min_position_votes']} "
                   "consensus votes agreeing with the assembled base. SHM required supported "
                   f"coverage of {shm['min_supported_subregion_fraction']} of each scored region's "
                   "reconstructed-germline occupied positions, excluding numbering gaps. "
                   "Whole V also required at least "
                   f"{shm['min_germline_occupied_v_positions']} occupied germline V positions.")
    methods.append("SHM coverage was tested independently on each retained contig representative "
                   "with the member's exact full sequence. For each region, a qualifying "
                   "source was selected by supported position count, original record support and "
                   "sequence ID. Whole V required one source to pass all five subregions and the "
                   "occupancy floor. Votes and partial coverage were never pooled across sources. "
                   "Each region's evidence source was recorded. SHM was unavailable when no "
                   "source qualified. Each eligible unique sequence contributed once to its "
                   "library-local clone median, and clone medians received equal weight.")
    methods.append("BASELINe selection (SHazaM calcBaseline) used the focused test over the IMGT V "
                   "CDR and FWR regions with the HH_S5F targeting model. The clone consensus "
                   "(SHazaM collapseClones) used a minimum base frequency of "
                   f"{params['selection']['collapse_min_freq']}. Each member, one distinct full "
                   "sequence within a library and locus, contributed once to this consensus "
                   "within its library and isotype group. The effective germline was built "
                   "separately with the pinned collapseClones mostCommon rule, with "
                   "deterministic ties and ambiguous bases excluded. "
                   "A consensus gap opposite an A, C, G or T germline base was set to N before "
                   "estimation. Incomplete codons were masked in that effective sequence before "
                   "observed mutations and expected opportunities were calculated, preserving "
                   "IMGT frame and numbering gaps. Non-covering members remained in the "
                   "consensus denominator. Partial sequences remained eligible. Regional usable "
                   "density counts, counts of clones with a covered position and 95% intervals "
                   "were reported without an additional evidence cutoff.")
    methods.append("Library pairs within one donor were compared at "
                   + ", ".join(params["overlap"]["loci"])
                   + ". Jaccard used group detection alone. Morisita-Horn used Horn's proportion "
                   "form over within-library read-support proportions. A pair with no support at "
                   "one endpoint carries its status and no similarity value.")
    pair_minima = "; ".join(
        f"{locus} = {params['overlap']['min_consensus_count'].get(locus, 0)}"
        for locus in params["overlap"]["loci"])
    methods.append("Bray-Curtis dissimilarity was half the sum of absolute support-proportion "
                   "differences over the complete union of receptor groups. "
                   "Differences in the released pair table compare endpoint 2 (q) with endpoint 1 "
                   "(p), whose order carries no chronology. Shared receptor groups were counted as "
                   "higher at endpoint 2, higher at endpoint 1 or equal by their unrounded "
                   "proportions, with no effect-size cutoff. "
                   "Each endpoint retains its full support denominator. Pair-support minima "
                   f"were {pair_minima}. Zero disables the positive minimum. Values below a "
                   "pair minimum remain released with a flag and are excluded from default "
                   "pair report summaries. Donor overlap means are released only when every "
                   "donor library has support at the locus, and they average every pair "
                   "without the pair-support filter. Report chronology uses supplied visit order only. "
                   "All differences describe recovered read support, not absolute cell abundance.")
    files = [f"{name} ({identity['name']}, SHA256 {identity['sha256']})"
             for name, identity in record["references"].items() if "sha256" in identity]
    if files:
        methods.append("Reference file identities: " + "; ".join(files) + ".")
    for identity in record["references"].values():
        if "imgt_snapshot" not in identity:
            continue
        text = (f"The {identity['name']} reference was curated from the IMGT/GENE-DB download dated "
                f"{identity['imgt_snapshot']['date']} in the analysis image, and the image's IgBLAST "
                "germline databases were checked to hold the same alleles and sequences.")
        auxiliary = identity["audit"]["auxiliary_data"]
        added = auxiliary["supplemented"]
        if added:
            text += (" The IgBLAST auxiliary data file was NCBI's file from the image with entries added for "
                     + ", ".join(f"{entry['allele']} ({'identical to' if entry['rule'] == 'identical_sequence' else 'from'}"
                                 f" {', '.join(entry['copied_from'])})" for entry in added)
                     + ". Each added entry was copied from an allele with an identical sequence, or from a "
                       "same-gene allele of equal length whose entry the file held and whose conserved codon "
                       f"and frame the allele preserves. These rules reproduced all "
                       f"{auxiliary['identical_pairs_checked']} identical-sequence pairs and all "
                       f"{auxiliary['sibling_pairs_checked']} same-gene pairs of equal length in the file.")
        absent = auxiliary["j_alleles_without_entry"]
        if absent:
            text += (" No auxiliary entry exists for J alleles "
                     + ", ".join(f"{entry['allele']} ({entry['functionality']})" for entry in absent)
                     + ". A record whose first-listed top J hit was one of these alleles received no "
                       "CDR3 from IgBLAST and was rejected by MakeDb.")
        methods.append(text)
    other_reference = sorted(library for library, identity in history["libraries"].items()
                             if not history["producers"][identity]["reference"]["same_as_annotation_reference"])
    if other_reference:
        methods.append(f"Assembly of {_libraries(len(other_reference))} read a reference file that differs "
                       "from the annotation reference: " + ", ".join(other_reference) + ".")
    analysis = record["images"][record["software"]["analysis"]["image"]]
    assembly = record["images"][record["software"]["configured_assembly"]["image"]]
    methods.append(f"The analysis image identity is {analysis['digest']} ({analysis['format']} form). "
                   f"The configured assembly image identity is {assembly['digest']} ({assembly['format']} form).")
    assembly_versions = [f"{entry['name']} {entry['version']}"
                         for entry in record["software"]["configured_assembly"]["tools"]
                         if entry["status"] == "measured"]
    if assembly_versions:
        methods.append("The queried assembly image contains " + ", ".join(assembly_versions) + ".")
    # A --post-trust4 run does not obtain the configured assembly image, so its
    # version query can fail while the reused assemblies keep their own record.
    unknown_assembly = [entry["name"] for entry in record["software"]["configured_assembly"]["tools"]
                        if entry["status"] != "measured"]
    if unknown_assembly:
        methods.append("The " + ", ".join(unknown_assembly)
                       + " version of the configured assembly image could not be queried.")
    configured = record["software"]["configured_assembly"]["image"]
    other_image = [library for library, producer_id in history["libraries"].items()
                   if history["producers"][producer_id]["image"] != configured]
    if other_image:
        methods.append(f"Assembly of {_libraries(len(other_image))} used an image other than the "
                       "configured assembly image. The assembly versions recorded above describe "
                       "those libraries.")
    methods.append(f"Assembly producer history is recorded for {_libraries(len(history['libraries']))} "
                   f"and unavailable for {len(history['unknown']['libraries'])}. "
                   "Configured assembly settings and queried image versions do not identify the "
                   "producers of reused outputs with unknown history.")
    return methods
