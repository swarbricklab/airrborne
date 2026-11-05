"""Assemble producer-owned population counts for one scope.

The audit carries library/donor counts. Estimate statuses and pair comparisons
remain in their metric tables. No biological predicate is applied here.
qc.tsv divides supplied support counts and sums the grouping reasons the
library AIRR exports record. qc_trail.tsv renames and reorders audit rows for
release.
"""

import csv
from decimal import Decimal, InvalidOperation
import gzip
from pathlib import Path
from typing import Iterable, Mapping, Sequence

COLUMNS = (
    "scope_id", "grain", "donor_id", "tissue_id", "library_id", "locus",
    "locus_basis", "analysis", "stage", "stage_index", "isotype_group", "v_region",
    "kind", "decision_status", "input_unit", "output_unit", "n_units_in",
    "n_units_out", "consensus_count_in", "consensus_count_out", "support_conserved",
)
# qc_audit.tsv, the scope's working audit, leads with the identifiers in
# IDENTITY_ORDER (scope_tables.py), as every collected table does. Producer
# shards keep COLUMNS, which qc_audit.R writes. Donor-level rows have blank
# library and tissue IDs.
SCOPE_COLUMNS = (
    "scope_id", "library_id", "tissue_id", "donor_id",
    *(column for column in COLUMNS
      if column not in ("scope_id", "library_id", "tissue_id", "donor_id")),
)
KEY = ("scope_id", "donor_id", "library_id", "locus", "locus_basis",
       "analysis", "stage", "isotype_group", "v_region")
COUNTS = ("n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out")
LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
# qc.tsv has one row per library and locus IgBLAST assigned. Each fraction's
# denominator is a support count in the same row: IgBLAST-annotated support for
# the first two fractions and common-QC support for the rest.
QC_COLUMNS = (
    "scope_id", "library_id", "tissue_id", "donor_id", "locus",
    "trust4_record_count", "trust4_consensus_count",
    "igblast_consensus_count", "nonproductive_fraction",
    "common_qc_consensus_count", "common_qc_fraction_of_annotated",
    "grouped_consensus_count", "grouping_eligible_fraction_of_common_qc",
    "v_call_unresolved_fraction", "j_call_unresolved_fraction",
    "junction_ineligible_fraction", "c_call_ineligible_fraction",
    "contig_representative_ineligible_fraction",
    "class_callable_fraction", "subclass_callable_fraction",
    "n_grouped_igh_members", "n_shm_measured_members_whole_v",
    "n_shm_eligible_members_whole_v", "shm_eligible_fraction_of_grouped_members_whole_v",
    "n_shm_eligible_members_any_subregion", "shm_eligible_fraction_of_grouped_members_any_subregion",
)
# The grouping failure reasons qc.tsv reports, named by their fraction column
# without the suffix. Every reason except V and J resolution has fixed tokens in
# grouping_status: igh_grouping.R and rearrangement_groups.R join a record's
# failures with ";", and grouping_constant_status() in gene_calls.R names the C
# failures. prepare_assembly_mapping() in assembly_backbone.R gives an otherwise
# eligible record whose contig's dominant variant fails a single status, which
# the export writes as ineligible_contig_representative.
GROUPING_REASONS = ("v_call_unresolved", "j_call_unresolved", "junction_ineligible",
                    "c_call_ineligible", "contig_representative_ineligible")
REASON_TOKENS = {
    "ineligible_junction": "junction_ineligible",
    "c_locus_conflict": "c_call_ineligible",
    "c_call_incomplete_candidates": "c_call_ineligible",
    "ineligible_contig_representative": "contig_representative_ineligible",
}
EXPORT_SUFFIX = "_rearrangements.tsv.gz"
# qc_trail.tsv is the release view of the library audit rows.
TRAIL_COLUMNS = (
    "scope_id", "library_id", "tissue_id", "donor_id", "locus", "stage_order", "stage",
    "analysis", "kind", "decision_status", "isotype_group", "v_region", "locus_source",
    "input_unit", "n_units_in", "output_unit", "n_units_out",
    "consensus_count_in", "consensus_count_out",
)
# Audit stage, released stage, analysis and kind, in the released order that
# stage_order numbers from 1. The kinds separate an exclusion, after which a
# record takes no further part, from eligibility for one analysis and from the
# masking of a constant-region call; docs/qc_matrix.md describes each stage.
# clone_observations at IGH and exact_groups at the other loci both become
# receptor_units, which counts the members that grouping combines.
TRAIL_STAGES = (
    ("trust4_contigs_with_cdr3_call", "trust4_contigs", "assembly", "checkpoint"),
    ("trust4_cdr3_calls", "trust4_cdr3_calls", "assembly", "checkpoint"),
    ("trust4_reported_airr", "trust4_records", "assembly", "checkpoint"),
    ("imputed_cdr3_filter", "observed_cdr3", "assembly", "exclusion"),
    ("igblast_locus_routing", "igblast_input", "annotation", "exclusion"),
    ("igblast_makedb", "igblast_annotated", "annotation", "exclusion"),
    ("igblast_locus_agreement", "locus_agreement", "annotation", "diagnostic"),
    ("qc_locus_match", "qc_locus_match", "common_qc", "criterion"),
    ("qc_d_locus_match", "qc_d_locus_match", "common_qc", "criterion"),
    ("qc_low_n", "qc_low_n", "common_qc", "criterion"),
    ("qc_productive", "qc_productive", "common_qc", "criterion"),
    ("qc_junction_mod3", "qc_junction_mod3", "common_qc", "criterion"),
    ("qc_positive_support", "qc_positive_support", "common_qc", "criterion"),
    ("qc_pass_total", "common_qc", "common_qc", "exclusion"),
    ("grouping_eligible", "grouping_eligible", "grouping", "eligibility"),
    ("assembly_backbone_collapse", "contig_collapse", "grouping", "grouping"),
    ("exact_sequence_collapse", "exact_sequence_collapse", "grouping", "grouping"),
    ("clone_observations", "receptor_units", "grouping", "grouping"),
    ("exact_groups", "receptor_units", "grouping", "grouping"),
    ("c_call_present", "c_call_present", "isotype", "criterion"),
    ("c_call_locus_match", "c_call_locus_match", "isotype", "criterion"),
    ("c_call_candidates_complete", "c_call_candidates_complete", "isotype", "criterion"),
    ("c_call_enough_evidence", "c_call_enough_evidence", "isotype", "criterion"),
    ("c_call_callable", "isotype_callable", "isotype", "masking"),
    ("isotype_class_callable", "class_callable", "isotype", "eligibility"),
    ("isotype_subclass_callable", "subclass_callable", "isotype", "eligibility"),
    ("shm_measured_records", "shm_measured", "shm", "eligibility"),
    ("shm_supported_positions", "shm_supported_positions", "shm", "eligibility"),
    ("shm_eligible_members", "shm_eligible_members", "shm", "eligibility"),
    ("selection_input", "selection_input", "selection", "eligibility"),
)
# Audit stages the trail leaves out because they repeat counts released
# elsewhere. provenance.json records the assembly history and clones.tsv the
# donor clones. shm.tsv counts the SHM-eligible clones, selection_pressure.tsv
# the selection clones and diversity.tsv the diversity repertoire.
# inference_representatives repeats the exact_sequence_collapse output, and
# selection_input_source_support counts the selection_input members again.
TRAIL_OMITTED = frozenset({
    "assembly_history", "inference_representatives", "donor_clones", "shm_eligible_clones",
    "selection_input_source_support", "selection_effective_clones", "selection_usable_pdfs",
    "diversity_repertoire"})
TRAIL_UNITS = {
    "contig": "contig", "cdr3_call": "cdr3_call", "rearrangement": "airr_record",
    "backbone_representative": "contig", "inference_representative": "member",
    "clone_member": "member", "selection_input_record": "member",
    "library_clone_observation": "clone", "clonotype": "clonotype",
    "record_position": "position",
}
LOCUS_SOURCES = {"native": "trust4", "annotated": "igblast", "not_applicable": ""}


def read_rows(path: Path, header: Sequence[str] | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"{path}: missing or repeated header")
        if header is not None and tuple(reader.fieldnames) != tuple(header):
            raise ValueError(f"{path}: QC audit header differs from the producer contract")
        rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ValueError(f"{path}: malformed table row")
    return rows


def number(value: str) -> Decimal | None:
    if value == "":
        return None
    try:
        result = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"Invalid audit number {value!r}") from error
    if not result.is_finite() or result < 0:
        raise ValueError(f"Invalid audit number {value!r}")
    return result


def validate_rows(rows: Sequence[Mapping[str, str]]) -> None:
    keys = set()
    for row in rows:
        key = tuple(row[field] for field in KEY)
        if key in keys:
            raise ValueError(f"Duplicate QC audit key: {key}")
        keys.add(key)
        if row["grain"] not in {"library", "donor"}:
            raise ValueError("QC audit grain must be library or donor")
        if row["locus"] not in {*LOCI, "all"} or row["locus_basis"] not in {
                "native", "annotated", "not_applicable"}:
            raise ValueError("Invalid QC audit locus context")
        if row["kind"] not in {"checkpoint", "criterion", "eligibility", "grouping", "diagnostic"}:
            raise ValueError("Invalid QC audit kind")
        if row["decision_status"] not in {
                "applied", "disabled", "not_applicable", "unknown_history", "unavailable"}:
            raise ValueError("Invalid QC audit decision status")
        values = {field: number(row[field]) for field in COUNTS}
        for field in ("stage_index", "n_units_in", "n_units_out"):
            value = number(row[field])
            if value is not None and value != value.to_integral_value():
                raise ValueError(f"Non-integer audit unit count: {field}")
        if row["kind"] == "checkpoint" and values["n_units_in"] != values["n_units_out"]:
            raise ValueError("QC audit checkpoint changed unit counts")
        if row["support_conserved"] not in {"TRUE", "FALSE", ""}:
            raise ValueError("Invalid support-conservation assertion")
        if row["support_conserved"] == "TRUE" and (
                values["consensus_count_in"] is None or
                values["consensus_count_in"] != values["consensus_count_out"]):
            raise ValueError("QC audit operation did not preserve defined support")
        if (row["input_unit"] == row["output_unit"] and
                row["kind"] in {"criterion", "eligibility", "grouping", "diagnostic"} and
                values["n_units_in"] is not None and values["n_units_out"] is not None and
                values["n_units_out"] > values["n_units_in"]):
            raise ValueError("QC audit subset exceeds its input population")


def contextualise(row: Mapping[str, str], manifest: Mapping[str, Mapping[str, str]],
                  scope_id: str) -> dict[str, str]:
    row = dict(row)
    if row["scope_id"] and row["scope_id"] != scope_id:
        raise ValueError("QC audit belongs to another scope")
    row["scope_id"] = scope_id
    if row["grain"] == "library":
        library = row["library_id"]
        if library not in manifest:
            raise ValueError("QC audit library is outside the scope")
        for field in ("donor_id", "tissue_id"):
            expected = manifest[library][field]
            if row[field] and row[field] != expected:
                raise ValueError(f"QC audit disagrees with manifest {field}")
            row[field] = expected
    elif (row["library_id"] or row["tissue_id"] or
          row["donor_id"] not in {item["donor_id"] for item in manifest.values()}):
        raise ValueError("QC audit donor row has invalid identity context")
    return row


def checkpoint(library: str, stage: str, count: str, support: str, *,
               locus: str, analysis: str, index: int, unit: str,
               isotype: str = "", region: str = "") -> dict[str, str]:
    row = dict.fromkeys(COLUMNS, "")
    row.update(library_id=library, grain="library", locus=locus, locus_basis="annotated",
               analysis=analysis, stage=stage, stage_index=str(index), kind="checkpoint",
               decision_status="applied", input_unit=unit, output_unit=unit,
               n_units_in=count, n_units_out=count,
               consensus_count_in=support, consensus_count_out=support,
               support_conserved="TRUE" if support else "", isotype_group=isotype, v_region=region)
    return row


def canonical_totals(rows: list[dict[str, str]]) -> None:
    """Add the canonical-locus total for the grouping stages that arrive split.

    IGH and non-IGH grouping are prepared by separate producers, so neither can
    write a seven-locus total on its own. Every other stage writes its own total,
    and where one is present it must equal the sum recomputed here. A total is
    added only when all seven loci are present, because a partial partition would
    understate the population.
    """
    groups = {}
    context = tuple(field for field in COLUMNS if field not in
                    {"locus", *COUNTS})
    for row in rows:
        if row["grain"] == "library":
            groups.setdefault(tuple(row[field] for field in context), []).append(row)
    additions = []
    for group in groups.values():
        by_locus = {row["locus"]: row for row in group}
        if len(by_locus) != len(group):
            raise ValueError("Repeated locus in QC audit partition")
        split_stage = group[0]["stage"] in {"grouping_eligible", "exact_sequence_collapse", "inference_representatives"}
        if split_stage and set(by_locus) - {"all"} != set(LOCI):
            raise ValueError("Grouping audit lacks a complete canonical-locus partition")
        if set(LOCI).issubset(by_locus):
            total = dict(group[0], locus="all")
            for field in COUNTS:
                values = [number(by_locus[locus][field]) for locus in LOCI]
                if all(value is None for value in values):
                    total[field] = ""
                elif any(value is None for value in values):
                    raise ValueError("Audit total mixes defined and unavailable counts")
                else:
                    total[field] = str(sum(values))
            if "all" in by_locus:
                if any(number(by_locus["all"][field]) != number(total[field]) for field in COUNTS):
                    raise ValueError("QC audit canonical-locus total disagrees")
            elif split_stage:
                additions.append(total)
    rows.extend(additions)


def fraction(numerator: Decimal | None, denominator: Decimal | None) -> str:
    """Divide two support counts exactly and keep 15 significant digits.

    The fraction is empty when either count is unavailable or the denominator is
    zero. Zero is an observed zero numerator over positive support.
    """
    if numerator is None or denominator is None or denominator == 0:
        return ""
    if not 0 <= numerator <= denominator:
        raise ValueError("A QC support fraction must lie between zero and one")
    return format(numerator / denominator, ".15g")


def grouping_reasons(status: str) -> frozenset[str]:
    """Name the GROUPING_REASONS one grouping_status records.

    An eligible record has none. Any v_call_ or j_call_ token, which names the
    V or J status that failed, is an unresolved call. An empty or unknown token
    stops the collection.
    """
    if status == "eligible":
        return frozenset()
    reasons = set()
    for token in status.split(";"):
        if token.startswith("v_call_") and len(token) > len("v_call_"):
            reasons.add("v_call_unresolved")
        elif token.startswith("j_call_") and len(token) > len("j_call_"):
            reasons.add("j_call_unresolved")
        elif token in REASON_TOKENS:
            reasons.add(REASON_TOKENS[token])
        else:
            raise ValueError(f"Unknown grouping_status reason {token!r}")
    return frozenset(reasons)


def export_support(path: Path) -> dict[str, dict[str, Decimal]]:
    """Sum one library export's read support by locus, eligibility and reason.

    The export holds the library's common-QC records at the locus IgBLAST
    assigned. Only locus, consensus_count and grouping_status are read, one row
    at a time. Each locus returns its common-QC support, its grouping-eligible
    support and the support of the records naming each grouping reason, which
    can overlap. A locus with no records has zero sums.
    """
    totals = {locus: dict.fromkeys(("common_qc", "grouping_eligible", *GROUPING_REASONS), Decimal(0))
              for locus in LOCI}
    reasons_of: dict[str, frozenset[str]] = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        columns = ("locus", "consensus_count", "grouping_status")
        if not header or any(header.count(column) != 1 for column in columns):
            raise ValueError(f"{path.name}: AIRR export lacks one {', '.join(columns)} column each")
        indices = [header.index(column) for column in columns]
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"{path.name}: malformed AIRR export row")
            locus, support, status = (row[index] for index in indices)
            if locus not in totals:
                raise ValueError(f"{path.name}: AIRR export record has locus {locus!r}")
            value = number(support)
            if value is None:
                raise ValueError(f"{path.name}: AIRR export record lacks read support")
            if status not in reasons_of:
                reasons_of[status] = grouping_reasons(status)
            sums = totals[locus]
            sums["common_qc"] += value
            if not reasons_of[status]:
                sums["grouping_eligible"] += value
            for reason in reasons_of[status]:
                sums[reason] += value
    return totals


def _populations(rows: Sequence[Mapping[str, str]]) -> dict[tuple[str, str, str], Mapping[str, str]]:
    """Index unstratified library rows at IgBLAST's locus by library, locus and stage."""
    return {(row["library_id"], row["locus"], row["stage"]): row for row in rows
            if row["grain"] == "library" and row["locus"] in LOCI and row["locus_basis"] == "annotated"
            and not row["isotype_group"] and not row["v_region"]}


def qc_table(rows: Sequence[Mapping[str, str]], exports: Mapping[str, Path],
             manifest: Mapping[str, Mapping[str, str]], scope_id: str) -> list[dict[str, str]]:
    """Build qc.tsv from validated audit rows and each library's AIRR export.

    The support counts come from the audit: the qc_productive input is the
    IgBLAST-annotated population, qc_pass_total the common-QC output and
    grouping_eligible the grouping-eligible output. The common-QC, grouping and
    isotype populations must each start from the one before, and each export's
    common-QC and grouping-eligible support must equal the audit's. The reason
    fractions come from the exports. Class and subclass fractions apply at IGH
    only and are empty elsewhere. Rows follow library ID, then LOCI.
    """
    if set(exports) != set(manifest):
        raise ValueError("The QC table requires one AIRR export for every manifest library")
    populations = _populations(rows)
    native = {(row["library_id"], row["locus"]): row for row in rows
              if row["stage"] == "trust4_reported_airr" and row["locus_basis"] == "native"}

    def population(library: str, locus: str, stage: str) -> Mapping[str, str]:
        if (library, locus, stage) not in populations:
            raise ValueError(f"QC table lacks the {stage} population of {library} {locus}")
        return populations[(library, locus, stage)]

    def support(row: Mapping[str, str], column: str) -> Decimal | None:
        return number(row[column])

    result = []
    for library in sorted(manifest):
        sums = export_support(exports[library])
        for locus in LOCI:
            productive = population(library, locus, "qc_productive")
            common = population(library, locus, "qc_pass_total")
            grouping = population(library, locus, "grouping_eligible")
            annotated = support(productive, "consensus_count_in")
            common_qc = support(common, "consensus_count_out")
            eligible = support(grouping, "consensus_count_out")
            if (support(common, "consensus_count_in") != annotated or
                    support(grouping, "consensus_count_in") != common_qc):
                raise ValueError(f"QC populations of {library} {locus} do not follow one another")
            if sums[locus]["common_qc"] != common_qc or sums[locus]["grouping_eligible"] != eligible:
                raise ValueError(f"AIRR export support of {library} {locus} differs from the audit")
            failed = support(productive, "consensus_count_out")
            row = {"scope_id": scope_id, "library_id": library,
                   "tissue_id": manifest[library]["tissue_id"],
                   "donor_id": manifest[library]["donor_id"], "locus": locus,
                   "trust4_record_count": native.get((library, locus), {}).get("n_units_out", ""),
                   "trust4_consensus_count": native.get((library, locus), {}).get("consensus_count_out", ""),
                   "igblast_consensus_count": productive["consensus_count_in"],
                   "nonproductive_fraction": fraction(
                       None if annotated is None or failed is None else annotated - failed, annotated),
                   "common_qc_consensus_count": common["consensus_count_out"],
                   "common_qc_fraction_of_annotated": fraction(common_qc, annotated),
                   "grouped_consensus_count": grouping["consensus_count_out"],
                   "grouping_eligible_fraction_of_common_qc": fraction(eligible, common_qc),
                   **{reason + "_fraction": fraction(sums[locus][reason], common_qc)
                      for reason in GROUPING_REASONS},
                   "class_callable_fraction": "", "subclass_callable_fraction": ""}
            if locus == "IGH":
                for level in ("class", "subclass"):
                    callable_row = population(library, locus, f"isotype_{level}_callable")
                    if support(callable_row, "consensus_count_in") != common_qc:
                        raise ValueError(f"Isotype {level} input of {library} differs from common QC")
                    row[f"{level}_callable_fraction"] = fraction(
                        support(callable_row, "consensus_count_out"), common_qc)
            row.update(dict.fromkeys(QC_COLUMNS[QC_COLUMNS.index("n_grouped_igh_members"):], ""))
            result.append(row)
    return result


def qc_trail(rows: Sequence[Mapping[str, str]],
             value_orders: Mapping[str, Sequence[str]]) -> list[dict[str, str]]:
    """Build the released QC trail from validated, contextualised audit rows.

    Donor rows and the TRAIL_OMITTED stages are left out, and every other stage
    must appear in TRAIL_STAGES. Units and locus bases take their released names.
    A receptor_units row reads its input members from the exact_sequence_collapse
    output of the same library and locus, whose support it must carry. Rows sort
    by library, by locus with all first and then LOCI, by stage_order, and by
    isotype group and V region in the orders `value_orders` lists.
    """
    orders = {name: order for order, name in
              enumerate(dict.fromkeys(name for _, name, _, _ in TRAIL_STAGES), start=1)}
    stages = {source: (orders[name], name, analysis, kind)
              for source, name, analysis, kind in TRAIL_STAGES}
    members = {(row["library_id"], row["locus"]): row for row in rows
               if row["grain"] == "library" and row["stage"] == "exact_sequence_collapse"}
    result = []
    for row in rows:
        if row["grain"] != "library" or row["stage"] in TRAIL_OMITTED:
            continue
        if row["stage"] not in stages:
            raise ValueError(f"QC trail has no released name for audit stage {row['stage']!r}")
        order, name, analysis, kind = stages[row["stage"]]
        incoming = (row["input_unit"], row["n_units_in"], row["consensus_count_in"])
        if name == "receptor_units":
            source = members.get((row["library_id"], row["locus"]))
            if source is None or number(source["consensus_count_out"]) != number(row["consensus_count_in"]):
                raise ValueError(f"Receptor units of {row['library_id']} {row['locus']} "
                                 "differ from their exact-sequence members")
            incoming = (source["output_unit"], source["n_units_out"], source["consensus_count_out"])
        units = (incoming[0], row["output_unit"])
        if any(unit not in TRAIL_UNITS for unit in units):
            raise ValueError(f"QC trail has no released unit for stage {row['stage']!r}")
        result.append({
            **{field: row[field] for field in ("scope_id", "library_id", "tissue_id", "donor_id", "locus")},
            "stage_order": str(order), "stage": name, "analysis": analysis, "kind": kind,
            "decision_status": row["decision_status"],
            "isotype_group": row["isotype_group"], "v_region": row["v_region"],
            "locus_source": LOCUS_SOURCES[row["locus_basis"]],
            "input_unit": TRAIL_UNITS[units[0]], "n_units_in": incoming[1],
            "output_unit": TRAIL_UNITS[units[1]], "n_units_out": row["n_units_out"],
            "consensus_count_in": incoming[2], "consensus_count_out": row["consensus_count_out"],
        })
    loci = {locus: rank for rank, locus in enumerate(("all", *LOCI))}
    ranks = {column: {"": 0, **{value: rank for rank, value in
                                enumerate(value_orders.get(column, ()), start=1)}}
             for column in ("isotype_group", "v_region")}
    unlisted = {row[column] for row in result for column in ranks if row[column] not in ranks[column]}
    if unlisted:
        raise ValueError(f"QC trail strata have no listed order: {sorted(unlisted)}")
    result.sort(key=lambda row: (row["library_id"], loci[row["locus"]], int(row["stage_order"]),
                                 ranks["isotype_group"][row["isotype_group"]],
                                 ranks["v_region"][row["v_region"]]))
    keys = {tuple(row[field] for field in ("library_id", "locus", "stage", "isotype_group", "v_region"))
            for row in result}
    if len(keys) != len(result):
        raise ValueError("QC trail repeats a library, locus, stage and stratum")
    return result


def collect_scope_qc_audit(manifest_file: Path, shards: Sequence[Path],
                           diversity_files: Sequence[Path], selection_file: Path,
                           output: Path, producer_libraries: Iterable[str] = (), *,
                           exports: Mapping[str, Path], qc_output: Path, trail_output: Path,
                           value_orders: Mapping[str, Sequence[str]],
                           coverage_file: Path | None = None) -> None:
    """Write qc_audit.tsv, qc.tsv and qc_trail.tsv for one scope.

    `shards` are the library and donor audit tables their producers wrote, in
    the COLUMNS header. Their rows gain scope and manifest context, and split
    grouping stages gain their canonical-locus totals. `diversity_files` add one
    diversity-population checkpoint per library and locus, and `selection_file`
    adds the selection clone populations. `producer_libraries` are the libraries
    with a recorded assembly producer. Every other manifest library gets an
    assembly-history row with status unknown_history. `exports` maps each
    manifest library to its AIRR export for qc_table(), and `value_orders` gives
    the trail's isotype group and V region orders for qc_trail(). Every table is
    built and validated before any is written.
    """
    manifest_rows = read_rows(manifest_file)
    manifest = {row["library_id"]: row for row in manifest_rows}
    scopes = {row["scope_id"] for row in manifest_rows}
    if (not manifest or len(scopes) != 1 or len(manifest) != len(manifest_rows) or
            any(not row[field] for row in manifest_rows for field in
                ("scope_id", "library_id", "donor_id", "tissue_id"))):
        raise ValueError("Audit requires one scope and unique manifest library identities")
    scope_id = next(iter(scopes))
    rows = [contextualise(row, manifest, scope_id)
            for path in shards for row in read_rows(path, COLUMNS)]
    canonical_totals(rows)
    assigned = {row["library_id"] for row in rows if row["stage"] == "clone_observations"}
    for path in diversity_files:
        for metric in read_rows(path):
            unit = "library_clone_observation" if metric["locus"] == "IGH" else "clonotype"
            row = checkpoint(metric["library_id"], "diversity_repertoire", metric["richness"],
                             metric["grouped_consensus_count"], locus=metric["locus"], analysis="diversity",
                             index=80, unit=unit)
            rows.append(contextualise(row, manifest, scope_id))
            # Assignment cannot name an absent library from an empty donor pool.
            # Its authoritative IGH metric supplies that observed empty population.
            if metric["locus"] == "IGH" and metric["library_id"] not in assigned:
                if number(metric["richness"]) != 0 or number(metric["grouped_consensus_count"]) != 0:
                    raise ValueError("Nonempty IGH repertoire lacks assignment audit")
                row.update(stage="clone_observations", analysis="grouping", stage_index="34")
                rows.append(contextualise(row, manifest, scope_id))
    for metric in read_rows(selection_file):
        for stage, column, region, index in (
                ("selection_effective_clones", "n_clone_consensus_sequences", "", 71),
                ("selection_usable_pdfs", "n_clone_consensus_sequences_used_cdr", "cdr", 72),
                ("selection_usable_pdfs", "n_clone_consensus_sequences_used_fwr", "fwr", 72)):
            row = checkpoint(metric["library_id"], stage, metric[column], "", locus="IGH",
                             analysis="selection", index=index, unit="effective_clone",
                             isotype=metric["isotype_group"], region=region)
            rows.append(contextualise(row, manifest, scope_id))
    known = set(producer_libraries)
    if not known.issubset(manifest):
        raise ValueError("Assembly producer identity is outside the scope")
    for library in manifest:
        row = checkpoint(library, "assembly_history", "", "", locus="all",
                         analysis="assembly", index=0, unit="library")
        row.update(locus_basis="not_applicable",
                   decision_status="applied" if library in known else "unknown_history")
        rows.append(contextualise(row, manifest, scope_id))
    validate_rows(rows)
    # One deterministic order for the working table. stage_index orders display
    # within a library; it does not describe a sequence of losses.
    rows.sort(key=lambda row: (row["donor_id"], row["library_id"], int(row["stage_index"]),
                              row["locus"], row["isotype_group"], row["v_region"]))
    qc = qc_table(rows, exports, manifest, scope_id)
    if coverage_file is not None:
        add_shm_coverage_counts(qc, coverage_file)
    trail = qc_trail(rows, value_orders)
    for path, columns, records in ((output, SCOPE_COLUMNS, rows), (qc_output, QC_COLUMNS, qc),
                                   (trail_output, TRAIL_COLUMNS, trail)):
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        try:
            with temporary.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(records)
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)


def add_shm_coverage_counts(qc: list[dict[str, str]], path: Path) -> None:
    """Count whole-V and any-subregion eligibility once per contiguous member."""
    counts = {row["library_id"]: [0, 0, 0, 0] for row in qc if row["locus"] == "IGH"}
    previous = None
    regions = set()
    eligible_any = False
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["library_id"], row["sequence_id"])
            if key != previous:
                if previous is not None:
                    if regions != {"whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3"}:
                        raise ValueError("SHM coverage requires six regions per member")
                    counts[previous[0]][3] += int(eligible_any)
                previous, regions, eligible_any = key, set(), False
            region = row["v_region"]
            if region in regions or key[0] not in counts:
                raise ValueError("SHM coverage repeats a region or names an unknown library")
            regions.add(region)
            eligible = row["shm_eligible"] == "TRUE"
            if region == "whole_v":
                counts[key[0]][0] += 1
                counts[key[0]][1] += int(row["measurement_status"] != "unmeasured")
                counts[key[0]][2] += int(eligible)
            else:
                eligible_any |= eligible
    if previous is not None:
        if regions != {"whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3"}:
            raise ValueError("SHM coverage requires six regions per member")
        counts[previous[0]][3] += int(eligible_any)
    for row in qc:
        if row["locus"] != "IGH":
            continue
        grouped, measured, whole, any_region = counts[row["library_id"]]
        row.update(n_grouped_igh_members=str(grouped), n_shm_measured_members_whole_v=str(measured),
                   n_shm_eligible_members_whole_v=str(whole),
                   shm_eligible_fraction_of_grouped_members_whole_v=str(whole / grouped) if grouped else "",
                   n_shm_eligible_members_any_subregion=str(any_region),
                   shm_eligible_fraction_of_grouped_members_any_subregion=str(any_region / grouped) if grouped else "")
