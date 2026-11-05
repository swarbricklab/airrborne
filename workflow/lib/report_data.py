"""Prepare one scope's authoritative tables for the static report.

Measurements are copied without recalculation. The headline must retain the
complete manifest. Large per-library tables stay in separate payloads so opening
the first chapter does not load every clone or rarefaction point. Gene-use
payloads are separate so group comparisons do not load clone details. The clone
SHM payload holds only the rows in_clone_shm_payload() names.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections.abc import Mapping
from pathlib import Path

from .release_bundle import release_path
from .report_metadata import load_library_metadata, load_tissue_metadata, join_scope_metadata

IDENTIFIERS = {
    "scope_id",
    "library_id",
    "tissue_id",
    "donor_id",
    "library_id_1",
    "library_id_2",
    "tissue_id_1",
    "tissue_id_2",
    "clone_id",
    "group_key",
}
LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
# The qc.tsv fractions, each over a support count in its own row.
QC_FRACTIONS = (
    "nonproductive_fraction", "common_qc_fraction_of_annotated",
    "grouping_eligible_fraction_of_common_qc",
    "v_call_unresolved_fraction", "j_call_unresolved_fraction", "junction_ineligible_fraction",
    "c_call_ineligible_fraction", "contig_representative_ineligible_fraction",
    "class_callable_fraction", "subclass_callable_fraction",
)
NUMERIC_COLUMNS = {
    "support_depth",
    "group_rank",
    "richness",
    "hill_shannon",
    "inverse_simpson",
    "rarefied_inverse_simpson",
    "consensus_count",
    "grouped_consensus_count",
    "support_fraction",
    "rarefaction_target",
    "min_consensus_count_for_diversity",
    "min_consensus_count_per_library",
    "shannon_entropy",
    "clonality",
    "top_1_support_fraction",
    "top_5_support_fraction",
    "top_10_support_fraction",
    "rarefied_richness",
    "rarefied_hill_shannon",
    "inverse_simpson_asymptote",
    "grouped_consensus_count_1",
    "grouped_consensus_count_2",
    "n_groups_1",
    "n_groups_2",
    "n_groups_shared",
    "jaccard",
    "morisita_horn",
    "shared_consensus_count_1",
    "shared_consensus_count_2",
    "shared_support_fraction_1",
    "shared_support_fraction_2",
    "bray_curtis_proportion",
    "n_shared_higher_in_2",
    "n_shared_higher_in_1",
    "n_shared_equal",
    "n_members_measured",
    "n_members_eligible",
    "n_members_zero_shm",
    "n_clone_consensus_sequences",
    "n_groups",
    "n_clones_eligible",
    "n_clones_zero_median_shm",
    "mean_clone_median_shm",
    "median_clone_median_shm",
    "q25_clone_median_shm",
    "q75_clone_median_shm",
    "iqr_clone_median_shm",
    "fraction_clones_zero_median_shm",
    "median_trust4_consensus_base_votes",
    "n_clone_consensus_sequences_used_cdr",
    "n_clone_consensus_sequences_used_fwr",
    "sigma_cdr",
    "sigma_cdr_ci_low",
    "sigma_cdr_ci_high",
    "sigma_fwr",
    "sigma_fwr_ci_low",
    "sigma_fwr_ci_high",
    "stage_order",
    "n_units_in",
    "n_units_out",
    "consensus_count_in",
    "consensus_count_out",
    "common_qc_consensus_count",
    "igblast_consensus_count",
    "grouped_consensus_count",
    *QC_FRACTIONS,
    "callable_consensus_count",
    "group_fraction",
    "median_shm_frequency",
    "eligible_consensus_count",
    "n_members",
    "n_members_eligible",
    "n_members_zero_shm",
    "n_members_positive_shm",
    "distance_bin_start",
    "distance_bin_end",
    "distance_bin_midpoint",
    "n_distances_in_bin",
    "observed_density",
    "near_component_density",
    "far_component_density",
    "clone_threshold",
    "n_donors",
    "n_donors_fitted",
    "n_libraries",
    "n_distances",
    "n_distances_fitted",
    "fitted_threshold",
    "near_weight",
    "near_shape",
    "near_scale",
    "far_weight",
    "far_mean",
    "far_sd",
    "absolute_log_likelihood",
    "fit_sensitivity",
    "n_unique_sequences",
}


SELECTION_CONTEXT_COLUMNS = {f"n_clone_consensus_sequences_with_positions_{region}" for region in ("cdr", "fwr")}
NUMERIC_COLUMNS |= SELECTION_CONTEXT_COLUMNS


def metric_definitions() -> list[dict]:
    metrics = []
    for locus in LOCI:
        metrics.append({
            "column": f"{locus.lower()}_yield_per_million_read_pairs",
            "label": f"{locus} receptor yield",
            "unit": "support per million input pairs",
            "source": "yield.tsv",
            "section": "signal",
        })
    for locus in ("igh", "trb"):
        for field, label, unit in (
            ("rarefied_inverse_simpson", "target-depth inverse Simpson", "effective groups"),
            ("inverse_simpson", "observed inverse Simpson", "effective groups"),
            ("top_5_support_fraction", "observed top-5 support share", "fraction"),
            ("top_10_support_fraction", "observed top-10 support share", "fraction"),
            ("grouped_consensus_count", "grouping-eligible support", "support units"),
            ("richness", "observed group count", "groups"),
        ):
            metrics.append({
                "column": f"{locus}_{field}",
                "label": f"{locus.upper()} {label}",
                "unit": unit,
                "source": "diversity.tsv",
                "reason_column": f"{locus}_rarefaction_status" if field == "rarefied_inverse_simpson" else f"{locus}_diversity_status",
                "section": "diversity",
            })
    for isotype in ("igm", "igd", "igg", "iga", "ige"):
        metrics.append({
            "column": f"class_{isotype}_fraction",
            "label": f"Ig{isotype[2:].upper()} class share",
            "unit": "class-callable support fraction",
            "source": "isotype_usage.tsv",
            "section": "antibody",
        })
    metrics.append({
        "column": "log2_igg_iga",
        "label": "Log2 IgG:IgA",
        "unit": "log2 support ratio",
        "source": "headline.tsv",
        "reason_column": "igg_iga_status",
        "section": "antibody",
    })
    for isotype, label in (("igh", "All IGH"), ("igm", "IgM"), ("igg", "IgG"), ("iga", "IgA")):
        for field, suffix, unit in (("mean_clone_median_shm", "mean clone-median whole-V SHM", "mutation frequency"), ("fraction_clones_zero_median_shm", "clones with zero median whole-V SHM", "eligible-clone fraction")):
            metrics.append({
                "column": f"{isotype}_{field}_whole_v",
                "label": f"{label} {suffix}",
                "unit": unit,
                "source": "shm.tsv",
                "section": "antibody",
            })
    metrics.append({
        "column": "input_read_pairs",
        "label": "Input read pairs",
        "unit": "read pairs",
        "source": "yield.tsv",
        "section": "signal",
    })
    for metric in metrics:
        metric["detail_source"] = metric["source"]
        metric["source"] = "headline.tsv" if metric["column"] in HEADLINE_FIELDS else metric["source"]
        column = metric["column"]
        metric["format"] = "percent" if any(token in column for token in ("fraction", "_shm")) else "number"
        metric["scale"] = "log" if any(token in column for token in ("yield_per_", "inverse_simpson", "consensus_count", "richness", "input_read_pairs")) else "linear"
        metric["bounded"] = "fraction" in column
    return metrics


# Library values the report shows beside the headline, copied from their topic
# tables onto each library row. The headline selects summaries; these stay in
# diversity.tsv, shm.tsv and isotype_usage.tsv.
CONTEXT_DIVERSITY_FIELDS = ("diversity_status", "inverse_simpson", "top_10_support_fraction",
                            "hill_shannon", "top_1_support_fraction")
CONTEXT_SHM_FIELDS = ("mean_clone_median_shm", "fraction_clones_zero_median_shm", "n_clones_eligible",
                      "n_clones_zero_median_shm", "n_members_eligible", "n_members_measured")
CONTEXT_ISOTYPE_CLASSES = ("IgM", "IgG", "IgA")
# Columns the report reads from the headline itself.
from workflow.lib.headline import HEADLINE_COLUMNS

HEADLINE_FIELDS = set(HEADLINE_COLUMNS) - {"scope_id", "library_id", "tissue_id", "donor_id"}


def add_library_context(scope_dir: Path, libraries: list[dict]) -> None:
    """Copy the report's non-headline library values from their topic tables.

    Each library row gains `<locus>_<field>` from diversity.tsv for IGH and TRB,
    `<group>_<field>_whole_v` from the whole-V rows of shm.tsv, and `<class>_consensus_count`,
    the IgM, IgG and IgA class support of isotype_usage.tsv. A value present in
    the headline keeps its headline copy.
    """
    by_id = {row["library_id"]: row for row in libraries}
    metrics = scope_dir / "metrics"

    def put(library_id, column, value):
        row = by_id.get(library_id)
        if row is None:
            raise ValueError(f"A topic table names an out-of-scope library: {library_id}")
        row.setdefault(column, value)

    for source in read_table(metrics / "yield.tsv"):
        put(source["library_id"], "yield_status", source["yield_status"])
    for source in read_table(metrics / "qc.tsv"):
        if source["locus"] in ("IGH", "TRB"):
            for field in ("common_qc_fraction_of_annotated", "grouping_eligible_fraction_of_common_qc"):
                put(source["library_id"], f"{source['locus'].lower()}_{field}", source[field])
    for source in read_table(metrics / "isotype_summary.tsv"):
        put(source["library_id"], "igg_iga_status", source["igg_iga_status"])
    for row in libraries:
        numerator, denominator = row["igh_n_shm_eligible_clones_whole_v"], row["igh_richness"]
        row["igh_shm_eligible_clone_fraction_whole_v"] = (
            numerator / denominator if numerator is not None and denominator else None)
    for source in read_table(metrics / "diversity.tsv"):
        if source["locus"] in ("IGH", "TRB"):
            for field in CONTEXT_DIVERSITY_FIELDS:
                put(source["library_id"], f"{source['locus'].lower()}_{field}", source[field])
    for source in read_table(metrics / "shm.tsv"):
        if source["v_region"] == "whole_v":
            for field in CONTEXT_SHM_FIELDS:
                put(source["library_id"], f"{('igh' if source['isotype_group'] == 'all_IGH' else source['isotype_group'].lower())}_{field}_whole_v", source[field])
    for source in read_table(metrics / "isotype_usage.tsv"):
        if source["isotype_level"] == "class" and source["isotype"] in CONTEXT_ISOTYPE_CLASSES:
            put(source["library_id"], f"{source['isotype'].lower()}_consensus_count",
                source["consensus_count"])


def read_table(path: Path, numeric: set[str] | None = None,
               categories: set[str] | None = None) -> list[dict]:
    """Read report rows with string identities and explicit missing numeric cells.

    Empty non-identifier cells become None. Literal NA remains a string in
    identifiers and declared categories. Numeric columns require finite values.
    """
    numeric = NUMERIC_COLUMNS | (numeric or set())
    categories = categories or set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"{path.name} has a missing or repeated header.")
        rows = []
        for line, source in enumerate(reader, start=2):
            row = {}
            for column, value in source.items():
                if column is None or value is None:
                    raise ValueError(f"{path.name} contains a malformed row.")
                # Identifiers retain literal NA and empty strings. Categories
                # retain literal NA but empty cells are missing. Other columns
                # treat both empty cells and NA as missing.
                if column not in IDENTIFIERS and (value == "" or (value == "NA" and column not in categories)):
                    value = None
                elif column in numeric:
                    try:
                        value = float(value)
                    except ValueError as error:
                        raise ValueError(
                            f"{path}:{line}: column {column} must contain a number or a missing value; "
                            f"found {value!r}."
                        ) from error
                    if not math.isfinite(value):
                        raise ValueError(f"{path}:{line}: column {column} must be finite or missing.")
                    if value.is_integer():
                        value = int(value)
                elif column in ("is_target", "is_endpoint"):
                    if value not in ("TRUE", "FALSE", "T", "F"):
                        raise ValueError(f"{path.name}.{column} must be a boolean.")
                    value = value in ("TRUE", "T")
                row[column] = value
            rows.append(row)
        return rows


def require_columns(path: Path, columns: list[str]) -> None:
    """Require a TSV header containing every named column, including empty tables."""
    with path.open(newline="", encoding="utf-8") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if not header:
        raise ValueError(f"{path} must contain a tab-separated header.")
    missing = set(columns) - set(header)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {', '.join(sorted(missing))}.")


def unique_rows(rows: list[dict], keys: list[str], source: str, nullable_keys: tuple[str, ...] = ()) -> dict:
    """Index rows by unique key tuples; only declared nullable keys may be missing."""
    result = {}
    for row in rows:
        key = tuple(row[column] for column in keys)
        missing_key = any(
            (value is None or value == "") and column not in nullable_keys
            for column, value in zip(keys, key)
        )
        if missing_key or key in result:
            raise ValueError(f"{source} requires unique, present {', '.join(keys)} keys.")
        result[key] = row
    return result


def in_clone_shm_payload(row: Mapping[str, object]) -> bool:
    """True for a clone_shm row the report reads.

    The library detail in report_sections.js plots whole V for all_IGH and lists
    whole V for IgM, IgG and IgA. The released table keeps every region and group,
    and check_release.py applies this filter before comparing it with the payload.
    """
    return row["v_region"] == "whole_v" and row["isotype_group"] in ("all_IGH", "IgM", "IgG", "IgA")


def pack(rows: list[dict]) -> dict:
    columns = list(rows[0]) if rows else []
    return {"columns": columns, "rows": [[row.get(column) for column in columns] for row in rows]}


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=True, allow_nan=False, separators=(",", ":")) + "\n", encoding="utf-8")


def _prepare_libraries(scope_dir: Path, metadata_sources: Mapping | None) -> tuple:
    """Validate the complete headline population and attach optional report context.

    Return the scope identity, library rows, metric and metadata descriptors,
    timepoint orders and default grouping. Metadata never changes source tables.
    """
    manifest_path = scope_dir / "manifest" / "scope_manifest.tsv"
    headline_path = scope_dir / "metrics" / "headline.tsv"
    require_columns(manifest_path, ["scope_id", "library_id", "tissue_id", "donor_id"])
    manifest = read_table(manifest_path)
    manifest_by_id = unique_rows(manifest, ["library_id"], manifest_path.name)
    scopes = {row["scope_id"] for row in manifest}
    if len(scopes) != 1:
        raise ValueError("A report requires one nonempty scope manifest.")
    scope_id = next(iter(scopes))
    metrics = metric_definitions()
    require_columns(headline_path, ["scope_id", "library_id", "tissue_id", "donor_id", *sorted(HEADLINE_FIELDS)])
    headline_numeric = {field for field in HEADLINE_FIELDS if not field.endswith("_status")}
    libraries = read_table(headline_path, headline_numeric)
    if any(column.startswith("meta_") for row in libraries for column in row):
        raise ValueError("The headline must contain scientific measurements only; rebuild it before preparing the report.")
    by_id = unique_rows(libraries, ["library_id"], headline_path.name)
    if set(by_id) != set(manifest_by_id):
        raise ValueError("The headline must retain exactly the scope manifest's libraries.")
    for key, row in by_id.items():
        if any(row[field] != manifest_by_id[key][field] for field in ("scope_id", "tissue_id", "donor_id")):
            raise ValueError(f"The headline identity differs from the manifest for {key[0]}.")
    add_library_context(scope_dir, libraries)
    for metric in metrics:
        if any(metric["column"] not in row for row in libraries):
            raise ValueError(f"No library value supplies the report metric {metric['column']}.")
    sources = metadata_sources or {}
    if set(sources) - {"tissue", "library"}:
        raise ValueError("Report metadata sources must be tissue or library tables")
    tissue = sources.get("tissue")
    library = sources.get("library")
    libraries, declarations = join_scope_metadata(
        libraries,
        load_tissue_metadata(Path(tissue["file"]), tissue.get("columns")) if tissue else None,
        load_library_metadata(Path(library["file"]), library.get("columns")) if library else None)
    for row in libraries:
        row["report_group"] = "All libraries"
        row["timepoint"] = row.get("meta_timepoint")
        row["timepoint_order"] = row.get("meta_timepoint_order")
    timepoint_orders = {row["timepoint"]: row["timepoint_order"] for row in libraries
                        if row["timepoint"] is not None and row["timepoint_order"] is not None}
    metadata = [dict(row) for row in declarations]
    categorical = [row for row in metadata if row["type"] == "category" and row["column"] != "meta_timepoint"]
    default_group = categorical[0]["column"] if categorical else "report_group"
    metadata.append({
        "column": "report_group",
        "label": "All libraries",
        "type": "category",
        "source_level": "library",
        "source": "scope manifest",
    })
    return scope_id, libraries, metrics, metadata, declarations, timepoint_orders, default_group


def prepare(scope_dir: Path, output_dir: Path, metadata_sources: Mapping | None = None,
            calibration_paths: list[Path] | None = None) -> dict:
    """Write one scope's report payloads and return the main payload.

    Reads the scope manifest and the headline and metric tables under
    `scope_dir`. `metadata_sources` optionally declares a `tissue` and a
    `library` source, each a `file` with an optional `columns` list.
    `calibration_paths` are the scope's threshold and donor calibration tables
    under `scope_dir`, which the payload names by their release paths.
    Writes report.json, one details/<sha256 of library_id>.json per library and
    one genes/<locus>.json per locus under `output_dir`. Stops when a table
    names another scope, an out-of-scope library or a conflicting identity.
    """
    (scope_id, libraries, metrics, metadata, declarations,
     timepoint_orders, default_group) = _prepare_libraries(scope_dir, metadata_sources)
    by_id = {(row["library_id"],): row for row in libraries}
    detail = {row["library_id"]: {} for row in libraries}

    def load(name: str, keys: list[str], columns: list[str], nullable_keys: tuple[str, ...] = ()) -> list[dict]:
        path = scope_dir / "metrics" / f"{name}.tsv"
        require_columns(path, columns)
        rows = read_table(path)
        unique_rows(rows, keys, path.name, nullable_keys)
        for row in rows:
            if "scope_id" in row and row["scope_id"] != scope_id:
                raise ValueError(f"{path.name} contains another scope.")
            if "library_id" in row and (row["library_id"],) not in by_id:
                raise ValueError(f"{path.name} contains an out-of-scope library.")
            if "library_id" in row:
                identity = by_id[(row["library_id"],)]
                if any(field in row and row[field] != identity[field] for field in ("donor_id", "tissue_id")):
                    raise ValueError(f"{path.name} contains a conflicting library identity.")
        return rows

    tables = {}
    for key, name, keys, columns in (
        ("shm", "shm", ["library_id", "isotype_group", "v_region"], [
            "library_id",
            "isotype_group",
            "v_region",
            "mean_clone_median_shm",
            "n_clones_eligible",
            "n_members_eligible",
            "n_members_measured",
            "fraction_clones_zero_median_shm",
        ]),
        ("selection", "selection_pressure", ["library_id", "isotype_group"], [
            "library_id",
            "isotype_group",
            "n_clone_consensus_sequences",
            "sigma_cdr",
            "sigma_fwr",
            "sigma_cdr_status",
            "sigma_fwr_status",
            "n_clone_consensus_sequences_used_cdr",
            "n_clone_consensus_sequences_used_fwr",
            "sigma_cdr_ci_low",
            "sigma_cdr_ci_high",
            "sigma_fwr_ci_low",
            "sigma_fwr_ci_high",
            *sorted(SELECTION_CONTEXT_COLUMNS),
        ]),
        ("isotype", "isotype_usage", ["library_id", "isotype_level", "isotype"], [
            "library_id",
            "isotype_level",
            "isotype",
            "consensus_count",
            "support_fraction",
            "callable_consensus_count",
        ]),
    ):
        tables[key] = pack(load(name, keys, columns))
    diversity = load("diversity", ["library_id", "locus"], [
        "library_id",
        "locus",
        "inverse_simpson",
        "rarefied_inverse_simpson",
        "grouped_consensus_count",
        "rarefaction_status",
        "diversity_status",
        "min_consensus_count_for_diversity",
    ])
    unique_rows(diversity, ["library_id", "locus"], "Diversity")
    tables["diversity"] = pack(diversity)
    # The report derives the Simpson difference, the shared-group fractions and
    # the support only at one endpoint from these columns and diversity.tsv.
    pairs = load("donor_overlap_pairs", ["locus", "library_id_1", "library_id_2"], [
        "library_id_1", "library_id_2", "locus", "pair_status", "min_consensus_count_per_library",
        "grouped_consensus_count_1", "grouped_consensus_count_2",
        "n_groups_shared", "shared_support_fraction_1", "shared_support_fraction_2",
        "jaccard", "morisita_horn", "bray_curtis_proportion",
        "n_shared_higher_in_1", "n_shared_higher_in_2", "n_shared_equal"])
    for pair in pairs:
        a, b = pair["library_id_1"], pair["library_id_2"]
        if (a,) not in by_id or (b,) not in by_id or a == b:
            raise ValueError("Overlap endpoints must be distinct libraries from one in-scope donor.")
        endpoint_a, endpoint_b = by_id[(a,)], by_id[(b,)]
        if endpoint_a["donor_id"] != endpoint_b["donor_id"]:
            raise ValueError("Overlap endpoints must be distinct libraries from one in-scope donor.")
        if (pair.get("donor_id") != endpoint_a["donor_id"]
                or pair.get("tissue_id_1") != endpoint_a["tissue_id"]
                or pair.get("tissue_id_2") != endpoint_b["tissue_id"]):
            raise ValueError("Overlap endpoint metadata differs from the scope manifest.")
        pair["pair_id"] = json.dumps([scope_id, pair["locus"], a, b], separators=(",", ":"))
        pair["pair_label"] = f"{pair['locus']}: {a} / {b}"
        pair["library_id_a"], pair["library_id_b"] = a, b
        endpoints = [by_id[(a,)], by_id[(b,)]]
        if all(row["timepoint"] is not None for row in endpoints):
            timepoints = [row["timepoint"] for row in endpoints]
            if all(label in timepoint_orders for label in timepoints):
                ordered = sorted(timepoints, key=timepoint_orders.get)
            else:
                ordered = sorted(timepoints)
            pair["timepoint_pair"] = json.dumps(ordered, separators=(",", ":"))
        else:
            pair["timepoint_pair"] = None
        # Choosing the horizontal endpoint depth is presentation, not a new metric.
        pair["smaller_endpoint_consensus_count"] = min(pair["grouped_consensus_count_1"],
                                                       pair["grouped_consensus_count_2"])
    for key, name, keys, columns in (
        ("curves", "rarefaction_curves", ["library_id", "locus", "support_depth"], [
            "library_id",
            "locus",
            "group_unit",
            "support_depth",
            "richness",
            "hill_shannon",
            "inverse_simpson",
            "is_target",
            "is_endpoint",
        ]),
        ("abundance", "group_abundance", ["library_id", "locus", "group_key"], ["library_id", "locus", "group_unit", "group_key", "consensus_count", "support_fraction"]),
        ("clone_shm", "clone_shm", ["library_id", "clone_id", "v_region", "isotype_group"], [
            "library_id",
            "clone_id",
            "v_region",
            "isotype_group",
            "consensus_count",
            "n_members_eligible",
            "n_members_zero_shm",
            "n_members_positive_shm",
            "median_shm_frequency",
            "mutation_state",
        ]),
    ):
        rows = load(name, keys, columns)
        if key == "clone_shm":
            rows = [row for row in rows if in_clone_shm_payload(row)]
        grouped = {library_id: [] for library_id in detail}
        for row in rows:
            grouped[row["library_id"]].append(row)
        for library_id, subset in grouped.items():
            detail[library_id][key] = pack(subset)
    # qc.tsv gives each library and locus its QC fractions; the trail gives each
    # stage's units and support. Both are released whole, for every library.
    tables["qc"] = pack(load("qc", ["library_id", "locus"], [
        "library_id", "locus", "igblast_consensus_count", "common_qc_consensus_count",
        "grouped_consensus_count", *QC_FRACTIONS]))
    trail_keys = ["library_id", "locus", "stage", "isotype_group", "v_region"]
    tables["qc_trail"] = pack(load("qc_trail", trail_keys, [
        *trail_keys, "stage_order", "analysis", "kind", "decision_status", "locus_source",
        "input_unit", "n_units_in", "output_unit", "n_units_out",
        "consensus_count_in", "consensus_count_out"], ("isotype_group", "v_region")))
    gene_usage = {locus: [] for locus in LOCI}
    keys = ["library_id", "locus", "gene_type", "gene_key"]
    for row in load("gene_usage", keys, [*keys, "group_unit", "gene_names", "gene_equivalence_group", "n_groups",
                                         "group_fraction", "consensus_count", "support_fraction"]):
        gene_usage[row["locus"]].append(row)
    gene_usage_files = {}
    for locus, rows in gene_usage.items():
        filename = f"genes/{locus}.json"
        gene_usage_files[locus] = filename
        write_json(output_dir / filename, pack(rows))
    detail_files = {}
    for library_id, payload in detail.items():
        filename = f"details/{hashlib.sha256(library_id.encode()).hexdigest()}.json"
        detail_files[library_id] = filename
        write_json(output_dir / filename, payload)
    calibration = {}
    for path in calibration_paths or []:
        require_columns(path, ["scope_id"])
        rows = read_table(path)
        if any(row["scope_id"] != scope_id for row in rows):
            raise ValueError(f"{path.name} contains another calibration scope.")
        relative = Path(path).resolve().relative_to(Path(scope_dir).resolve()).as_posix()
        calibration[release_path(relative)] = pack(rows)
    result = {
        "scope_id": scope_id,
        "title": f"{scope_id} repertoire report",
        "libraries": pack(libraries),
        "pairs": pack(pairs),
        "tables": tables,
        "metrics": metrics,
        "metadata": metadata,
        "has_metadata": bool(declarations),
        "detail_files": detail_files,
        "gene_usage_files": gene_usage_files,
        "calibration": calibration,
        "default_context": {
            "chapter": "signal",
            "x_metric": "igh_yield_per_million_read_pairs",
            "y_metric": "igh_rarefied_inverse_simpson",
            "group_field": default_group,
            "colour_field": default_group,
            "locus": "IGH",
            "facet": any(row["timepoint"] is not None for row in libraries),
            "trajectories": False,
        },
    }
    write_json(output_dir / "report.json", result)
    return result
