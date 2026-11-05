#!/usr/bin/env python3
"""Check a standalone scope folder against its shipped observations.

check_tables.py separately checks dictionary headers and scientific types.
These checks read the actual handover: identities, supplied metadata, browser
payloads, AIRR links and transferred grouping support. The official AIRR 2.0
environment is required; Plotly and a browser are not. Optional source results
also verify every decompressed annotation cell and analysis-member addition.
"""

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
import gzip
import hashlib
from html.parser import HTMLParser
from itertools import groupby
import json
import math
from operator import itemgetter
from pathlib import Path, PurePosixPath
import re
import statistics
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from workflow.lib.airr_exports import (EXPORT_DROPPED, EXPORT_NAMES, EXPORT_VALUES, GERMLINE_FIELDS,
                                       order_export_columns)
from workflow.lib.alignment_export import ALIGNMENT_FIELDS
from workflow.lib.headline import HEADLINE_COLUMNS
from workflow.lib.release_bundle import release_path
from workflow.lib.repertoire_metadata import SEMANTIC_COUNTS, validate_scope_repertoires
from workflow.lib.release_methods import FIXED_METHODS, read_provenance, read_methods_inputs, methods_paragraphs

IDENTIFIERS = ("scope_id", "library_id", "donor_id", "tissue_id")
LOCI = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
INTERNAL_FILES = {".complete.json", ".snakemake_timestamp"}
RELEASE_FILES = "manifest/release_files.tsv"
TABLES = {
    "shm": (("shm",), ("library_id", "isotype_group", "v_region")),
    "selection": (("selection_pressure",), ("library_id", "isotype_group")),
    "isotype": (("isotype_usage",), ("library_id", "isotype_level", "isotype")),
    "diversity": (("diversity",), ("library_id", "locus")),
    "qc": (("qc",), ("library_id", "locus")),
    "qc_trail": (("qc_trail",), ("library_id", "locus", "stage", "isotype_group", "v_region")),
    "curves": (("rarefaction_curves",), ("library_id", "locus", "support_depth")),
    "abundance": (("group_abundance",), ("library_id", "locus", "group_key")),
    "clone_shm": (("clone_shm",), ("library_id", "clone_id", "v_region", "isotype_group")),
    "gene_usage": (("gene_usage",), ("library_id", "locus", "gene_type", "gene_key")),
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def read_tsv(path):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        require(columns and len(columns) == len(set(columns)), f"Invalid header: {path}")
        rows = list(reader)
    require(all(None not in row and None not in row.values() for row in rows), f"Malformed row: {path}")
    return columns, rows


def rows_by_library(path, libraries):
    """Yield each library's rows from a table sorted by library, holding one library at a time.

    `libraries` are visited in sorted order, and a library with no rows yields an
    empty list. Stops when the table holds rows that order does not reach.
    """
    with path.open(newline="", encoding="utf-8") as handle:
        groups = groupby(csv.DictReader(handle, delimiter="\t"), key=itemgetter("library_id"))
        pending = next(groups, None)
        for library_id in sorted(libraries):
            rows = []
            if pending is not None and pending[0] == library_id:
                rows = list(pending[1])
                pending = next(groups, None)
            yield library_id, rows
        require(pending is None, f"{path.name} is unsorted or names another library: {pending and pending[0]}")


def keyed(rows, keys, name):
    result = {}
    for row in rows:
        key = tuple("" if row[column] is None else str(row[column]) for column in keys)
        require(key not in result, f"Duplicate {keys} in {name}: {key}")
        result[key] = row
    return result


def unpack(table):
    columns, rows = table["columns"], table["rows"]
    require(len(columns) == len(set(columns)) and all(len(row) == len(columns) for row in rows),
            "Malformed packed report table")
    return [dict(zip(columns, row)) for row in rows]


def local_file(root, relative):
    value = str(relative)
    parsed = urlsplit(value)
    require(not parsed.scheme and not parsed.netloc and not parsed.query, f"Nonlocal artifact link: {value}")
    path = PurePosixPath(unquote(parsed.path))
    require(path.parts and not path.is_absolute() and ".." not in path.parts and "\\" not in value,
            f"Artifact path escapes its folder: {value}")
    result = root.joinpath(*path.parts)
    require(result.resolve().is_relative_to(root.resolve()) and result.is_file(), f"Missing local artifact: {value}")
    return result


def same_cell(value, source, context):
    if value is None:
        require(source in ("", "NA"), f"Report lost an available value: {context}")
    elif isinstance(value, bool):
        require(source in (("T", "TRUE") if value else ("F", "FALSE")), f"Boolean differs: {context}")
    elif isinstance(value, (int, float)):
        require(math.isfinite(value) and source != "" and float(source) == value, f"Numeric value differs: {context}")
    else:
        require(value == source, f"String value differs: {context}")


def check_identity(row, manifest, name):
    allowed = {column: {record[column] for record in manifest.values()} for column in IDENTIFIERS}
    for column, values in allowed.items():
        for field in (column, column + "_1", column + "_2"):
            if field in row:
                require(isinstance(row[field], str) and row[field] in values, f"Out-of-scope {field} in {name}: {row[field]!r}")
    if "library_id" in row:
        identity = manifest[row["library_id"]]
        require(all(row[column] == identity[column] for column in IDENTIFIERS if column in row),
                f"Conflicting library identity in {name}")


def read_script(path, library_id=None, namespace="REPORT_DETAIL"):
    text = path.read_text(encoding="utf-8")
    if library_id is None:
        prefix = "window.REPORT_DATA = "
        require(text.startswith(prefix) and text.rstrip().endswith(";"), f"Unexpected report script: {path}")
        return json.loads(text[len(prefix):].rstrip().removesuffix(";"))
    match = re.fullmatch(r'window\.' + re.escape(namespace) + r'\[("(?:[^"\\]|\\.)*")\] = (.*);\s*', text, re.S)
    require(match is not None and json.loads(match.group(1)) == library_id, f"Detail script identity differs: {path}")
    return json.loads(match.group(2))


class AssetLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        require(tag != "base" and attrs.get("type") != "module", "Report requires local classic scripts")
        if "src" in attrs:
            self.links.append(attrs["src"])
        if tag == "link" and "href" in attrs:
            self.links.append(attrs["href"])


def check_report(bundle, manifest, headline, provenance, tables):
    root = bundle / "report"
    # The detail payloads hold the clone SHM rows the report reads, not the whole table.
    tables = {**tables, "clone_shm": [row for row in tables["clone_shm"] if row["v_region"] == "whole_v"
                                      and row["isotype_group"] in REPORT_CLONE_SHM_GROUPS]}
    parser = AssetLinks()
    parser.feed((root / "index.html").read_text())
    for link in parser.links:
        local_file(root, link)
    for stylesheet in root.rglob("*.css"):
        for link in re.findall(r"url\(\s*['\"]?([^)'\"]+)['\"]?\s*\)", stylesheet.read_text()):
            if not link.startswith("data:"):
                local_file(stylesheet.parent, link)
    data = read_script(root / "data.js")
    require(data["scope_id"] == next(iter(manifest.values()))["scope_id"], "Report has another scope")
    libraries = keyed(unpack(data["libraries"]), ("library_id",), "Report libraries")
    require({key[0] for key in libraries} == set(manifest), "Report libraries differ from the manifest")
    for row in libraries.values():
        check_identity(row, manifest, "Report libraries")
        for column, value in headline[row["library_id"]].items():
            same_cell(row[column], value, f"headline/{row['library_id']}/{column}")
    require(data["provenance"] == provenance, "Report provenance differs from the shipped record")
    for name, packed in data["calibration"].items():
        source_columns, source = read_tsv(local_file(bundle, name))
        observed = unpack(packed)
        require(len(source) == len(observed), f"Report calibration count differs: {name}")
        if source:
            require(packed["columns"] == source_columns, f"Report calibration columns differ: {name}")
        for before, after in zip(source, observed):
            check_identity(after, manifest, "Report calibration " + name)
            for column, value in after.items():
                same_cell(value, before[column], f"{name}/{column}")

    def compare(name, packed, library_id=None, locus=None):
        source, keys = tables[name], TABLES[name][1]
        if library_id is not None:
            source = [row for row in source if row["library_id"] == library_id]
        if locus is not None:
            source = [row for row in source if row["locus"] == locus]
        require(packed["columns"] == (list(source[0]) if source else []),
                f"Report {name} columns differ from its source table")
        expected = keyed(source, keys, name + " source")
        observed = keyed(unpack(packed), keys, name + " report")
        require(set(expected) == set(observed), f"Report {name} omits or adds source rows")
        for key, row in observed.items():
            check_identity(row, manifest, "Report " + name)
            for column, value in row.items():
                same_cell(value, expected[key][column], f"{name}/{key}/{column}")

    for name, packed in data["tables"].items():
        compare(name, packed)
    _, overlaps = read_tsv(bundle / "metrics/donor_overlap_pairs.tsv")
    expected_pairs = keyed(overlaps, ("locus", "library_id_1", "library_id_2"), "Overlap")
    pairs = unpack(data["pairs"])
    require(len(pairs) == len(expected_pairs), "Report overlap count differs")
    keyed(pairs, ("pair_id",), "Report overlap")
    for pair in pairs:
        check_identity(pair, manifest, "Report overlap")
        key = (pair["locus"], pair["library_id_a"], pair["library_id_b"])
        require(key in expected_pairs, "Report has an unshipped overlap pair")
        require(json.loads(pair["pair_id"]) == [data["scope_id"], *key], "Pair identity differs from its endpoints")
        for column, value in expected_pairs[key].items():
            same_cell(pair[column], value, "Overlap " + column)
    require(set(data["detail_files"]) == set(manifest), "Detail index differs from the scope libraries")
    paths = set()
    for library_id, relative in data["detail_files"].items():
        path = local_file(root, relative)
        require(path not in paths, "Libraries share a detail file")
        paths.add(path)
        detail = read_script(path, library_id)
        require(set(detail) == {"curves", "abundance", "clone_shm"}, "Unexpected detail table family")
        for name, packed in detail.items():
            compare(name, packed, library_id)
    compare("curves", read_script(local_file(root, data["curve_file"]), "cohort", "REPORT_CURVES"))
    require(set(data["gene_usage_files"]) == set(LOCI), "Gene-use index differs from the receptor loci")
    gene_paths = set()
    for locus, relative in data["gene_usage_files"].items():
        path = local_file(root, relative)
        require(path not in gene_paths, "Loci share a gene-use file")
        gene_paths.add(path)
        compare("gene_usage", read_script(path, locus, "REPORT_GENES"), locus=locus)
    return len(paths)


# Exported AIRR extensions use the T/F encoding of the standard AIRR booleans.
AIRR_FLAGS = ("grouping_eligible", "is_contig_representative", "is_member_representative",
              "c_call_present", "c_call_locus_match", "c_call_enough_evidence", "c_call_callable")


def boolean(value):
    require(value in ("TRUE", "FALSE", "T", "F"), f"Invalid grouping flag: {value}")
    return value in ("TRUE", "T")


def check_grouping(rows, library_id, trail):
    by_id = keyed(rows, ("sequence_id",), "Exported AIRR")
    transferred = {}
    sequence_keys = set()
    for row in rows:
        for flag in AIRR_FLAGS:
            require(row.get(flag, "") in ("T", "F", ""), f"AIRR flag {flag} is not T or F: {row.get(flag)!r}")
        eligible, representative = boolean(row["grouping_eligible"]), boolean(row["is_member_representative"])
        contig_representative = boolean(row["is_contig_representative"])
        contig_link = row["contig_representative_sequence_id"]
        link = row["member_representative_sequence_id"]
        require(row["contig_id"].startswith(library_id + "~"), "contig_id is not library-qualified")
        if eligible:
            require((contig_link,) in by_id, "Eligible AIRR row has an absent contig representative")
            contig = by_id[(contig_link,)]
            require(boolean(contig["grouping_eligible"]) and boolean(contig["is_contig_representative"]),
                    "Contig link does not name an eligible contig representative")
            require(contig_representative == (row["sequence_id"] == contig_link),
                    "Contig representative flag differs from its link")
            require((link,) in by_id, "Eligible AIRR row has an absent representative")
            target = by_id[(link,)]
            require(boolean(target["is_member_representative"]) and
                    boolean(target["is_contig_representative"]) and
                    target["locus"] == contig["locus"] == row["locus"], "Representative locus or flag differs")
            require(representative == (row["sequence_id"] == link), "Representative flag differs from its link")
            if row["locus"] in ("IGH", "IGK", "IGL"):
                require(contig["contig_id"] == row["contig_id"], "BCR contig representative crosses contigs")
            else:
                require(contig_link == row["sequence_id"], "TCR contig alternative was collapsed")
            require(contig["member_representative_sequence_id"] == link and
                    contig["sequence"] == target["sequence"],
                    "Member link does not preserve its contig representative's exact sequence")
            transferred[link] = transferred.get(link, 0) + float(row["consensus_count"])
        else:
            require(not representative and not contig_representative and link == contig_link == "",
                    "Ineligible AIRR row contributes a representative link")
        require(representative or float(row["member_consensus_count"]) == 0, "Nonrepresentative carries inference support")
        if representative:
            key = (row["locus"], row["sequence"])
            require(key not in sequence_keys, "Inference representatives repeat an exact sequence within a library/locus")
            sequence_keys.add(key)
    for link, support in transferred.items():
        require(float(by_id[(link,)]["member_consensus_count"]) == support, "Representative support differs from eligible linked members")
    stages = keyed([row for row in trail if row["library_id"] == library_id
                    and not row["isotype_group"] and not row["v_region"]],
                   ("locus", "stage"), "QC trail")
    for locus in (*LOCI, "all"):
        source = [row for row in rows if locus == "all" or row["locus"] == locus]
        eligible = [row for row in source if boolean(row["grouping_eligible"])]
        backbones = [row for row in eligible if boolean(row["is_contig_representative"])]
        representatives = [row for row in eligible if boolean(row["is_member_representative"])]
        incoming = sum(float(row["consensus_count"]) for row in source)
        retained = sum(float(row["consensus_count"]) for row in eligible)
        expected_stages = [("grouping_eligible", 15, len(source), len(eligible), incoming, retained),
                           ("exact_sequence_collapse", 17, len(backbones), len(representatives), retained, retained)]
        if locus in ("IGH", "IGK", "IGL"):
            expected_stages.append(("contig_collapse", 16, len(eligible), len(backbones), retained, retained))
        # Members join a clone at IGH and a clonotype elsewhere, and clone_id names either.
        if locus != "all":
            groups = {row["clone_id"] for row in representatives}
            expected_stages.append(("receptor_units", 18, len(representatives), len(groups), retained, retained))
        for stage, order, *expected in expected_stages:
            require((locus, stage) in stages, f"Missing grouping trail row: {library_id}/{locus}/{stage}")
            row = stages[(locus, stage)]
            require(int(row["stage_order"]) == order, "Grouping trail stage order differs")
            observed = [float(row[column]) for column in ("n_units_in", "n_units_out", "consensus_count_in", "consensus_count_out")]
            require(observed == expected, f"Grouping trail differs from exported AIRR: {library_id}/{locus}/{stage}")


# qc.tsv names each grouping reason by the grouping_status tokens it sums.
QC_REASONS = ("v_call_unresolved", "j_call_unresolved", "junction_ineligible",
              "c_call_ineligible", "contig_representative_ineligible")
QC_REASON_TOKENS = {"ineligible_junction": "junction_ineligible", "c_locus_conflict": "c_call_ineligible",
                    "c_call_incomplete_candidates": "c_call_ineligible",
                    "ineligible_contig_representative": "contig_representative_ineligible"}


def quotient(numerator, denominator):
    """The exact quotient at the 15 significant digits qc.tsv keeps, or None at a zero denominator.

    scope_qc_audit.py divides exact decimals. A float quotient can round the
    15th significant digit the other way, as 322/349 does.
    """
    numerator, denominator = Decimal(numerator), Decimal(denominator)
    return None if denominator == 0 else float(format(numerator / denominator, ".15g"))


def check_qc_table(qc_rows, trail, libraries):
    """qc.tsv must carry the trail's support counts and their exact quotients."""
    qc = keyed(qc_rows, ("library_id", "locus"), "QC")
    require(set(qc) == {(library, locus) for library in libraries for locus in LOCI},
            "QC table lacks the complete library/locus grid")
    stages = keyed([row for row in trail if row["locus"] in LOCI and row["locus_source"] == "igblast"
                    and not row["isotype_group"] and not row["v_region"]],
                   ("library_id", "locus", "stage"), "QC trail")
    for (library, locus), row in qc.items():
        productive = stages[(library, locus, "qc_productive")]
        common = stages[(library, locus, "common_qc")]["consensus_count_out"]
        eligible = stages[(library, locus, "grouping_eligible")]["consensus_count_out"]
        annotated = productive["consensus_count_in"]
        for column, value in (("igblast_consensus_count", annotated), ("common_qc_consensus_count", common),
                              ("grouped_consensus_count", eligible)):
            require(row[column] == value, f"QC count differs from the trail: {library}/{locus}/{column}")
        fractions = {"nonproductive_fraction": (Decimal(annotated) - Decimal(productive["consensus_count_out"]), annotated),
                     "common_qc_fraction_of_annotated": (common, annotated),
                     "grouping_eligible_fraction_of_common_qc": (eligible, common)}
        for level in ("class", "subclass"):
            fractions[f"{level}_callable_fraction"] = (
                (stages[(library, locus, f"{level}_callable")]["consensus_count_out"], common)
                if locus == "IGH" else (0, 0))
        for column, (numerator, denominator) in fractions.items():
            same_cell(quotient(numerator, denominator), row[column], f"QC fraction/{library}/{locus}/{column}")
    return qc


def check_qc_reasons(rows, library_id, qc):
    """Rebuild a library's common-QC support and grouping-reason fractions from its AIRR export."""
    sums = {locus: {name: Decimal(0) for name in ("common_qc", "eligible", *QC_REASONS)} for locus in LOCI}
    for row in rows:
        status = row["grouping_status"]
        require((status == "eligible") == boolean(row["grouping_eligible"]),
                f"grouping_status disagrees with grouping_eligible: {row['sequence_id']}")
        tokens = [] if status == "eligible" else status.split(";")
        reasons = set()
        for token in tokens:
            reason = ("v_call_unresolved" if token.startswith("v_call_") else
                      "j_call_unresolved" if token.startswith("j_call_") else QC_REASON_TOKENS.get(token))
            require(reason is not None, f"Unknown grouping_status reason: {token!r}")
            reasons.add(reason)
        support = Decimal(row["consensus_count"])
        locus = sums[row["locus"]]
        locus["common_qc"] += support
        locus["eligible"] += support if not reasons else 0
        for reason in reasons:
            locus[reason] += support
    for name, totals in sums.items():
        row = qc[(library_id, name)]
        require(Decimal(row["common_qc_consensus_count"]) == totals["common_qc"] and
                Decimal(row["grouped_consensus_count"]) == totals["eligible"],
                f"QC support differs from exported AIRR: {library_id}/{name}")
        for reason in QC_REASONS:
            same_cell(quotient(totals[reason], totals["common_qc"]), row[reason + "_fraction"],
                      f"QC reason fraction/{library_id}/{name}/{reason}")


def check_vj_pair_marginals(pairs, usage):
    """Summing vj_pair_usage over one segment must give gene_usage for the other.

    Counts and support compare exactly. The fractions compare within 1e-9, since
    each file keeps its own rounded quotients. Both tables must name the same
    entities at each library and locus, with one label and one group_unit.
    """
    marginals = {}
    for row in pairs:
        for gene_type, key, names in (("V", row["v_gene_key"], row["v_gene_names"]),
                                      ("J", row["j_gene_key"], row["j_gene_names"])):
            total = marginals.setdefault((row["library_id"], row["locus"], gene_type, key), dict(
                gene_names=names, group_unit=row["group_unit"], n_groups=Decimal(0),
                consensus_count=Decimal(0), group_fraction=0.0, support_fraction=0.0))
            require(total["gene_names"] == names and total["group_unit"] == row["group_unit"],
                    f"V and J pair usage names one entity two ways: {row['library_id']}/{row['locus']}/{key}")
            for column in ("n_groups", "consensus_count"):
                total[column] += Decimal(row[column])
            for column in ("group_fraction", "support_fraction"):
                total[column] += float(row[column])
    expected = keyed(usage, ("library_id", "locus", "gene_type", "gene_key"), "Gene usage")
    require(set(marginals) == set(expected), "V and J pair usage names other entities than gene usage")
    for key, row in expected.items():
        total = marginals[key]
        for column in ("gene_names", "group_unit"):
            require(total[column] == row[column], f"V and J pair marginal differs from gene usage: {key}/{column}")
        for column in ("n_groups", "consensus_count"):
            require(total[column] == Decimal(row[column]), f"V and J pair marginal differs from gene usage: {key}/{column}")
        for column in ("group_fraction", "support_fraction"):
            require(math.isclose(total[column], float(row[column]), abs_tol=1e-9),
                    f"V and J pair marginal differs from gene usage: {key}/{column}")


# The columns of donor_fits.tsv that only a fitted donor holds.
FIT_PARAMETERS = ("near_weight", "near_shape", "near_scale", "far_weight", "far_mean", "far_sd",
                  "absolute_log_likelihood", "fit_sensitivity")
COMPONENT_TOLERANCE = 1e-9


def gamma_density(x, shape, scale):
    """The gamma probability density at `x` above zero, through its logarithm."""
    return math.exp((shape - 1) * math.log(x) - x / scale - math.lgamma(shape) - shape * math.log(scale))


def normal_density(x, mean, sd):
    """The normal probability density at `x`."""
    return math.exp(-0.5 * ((x - mean) / sd) ** 2) / (sd * math.sqrt(2 * math.pi))


def check_threshold_components(fits, histograms):
    """Each fitted donor's mixture must reproduce its histogram's component densities.

    `fits` are donor_fits.tsv rows and `histograms` donor_distance_histograms.tsv
    rows. At each bin midpoint the near component is near_weight times a gamma
    density and the far component far_weight times a normal density. Each must
    match its released cell within COMPONENT_TOLERANCE, relative to the larger
    value. A subnormal double holds too few bits for that comparison, so two values
    also agree when they differ by less than the smallest normal double. A donor
    with no fit holds no parameters and no component densities.
    """
    bins = defaultdict(list)
    for row in histograms:
        bins[row["donor_id"]].append(row)
    require(set(bins) == {row["donor_id"] for row in fits}, "Donor fits and distance histograms name other donors")
    for fit in fits:
        donor = fit["donor_id"]
        if fit["fit_status"] != "fitted":
            require(not any(fit[column] for column in FIT_PARAMETERS),
                    f"A donor without a fit holds mixture parameters: {donor}")
            require(not any(row["near_component_density"] or row["far_component_density"] for row in bins[donor]),
                    f"A donor without a fit holds component densities: {donor}")
            continue
        require(all(fit[column] for column in FIT_PARAMETERS), f"A fitted donor lacks mixture parameters: {donor}")
        near_weight, near_shape, near_scale, far_weight, far_mean, far_sd = (
            float(fit[column]) for column in FIT_PARAMETERS[:6])
        for row in bins[donor]:
            midpoint = float(row["distance_bin_midpoint"])
            for column, rebuilt in (("near_component_density", near_weight * gamma_density(midpoint, near_shape, near_scale)),
                                    ("far_component_density", far_weight * normal_density(midpoint, far_mean, far_sd))):
                require(row[column] != "" and math.isclose(float(row[column]), rebuilt, rel_tol=COMPONENT_TOLERANCE,
                                                           abs_tol=sys.float_info.min),
                        f"Donor mixture does not reproduce its histogram: {donor}/{midpoint}/{column}: "
                        f"released {row[column]!r}, rebuilt {rebuilt!r}")


def check_clonotype_members(members, groups, rows, library_id):
    """Reconstruct non-IGH group support through released representative links."""
    representatives = {row["sequence_id"]: row for row in rows
                       if row["locus"] != "IGH" and boolean(row["is_member_representative"])}
    mapped = keyed(members, ("sequence_id",), "Clonotype membership")
    require({key[0] for key in mapped} == set(representatives),
            "Clonotype membership must cover each non-IGH inference representative once")
    by_group = {key[0]: row for key, row in keyed(groups, ("clonotype_id",), "Clonotypes").items()}
    support = {}
    for member in members:
        source = representatives[member["sequence_id"]]
        group = by_group.get(member["clonotype_id"])
        require(group is not None, "Clonotype member names an absent group")
        require(member["locus"] == source["locus"] == group["locus"], "Clonotype member locus differs from its source")
        expected_key = "~".join((source["junction"].upper(), source["v_gene_key"], source["j_gene_key"]))
        require(group["clonotype_key"] == expected_key, "Clonotype key differs from its member's source")
        require(member["member_consensus_count"] == source["member_consensus_count"],
                "Clonotype member support differs from transferred representative support")
        support[member["clonotype_id"]] = support.get(member["clonotype_id"], 0) + int(member["member_consensus_count"])
    require(set(support) == set(by_group), "Clonotype groups differ from membership")
    for clonotype_id, total in support.items():
        require(total == int(by_group[clonotype_id]["consensus_count"]), "Clonotype support differs from its members")
    for clonotype_id, group in by_group.items():
        representative = mapped.get((group["representative_sequence_id"],))
        require(representative is not None and representative["clonotype_id"] == clonotype_id,
                "Clonotype representative is not one of its members")
        require(representatives[group["representative_sequence_id"]]["junction_length"] == group["junction_length"],
                "Clonotype junction length differs from its representative's AIRR record")


def check_shared_representatives(shared, rows, library_id):
    """Each shared group's representative ID names the AIRR record of its junction.

    `shared` holds the shared_groups rows whose representative comes from
    `library_id`, and `rows` is that library's AIRR export.
    """
    records = {row["sequence_id"]: row for row in rows}
    for row in shared:
        record = records.get(row["representative_sequence_id"])
        require(record is not None and boolean(record["is_member_representative"]),
                f"Shared-clone representative is not a member in its library: {row['representative_sequence_id']}")
        require(record["locus"] == row["locus"] and record["junction"] == row["representative_junction"],
                "Shared-clone descriptors differ from their representative record")


# MakeDb's V-domain regions in sequence order. Together they spell the IMGT-gapped
# sequence_alignment, and the export leaves empty a region with no sequenced
# base, that is one made only of '.' and '-'.
V_DOMAIN = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]


def check_source(results, library_id, exported_columns, exported):
    root = results / "libraries" / library_id
    source_columns, source = read_tsv(root / "airr" / f"{library_id}_airr_qc_pass.tsv")
    member_columns, members = read_tsv(root / "metrics" / f"{library_id}_igh_grouping.tsv")
    nonigh_columns, nonigh = read_tsv(root / "metrics" / f"{library_id}_nonigh_grouping.tsv")
    require(member_columns == nonigh_columns, "Grouping branch headers differ")
    members += nonigh
    members = keyed(members, ("sequence_id",), "Source analysis members")
    extra = [column for column in member_columns if column not in ("library_id", "sequence_id", "locus")]
    _, makedb = read_tsv(root / "airr" / f"{library_id}_makedb_pass.tsv.gz")
    makedb = keyed(makedb, ("sequence_id",), "Source MakeDb records")
    kept = [EXPORT_NAMES.get(column, column) for column in source_columns if column not in EXPORT_DROPPED]
    call_columns, calls = read_tsv(root / "metrics" / f"{library_id}_isotype_calls.tsv")
    call_fields = ["c_call_present", "c_call_locus_match",
                   "c_call_enough_evidence", "c_call_callable", "isotype_class", "isotype_subclass"]
    regions = [column for column in V_DOMAIN if column != "cdr3"]
    calls = keyed(calls, ("sequence_id",), "Source qualified calls")
    require(set(calls) == {(row["sequence_id"],) for row in source if row["locus"] == "IGH"},
            "Qualification does not cover the source IGH population")
    declared = [*kept, "source_sequence_id", "repertoire_id", "sample_processing_id", "data_processing_id",
                *(EXPORT_NAMES.get(column, column) for column in extra), *call_fields,
                *GERMLINE_FIELDS, *regions, *ALIGNMENT_FIELDS, "changeo_rev_comp", "clone_id"]
    require(exported_columns == order_export_columns(declared),
            "Decompressed AIRR columns differ from source and declared additions")
    require(len(source) == len(exported) == len(members), "Decompressed AIRR omits source records")
    for before, after in zip(source, exported):
        identity = before["sequence_id"]
        require(after["source_sequence_id"] == identity, "AIRR source row order or identity differs")
        for column, value in before.items():
            if column in EXPORT_DROPPED:
                continue
            if column == "rev_comp":
                require(after["changeo_rev_comp"] == value and after["rev_comp"] == "F",
                        "Stored query orientation or historical flag differs")
                continue
            if column in ("v_germline_start", "v_germline_end"):
                require(after["changeo_" + column] == value, "Native V coordinates were not retained")
                continue
            if column == "d_cigar" and after["d_alignment_status"] != "available":
                require(after[column] == "", "Unsupported D CIGAR must be empty")
                continue
            expected = library_id + "~" + value if column == "sequence_id" else value
            require(after[EXPORT_NAMES.get(column, column)] == expected,
                    f"Decompressed source cell differs: {library_id}/{identity}/{column}")
        for column in extra:
            expected = members[(identity,)][column]
            exported = EXPORT_NAMES.get(column, column)
            if expected and column in ("assembly_backbone_id", "backbone_representative_sequence_id",
                                        "inference_representative_sequence_id"):
                expected = library_id + "~" + expected
            if exported in AIRR_FLAGS:
                expected = {"TRUE": "T", "FALSE": "F"}.get(expected, expected)
            expected = EXPORT_VALUES.get(exported, {}).get(expected, expected)
            require(after[exported] == expected, f"Exported analysis-member field differs: {library_id}/{identity}/{column}")
        for column in call_fields:
            expected = calls[(identity,)][column] if before["locus"] == "IGH" else ""
            if column in AIRR_FLAGS:
                expected = {"TRUE": "T", "FALSE": "F"}.get(expected, expected)
            require(after[column] == expected, f"Exported qualified-call field differs: {library_id}/{identity}/{column}")
        if before["locus"] == "IGH":
            call = calls[(identity,)]
            require(call["c_call"] == before["c_call"], "Qualification lost its original call")
            # The export leaves out qualified_c_call, which is c_call on a callable record.
            require(call["qualified_c_call"] == (before["c_call"] if call["c_call_callable"] == "TRUE" else ""),
                    f"Qualified call is not c_call on a callable record: {library_id}/{identity}")
        source_regions = makedb[(identity,)]
        require(after["changeo_d_cigar"] == source_regions["d_cigar"], "Native D CIGAR differs")
        for column in (*GERMLINE_FIELDS, "cdr3"):
            published = "changeo_" + column if column in ("d_germline_start", "d_germline_end") else column
            require(after[published] == source_regions[column], f"Exported MakeDb field differs: {library_id}/{identity}/{column}")
        require("".join(source_regions[column] for column in V_DOMAIN) == after["sequence_alignment"],
                f"MakeDb V-domain regions do not spell sequence_alignment: {library_id}/{identity}")
        # MakeDb's region lengths locate each region in sequence_alignment. A
        # nonempty exported region is that slice, and an empty one is a slice
        # with no sequenced base.
        start = 0
        for column in V_DOMAIN:
            piece = after["sequence_alignment"][start:start + len(source_regions[column])]
            start += len(source_regions[column])
            if column == "cdr3":
                continue
            require(after[column] == piece if after[column] else set(piece) <= {".", "-"},
                    f"Exported V-domain region differs from its alignment slice: {library_id}/{identity}/{column}")


def check_member_assignments(rows, clone_members, clonotype_members):
    """clone_id names each member's scope clone or clonotype, and nothing else."""
    assigned = {member["sequence_id"]: member["clone_id"] for member in clone_members}
    assigned.update({member["sequence_id"]: member["clonotype_id"] for member in clonotype_members})
    for row in rows:
        require(row["clone_id"] == assigned.get(row["sequence_id"], ""),
                f"AIRR clone_id differs from the scope members: {row['sequence_id']}")
    require(set(assigned) <= {row["sequence_id"] for row in rows}, "A scope member is absent from its AIRR export")


def check_released_member_shm(members, member_shm, records):
    """Released members are inference representatives with SHM only where eligible.

    `members` are one library's clone_members rows and `member_shm` its
    clone_member_shm rows. Each member and region appears at most once. A row's
    mutation count, rate and vote median must be empty where the member is not
    SHM-eligible in that region.
    """
    sources = keyed(records, ("sequence_id",), "Exported SHM sources")
    for row in (*members, *member_shm):
        identity = row["sequence_id"]
        require((identity,) in sources, "Clone member is absent from exported AIRR")
        require(boolean(sources[(identity,)]["is_member_representative"]),
                "Clone member is not an inference representative")
    keyed(member_shm, ("sequence_id", "v_region"), "Clone member SHM")
    for row in member_shm:
        if not boolean(row["shm_eligible"]):
            require(row["n_mutations"] == row["shm_frequency"] == row["median_trust4_consensus_base_votes"] == "",
                    "Released member SHM outside eligibility must be empty")


def check_member_coverage(members, sparse, coverage, records, qc_rows):
    """Reconcile complete coverage against sparse measurements and released evidence."""
    source = {row["sequence_id"]: row for row in records}
    member = {row["sequence_id"]: row for row in members}
    measured = keyed(sparse, ("sequence_id", "v_region"), "Measured member regions")
    audit = keyed(coverage, ("sequence_id", "v_region"), "Complete member regions")
    require(set(audit) == {(identity, region) for identity in member for region in SHM_REGIONS},
            "Coverage audit must contain all six regions of every IGH member")
    whole_measured, whole_eligible, any_eligible = set(), set(), set()
    for key, row in audit.items():
        identity, region = key
        record = source[identity]
        require(row["clone_id"] == member[identity]["clone_id"], "Coverage clone differs from membership")
        eligible = boolean(row["shm_eligible"])
        expected = "unmeasured" if key not in measured else ("eligible" if eligible else "measured_ineligible")
        require(row["measurement_status"] == expected, "Coverage state disagrees with sparse measurement")
        require(not eligible or key in measured, "Eligible region has no sparse measurement")
        if key in measured:
            require(eligible == boolean(measured[key]["shm_eligible"]), "Coverage eligibility differs from sparse SHM")
        require(row["germline_reconstruction_status"] == record["germline_reconstruction_status"],
                "Coverage reconstruction status differs from AIRR evidence")
        evidence = row["evidence_sequence_id"]
        if eligible:
            require(evidence in source, "Coverage evidence is absent from AIRR")
            observed = source[evidence]
            require(boolean(observed["is_contig_representative"]) and
                    observed["member_representative_sequence_id"] == identity and
                    observed["sequence"] == record["sequence"], "Coverage evidence does not support its member")
        else:
            require(evidence == "", "Ineligible region names qualifying evidence")
        if region == "whole_v":
            if key in measured: whole_measured.add(identity)
            if eligible: whole_eligible.add(identity)
        elif eligible:
            any_eligible.add(identity)
    count = len(member)
    expected = dict(n_grouped_igh_members=count, n_shm_measured_members_whole_v=len(whole_measured),
        n_shm_eligible_members_whole_v=len(whole_eligible),
        shm_eligible_fraction_of_grouped_members_whole_v=len(whole_eligible) / count if count else None,
        n_shm_eligible_members_any_subregion=len(any_eligible),
        shm_eligible_fraction_of_grouped_members_any_subregion=len(any_eligible) / count if count else None)
    for row in qc_rows:
        if row["locus"] == "IGH":
            check_rebuilt_cells(row, expected, "Complete SHM coverage")
        else:
            require(all(row[field] == "" for field in expected), "IGH coverage populated for another locus")


def check_clone_representatives(clones, members_by_library, shared):
    """Each clone representative is one of its members, consistent with shared detail."""
    members = {row["sequence_id"]: row for rows in members_by_library.values() for row in rows}
    identities = {}
    for clone in clones:
        identity = clone["representative_sequence_id"]
        require(identity in members, "Clone representative is not a released member")
        row = members[identity]
        require(all(clone[field] == row[field] for field in ("scope_id", "donor_id", "clone_id")) and
                clone["representative_library_id"] == row["library_id"], "Clone representative crosses its unit")
        identities[clone["clone_id"]] = identity
    for row in shared:
        if row["locus"] == "IGH":
            require(identities[row["group_key"]] == row["representative_sequence_id"],
                    "Shared-group and clone representatives disagree")


# The SHM strata of shm_summaries.R. Every member enters all_IGH, class groups
# read isotype_class and subclass groups read isotype_subclass.
SHM_REGIONS = ("whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
SHM_CLASS_GROUPS = ("IgM", "IgG", "IgA")
SHM_SUBCLASS_GROUPS = ("IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2")
SHM_GROUPS = ("all_IGH", *SHM_CLASS_GROUPS, *SHM_SUBCLASS_GROUPS)
REPORT_CLONE_SHM_GROUPS = ("all_IGH", *SHM_CLASS_GROUPS)
REBUILD_TOLERANCE = 1e-9


def type7_quantile(values, probability):
    """The quantile of nonempty `values` under R's default type 7."""
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def shm_groups(member):
    """The SHM strata one clone_members row enters."""
    groups = ["all_IGH"]
    if member["isotype_class"] in SHM_CLASS_GROUPS:
        groups.append(member["isotype_class"])
    if member["isotype_subclass"] in SHM_SUBCLASS_GROUPS:
        groups.append(member["isotype_subclass"])
    return groups


def check_rebuilt_cells(row, expected, context):
    """Each released cell must hold its rebuilt value within tolerance, or be empty where there is none."""
    for column, value in expected.items():
        observed = row[column]
        if value is None:
            matched = observed == ""
        elif isinstance(value, str):
            matched = observed == value
        else:
            matched = observed != "" and abs(float(observed) - value) <= REBUILD_TOLERANCE
        require(matched, f"Released cell differs from its rebuild from member rows: {context}/{column}: "
                         f"released {observed!r}, rebuilt {value!r}")


def check_shm_reconstruction(members, member_shm, shm_rows, clone_shm_rows, library_id):
    """Rebuild one library's shm.tsv and clone_shm.tsv rows from its released member rows.

    `members` and `member_shm` are the library's clone_members and clone_member_shm
    rows, the only inputs of the rebuild. `shm_rows` and `clone_shm_rows` are its
    rows of the two summaries. Each region and group of shm.tsv counts its member
    SHM rows as measured, and its eligible rows give every other column. Each clone
    gives one median, so the mean, median and quantiles weight clones equally.
    clone_shm.tsv must hold exactly the clones with an eligible member, and its
    n_members and consensus_count cover all of the subset's members.
    """
    by_id = {member["sequence_id"]: member for member in members}
    subsets = defaultdict(lambda: [0, 0])
    for member in members:
        for group in shm_groups(member):
            subsets[(member["clone_id"], group)][0] += 1
            subsets[(member["clone_id"], group)][1] += int(member["member_consensus_count"])
    cells = defaultdict(list)
    for row in member_shm:
        member = by_id.get(row["sequence_id"])
        require(member is not None and member["clone_id"] == row["clone_id"],
                f"Member SHM row names no clone member: {library_id}/{row['sequence_id']}")
        require(row["v_region"] in SHM_REGIONS, f"Unknown member SHM region: {row['v_region']}")
        for group in shm_groups(member):
            cells[(group, row["v_region"])].append((row, int(member["member_consensus_count"])))
    summaries = keyed(shm_rows, ("isotype_group", "v_region"), f"SHM of {library_id}")
    require(set(summaries) == {(group, region) for group in SHM_GROUPS for region in SHM_REGIONS},
            f"SHM lacks the complete group and region grid: {library_id}")
    clone_summaries = keyed(clone_shm_rows, ("clone_id", "isotype_group", "v_region"),
                            f"Clone SHM of {library_id}")
    rebuilt_clones = set()
    for (group, region), summary in summaries.items():
        measured = cells.get((group, region), [])
        by_clone = defaultdict(list)
        for row, support in measured:
            if boolean(row["shm_eligible"]):
                by_clone[row["clone_id"]].append((float(row["n_mutations"]), float(row["shm_frequency"]), support))
        votes = [float(row["median_trust4_consensus_base_votes"]) for row, _ in measured
                 if boolean(row["shm_eligible"]) and row["median_trust4_consensus_base_votes"]]
        medians, zero_medians = [], 0
        for clone_id, clone in by_clone.items():
            mutations = [item[0] for item in clone]
            median = statistics.median(item[1] for item in clone)
            medians.append(median)
            zero_medians += statistics.median(mutations) == 0
            zero, positive = sum(value == 0 for value in mutations), sum(value > 0 for value in mutations)
            key = (clone_id, group, region)
            rebuilt_clones.add(key)
            require(key in clone_summaries, f"Clone SHM lacks a clone with an eligible member: {library_id}/{key}")
            check_rebuilt_cells(clone_summaries[key], {
                "n_members": subsets[(clone_id, group)][0], "consensus_count": subsets[(clone_id, group)][1],
                "n_members_eligible": len(clone), "n_members_zero_shm": zero, "n_members_positive_shm": positive,
                "median_shm_frequency": median,
                "mutation_state": ("all_zero" if zero == len(clone) else
                                   "all_positive" if positive == len(clone) else "mixed"),
                "eligible_consensus_count": sum(item[2] for item in clone)}, (library_id, *key))
        q25 = type7_quantile(medians, 0.25) if medians else None
        q75 = type7_quantile(medians, 0.75) if medians else None
        check_rebuilt_cells(summary, {
            "n_members_measured": len(measured),
            "n_members_eligible": sum(len(clone) for clone in by_clone.values()),
            "n_members_zero_shm": sum(item[0] == 0 for clone in by_clone.values() for item in clone),
            "n_clones_eligible": len(medians), "n_clones_zero_median_shm": zero_medians,
            "fraction_clones_zero_median_shm": zero_medians / len(medians) if medians else None,
            "mean_clone_median_shm": statistics.fmean(medians) if medians else None,
            "median_clone_median_shm": statistics.median(medians) if medians else None,
            "q25_clone_median_shm": q25, "q75_clone_median_shm": q75,
            "iqr_clone_median_shm": q75 - q25 if medians else None,
            "median_trust4_consensus_base_votes": statistics.median(votes) if votes else None}, (library_id, group, region))
    require(set(clone_summaries) == rebuilt_clones, f"Clone SHM holds a row with no eligible member: {library_id}")


def check_clone_isotype(rows, members, library_id):
    """One library's clone_isotype rows must sum to its clone_members rows within each clone and level.

    Each member resolving at a level enters one row there, through isotype_class
    or isotype_subclass, and a clone with no such member has no row at that level.
    """
    rebuilt = {}
    callable_support = defaultdict(int)
    for member in members:
        support = int(member["member_consensus_count"])
        for level, column in (("class", "isotype_class"), ("subclass", "isotype_subclass")):
            if member[column]:
                n_members, consensus_count = rebuilt.get((member["clone_id"], level, member[column]), (0, 0))
                rebuilt[(member["clone_id"], level, member[column])] = (n_members + 1, consensus_count + support)
                callable_support[(member["clone_id"], level)] += support
    released = keyed(rows, ("clone_id", "isotype_level", "isotype"), f"Clone isotype of {library_id}")
    require(set(released) == set(rebuilt), f"Clone isotype rows differ from the resolved members: {library_id}")
    for key, row in released.items():
        n_members, consensus_count = rebuilt[key]
        total = callable_support[key[:2]]
        check_rebuilt_cells(row, {"n_members": n_members, "consensus_count": consensus_count,
                                  "callable_consensus_count": total, "support_fraction": consensus_count / total},
                            (library_id, *key))


def check_shm_evidence(members, records):
    """Evidence links must name an identical retained contig representative in the same unit.

    `members` are rows of the donor member shards, which keep the evidence IDs
    the release leaves out.
    """
    sources = keyed(records, ("sequence_id",), "Exported SHM sources")
    for member in members:
        identity = member["sequence_id"]
        require((identity,) in sources, "Clone member is absent from exported AIRR")
        canonical = sources[(identity,)]
        require(boolean(canonical["is_member_representative"]), "Clone member is not an inference representative")
        for region in ("whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3"):
            evidence_id = member[f"shm_evidence_sequence_id_{region}"]
            eligible = boolean(member[f"shm_eligible_{region}"])
            require(bool(evidence_id) == eligible, "SHM evidence link disagrees with eligibility")
            if not eligible:
                continue
            require((evidence_id,) in sources, "SHM evidence source is absent from exported AIRR")
            evidence = sources[(evidence_id,)]
            require(boolean(evidence["is_contig_representative"]) and
                    evidence["member_representative_sequence_id"] == identity and
                    evidence["sequence"] == canonical["sequence"] and
                    evidence["library_id"] == member["library_id"] and evidence["locus"] == "IGH",
                    "SHM evidence must come from an identical retained contig representative")


def check_semantic_counts(validation):
    """The AIRR receipt holds integer counts at every locus, which cover its linked records once."""
    by_locus = validation["semantic_checks"]["by_locus"]
    require(set(by_locus) == set(LOCI), "AIRR validation lacks semantic counts for a receptor locus")
    for locus, counts in by_locus.items():
        require(set(counts) == set(SEMANTIC_COUNTS) and all(type(value) is int for value in counts.values()),
                f"AIRR validation semantic counts are incomplete at {locus}")
    require(sum(counts["n_records"] for counts in by_locus.values()) == sum(validation["linked_records"].values()),
            "AIRR validation n_records does not sum to the linked records")


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def check_folder(bundle):
    """Compare the delivered files with their manifest and receipt, and count them.

    The count leaves out the publisher and Snakemake receipts, as the receipt's
    own count does, and includes the manifest, which lists every other file.
    """
    files = set()
    for path in bundle.rglob("*"):
        if path.is_file() and path.name not in INTERNAL_FILES:
            files.add(path.relative_to(bundle).as_posix())
    require(files, "Bundle contains no delivered files")
    columns, rows = read_tsv(bundle / RELEASE_FILES)
    require(columns == ["path", "n_bytes", "sha256"], "Release file manifest has other columns")
    listed = [row["path"] for row in rows]
    require(listed == sorted(set(listed)), "Release file manifest must be sorted by path without repeats")
    expected = files - {RELEASE_FILES}
    require(set(listed) == expected, f"Release file manifest lists absent files {sorted(set(listed) - expected)} "
                                     f"and omits delivered files {sorted(expected - set(listed))}")
    for row in rows:
        path = bundle / PurePosixPath(row["path"])
        require(row["n_bytes"] == str(path.stat().st_size), f"Released file size differs from its manifest: {row['path']}")
        require(row["sha256"] == file_sha256(path), f"Released file hash differs from its manifest: {row['path']}")
    receipt = json.loads((bundle / ".complete.json").read_text())
    require(receipt["release_files_sha256"] == file_sha256(bundle / RELEASE_FILES),
            "Release receipt hash differs from the file manifest")
    require(receipt["files"] == len(files), "Release receipt file count differs from the folder")
    return len(files)


def validate_artifacts(bundle, results=None):
    bundle = Path(bundle).resolve()
    _, manifest_rows = read_tsv(bundle / "manifest/scope_manifest.tsv")
    manifest = {key[0]: row for key, row in keyed(manifest_rows, ("library_id",), "Manifest").items()}
    require(manifest and len({row["scope_id"] for row in manifest.values()}) == 1, "Manifest must contain one nonempty scope")
    scope_id = next(iter(manifest.values()))["scope_id"]
    for path in bundle.rglob("*.tsv"):
        if path.name == "scope_assignment.tsv":
            require(scope_id == "cohort", "Project release contains the whole-cohort assignment register")
            continue
        _, rows = read_tsv(path)
        for row in rows:
            check_identity(row, manifest, str(path.relative_to(bundle)))
    columns, headline_rows = read_tsv(bundle / "metrics/headline.tsv")
    require(columns == list(HEADLINE_COLUMNS), "Headline header differs from the 41-column contract")
    headline = {key[0]: row for key, row in keyed(headline_rows, ("library_id",), "Headline").items()}
    require(set(headline) == set(manifest), "Headline differs from complete manifest")
    require(not any(column.startswith("meta_") for column in columns), "Headline contains report metadata")
    require(not (bundle / "metadata_descriptors.json").exists(), "Metadata descriptors belong only to the report")
    provenance = read_provenance(bundle / "provenance.json", scope_id, set(manifest))
    for relative in provenance["files"].values():
        local_file(bundle, relative)
    config, method_libraries, threshold = read_methods_inputs(
        bundle / "effective_config.yaml", bundle / "manifest/scope_manifest.tsv",
        bundle / release_path(f"clones/{scope_id}_clone_threshold.tsv"))
    require(config["fixed_methods"] == FIXED_METHODS, "effective_config.yaml fixed_methods differ from the code's")
    methods = methods_paragraphs(provenance, config, method_libraries, threshold)
    payload = json.loads((bundle / "report/data.js").read_text().split(" = ", 1)[1].rstrip(";\n"))
    require(payload["methods"] == methods, "Report methods differ from the explicit settings and calibration")
    readme = (bundle / "README.md").read_text()
    require(all(paragraph in readme for paragraph in methods), "README and report methods differ")
    tables = {name: [row for table in source_names for row in read_tsv(bundle / release_path(f"metrics/{table}.tsv"))[1]]
              for name, (source_names, _) in TABLES.items()}
    require(all(row["mutation_state"] != "no_evaluable_members" for row in tables["clone_shm"]),
            "Released clone SHM keeps a row with no evaluable member")
    for library, row in headline.items():
        shm = keyed([source for source in tables["shm"] if source["library_id"] == library
                     and source["isotype_group"] == "all_IGH" and source["v_region"] == "whole_v"],
                    ("library_id",), "Whole-V SHM")[(library,)]
        eligible, total = shm["n_members_eligible"], shm["n_members_measured"]
        expected = float(eligible) / float(total) if total and eligible and float(total) > 0 else None
        same_cell(expected, row["igh_shm_eligible_fraction_of_measured_members_whole_v"],
                  f"SHM measured-member fraction/{library}")
    check_vj_pair_marginals(read_tsv(bundle / release_path("metrics/vj_pair_usage.tsv"))[1], tables["gene_usage"])
    check_threshold_components(read_tsv(bundle / release_path(f"clones/{scope_id}_donor_clone_thresholds.tsv"))[1],
                               read_tsv(bundle / release_path(f"clones/{scope_id}_clone_distance_diagnostics.tsv"))[1])
    qc = check_qc_table(tables["qc"], tables["qc_trail"], manifest)
    details = check_report(bundle, manifest, headline, provenance, tables)
    repertoire = json.loads((bundle / "airr/Repertoire.json").read_text())
    require(repertoire["scope_id"] == scope_id and isinstance(repertoire.get("processing_identity"), str) and repertoire["processing_identity"],
            "Repertoire needs its scope and processing identity")
    repertoires = keyed(repertoire["Repertoire"], ("repertoire_name",), "Repertoires")
    require({key[0] for key in repertoires} == set(manifest), "Repertoires differ from the scope")
    paths = {}
    airr_records = 0
    members_by_library = {library: [] for library in manifest}
    for member in read_tsv(bundle / release_path("metrics/clone_members.tsv"))[1]:
        members_by_library[member["library_id"]].append(member)
    # The donor member shards keep the SHM evidence IDs the release leaves out.
    shards_by_library = {library: [] for library in manifest}
    if results is not None:
        for path in sorted((Path(results) / "scopes" / scope_id / "intermediate").glob("*/*_clone_members.tsv")):
            for member in read_tsv(path)[1]:
                shards_by_library[member["library_id"]].append(member)
    clonotype_members = read_tsv(bundle / release_path("metrics/clonotype_members.tsv"))[1]
    clonotypes = read_tsv(bundle / release_path("metrics/clonotypes.tsv"))[1]
    keyed(clonotypes, ("clonotype_id",), "Clonotypes across the release")
    shared_by_library = {library: [] for library in manifest}
    shared = read_tsv(bundle / release_path("metrics/shared_groups.tsv"))[1]
    check_clone_representatives(read_tsv(bundle / release_path("metrics/clones.tsv"))[1], members_by_library, shared)
    for row in shared:
        sources = [library for library, identity in manifest.items()
                   if identity["donor_id"] == row["donor_id"]
                   and row["representative_sequence_id"].startswith(library + "~")]
        require(sources, "Shared-clone representative comes from outside the donor's libraries")
        shared_by_library[max(sources, key=len)].append(row)
    shm_by_library, clone_shm_by_library = defaultdict(list), defaultdict(list)
    for row in tables["shm"]:
        shm_by_library[row["library_id"]].append(row)
    for row in tables["clone_shm"]:
        clone_shm_by_library[row["library_id"]].append(row)
    # clone_member_shm.tsv holds several rows per member, so only one library's
    # rows stay in memory.
    coverage_groups = rows_by_library(bundle / release_path("metrics/clone_member_shm_coverage.tsv"), manifest)
    for library_id, member_shm in rows_by_library(bundle / release_path("metrics/clone_member_shm.tsv"), manifest):
        coverage_library, coverage = next(coverage_groups)
        require(coverage_library == library_id, "Coverage library order differs")
        record = repertoires[(library_id,)]
        paths[library_id] = record["data_processing"][0]["data_processing_files"][0]
        require("supplied_metadata" not in record, "AIRR export contains report metadata")
        fields, rows = read_tsv(local_file(bundle, paths[library_id]))
        check_grouping(rows, library_id, tables["qc_trail"])
        check_qc_reasons(rows, library_id, qc)
        check_released_member_shm(members_by_library[library_id], member_shm, rows)
        check_member_coverage(members_by_library[library_id], member_shm, coverage, rows,
                              [row for row in tables["qc"] if row["library_id"] == library_id])
        check_shm_reconstruction(members_by_library[library_id], member_shm, shm_by_library[library_id],
                                 clone_shm_by_library[library_id], library_id)
        if results is not None:
            require({member["sequence_id"] for member in shards_by_library[library_id]} ==
                    {member["sequence_id"] for member in members_by_library[library_id]},
                    "Released clone members differ from their donor shards")
            check_shm_evidence(shards_by_library[library_id], rows)
        check_clonotype_members([row for row in clonotype_members if row["library_id"] == library_id],
                                [row for row in clonotypes if row["library_id"] == library_id],
                                rows, library_id)
        check_member_assignments(rows, members_by_library[library_id],
                                 [row for row in clonotype_members if row["library_id"] == library_id])
        check_shared_representatives(shared_by_library[library_id], rows, library_id)
        if results is not None:
            check_source(Path(results), library_id, fields, rows)
        airr_records += len(rows)
    for library_id, rows in rows_by_library(bundle / release_path("metrics/clone_isotype.tsv"), manifest):
        check_clone_isotype(rows, members_by_library[library_id], library_id)
    require(set(paths.values()) == {path.relative_to(bundle).as_posix() for path in (bundle / "airr").glob("*.tsv.gz")}, "Bundle contains unlinked AIRR annotation")
    observed = validate_scope_repertoires(repertoire, bundle, manifest_rows, paths, repertoire["processing_identity"],
                                          member_assignments=True)
    recorded = json.loads((bundle / "airr/validation.json").read_text())
    require(recorded["cigar_reference_spans_and_bases_valid"] is True, "Staged reference validation was not performed")
    # The reference snapshot is external to the standalone bundle. The publisher
    # checks its bases; this check verifies that receipt and repeats all local checks.
    observed["cigar_reference_spans_and_bases_valid"] = recorded["cigar_reference_spans_and_bases_valid"]
    require(observed == recorded, "Recorded AIRR validation differs from the actual files")
    check_semantic_counts(observed)
    return {"scope_id": scope_id, "libraries": len(manifest), "headline_columns": len(columns),
            "airr_records": airr_records, "detail_payloads": details, "files": check_folder(bundle),
            "source_checked": results is not None, "airr_structural_valid": observed["structural_valid"],
            "airr_strict_metadata_valid": observed["strict_valid"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--results", type=Path, help="Optional results root containing libraries/ source products")
    args = parser.parse_args()
    print(json.dumps(validate_artifacts(args.bundle, args.results), indent=2, sort_keys=True))
