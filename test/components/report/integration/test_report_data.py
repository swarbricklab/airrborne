#!/usr/bin/env python3
"""Check report preparation against small authoritative-table contracts.

The fixture measurements are supplied values. These tests check preservation,
identity, missingness and portable rendering, not the numerical estimators.

Run: python test/run_tests.py --component report --level integration
"""

from __future__ import annotations

import csv
import copy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from test.support.release_context import write_release_context
from test.components.overlap.integration.reference_calculations import pairwise
from workflow.lib.report_data import LOCI, metric_definitions, prepare, read_table
from workflow.lib.headline import HEADLINE_COLUMNS, build_headline, write_headline
from workflow.lib.release_bundle import release_provenance


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if value is None else "TRUE" if value is True else "FALSE" if value is False else value for key, value in row.items()})


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read fixture cells as text before deliberately editing an input table."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unpack(table: dict) -> list[dict]:
    return [dict(zip(table["columns"], row)) for row in table["rows"]]


def create_fixture(scope: Path, with_metadata: bool = True) -> dict:
    identities = [
        {"scope_id": "project", "library_id": "0460", "donor_id": "d1", "tissue_id": "t1"},
        {"scope_id": "project", "library_id": "NA", "donor_id": "d1", "tissue_id": "t1"},
        {"scope_id": "project", "library_id": "third", "donor_id": "d1", "tissue_id": "t3"},
        {"scope_id": "project", "library_id": "other", "donor_id": "d2", "tissue_id": "t4"},
    ]
    write_tsv(scope / "manifest" / "scope_manifest.tsv", identities)
    headline = []
    for index, identity in enumerate(identities):
        row = {**identity, **{metric["column"]: None for metric in metric_definitions()}}
        for locus_index, locus in enumerate(LOCI):
            row[f"{locus.lower()}_yield_per_million_read_pairs"] = 0 if index % 9 == 0 else 10 ** ((index * 17 + locus_index * 11) % 70 / 20 - 2)
        for locus in ("igh", "trb"):
            row.update({f"{locus}_rarefaction_status": "below_target" if index == 0 else "available", f"{locus}_diversity_status": "available", f"{locus}_rarefaction_target": 5, f"{locus}_grouped_consensus_count": 12, f"{locus}_richness": 3, f"{locus}_common_qc_fraction_of_annotated": 6 / 7, f"{locus}_grouping_eligible_fraction_of_common_qc": 1})
            if index == 0:
                row[f"{locus}_rarefied_inverse_simpson"] = None
        row.update(input_read_pairs=1000000, yield_status="available", igg_iga_status="available", igm_consensus_count=3, igg_consensus_count=5, iga_consensus_count=2, class_callable_consensus_count=10)
        for isotype in ("igh", "igm", "igg", "iga"):
            row[f"{isotype}_n_shm_eligible_clones_whole_v"] = 2
        if with_metadata:
            row.update(meta_group="A" if index < 3 else None, meta_score=0 if index < 2 else None, meta_timepoint="T1" if index < 2 else "T3", meta_timepoint_order=1 if index < 2 else 3)
            row["meta_subtype"] = ["Subtype A", "Subtype A", "Subtype B", None][index % 4]
        headline.append(row)
    shm = []
    summary = []
    selection = []
    isotype = []
    diversity = []
    curves = []
    abundance = []
    qc = []
    trail = []
    clone_shm = []
    gene_usage = []
    for index, identity in enumerate(identities):
        for group in ("all_IGH", "IgM", "IgG", "IgA"):
            for region in ("whole_v", "cdr2", "fwr3"):
                missing = index % 11 == 0
                shm.append({**identity, "isotype_group": group, "v_region": region, "n_members_measured": 8, "n_members_eligible": 0 if missing else 4, "n_clones_eligible": 0 if missing else 2, "n_clones_zero_median_shm": 0 if missing else index % 3, "mean_clone_median_shm": None if missing else (index % 13) / 100, "fraction_clones_zero_median_shm": None if missing else (index % 3) / 2})
            selection.append({**identity, "isotype_group": group, "n_clone_consensus_sequences": 3,
                              "n_clone_consensus_sequences_with_positions_cdr": 2, "n_clone_consensus_sequences_used_cdr": 0 if index == 0 else 2,
                              "sigma_cdr_status": "no_usable_mutation_pdfs" if index == 0 else "estimated",
                              "sigma_cdr": None if index == 0 else 0.4, "sigma_cdr_ci_low": None if index == 0 else 0.1,
                              "sigma_cdr_ci_high": None if index == 0 else 0.6,
                              "n_clone_consensus_sequences_with_positions_fwr": 2, "n_clone_consensus_sequences_used_fwr": 2,
                              "sigma_fwr_status": "estimated", "sigma_fwr": -0.4,
                              "sigma_fwr_ci_low": -0.6, "sigma_fwr_ci_high": -0.1})
        for resolution, classes in (("class", ("IgM", "IgD", "IgG", "IgA", "IgE")), ("subclass", ("IgM", "IgD", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "IgE"))):
            weights = {"IgM": (index * 7) % 21, "IgD": index % 3, "IgG": (index * 11) % 31, "IgA": 0 if index % 6 == 0 else (index * 3) % 17 + 1, "IgE": 0}
            callable_count = sum(weights.values())
            for name in classes:
                parts = 4 if name.startswith("IgG") else 2
                count = weights[name] if name in weights else weights[name[:3]] // parts + (int(name[-1]) <= weights[name[:3]] % parts)
                isotype.append({**identity, "isotype_level": resolution, "isotype": name, "consensus_count": count, "support_fraction": count / callable_count if callable_count else None, "callable_consensus_count": callable_count})
        for locus in LOCI:
            diversity.append({**identity, "locus": locus, "group_unit": "clone" if locus == "IGH" else "clonotype", "grouped_consensus_count": 12, "min_consensus_count_for_diversity": 0, "diversity_status": "available", "richness": 3, "shannon_entropy": 1, "hill_shannon": 2.7, "inverse_simpson": 2.5, "clonality": 0.4, "top_1_support_fraction": 0.6, "top_5_support_fraction": 1, "top_10_support_fraction": 1, "rarefaction_target": 5, "rarefaction_status": "below_target" if index == 0 else "available", "rarefied_richness": None if index == 0 else 2.9, "rarefied_hill_shannon": None if index == 0 else 2.2, "rarefied_inverse_simpson": None if index == 0 else 2, "inverse_simpson_asymptote": 5, "asymptote_status": "available"})
            for depth in (1, 5, 12):
                curves.append({**identity, "locus": locus, "group_unit": "clone" if locus == "IGH" else "clonotype", "support_depth": depth, "richness": 1 if depth == 1 else 3, "hill_shannon": 1 if depth == 1 else 2.7, "inverse_simpson": 1 if depth == 1 else 2.5, "is_target": depth == 5, "is_endpoint": depth == 12})
            for clone, support in (("shared", 7), (f"local_{index}", 3), ("other_shared", 2)):
                abundance.append({**identity, "locus": locus, "group_unit": "clone" if locus == "IGH" else "clonotype", "group_key": clone, "group_rank": 1 if support == 7 else 2 if support == 3 else 3, "consensus_count": support, "support_fraction": support / 12})
            qc.append({**identity, "locus": locus, "igblast_consensus_count": 14, "nonproductive_fraction": 1 / 7,
                       "common_qc_consensus_count": 12, "common_qc_fraction_of_annotated": 6 / 7,
                       "grouped_consensus_count": 12, "grouping_eligible_fraction_of_common_qc": 1,
                       **{f"{reason}_fraction": 0 for reason in (
                           "v_call_unresolved", "j_call_unresolved", "junction_ineligible", "c_call_ineligible",
                           "contig_representative_ineligible")},
                       "class_callable_fraction": 10 / 12 if locus == "IGH" else None,
                       "subclass_callable_fraction": 10 / 12 if locus == "IGH" else None})
            for order, stage, analysis, kind, region in ((14, "common_qc", "common_qc", "exclusion", ""),
                                                         (28, "shm_eligible_members", "shm", "eligibility", "whole_v")):
                trail.append({**identity, "locus": locus, "stage_order": order, "stage": stage, "analysis": analysis,
                              "kind": kind, "decision_status": "applied", "isotype_group": "all_IGH" if region else "",
                              "v_region": region, "locus_source": "igblast", "input_unit": "airr_record",
                              "n_units_in": 4, "output_unit": "airr_record", "n_units_out": 3,
                              "consensus_count_in": 14, "consensus_count_out": 12})
            for gene in ("V", "J"):
                for call, count in ((1, index % 4), (2, 3 - index % 4)):
                    if count:
                        gene_usage.append({**identity, "locus": locus, "group_unit": "clone" if locus == "IGH" else "clonotype", "gene_type": gene, "gene_key": f"{locus}{gene}{call}", "gene_names": f"{locus}{gene}{call}", "gene_equivalence_group": f"{locus}{gene}{call}", "identical_alleles": "", "n_groups": count, "group_fraction": count / 3, "consensus_count": 10 * count, "support_fraction": count / 3})
        # Released clone SHM keeps rows with an eligible member, so the "shared"
        # clone has abundance and no whole-V SHM.
        for clone, count, median, state in ((f"local_{index}", 3, 0, "mixed"), ("other_shared", 2, 0.2, "all_positive")):
            # The cdr2 and IgG1 rows stay out of the report's payload.
            for region, group in (("whole_v", "all_IGH"), ("whole_v", "IgM"), ("whole_v", "IgG"), ("whole_v", "IgA"),
                                  ("cdr2", "all_IGH"), ("whole_v", "IgG1")):
                clone_shm.append({**identity, "clone_id": clone, "v_region": region, "isotype_group": group, "consensus_count": count, "n_members": 3, "n_members_eligible": 0 if median is None else 3, "n_members_zero_shm": 3 if state == "all_zero" else 2 if median == 0 else 0, "n_members_positive_shm": 0 if state == "all_zero" else 1 if median == 0 else 0 if median is None else 3, "median_shm_frequency": median, "mutation_state": state, "eligible_consensus_count": count})
    pairs = []
    # Keep the headline and detail views on the same supplied measurements.
    # These constructed tables exercise display contracts, not estimator accuracy.
    for index, row in enumerate(headline):
        library_id = row["library_id"]
        for locus in ("igh", "trb"):
            detail_row = next(item for item in diversity if item["library_id"] == library_id and item["locus"] == locus.upper())
            for field in ("rarefied_inverse_simpson", "top_5_support_fraction", "richness", "rarefaction_target", "rarefaction_status", "diversity_status"):
                row[f"{locus}_{field}"] = detail_row[field]
            row[f"{locus}_grouped_consensus_count"] = detail_row["grouped_consensus_count"]
        classes = [item for item in isotype if item["library_id"] == library_id and item["isotype_level"] == "class"]
        for item in classes:
            row[f"class_{item['isotype'].lower()}_fraction"] = item["support_fraction"]
            if item["isotype"] in ("IgM", "IgG", "IgA"):
                row[f"{item['isotype'].lower()}_consensus_count"] = item["consensus_count"]
        row["class_callable_consensus_count"] = classes[0]["callable_consensus_count"]
        igg, iga = row["igg_consensus_count"], row["iga_consensus_count"]
        row["log2_igg_iga"] = math.log2(igg / iga) if igg and iga else None
        row["igg_iga_status"] = "available" if igg and iga else "no_class_callable_evidence" if not row["class_callable_consensus_count"] else "zero_iga" if not iga else "zero_igg"
        summary.append({**{key: row[key] for key in ("scope_id", "library_id", "tissue_id", "donor_id", "igm_consensus_count", "igg_consensus_count", "iga_consensus_count", "class_callable_consensus_count", "log2_igg_iga", "igg_iga_status")},
                        "common_qc_consensus_count": row["class_callable_consensus_count"] + 2})
        for group in ("all_IGH", "IgM", "IgG", "IgA"):
            item = next(item for item in shm if item["library_id"] == library_id and item["isotype_group"] == group and item["v_region"] == "whole_v")
            for field, name in (("mean_clone_median_shm", "mean_clone_median_shm"),
                                ("fraction_clones_zero_median_shm", "fraction_clones_zero_median_shm"),
                                ("n_clones_eligible", "n_shm_eligible_clones")):
                row[f"{('igh' if group == 'all_IGH' else group.lower())}_{name}_whole_v"] = item[field]
        eligible, measured = row["igh_n_shm_eligible_clones_whole_v"], row["igh_richness"]
        row["igh_shm_eligible_clone_fraction_whole_v"] = eligible / measured if measured else None
        item = next(item for item in shm if item["library_id"] == library_id and item["isotype_group"] == "all_IGH" and item["v_region"] == "whole_v")
        row["igh_shm_eligible_fraction_of_measured_members_whole_v"] = item["n_members_eligible"] / item["n_members_measured"]
    write_tsv(scope / "metrics" / "headline.tsv",
              [{key: row[key] for key in HEADLINE_COLUMNS} for row in headline])
    pair_indices = [(0, 1), (1, 2)]
    for a, b in pair_indices:
        for locus in LOCI:
            pairs.append({"scope_id": "project", "donor_id": identities[a]["donor_id"], "locus": locus, "group_unit": "clone" if locus == "IGH" else "clonotype", "library_id_1": identities[a]["library_id"], "tissue_id_1": identities[a]["tissue_id"], "library_id_2": identities[b]["library_id"], "tissue_id_2": identities[b]["tissue_id"], "n_groups_1": 3, "n_groups_2": 3, "grouped_consensus_count_1": 12, "grouped_consensus_count_2": 12, "n_groups_shared": 2, "jaccard": 0.5, "morisita_horn": None if (a + b) % 19 == 0 else ((a * 13 + b * 7) % 101) / 100, "shared_consensus_count_1": 9, "shared_consensus_count_2": 9, "shared_support_fraction_1": 0.75, "shared_support_fraction_2": 0.75})
            vectors = [{row["group_key"]: row["consensus_count"] for row in abundance
                        if row["library_id"] == identities[index]["library_id"] and row["locus"] == locus}
                       for index in (a, b)]
            measurements = pairwise(*vectors)
            measurements = {key: value for key, value in measurements.items() if key not in pairs[-1]}
            pairs[-1].update(measurements)
    for name, rows in {
        "shm": shm, "selection_pressure": selection, "isotype_usage": isotype, "isotype_summary": summary,
        "diversity": diversity,
        "yield": [{key: row[key] for key in ("scope_id", "library_id", "tissue_id", "donor_id", "input_read_pairs", "yield_status")} for row in headline],
        "rarefaction_curves": curves, "group_abundance": abundance, "qc": qc, "qc_trail": trail,
        "clone_shm": clone_shm, "donor_overlap_pairs": pairs,
        "gene_usage": gene_usage,
    }.items():
        write_tsv(scope / "metrics" / f"{name}.tsv", rows)
    if not with_metadata:
        return {}
    columns = ["group", "score", "timepoint", "timepoint_order", "subtype"]
    tissues = {row["tissue_id"]: {"tissue_id": row["tissue_id"],
               **{column: row["meta_" + column] for column in columns}} for row in headline}
    path = scope / "tissues.tsv"
    write_tsv(path, list(tissues.values()))
    return {"tissue": {"file": str(path), "columns": columns}}


class ReportContracts(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.scope = self.root / "scope"
        self.metadata = create_fixture(self.scope)

    def tearDown(self) -> None:
        self.temporary.cleanup()


    def test_declared_categories_keep_na_and_mark_empty_missing(self) -> None:
        path = self.root / "categories.tsv"
        path.write_text("library_id\tcategory\nNA\tNA\n0460\t\n")
        self.assertEqual(read_table(path, categories={"category"}), [
            {"library_id": "NA", "category": "NA"},
            {"library_id": "0460", "category": None},
        ])

    def test_preserve_values_identifiers_and_complete_populations(self) -> None:
        report = prepare(self.scope, self.root / "prepared", self.metadata)
        libraries = {row["library_id"]: row for row in unpack(report["libraries"])}
        self.assertEqual(set(libraries), {"0460", "NA", "third", "other"})
        self.assertIsNone(libraries["0460"]["igh_rarefied_inverse_simpson"])
        self.assertEqual(libraries["NA"]["igh_mean_clone_median_shm_whole_v"], 0.01)
        self.assertEqual(libraries["0460"]["meta_score"], 0)
        self.assertIsNone(libraries["third"]["meta_score"])
        self.assertEqual(libraries["third"]["timepoint_order"], 3)
        payload = json.loads((self.root / "prepared" / report["detail_files"]["0460"]).read_text())
        clone_shm = unpack(payload["clone_shm"])
        self.assertEqual(len(clone_shm), 8)
        self.assertEqual({row["v_region"] for row in clone_shm}, {"whole_v"})
        self.assertEqual({row["isotype_group"] for row in clone_shm}, {"all_IGH", "IgM", "IgG", "IgA"})
        # The released table keeps the rows the report does not read.
        self.assertEqual(len(read_tsv(self.scope / "metrics/clone_shm.tsv")), 48)
        self.assertEqual(len(unpack(payload["abundance"])), 21)
        self.assertNotIn("sequence", payload["abundance"]["columns"])
        self.assertEqual(set(payload), {"curves", "abundance", "clone_shm"})

    def test_qc_tables_keep_every_row_and_unavailable_fractions(self) -> None:
        report = prepare(self.scope, self.root / "prepared", self.metadata)
        qc = {(row["library_id"], row["locus"]): row for row in unpack(report["tables"]["qc"])}
        self.assertEqual(len(qc), 4 * len(LOCI))
        self.assertEqual(qc[("NA", "IGH")]["nonproductive_fraction"], 1 / 7)
        self.assertEqual(qc[("NA", "IGH")]["grouping_eligible_fraction_of_common_qc"], 1)
        self.assertIsNone(qc[("NA", "TRB")]["class_callable_fraction"])
        trail = unpack(report["tables"]["qc_trail"])
        self.assertEqual(len(trail), 2 * 4 * len(LOCI))
        stratified = next(row for row in trail if row["stage"] == "shm_eligible_members")
        self.assertEqual((stratified["stage_order"], stratified["isotype_group"], stratified["v_region"]),
                         (28, "all_IGH", "whole_v"))
        unstratified = next(row for row in trail if row["stage"] == "common_qc")
        self.assertIsNone(unstratified["isotype_group"])
        path = self.scope / "metrics" / "qc_trail.tsv"
        rows = read_tsv(path)
        write_tsv(path, rows + rows[:1])
        with self.assertRaisesRegex(ValueError, "qc_trail.tsv requires unique"):
            prepare(self.scope, self.root / "prepared", self.metadata)


    def test_comparison_labels_follow_supplied_timepoint_order(self) -> None:
        path = self.scope / "ordered_visits.csv"
        path.write_text("tissue_id,timepoint,timepoint_order\nt1,Z baseline,1\nt3,A follow-up,2\n")
        report = prepare(self.scope, self.root / "prepared", {"tissue": {"file": str(path)}})
        pairs = [row for row in unpack(report["pairs"]) if row["library_id_2"] == "third"]
        self.assertTrue(pairs)
        for pair in pairs:
            self.assertEqual(json.loads(pair["timepoint_pair"]), ["Z baseline", "A follow-up"])
            self.assertEqual(pair["library_id_a"], pair["library_id_1"])
            self.assertEqual(pair["library_id_b"], pair["library_id_2"])

    def test_selection_context_and_uncertainty_remain_numeric_without_a_count_cutoff(self) -> None:
        report = prepare(self.scope, self.root / "prepared", self.metadata)
        selection = unpack(report["tables"]["selection"])
        self.assertTrue(selection)
        for row in selection:
            self.assertEqual(row["n_clone_consensus_sequences"], 3)
            self.assertEqual(row["n_clone_consensus_sequences_used_fwr"], 2)
            self.assertEqual(row["sigma_fwr"], -0.4)
            self.assertEqual(row["sigma_fwr_ci_low"], -0.6)
            self.assertEqual(row["sigma_fwr_status"], "estimated")
            self.assertEqual(row["n_clone_consensus_sequences_with_positions_fwr"], 2)

    def test_comparison_label_order_does_not_fill_missing_endpoint_order(self) -> None:
        path = self.scope / "ordered_visits.csv"
        path.write_text("tissue_id,timepoint,timepoint_order\nt1,Z baseline,\nt3,A follow-up,2\nt4,Z baseline,1\n")
        report = prepare(self.scope, self.root / "prepared", {"tissue": {"file": str(path)}})
        libraries = {row["library_id"]: row for row in unpack(report["libraries"])}
        self.assertIsNone(libraries["0460"]["timepoint_order"])
        pairs = [row for row in unpack(report["pairs"]) if row["library_id_2"] == "third"]
        self.assertTrue(pairs)
        for pair in pairs:
            self.assertEqual(json.loads(pair["timepoint_pair"]), ["Z baseline", "A follow-up"])
            self.assertIsNone(libraries[pair["library_id_1"]]["timepoint_order"])
            self.assertEqual(pair["library_id_a"], pair["library_id_1"])
            self.assertEqual(pair["library_id_b"], pair["library_id_2"])

    def test_pair_minimum_and_flag_do_not_replace_measurements_or_endpoint_depth(self) -> None:
        path = self.scope / "metrics" / "donor_overlap_pairs.tsv"
        rows = read_tsv(path)
        for row in rows:
            row.update(min_consensus_count_per_library=100, pair_status="below_minimum_support")
        write_tsv(path, rows)
        report = prepare(self.scope, self.root / "prepared", self.metadata)
        for source, pair in zip(rows, unpack(report["pairs"])):
            self.assertEqual(pair["min_consensus_count_per_library"], 100)
            self.assertEqual(pair["smaller_endpoint_consensus_count"], 12)
            self.assertEqual(pair["pair_status"], "below_minimum_support")
            for column in ("bray_curtis_proportion", "n_shared_higher_in_2", "n_shared_higher_in_1", "n_shared_equal",
                           "shared_support_fraction_1", "shared_support_fraction_2"):
                self.assertEqual(pair[column], float(source[column]))
            # The report derives the dropped measurements in the browser.
            self.assertNotIn("simpson_concentration_difference", pair)
            self.assertEqual(pair["library_id_a"], source["library_id_1"])
            self.assertEqual(pair["library_id_b"], source["library_id_2"])


    def test_missing_and_ambiguous_headline_rows_fail(self) -> None:
        path = self.scope / "metrics" / "headline.tsv"
        rows = read_tsv(path)
        write_tsv(path, rows[1:])
        with self.assertRaisesRegex(ValueError, "exactly"):
            prepare(self.scope, self.root / "prepared", self.metadata)
        write_tsv(path, rows + [rows[0]])
        with self.assertRaisesRegex(ValueError, "unique"):
            prepare(self.scope, self.root / "prepared", self.metadata)

    def test_out_of_scope_overlap_endpoint_fails(self) -> None:
        path = self.scope / "metrics" / "donor_overlap_pairs.tsv"
        rows = read_tsv(path)
        rows[0]["library_id_2"] = "unpublished"
        write_tsv(path, rows)
        with self.assertRaisesRegex(ValueError, "endpoints"):
            prepare(self.scope, self.root / "prepared", self.metadata)

    def test_report_without_metadata_preserves_the_scientific_headline(self) -> None:
        headline = self.scope / "metrics/headline.tsv"
        before = headline.read_bytes()
        report = prepare(self.scope, self.root / "prepared", self.metadata)
        self.assertTrue(report["has_metadata"])
        self.assertEqual(headline.read_bytes(), before)
        report = prepare(self.scope, self.root / "prepared")
        self.assertFalse(report["has_metadata"])
        self.assertEqual(report["default_context"]["group_field"], "report_group")
        self.assertFalse(any(column.startswith("meta_") for column in report["libraries"]["columns"]))
        self.assertEqual(headline.read_bytes(), before)

    def test_scope_sources_select_fields_and_keep_library_annotations_separate(self) -> None:
        sources = copy.deepcopy(self.metadata)
        path = self.scope / "libraries.tsv"
        write_tsv(path, [{"library_id": "0460", "batch": "01", "private": "Do not expose"},
                         {"library_id": "NA", "batch": "02", "private": "Do not expose"},
                         {"library_id": "outside", "batch": "03", "private": "Do not expose"}])
        sources["library"] = {"file": str(path), "columns": ["batch"]}
        sources["tissue"]["columns"] = ["score"]
        report = prepare(self.scope, self.root / "prepared", sources)
        rows = unpack(report["libraries"])
        self.assertEqual([row["meta_batch"] for row in rows], ["01", "02", None, None])
        self.assertEqual(rows[0]["meta_score"], 0)
        self.assertNotIn("meta_group", rows[0])
        self.assertNotIn("Do not expose", json.dumps(report))
        second = prepare(self.scope, self.root / "second", self.metadata)
        self.assertNotIn("meta_batch", unpack(second["libraries"])[0])
        self.assertIn("meta_group", unpack(second["libraries"])[0])

    def test_metadata_value_edits_change_only_the_prepared_report(self) -> None:
        measurements = {path: path.read_bytes() for path in (self.scope / "metrics").glob("*.tsv")}
        first = prepare(self.scope, self.root / "prepared", self.metadata)
        path = Path(self.metadata["tissue"]["file"])
        rows = read_tsv(path)
        rows[0]["score"] = "-0.25"
        write_tsv(path, rows)
        second = prepare(self.scope, self.root / "prepared", self.metadata)
        self.assertEqual(unpack(first["libraries"])[0]["meta_score"], 0)
        self.assertEqual(unpack(second["libraries"])[0]["meta_score"], -0.25)
        self.assertEqual(measurements, {path: path.read_bytes() for path in measurements})

    def test_file_only_metadata_preserves_original_headers_without_definitions(self) -> None:
        path = self.scope / "project.csv"
        path.write_text('tissue_id,Score (%),Treatment arm\nt1,-0.25,A\nt2,0,B\n')
        report = prepare(self.scope, self.root / "prepared", {"tissue": {"file": str(path)}})
        fields = [field for field in report["metadata"] if field["column"] != "report_group"]
        self.assertEqual(fields, [
            {"column": "meta_Score (%)", "label": "Score (%)", "type": "number", "source_level": "tissue"},
            {"column": "meta_Treatment arm", "label": "Treatment arm", "type": "category", "source_level": "tissue"}])
        self.assertEqual(unpack(report["libraries"])[0]["meta_Score (%)"], -0.25)

    def test_conflicting_timepoints_fail_before_report_files_are_written(self) -> None:
        path = self.scope / "visits.csv"
        sources = {"tissue": {"file": str(path)}}
        scientific_files = {file: file.read_bytes() for file in self.scope.rglob("*.tsv")}
        for second, message in (("T1,2", "conflicting timepoint_order"), ("T2,1", "different timepoint labels")):
            with self.subTest(second=second):
                path.write_text("tissue_id,timepoint,timepoint_order\nt1,T1,1\nt3," + second + "\n")
                with self.assertRaisesRegex(ValueError, message):
                    prepare(self.scope, self.root / "prepared", sources)
                self.assertFalse((self.root / "prepared").exists())
                self.assertEqual(scientific_files, {file: file.read_bytes() for file in scientific_files})


    def test_category_na_is_literal_and_empty_category_is_missing(self) -> None:
        path = Path(self.metadata["tissue"]["file"])
        rows = read_tsv(path)
        rows[0]["group"] = "NA"
        rows[1]["group"] = ""
        write_tsv(path, rows)
        report = prepare(self.scope, self.root / "prepared", self.metadata)
        libraries = unpack(report["libraries"])
        self.assertEqual(libraries[0]["meta_group"], "NA")
        self.assertEqual(libraries[1]["meta_group"], "NA")
        self.assertIsNone(libraries[2]["meta_group"])

    def test_authoritative_headline_producer_feeds_report_without_aliases(self) -> None:
        scope = self.root / "from-headline-producer"
        create_fixture(scope, with_metadata=False)
        metrics = scope / "metrics"
        source = read_tsv(metrics / "headline.tsv")
        write_tsv(metrics / "yield.tsv", [
            {key: ("available" if key == "yield_status" else row[key]) for key in ("library_id", "tissue_id", "donor_id", "input_read_pairs", "yield_status", *[f"{locus.lower()}_yield_per_million_read_pairs" for locus in LOCI])}
            for row in source])
        rows = build_headline(
            scope / "manifest" / "scope_manifest.tsv", metrics / "yield.tsv", metrics / "qc.tsv",
            metrics / "diversity.tsv", metrics / "isotype_usage.tsv", metrics / "isotype_summary.tsv", metrics / "shm.tsv")
        write_headline(rows, metrics / "headline.tsv")
        report = prepare(scope, self.root / "prepared")
        libraries = unpack(report["libraries"])
        self.assertEqual([row["library_id"] for row in libraries], ["0460", "NA", "third", "other"])
        self.assertIsNone(libraries[0]["igh_mean_clone_median_shm_whole_v"])
        self.assertEqual(libraries[1]["igh_mean_clone_median_shm_whole_v"], 0.01)
        self.assertEqual(libraries[1]["igh_rarefied_inverse_simpson"], 2)
        self.assertEqual(libraries[1]["trb_common_qc_fraction_of_annotated"], 6 / 7)
        self.assertEqual(libraries[1]["igh_grouping_eligible_fraction_of_common_qc"], 1)
        self.assertEqual(libraries[1]["yield_status"], "available")
        # Class support for the IgG:IgA hover comes from isotype_usage.tsv.
        self.assertEqual(libraries[0]["igg_consensus_count"], 0)
        self.assertNotIn("class_callable_fraction", libraries[0])
        self.assertFalse(report["has_metadata"])

    def test_explicit_calibration_is_preserved(self) -> None:
        threshold = self.scope / "clones" / "project_clone_threshold.tsv"
        write_tsv(threshold, [{"scope_id": "project", "clone_threshold": 0.1, "threshold_source": "fallback", "n_donors": 1, "n_donors_fitted": 0}])
        report = prepare(self.scope, self.root / "prepared", self.metadata, [threshold])
        self.assertNotIn("provenance", report)
        # The payload names the table where the release places it.
        self.assertEqual(list(report["calibration"]), ["igh_clone_threshold/scope_threshold.tsv"])
        self.assertEqual(unpack(report["calibration"]["igh_clone_threshold/scope_threshold.tsv"])[0]["clone_threshold"], 0.1)
        write_tsv(threshold, [{"scope_id": "another_project", "clone_threshold": 0.1}])
        with self.assertRaisesRegex(ValueError, "calibration scope"):
            prepare(self.scope, self.root / "prepared", self.metadata, [threshold])

    def test_donor_fit_and_histogram_measurements_stay_numeric_and_missing(self) -> None:
        fits = self.scope / "clones" / "project_donor_clone_thresholds.tsv"
        histograms = self.scope / "clones" / "project_clone_distance_diagnostics.tsv"
        parameters = ("near_weight", "near_shape", "near_scale", "far_weight", "far_mean", "far_sd",
                      "absolute_log_likelihood", "fit_sensitivity")
        write_tsv(fits, [{"scope_id": "project", "donor_id": "d1", "n_distances": 600, "n_distances_fitted": 500,
                          "fit_status": "fitted", "fitted_threshold": 0.2, **dict.fromkeys(parameters, 0.25)},
                         {"scope_id": "project", "donor_id": "d2", "n_distances": 4, "n_distances_fitted": None,
                          "fit_status": "too_few_distances", "fitted_threshold": None,
                          **dict.fromkeys(parameters, None)}])
        write_tsv(histograms, [{"scope_id": "project", "donor_id": "d1", "distance_bin_midpoint": 0.005,
                                "observed_density": 1.5, "near_component_density": 0.5, "far_component_density": None}])
        report = prepare(self.scope, self.root / "prepared", self.metadata, [fits, histograms])
        fitted, unfitted = unpack(report["calibration"]["igh_clone_threshold/donor_fits.tsv"])
        self.assertEqual(fitted["n_distances_fitted"], 500)
        self.assertEqual({fitted[column] for column in parameters}, {0.25})
        self.assertEqual(fitted["fit_status"], "fitted")
        self.assertIsNone(unfitted["n_distances_fitted"])
        self.assertEqual({unfitted[column] for column in parameters}, {None})
        histogram, = unpack(report["calibration"]["igh_clone_threshold/donor_distance_histograms.tsv"])
        self.assertEqual(histogram["observed_density"], 1.5)
        self.assertIsNone(histogram["far_component_density"])

    def test_gene_usage_keeps_each_locus_unit(self) -> None:
        report = prepare(self.scope, self.root / "prepared", self.metadata)
        for locus, unit in (("IGH", "clone"), ("TRB", "clonotype")):
            detail = json.loads((self.root / "prepared" / report["gene_usage_files"][locus]).read_text())
            self.assertEqual({row["group_unit"] for row in unpack(detail)}, {unit})

    def test_gene_usage_rejects_a_repeated_gene(self) -> None:
        path = self.scope / "metrics" / "gene_usage.tsv"
        rows = read_tsv(path)
        rows.append(dict(rows[0]))
        write_tsv(path, rows)
        with self.assertRaises(ValueError):
            prepare(self.scope, self.root / "prepared", self.metadata)


    def test_renderer_methods_use_explicit_inputs_and_preserve_measurements(self) -> None:
        from render_scope_report import render
        prepared, rendered = self.root / "prepared", self.root / "report"
        current, provenance, config, threshold = write_release_context(self.scope)
        prepare(self.scope, prepared, self.metadata)
        before = {str(path): path.read_bytes() for path in prepared.rglob("*") if path.is_file()}
        arguments = (prepared, rendered, provenance, config, self.scope / "manifest/scope_manifest.tsv", threshold)
        render(*arguments)
        data = json.loads((rendered / "data.js").read_text().split(" = ", 1)[1].rstrip(";\n"))
        # The report carries the record as the release does, with release paths.
        self.assertEqual(data["provenance"], release_provenance(current))
        self.assertEqual(data["provenance"]["files"]["threshold"], "igh_clone_threshold/scope_threshold.tsv")
        self.assertIn("threshold 0.1 (fallback)", " ".join(data["methods"]))
        self.assertIn("unavailable for 4", " ".join(data["methods"]))
        self.assertEqual(before, {str(path): path.read_bytes() for path in prepared.rglob("*") if path.is_file()})
        self.assertEqual(data["libraries"], json.loads((prepared / "report.json").read_text())["libraries"])
        provenance.write_text(json.dumps({**current, "scope_id": "other"}))
        with self.assertRaisesRegex(ValueError, "another scope"):
            render(*arguments)
        provenance.write_text(json.dumps(current))
        threshold.unlink()
        with self.assertRaises(FileNotFoundError):
            render(*arguments)


if __name__ == "__main__":
    unittest.main()
