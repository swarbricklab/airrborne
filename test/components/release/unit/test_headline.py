"""The headline preserves detailed values and complete scope observations.

Run: python test/run_tests.py --component release --level unit
"""

import copy
import csv
from pathlib import Path
import tempfile
import unittest

from workflow.lib.headline import HEADLINE_COLUMNS, build_headline, write_headline


ROOT = Path(__file__).resolve().parents[4]


def write_table(path, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class HeadlineTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.paths = {name: self.directory / (name + ".tsv") for name in
                      ("manifest", "yields", "qc", "diversity", "isotype", "isotype_summary", "shm")}
        self.manifest = [
            {"scope_id": "NA", "library_id": "0460", "donor_id": "001", "tissue_id": "T1"},
            {"scope_id": "NA", "library_id": "NA", "donor_id": "001", "tissue_id": "T1"},
            {"scope_id": "NA", "library_id": "empty", "donor_id": "002", "tissue_id": "T2"},
        ]
        self.tables = {name: [] for name in self.paths}
        self.tables["manifest"] = self.manifest
        for index, identity in enumerate(self.manifest):
            empty = index == 2
            # Collected library tables have no scope_id and append tissue/donor.
            collected = {key: identity[key] for key in ("library_id", "tissue_id", "donor_id")}
            yields = {locus + "_yield_per_million_read_pairs": "0" if empty else str((index + 1) * value)
                      for value, locus in enumerate(("igh", "igk", "igl", "tra", "trb", "trd", "trg"), start=1)}
            self.tables["yields"].append({**yields, "input_read_pairs": "1000000",
                                          "yield_status": "available", **collected})
            for locus in ("IGH", "TRB"):
                self.tables["qc"].append({**identity, "locus": locus,
                                          "common_qc_fraction_of_annotated": "" if empty else "0.857142857142857",
                                          "grouping_eligible_fraction_of_common_qc": "" if empty else "0.75"})
                self.tables["diversity"].append({
                    "locus": locus, "rarefied_inverse_simpson": "NA" if index != 1 else "2.125000000000003",
                    "inverse_simpson": "NA" if empty else "3.400000000000001",
                    "top_5_support_fraction": "NA" if empty else "0.950000000000001",
                    "top_10_support_fraction": "NA" if empty else "1",
                    "grouped_consensus_count": "0" if empty else "20",
                    "richness": "0" if empty else "7",
                    "rarefaction_target": "1000" if locus == "IGH" else "400",
                    "rarefaction_status": "empty_repertoire" if empty else "available" if index == 1 else "below_target",
                    "diversity_status": "empty_repertoire" if empty else "available",
                    "hill_shannon": "NA" if empty else "4.5", "top_1_support_fraction": "NA" if empty else "0.4",
                    **identity,
                })
            for group, fraction in (("IgM", "0.2"), ("IgD", "0"), ("IgG", "0.6"), ("IgA", "0.2"), ("IgE", "0")):
                self.tables["isotype"].append({"isotype_level": "class", "isotype": group,
                                               "support_fraction": "NA" if empty else fraction, **collected})
            self.tables["isotype_summary"].append({
                "common_qc_consensus_count": "0" if empty else "200",
                "class_callable_consensus_count": "0" if empty else "100",
                "igm_consensus_count": "0" if empty else "20", "igg_consensus_count": "0" if empty else "60",
                "iga_consensus_count": "0" if empty else "20", "log2_igg_iga": "NA" if empty else "1.58496250072116",
                "igg_iga_status": "no_class_callable_evidence" if empty else "available", **collected,
            })
            for group in ("all_IGH", "IgM", "IgG", "IgA"):
                for region in ("whole_v", "cdr2"):
                    self.tables["shm"].append({**identity, "isotype_group": group, "v_region": region,
                                               "mean_clone_median_shm": "NA" if empty else "0" if index == 0 else "0.123456789012345",
                                               "fraction_clones_zero_median_shm": "NA" if empty else "1" if index == 0 else "0.25",
                                               "n_clones_eligible": "0" if empty else "4", "n_clones_zero_median_shm": "0" if empty else "4" if index == 0 else "1",
                                               "n_members_eligible": "0" if empty else "6", "n_members_measured": "0" if empty else "9"})
        self.save()

    def save(self):
        for name, rows in self.tables.items():
            write_table(self.paths[name], rows)

    def test_exact_source_values_and_curated_column_order(self):
        rows = build_headline(**self.paths)
        self.assertEqual(list(rows[0]), list(HEADLINE_COLUMNS))
        locus_block = ("grouped_consensus_count", "richness", "top_5_support_fraction",
                       "rarefaction_target", "rarefaction_status", "rarefied_inverse_simpson")
        self.assertEqual(len(HEADLINE_COLUMNS), 41)
        self.assertEqual(list(HEADLINE_COLUMNS[:24]), [
            "scope_id", "library_id", "tissue_id", "donor_id", "input_read_pairs",
            *(locus + "_yield_per_million_read_pairs" for locus in ("igh", "igk", "igl", "tra", "trb", "trd", "trg")),
            *("igh_" + field for field in locus_block), *("trb_" + field for field in locus_block)])
        for dropped in ("yield_status", "igg_iga_status", "igh_common_qc_fraction_of_annotated",
                        "trb_grouping_eligible_fraction_of_common_qc", "igh_diversity_status",
                        "igh_shm_eligible_clone_fraction_whole_v"):
            self.assertNotIn(dropped, rows[0])
        self.assertEqual([tuple(row[key] for key in ("scope_id", "library_id", "donor_id", "tissue_id")) for row in rows],
                         [("NA", "0460", "001", "T1"), ("NA", "NA", "001", "T1"), ("NA", "empty", "002", "T2")])
        self.assertEqual(rows[1]["igh_rarefied_inverse_simpson"], "2.125000000000003")
        self.assertEqual(rows[0]["igh_top_5_support_fraction"], "0.950000000000001")
        self.assertEqual(rows[1]["igg_mean_clone_median_shm_whole_v"], "0.123456789012345")
        self.assertEqual(rows[0]["log2_igg_iga"], "1.58496250072116")
        self.assertEqual(rows[0]["class_callable_consensus_count"], "100")
        self.assertEqual(rows[0]["igh_grouped_consensus_count"], "20")
        self.assertEqual(rows[0]["igg_n_shm_eligible_clones_whole_v"], "4")
        for absent in ("trb_inverse_simpson", "igh_top_10_support_fraction", "class_callable_fraction",
                       "common_qc_consensus_count", "igm_consensus_count", "igh_n_members_eligible_whole_v",
                       "igg_fraction_clones_zero_median_shm_whole_v"):
            self.assertNotIn(absent, rows[0])
        self.assertEqual(rows[0]["igh_shm_eligible_fraction_of_measured_members_whole_v"], 6 / 9)
        self.assertNotIn("igh_hill_shannon", rows[0])
        self.assertNotIn("igh_top_1_support_fraction", rows[0])
        self.assertEqual([row["igh_yield_per_million_read_pairs"] for row in rows],
                         ["1", "2", "0"])

    def test_unavailable_metrics_and_measured_zeros_remain_distinct(self):
        rows = build_headline(**self.paths)
        observed, _, empty = rows
        self.assertIsNone(observed["igh_rarefied_inverse_simpson"])
        self.assertEqual(observed["igh_rarefaction_status"], "below_target")
        self.assertEqual(observed["igh_mean_clone_median_shm_whole_v"], "0")
        self.assertEqual(observed["class_igd_fraction"], "0")
        self.assertEqual(empty["igh_yield_per_million_read_pairs"], "0")
        self.assertEqual(empty["igh_grouped_consensus_count"], "0")
        for field in ("igh_rarefied_inverse_simpson", "class_igd_fraction", "igh_mean_clone_median_shm_whole_v", "igh_fraction_clones_zero_median_shm_whole_v", "igh_shm_eligible_fraction_of_measured_members_whole_v"):
            self.assertIsNone(empty[field])

    def test_shm_fractions_preserve_zero_eligibility_and_unavailable_counts(self):
        source = self.tables["shm"][0]
        source.update(n_members_eligible="0", n_clones_eligible="0")
        self.save()
        row = build_headline(**self.paths)[0]
        self.assertEqual(row["igh_shm_eligible_fraction_of_measured_members_whole_v"], 0)
        for missing in ("n_members_eligible", "n_members_measured"):
            original = source[missing]
            source[missing] = ""
            self.save()
            row = build_headline(**self.paths)[0]
            self.assertIsNone(row["igh_shm_eligible_fraction_of_measured_members_whole_v"])
            source[missing] = original
        source.update(n_members_eligible="10")
        self.save()
        with self.assertRaisesRegex(ValueError, "SHM member eligibility"):
            build_headline(**self.paths)

    def test_missing_or_duplicate_required_detail_rows_fail(self):
        for name in self.paths:
            if name == "manifest":
                continue
            original = copy.deepcopy(self.tables[name])
            with self.subTest(table=name, defect="missing"):
                write_table(self.paths[name], original[1:])
                with self.assertRaisesRegex(ValueError, "missing required key"):
                    build_headline(**self.paths)
            with self.subTest(table=name, defect="duplicate"):
                write_table(self.paths[name], original + [original[0]])
                with self.assertRaisesRegex(ValueError, "duplicate key"):
                    build_headline(**self.paths)
            write_table(self.paths[name], original)

    def test_conflicting_source_scope_tissue_or_donor_fails(self):
        original = copy.deepcopy(self.tables["shm"])
        for column in ("scope_id", "tissue_id", "donor_id"):
            with self.subTest(column=column):
                rows = copy.deepcopy(original)
                rows[0][column] = "wrong"
                write_table(self.paths["shm"], rows)
                with self.assertRaisesRegex(ValueError, "conflicting " + column):
                    build_headline(**self.paths)
        write_table(self.paths["shm"], original)
        with self.assertRaisesRegex(ValueError, "matching the manifest"):
            build_headline(**self.paths, scope_id="wrong")

    def test_missing_numeric_and_malformed_source_cells_fail(self):
        rows = copy.deepcopy(self.tables["diversity"])
        for row in rows:
            del row["rarefaction_status"]
        write_table(self.paths["diversity"], rows)
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            build_headline(**self.paths)
        rows = copy.deepcopy(self.tables["diversity"])
        rows[0]["rarefied_inverse_simpson"] = "Inf"
        write_table(self.paths["diversity"], rows)
        with self.assertRaisesRegex(ValueError, "must be finite"):
            build_headline(**self.paths)

    def test_a_library_without_input_read_pairs_keeps_its_status_and_unavailable_yield(self):
        self.tables["yields"][0].update(yield_status="no_input_read_pairs", input_read_pairs="0",
                                        igh_yield_per_million_read_pairs="")
        self.save()
        row = build_headline(**self.paths)[0]
        self.assertEqual(row["input_read_pairs"], "0")
        self.assertNotIn("yield_status", row)
        self.assertIsNone(row["igh_yield_per_million_read_pairs"])
        self.tables["yields"][0]["yield_status"] = ""
        self.save()
        with self.assertRaisesRegex(ValueError, "yield_status requires an availability reason"):
            build_headline(**self.paths)

    def test_input_read_pairs_come_from_the_yield_table(self):
        self.tables["yields"][0]["input_read_pairs"] = "2"
        self.save()
        self.assertEqual(build_headline(**self.paths)[0]["input_read_pairs"], "2")
        for row in self.tables["yields"]:
            del row["input_read_pairs"]
        self.save()
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            build_headline(**self.paths)

    def test_scientific_sheet_is_atomic_and_rejects_extra_columns(self):
        rows = build_headline(**self.paths)
        output = self.directory / "headline.tsv"
        write_headline(rows, output)
        self.assertEqual(list(read_table(output)[0]), list(HEADLINE_COLUMNS))
        self.assertEqual(read_table(output)[2]["igh_rarefied_inverse_simpson"], "")
        released = self.directory / "released_headline.tsv"
        released.hardlink_to(output)
        previous = released.read_bytes()
        rows[0]["igh_rarefied_inverse_simpson"] = "2"
        write_headline(rows, output)
        self.assertEqual(released.read_bytes(), previous)
        self.assertEqual(read_table(output)[0]["igh_rarefied_inverse_simpson"], "2")
        rows[0]["meta_group"] = "Private annotation"
        with self.assertRaisesRegex(ValueError, "fields not in fieldnames"):
            write_headline(rows, output)
        self.assertNotIn("Private annotation", output.read_text())


if __name__ == "__main__":
    unittest.main()
