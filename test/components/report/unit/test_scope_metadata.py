"""Display metadata must preserve identities, observations and unavailable values.

Run: python test/run_tests.py --component report --level unit
"""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from workflow.lib.report_metadata import (
    join_scope_metadata,
    load_library_metadata,
    load_tissue_metadata,
)


class ScopeMetadataTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.libraries = [
            {"library_id": "0460", "tissue_id": "T1", "donor_id": "001", "metric": 0},
            {"library_id": "NA", "tissue_id": "T1", "donor_id": "001", "metric": None},
            {"library_id": "L3", "tissue_id": "T2", "donor_id": "001", "metric": 2.5},
            {"library_id": "L4", "tissue_id": "T3", "donor_id": "002", "metric": 4},
        ]

    def write(self, text, name="metadata.csv"):
        path = self.root / name
        path.write_text(text)
        return path

    def test_tissue_join_preserves_repeated_libraries_and_missing_tissues(self):
        metadata = load_tissue_metadata(self.write("tissue_id,response\nT1,NA\nT2,\n"),
                                        ["response"])
        before = copy.deepcopy(self.libraries)
        rows, descriptions = join_scope_metadata(self.libraries, metadata)
        self.assertEqual(len(rows), len(before))
        self.assertEqual([r["library_id"] for r in rows], ["0460", "NA", "L3", "L4"])
        self.assertEqual([r["meta_response"] for r in rows], ["NA", "NA", None, None])
        self.assertEqual([r["metric"] for r in rows], [0, None, 2.5, 4])
        self.assertEqual(self.libraries, before)
        self.assertEqual(descriptions, [{"column": "meta_response", "type": "category",
                                        "label": "response", "source_level": "tissue"}])
        json.dumps({"rows": rows, "descriptors": descriptions}, allow_nan=False)

    def test_numeric_missing_values_are_not_zero(self):
        path = self.write("tissue_id,purity\nT1,NA\nT2,\nT3,0\n")
        metadata = load_tissue_metadata(path, ["purity"])
        rows, descriptions = join_scope_metadata(self.libraries, metadata)
        self.assertEqual([r["meta_purity"] for r in rows], [None, None, None, 0.0])
        self.assertEqual(descriptions[0]["type"], "number")

    def test_shared_tissue_file_exposes_only_the_selected_columns(self):
        path = self.write("tissue_id,label,score,private\nT1,A,3,unreleased\nT2,B,4,unreleased\n")
        first = load_tissue_metadata(path, ["label"])
        second = load_tissue_metadata(path, ["score"])
        a, _ = join_scope_metadata(self.libraries, first)
        b, _ = join_scope_metadata(self.libraries, second)
        self.assertIn("meta_label", a[0])
        self.assertNotIn("meta_score", a[0])
        self.assertEqual(b[0]["meta_score"], 3)
        self.assertNotIn("meta_label", b[0])
        self.assertNotIn("unreleased", json.dumps([a, b]))

    def test_csv_and_tsv_preserve_literal_identifier_na(self):
        libraries = [{"library_id": "NA", "tissue_id": "0460", "donor_id": "NA"}]
        for suffix, delimiter in [("csv", ","), ("tsv", "\t")]:
            with self.subTest(suffix=suffix):
                text = delimiter.join(["tissue_id", "donor_id", "label"]) + "\n"
                text += delimiter.join(["0460", "NA", "0460"]) + "\n"
                metadata = load_tissue_metadata(self.write(text, "metadata." + suffix),
                                                ["label"])
                rows, _ = join_scope_metadata(libraries, metadata)
                self.assertEqual(rows[0], {**libraries[0], "meta_label": "0460"})

    def test_extra_tissues_are_validated_but_never_exported(self):
        metadata = load_tissue_metadata(self.write("tissue_id,label\nT1,inside\noutside,private\n"),
                                        ["label"])
        rows, _ = join_scope_metadata(self.libraries, metadata)
        self.assertNotIn("private", json.dumps(rows))
        with self.assertRaisesRegex(ValueError, "duplicate tissue_id"):
            load_tissue_metadata(self.write("tissue_id,label\noutside,a\noutside,b\n"),
                                  ["label"])
        with self.assertRaisesRegex(ValueError, "invalid numeric"):
            load_tissue_metadata(self.write("tissue_id,timepoint_order\noutside,unknown\n"))

    def test_tissue_annotations_can_vary_within_a_donor(self):
        metadata = load_tissue_metadata(self.write("tissue_id,diagnosis\nT1,A\nT2,B\n"))
        rows, _ = join_scope_metadata(self.libraries, metadata)
        self.assertEqual([row["meta_diagnosis"] for row in rows], ["A", "A", "B", None])

    def test_optional_donor_ids_check_only_supplied_values(self):
        metadata = load_tissue_metadata(self.write("tissue_id,donor_id,label\nT1,,A\nT2,001,B\n"),
                                        ["label"])
        rows, _ = join_scope_metadata(self.libraries, metadata)
        self.assertEqual(rows[0]["donor_id"], "001")
        metadata = load_tissue_metadata(self.write("tissue_id,donor_id,label\nT1,1,A\n"),
                                        ["label"])
        with self.assertRaisesRegex(ValueError, "supplied donor_id"):
            join_scope_metadata(self.libraries, metadata)

    def test_library_context_does_not_spread_to_other_libraries(self):
        path = self.write("library_id,tissue_id,donor_id,batch,path,unexposed\n"
                          "0460,T1,001,A,fq/one,private\nNA,T1,001,B,fq/two,private\n",
                          "samplesheet.csv")
        metadata = load_library_metadata(path, ["batch"])
        rows, descriptions = join_scope_metadata(self.libraries, library_metadata=metadata)
        self.assertEqual([r["meta_batch"] for r in rows], ["A", "B", None, None])
        self.assertEqual(descriptions[0]["source_level"], "library")
        self.assertNotIn("private", json.dumps(rows))
        self.assertNotIn("meta_path", rows[0])

    def test_library_context_checks_tissue_and_donor_links(self):
        for column in ("tissue_id", "donor_id"):
            with self.subTest(column=column):
                path = self.write(f"library_id,{column},batch\n0460,wrong,A\n", "samplesheet.csv")
                metadata = load_library_metadata(path, ["batch"])
                with self.assertRaisesRegex(ValueError, f"supplied {column}"):
                    join_scope_metadata(self.libraries, library_metadata=metadata)

    def test_timepoint_fields_do_not_reorder_or_fill_missing_visits(self):
        path = self.write("tissue_id,timepoint,timepoint_order\nT1,T3,3\nT2,T1,1\nT3,NA,\n")
        metadata = load_tissue_metadata(path, ["timepoint", "timepoint_order"])
        rows, _ = join_scope_metadata(self.libraries, metadata)
        self.assertEqual([r["library_id"] for r in rows], ["0460", "NA", "L3", "L4"])
        self.assertEqual([r["meta_timepoint_order"] for r in rows], [3, 3, 1, None])
        self.assertEqual(rows[-1]["meta_timepoint"], "NA")

    def test_multiple_tissues_can_share_a_visit_order(self):
        path = self.write("tissue_id,timepoint_order\nT1,1\nT2,1\n")
        metadata = load_tissue_metadata(path, ["timepoint_order"])
        rows, _ = join_scope_metadata(self.libraries, metadata)
        self.assertEqual([r["meta_timepoint_order"] for r in rows], [1, 1, 1, None])

    def test_numeric_visit_labels_remain_categories_and_order_is_numeric(self):
        metadata = load_tissue_metadata(self.write("tissue_id,timepoint,timepoint_order\nT1,01,1\nT2,2,2\n"))
        rows, _ = join_scope_metadata(self.libraries, metadata)
        self.assertEqual(rows[0]["meta_timepoint"], "01")
        self.assertEqual(rows[0]["meta_timepoint_order"], 1)
        with self.assertRaisesRegex(ValueError, "tissue-level"):
            load_library_metadata(self.write("library_id,timepoint\n0460,1\n"))

    def test_timepoint_label_requires_one_order_across_donors(self):
        metadata = load_tissue_metadata(self.write(
            "tissue_id,timepoint,timepoint_order\nT1,Visit A,1\nT3,Visit A,2\n"))
        with self.assertRaisesRegex(ValueError, "Timepoint 'Visit A'.*1 .*T1.*2 .*T3"):
            join_scope_metadata(self.libraries, metadata)

    def test_distinct_timepoint_labels_require_distinct_orders(self):
        metadata = load_tissue_metadata(self.write(
            "tissue_id,timepoint,timepoint_order\nT1,Visit A,0\nT2,Visit B,0.0\n"))
        with self.assertRaisesRegex(ValueError, "timepoint_order 0.*'Visit A'.*T1.*'Visit B'.*T2"):
            join_scope_metadata(self.libraries, metadata)

    def test_shared_visits_and_incomplete_timepoints_preserve_all_libraries(self):
        for last_visit in ("Visit A,", ",2", ","):
            with self.subTest(last_visit=last_visit):
                metadata = load_tissue_metadata(self.write(
                    "tissue_id,timepoint,timepoint_order\nT1,Visit A,0\nT2,Visit A,0.0\nT3," + last_visit + "\n"))
                rows, _ = join_scope_metadata(self.libraries, metadata)
                self.assertEqual([row["library_id"] for row in rows], ["0460", "NA", "L3", "L4"])
                self.assertEqual([row["meta_timepoint_order"] for row in rows[:3]], [0, 0, 0])
                self.assertEqual(rows[-1]["meta_timepoint"], "Visit A" if last_visit.startswith("Visit") else None)
                self.assertEqual(rows[-1]["meta_timepoint_order"], 2 if last_visit == ",2" else None)

    def test_timepoint_conflicts_apply_only_to_selected_scope_and_columns(self):
        path = self.write("tissue_id,timepoint,timepoint_order\nT1,Visit A,1\nT2,Visit A,2\n")
        metadata = load_tissue_metadata(path)
        for libraries in (self.libraries[:2], self.libraries[2:]):
            rows, _ = join_scope_metadata(libraries, metadata)
            self.assertEqual(len(rows), len(libraries))
        for columns in (["timepoint"], ["timepoint_order"]):
            rows, _ = join_scope_metadata(self.libraries, load_tissue_metadata(path, columns))
            self.assertEqual(len(rows), len(self.libraries))

    def test_visit_order_does_not_infer_treatment_phase_or_require_consecutive_integers(self):
        metadata = load_tissue_metadata(self.write(
            "tissue_id,timepoint,timepoint_order,treatment_phase\n"
            "T1,3,-2,midtreatment\nT2,1,0.5,midtreatment\nT3,NA,4,pretreatment\n"))
        rows, _ = join_scope_metadata(self.libraries, metadata)
        self.assertEqual([row["meta_timepoint"] for row in rows], ["3", "3", "1", "NA"])
        self.assertEqual([row["meta_timepoint_order"] for row in rows], [-2, -2, 0.5, 4])
        self.assertEqual([row["meta_treatment_phase"] for row in rows], ["midtreatment"] * 3 + ["pretreatment"])

    def test_file_alone_preserves_headers_and_infers_types(self):
        path = self.write('tissue_id,donor_id,Score (%),Treatment arm,Code,Empty\nT1,001,-2.5,A,01,\nT2,001,0,B,02,\nT3,002,NA,C,03,\n')
        metadata = load_tissue_metadata(path)
        rows, descriptors = join_scope_metadata(self.libraries, metadata)
        self.assertEqual([item["label"] for item in descriptors], ["Score (%)", "Treatment arm", "Code", "Empty"])
        self.assertEqual([item["type"] for item in descriptors], ["number", "category", "category", "category"])
        self.assertEqual([row["meta_Score (%)"] for row in rows], [-2.5, -2.5, 0, None])
        self.assertEqual(rows[0]["meta_Code"], "01")
        self.assertIsNone(rows[0]["meta_Empty"])
        self.assertNotIn("meta_donor_id", rows[0])
        self.assertTrue(all(set(item) == {"column", "label", "type", "source_level"} for item in descriptors))

    def test_mixed_columns_remain_text_and_type_can_change_after_edit(self):
        path = self.write("tissue_id,score\nT1,1\nT2,2\n")
        self.assertEqual(load_tissue_metadata(path).descriptors[0]["type"], "number")
        path.write_text("tissue_id,score\nT1,1\nT2,unknown\n")
        metadata = load_tissue_metadata(path)
        self.assertEqual(metadata.descriptors[0]["type"], "category")
        self.assertEqual(metadata.rows["T1"]["meta_score"], "1")

    def test_columns_selection_rejects_definitions_and_duplicate_names(self):
        path = self.write("tissue_id,score\nT1,1\n")
        for columns in ([], ["score", "score"], [""], "score", {"score": {"type": "number"}}):
            with self.subTest(columns=columns), self.assertRaisesRegex(ValueError, "list"):
                load_tissue_metadata(path, columns)

    def test_duplicate_headers_malformed_rows_and_missing_declarations_fail(self):
        for text, message in [("tissue_id,label,label\nT1,A,B\n", "duplicate column"),
                               ("tissue_id,label\nT1,A\nT1,B\n", "duplicate tissue_id"),
                               ("tissue_id,label\nT1,A,B\n", "row width"),
                               ("tissue_id,label\n,A\n", "empty tissue_id"),
                               ("tissue_id\nT1\n", "missing selected")]:
            with self.subTest(text=text):
                with self.assertRaisesRegex(ValueError, message):
                    load_tissue_metadata(self.write(text), ["label"])

    def test_invalid_visit_order_is_rejected(self):
        for value in ("unknown", "NaN", "nan", "inf", "-inf"):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    load_tissue_metadata(self.write(f"tissue_id,timepoint_order\nT1,{value}\n"),
                                          ["timepoint_order"])

    def test_metadata_names_and_namespaces_cannot_collide(self):
        with self.assertRaisesRegex(ValueError, "identifiers"):
            load_tissue_metadata(self.write("tissue_id,donor_id\nT1,001\n"), ["donor_id"])
        tissue = load_tissue_metadata(self.write("tissue_id,batch\nT1,A\n"), ["batch"])
        library = load_library_metadata(self.write("library_id,batch\n0460,B\n", "samplesheet.csv"),
                                         ["batch"])
        with self.assertRaisesRegex(ValueError, "collide"):
            join_scope_metadata(self.libraries, tissue, library)
        with self.assertRaisesRegex(ValueError, "overwrite"):
            join_scope_metadata([{**self.libraries[0], "meta_batch": "existing"}], tissue)

    def test_scope_crosswalk_must_be_unambiguous(self):
        with self.assertRaisesRegex(ValueError, "Duplicate library_id"):
            join_scope_metadata([self.libraries[0], self.libraries[0]])
        with self.assertRaisesRegex(ValueError, "more than one donor"):
            join_scope_metadata([self.libraries[0], {**self.libraries[1], "donor_id": "002"}])
        with self.assertRaisesRegex(ValueError, "nonempty string"):
            join_scope_metadata([{**self.libraries[0], "library_id": 460}])

    def test_absent_metadata_adds_no_columns_or_controls(self):
        rows, descriptions = join_scope_metadata(self.libraries)
        self.assertEqual(rows, self.libraries)
        self.assertEqual(descriptions, [])
        self.assertIsNot(rows[0], self.libraries[0])
        self.assertEqual(join_scope_metadata([]), ([], []))

    def test_scope_specific_labels_can_differ_for_same_library(self):
        first = load_tissue_metadata(self.write("tissue_id,label\nT1,A\n", "first.csv"),
                                     ["label"])
        second = load_tissue_metadata(self.write("tissue_id,label\nT1,B\n", "second.csv"),
                                      ["label"])
        self.assertEqual(join_scope_metadata(self.libraries, first)[0][0]["meta_label"], "A")
        self.assertEqual(join_scope_metadata(self.libraries, second)[0][0]["meta_label"], "B")


if __name__ == "__main__":
    unittest.main()
