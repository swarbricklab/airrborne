"""Check reviewed values, exact identifiers, missing values and numerical tolerances.
Constructed cases anchor the independent system-result calculations.

Run: python test/run_tests.py --component inputs --level unit
"""

import unittest
import csv
import json
from pathlib import Path
import tempfile
from test.support.expectations import matches_expected, verify_expected_values


class CommittedExpectationTests(unittest.TestCase):
    """The comparison rules the committed expectations in test/ rely on."""

    def write_example(self, directory: Path, table: str, selector: dict,
                      column: str, expected: str, kind: str) -> Path:
        examples = directory / "examples.tsv"
        with examples.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["table", "selector", "column", "expected", "type",
                             "review_basis"])
            writer.writerow([table, json.dumps(selector), column, expected, kind,
                             "Constructed case."])
        return examples

    def test_identifier_zero_padding_is_significant(self):
        self.assertFalse(matches_expected("1", "001", "string"))
        with tempfile.TemporaryDirectory() as workdir:
            root = Path(workdir)
            (root / "counts.tsv").write_text("library_id\tsupport\n001\t2\n1\t3\n")
            examples = self.write_example(root, "counts.tsv", {"library_id": "001"},
                                          "support", "2", "integer")
            self.assertEqual(verify_expected_values(root, examples), 1)
            examples = self.write_example(root, "counts.tsv", {"library_id": "001"},
                                          "support", "3", "integer")
            self.assertRaises(AssertionError, verify_expected_values, root, examples)

    def test_missing_tokens_and_zero_are_distinct(self):
        for observed, expected in (("", "NA"), ("NaN", "nan"), ("NA", "0"), ("", "0")):
            for kind in ("string", "integer", "number"):
                with self.subTest(observed=observed, expected=expected, kind=kind):
                    self.assertFalse(matches_expected(observed, expected, kind))
        self.assertTrue(matches_expected("NA", "NA", "number"))

    def test_numeric_tolerance_is_explicit(self):
        self.assertFalse(matches_expected("1.000000001", "1", "string"))
        self.assertTrue(matches_expected("1.000000001", "1", "number"))
        self.assertFalse(matches_expected("1.1", "1", "number"))

    def test_integer_counts_are_exact_above_float_precision(self):
        self.assertTrue(matches_expected("2.0", "2", "integer"))
        self.assertFalse(matches_expected("9007199254740993", "9007199254740992",
                                          "integer"))

    def test_expected_examples_match_exactly_one_row(self):
        with tempfile.TemporaryDirectory() as workdir:
            root = Path(workdir)
            table = root / "counts.tsv"
            examples = self.write_example(root, "counts.tsv", {"library_id": "001"},
                                          "support", "2", "integer")
            table.write_text("library_id\tsupport\n001\t2\n001\t2\n")
            self.assertRaises(AssertionError, verify_expected_values, root, examples)
            table.write_text("library_id\tsupport\n002\t3\n")
            self.assertRaises(AssertionError, verify_expected_values, root, examples)

    def test_selected_scope_requires_its_own_expectations(self):
        with tempfile.TemporaryDirectory() as workdir:
            root = Path(workdir)
            table = root / "scopes/project_a/counts.tsv"
            table.parent.mkdir(parents=True)
            table.write_text("library_id\tsupport\n001\t2\n")
            examples = self.write_example(root, "scopes/project_a/counts.tsv",
                                          {"library_id": "001"}, "support", "2", "integer")
            self.assertEqual(verify_expected_values(root, examples, "project_a", {"001"}), 1)
            self.assertRaisesRegex(AssertionError, "No committed expectations",
                                   verify_expected_values, root, examples, "project_b", {"002"})
