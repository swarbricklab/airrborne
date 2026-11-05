"""Donor collection preserves source populations, identifiers, missing values and key order.

Run: python test/run_tests.py --component scopes --level unit
"""

from pathlib import Path
import tempfile
import unittest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib import scope_tables as MODULE


class DonorCollectionTest(unittest.TestCase):
    def test_merge_preserves_identifier_and_missing_strings(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first, second, empty, output = [root / name for name in ("a", "b", "c", "out")]
            first.write_text("library_id\tvalue\n0460\tNA\nz\t\n")
            second.write_text("library_id\tvalue\nNA\t0\n")
            empty.write_text("library_id\tvalue\n")
            MODULE.collect_donor_tables([first, second, empty], output, ["library_id"])
            self.assertEqual(output.read_text(),
                             "library_id\tvalue\n0460\tNA\nNA\t0\nz\t\n")

    def test_numeric_support_order_is_projected_without_internal_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "shard", Path(directory) / "out"
            source.write_text("library_id\tclone_id\tsort_group_support\n"
                              "L1\tc10\t-10\nL1\tc2\t-2\nL1\tc1\t-1\n")
            MODULE.collect_donor_tables([source], output,
                ["library_id", "sort_group_support", "clone_id"])
            self.assertEqual(output.read_text(), "library_id\tclone_id\nL1\tc10\nL1\tc2\nL1\tc1\n")

    def test_duplicate_source_identifiers_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.tsv"
            source.write_text("library_id\tvalue\n0460\t1\n0460\t2\n")
            with self.assertRaises(ValueError):
                MODULE.collect_donor_tables([source], root / "out.tsv", ["library_id"])

    def test_merge_follows_declared_value_order(self):
        keys = ["library_id", "isotype_group", "v_region"]
        orders = {"isotype_group": ["all_IGH", "IgM", "IgG", "IgA", "IgA1"],
                  "v_region": ["whole_v", "fwr1", "cdr1"]}
        header = "library_id\tisotype_group\tv_region\tvalue\n"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first, second, output = [root / name for name in ("a", "b", "out")]
            first.write_text(header + "L1\tall_IGH\twhole_v\t1\nL1\tall_IGH\tcdr1\t2\n"
                                      "L1\tIgA\twhole_v\t3\nL1\tIgA1\twhole_v\t4\n")
            second.write_text(header + "L0\tIgM\tfwr1\t5\nL2\tall_IGH\tfwr1\t6\n")
            MODULE.collect_donor_tables([first, second], output, keys, orders)
            # all_IGH precedes the classes, and whole_v precedes its subregions,
            # although text order would place both last.
            self.assertEqual(output.read_text(), header + "L0\tIgM\tfwr1\t5\n"
                             "L1\tall_IGH\twhole_v\t1\nL1\tall_IGH\tcdr1\t2\n"
                             "L1\tIgA\twhole_v\t3\nL1\tIgA1\twhole_v\t4\n"
                             "L2\tall_IGH\tfwr1\t6\n")
            collected = output.read_text()
            lexical = root / "lexical"
            lexical.write_text(header + "L1\tIgA\twhole_v\t3\nL1\tall_IGH\twhole_v\t1\n")
            with self.assertRaisesRegex(ValueError, "Unsorted"):
                MODULE.collect_donor_tables([lexical], output, keys, orders)
            undeclared = root / "undeclared"
            undeclared.write_text(header + "L1\tIgE\twhole_v\t1\n")
            with self.assertRaisesRegex(ValueError, "declared order omits"):
                MODULE.collect_donor_tables([undeclared], output, keys, orders)
            self.assertEqual(output.read_text(), collected)

    def test_gated_measurements_are_empty_outside_eligibility(self):
        header = "library_id\tsequence_id\trate\tcount\teligible\tother\n"
        shard = (header + "L1\ts1\t0.1\t3\tTRUE\t0\n"
                          "L1\ts2\t0.2\t0\tFALSE\t0\n"
                          "L1\ts3\t\t\tFALSE\tNA\n")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "shard", root / "out"
            source.write_text(shard)
            MODULE.collect_donor_tables([source], output, ["library_id", "sequence_id"],
                                        gated_columns={"rate": "eligible", "count": "eligible"})
            # A measured zero outside eligibility is emptied too; other columns
            # and the donor shard keep their values.
            self.assertEqual(output.read_text(), header + "L1\ts1\t0.1\t3\tTRUE\t0\n"
                             "L1\ts2\t\t\tFALSE\t0\nL1\ts3\t\t\tFALSE\tNA\n")
            self.assertEqual(source.read_text(), shard)
            with self.assertRaisesRegex(ValueError, "gated and eligibility columns"):
                MODULE.collect_donor_tables([source], output, ["library_id", "sequence_id"],
                                            gated_columns={"rate": "absent"})

    def test_projection_renames_and_dropped_rows_leave_the_shard_complete(self):
        header = "sequence_id\tlibrary_id\tv_gene\tsequence\tstate\n"
        shard = (header + "s1\tL1\tIGHV1-2\tACGT\tscored\n"
                          "s2\tL1\tIGHV3-23\tTTTT\tnone\n"
                          "s0\tL2\tIGHV1-2\tGGGG\tscored\n")
        keys = ["library_id", "sequence_id"]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "shard", root / "out"
            source.write_text(shard)
            # The output takes its columns and names from the projection and
            # leads with the identifiers; a dropped row still has to be in order.
            MODULE.collect_donor_tables([source], output, keys,
                                        output_columns={"sequence_id": "sequence_id",
                                                        "v_gene_key": "v_gene",
                                                        "library_id": "library_id"},
                                        dropped_values={"state": ["none"]})
            self.assertEqual(output.read_text(), "library_id\tsequence_id\tv_gene_key\n"
                             "L1\ts1\tIGHV1-2\nL2\ts0\tIGHV1-2\n")
            self.assertEqual(source.read_text(), shard)
            unsorted = root / "unsorted"
            unsorted.write_text(header + "s2\tL1\tIGHV3-23\tTTTT\tnone\n"
                                         "s1\tL1\tIGHV1-2\tACGT\tscored\n")
            with self.assertRaisesRegex(ValueError, "Unsorted"):
                MODULE.collect_donor_tables([unsorted], output, keys,
                                            dropped_values={"state": ["none"]})
            for arguments, message in (({"output_columns": {"v_gene_key": "absent"}}, "every column"),
                                       ({"dropped_values": {"absent": ["none"]}}, "drop rows")):
                with self.subTest(arguments=arguments), self.assertRaisesRegex(ValueError, message):
                    MODULE.collect_donor_tables([source], output, keys, **arguments)


if __name__ == "__main__":
    unittest.main()
