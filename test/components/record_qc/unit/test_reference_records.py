"""The fixture oracle must reject plausible missing or misreported QC evidence.

Run: python test/run_tests.py --component record_qc --level unit
"""

import copy
import csv
from pathlib import Path
import unittest

from test.components.record_qc.integration.reference_records import check_audit, check_clone_partitions, expected_early_rows
from test.support.reference_tables import read_table


class EarlyAuditTests(unittest.TestCase):
    def setUp(self):
        root = Path(__file__).resolve().parents[3]
        records = [row for row in read_table(root / "data/expected/expected_records.tsv")
                   if row["library_id"] == "SYN_D2_T2"]
        with (root / "data/post_trust4/libraries/SYN_D2_T2/trust4/SYN_D2_T2_cdr3.out").open() as handle:
            cdr3 = list(csv.reader(handle, delimiter="\t"))
        self.expected = expected_early_rows(records, cdr3, .1)
        self.actual = [{key: "" if value is None else str(value).upper() if isinstance(value, bool)
                       else str(value) for key, value in row.items()} for row in self.expected]

    def test_productivity_loss_has_original_igh_denominator(self):
        row, = [row for row in self.expected if row["stage"] == "qc_productive" and row["locus"] == "IGH"]
        self.assertEqual((row["n_units_in"], row["n_units_out"], row["consensus_count_in"],
                          row["consensus_count_out"]), (3, 2, 1710, 1536))
        self.assertEqual(check_audit(self.actual, self.expected, "D2/T2"), 91)

    def test_deleting_a_stage_or_zero_locus_is_detected(self):
        for remove in (lambda row: row["stage"] == "igblast_makedb",
                       lambda row: row["locus"] == "IGK"):
            with self.assertRaisesRegex(AssertionError, "early audit stages/loci"):
                check_audit([row for row in self.actual if not remove(row)], self.expected, "deleted")

    def test_counts_status_units_and_diagnostic_semantics_are_checked(self):
        for stage, column, value in (
                ("qc_productive", "n_units_out", "3"),
                ("qc_productive", "consensus_count_out", "1710"),
                ("qc_low_n", "decision_status", "disabled"),
                ("igblast_makedb", "input_unit", "contig"),
                ("igblast_locus_agreement", "kind", "eligibility"),
                ("igblast_makedb", "locus_basis", "annotated"),
                ("trust4_reported_airr", "support_conserved", "FALSE")):
            with self.subTest(stage=stage, column=column):
                rows = copy.deepcopy(self.actual)
                row, = [row for row in rows if row["stage"] == stage and row["locus"] == "IGH"]
                row[column] = value
                with self.assertRaisesRegex(AssertionError, column):
                    check_audit(rows, self.expected, "corrupted")

    def test_disabled_ceiling_and_empty_library_remain_explicit(self):
        rows = expected_early_rows([], [], 1)
        self.assertEqual(len(rows), 91)
        self.assertTrue(all(row["n_units_in"] == row["n_units_out"] == 0 for row in rows))
        self.assertEqual({row["decision_status"] for row in rows if row["stage"] == "qc_low_n"}, {"disabled"})
        self.assertTrue(all(row["consensus_count_in"] is None for row in rows if row["analysis"] == "assembly"))


class ClonePartitionTests(unittest.TestCase):
    def test_labels_can_change_but_merging_splitting_or_double_support_fail(self):
        expected = {"early": ("family", 7), "late": ("family", 11),
                    "distinct": ("distinct", 5), "other_donor": ("distinct", 13)}
        members = {key: dict(donor_id="D2" if key == "other_donor" else "D1",
                            clone_id="123" if key == "distinct" else "045",
                            member_consensus_count=str(value[1])) for key, value in expected.items()}
        check_clone_partitions(members, expected, "valid renumbering")
        for key, column, value in (("late", "clone_id", "999"),
                                   ("distinct", "clone_id", "045"),
                                   ("early", "member_consensus_count", "14")):
            with self.subTest(key=key, column=column):
                altered = copy.deepcopy(members)
                altered[key][column] = value
                with self.assertRaises(AssertionError):
                    check_clone_partitions(altered, expected, "corrupted")
        with self.assertRaisesRegex(AssertionError, "member set"):
            check_clone_partitions({key: row for key, row in members.items() if key != "late"}, expected, "missing")
