"""Reconcile complete coverage with all-member QC denominators independently."""
import csv
from pathlib import Path
import tempfile
import unittest

from workflow.lib.scope_qc_audit import add_shm_coverage_counts


class CoverageAuditTests(unittest.TestCase):
    def test_any_region_counts_a_member_once_and_missing_is_not_zero_shm(self):
        regions = ("whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "coverage.tsv"
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("library_id", "sequence_id", "v_region", "measurement_status", "shm_eligible"))
                for member in range(3):
                    for region in regions:
                        eligible = member == 0 or member == 1 and region in ("cdr1", "cdr2")
                        state = "unmeasured" if member == 2 else "eligible" if eligible else "measured_ineligible"
                        writer.writerow(("L1", str(member), region, state, "TRUE" if eligible else "FALSE"))
            qc = [{"library_id": "L1", "locus": "IGH"}, {"library_id": "empty", "locus": "IGH"},
                  {"library_id": "L1", "locus": "TRB"}]
            add_shm_coverage_counts(qc, path)
            self.assertEqual(qc[0]["n_grouped_igh_members"], "3")
            self.assertEqual(qc[0]["n_shm_measured_members_whole_v"], "2")
            self.assertEqual(qc[0]["n_shm_eligible_members_whole_v"], "1")
            self.assertEqual(qc[0]["n_shm_eligible_members_any_subregion"], "2")
            self.assertEqual(float(qc[0]["shm_eligible_fraction_of_grouped_members_whole_v"]), 1 / 3)
            self.assertEqual(qc[1]["shm_eligible_fraction_of_grouped_members_whole_v"], "")
            self.assertEqual(qc[2], {"library_id": "L1", "locus": "TRB"})
            path.write_text(path.read_text().rsplit("\n", 2)[0] + "\n")
            with self.assertRaisesRegex(ValueError, "six regions"):
                add_shm_coverage_counts(qc, path)
