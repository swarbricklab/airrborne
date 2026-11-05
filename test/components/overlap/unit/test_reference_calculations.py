"""Check overlap denominators and unavailable endpoints independently.
Constructed cases anchor the independent system-result calculations.

Run: python test/run_tests.py --component overlap --level unit
"""

import unittest
import csv
import json
from pathlib import Path
import tempfile
from itertools import combinations
from test.components.overlap.integration.reference_calculations import pairwise, verify_sharing


class SharingPopulationTests(unittest.TestCase):
    def test_missing_pairs_and_zero_support_rows_fail(self):
        with tempfile.TemporaryDirectory() as workdir:
            root = Path(workdir)
            (root / "metrics").mkdir()
            (root / "manifest").mkdir()
            (root / "manifest/scope_manifest.tsv").write_text(
                "donor_id\tlibrary_id\ttissue_id\nD1\tL1\tT2\nD1\tL2\tT1\nD1\tL3\tT3\n")
            (root / "metrics/group_abundance.tsv").write_text(
                "donor_id\tlocus\tlibrary_id\tgroup_key\tconsensus_count\n"
                "D1\tIGH\tL1\tshared\t2\nD1\tIGH\tL1\tprivate\t8\n"
                "D1\tIGH\tL2\tshared\t3\n")
            # Tissue order reverses the first two library IDs. Their unequal
            # shared fractions expose a numerical check using the wrong end.
            vectors = {"L2": {"shared": 3}, "L1": {"shared": 2, "private": 8}, "L3": {}}
            pairs = [dict(donor_id="D1", locus="IGH", group_unit="clone", library_id_1=a, library_id_2=b,
                          **pairwise(vectors[a], vectors[b]))
                     for a, b in combinations(vectors, 2)]
            shared = [dict(donor_id="D1", locus="IGH", group_key="shared", library_id=library,
                           n_libraries_in_donor=3, n_libraries_detected=2,
                           consensus_count=vector.get("shared", 0),
                           support_fraction=vector["shared"] / sum(vector.values()) if vector else None)
                      for library, vector in sorted(vectors.items())]

            def write(name: str, rows: list[dict]) -> None:
                with (root / "metrics" / name).open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                    writer.writeheader()
                    writer.writerows(rows)

            params = {"overlap": {"loci": ["IGH"], "min_consensus_count": {}}}
            write("donor_overlap_pairs.tsv", pairs)
            write("shared_groups.tsv", shared)
            summary = dict(donor_id="D1", locus="IGH", group_unit="clone", n_libraries_in_donor=3,
                           n_libraries_nonempty=2, n_pairs_total=3, n_pairs_available=1, n_groups_union=2,
                           n_groups_shared_2plus=1, n_groups_shared_all=0)
            write("donor_overlap_summary.tsv", [summary])
            self.assertEqual(verify_sharing(root, params), (3, 1))
            for column, value in (("n_pairs_total", 2), ("n_pairs_available", 3), ("n_groups_union", 3),
                                  ("n_groups_shared_2plus", 2), ("group_unit", "clonotype")):
                write("donor_overlap_summary.tsv", [{**summary, column: value}])
                self.assertRaisesRegex(AssertionError, column, verify_sharing, root, params)
            write("donor_overlap_summary.tsv", [summary])
            self.assertEqual(pairs[0]["shared_support_fraction_1"], 1)
            self.assertEqual(pairs[0]["shared_support_fraction_2"], .2)
            incorrect = {**pairs[0], "shared_support_fraction_1": .2,
                         "shared_support_fraction_2": 1}
            write("donor_overlap_pairs.tsv", [incorrect, *pairs[1:]])
            self.assertRaises(AssertionError, verify_sharing, root, params)
            for column in ("grouped_consensus_count_1", "grouped_consensus_count_2"):
                incorrect = {**pairs[0], column: pairs[0][column] + 1}
                write("donor_overlap_pairs.tsv", [incorrect, *pairs[1:]])
                self.assertRaisesRegex(AssertionError, column, verify_sharing, root, params)
            for duplicate in (pairs[0], {**pairs[0], "library_id_1": "L1", "library_id_2": "L2"}):
                write("donor_overlap_pairs.tsv", [*pairs, duplicate])
                self.assertRaises(AssertionError, verify_sharing, root, params)
            write("donor_overlap_pairs.tsv", pairs[:1])
            self.assertRaisesRegex(AssertionError, "comparison set", verify_sharing, root, params)
            write("donor_overlap_pairs.tsv", pairs)
            write("shared_groups.tsv", shared[:2])
            self.assertRaisesRegex(AssertionError, "clone-by-library", verify_sharing, root, params)
            write("shared_groups.tsv", shared[::-1])
            self.assertRaisesRegex(AssertionError, "clone-support order", verify_sharing, root, params)

from test.components.overlap.integration.reference_calculations import pairwise


class AnalyticalReferenceTests(unittest.TestCase):
    def test_sharing_uses_all_support_and_preserves_unavailable_endpoint(self):
            row = pairwise({"shared": 2, "private": 98}, {"shared": 5})
            self.assertEqual(row["shared_support_fraction_1"], .02)
            self.assertEqual(row["shared_support_fraction_2"], 1)
            row = pairwise({}, {"private": 5})
            self.assertIsNone(row["shared_support_fraction_1"])
            self.assertEqual(row["shared_support_fraction_2"], 0)
            self.assertIsNone(row["jaccard"])
