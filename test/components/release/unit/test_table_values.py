"""Preserve release identities, missing values and scientific table types.

Run: python test/run_tests.py --component release --level unit
"""

import csv
import gzip
from pathlib import Path
import tempfile
import unittest

from test.components.release.integration.check_tables import (
    Reference, TableContract, check_header, check_references, check_values, load_dictionary,
    read_references)


class ReleaseContractTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "rearrangements.tsv.gz"
        self.contract = TableContract(
            columns=["library_id", "count", "fraction", "status"],
            options={"status": {"available", "empty"}}, primary_key=["library_id"],
            types={"library_id": "string", "count": "integer", "fraction": "numeric", "status": "categorical"})

    def write(self, rows):
        with gzip.open(self.path, "wt", newline="") as handle:
            handle.write("library_id\tcount\tfraction\tstatus\n" + rows)

    def test_gzip_preserves_string_keys_zero_and_empty_measurements(self):
        self.write("0460\t2\t0\tavailable\nNA\t0\t\tempty\n")
        self.assertEqual(check_header(self.path, self.contract.columns, "fixture"), [])
        self.assertEqual(check_values(self.path, self.contract, self.contract.columns), [])

    def test_duplicate_keys_nonfinite_numbers_and_noninteger_counts_fail(self):
        self.write("0460\t1.5\tinf\tavailable\n0460\t2\tNaN\tunknown\n")
        errors = "\n".join(check_values(self.path, self.contract, self.contract.columns))
        for expected in ("primary key", "invalid integer", "invalid numeric", "vocabulary"):
            self.assertIn(expected, errors)


    def test_empty_cell_in_a_required_column_fails(self):
        self.write("0460\t\t0\tavailable\n")
        self.contract.required = {"library_id", "count"}
        errors = "\n".join(check_values(self.path, self.contract, self.contract.columns))
        self.assertIn("count is declared not nullable but has 1 empty cell(s)", errors)

    def test_reference_alternatives_read_in_step(self):
        row = dict(table="abundance", column="group_key", references_table="clones|clonotypes",
                   references_column="clone_id|clonotype_key",
                   references_when="group_unit=clone|group_unit!=clone")
        self.assertEqual(read_references(row), [
            Reference("group_key", "clones", "clone_id", ("group_unit", "=", "clone")),
            Reference("group_key", "clonotypes", "clonotype_key", ("group_unit", "!=", "clone"))])
        self.assertEqual(read_references(dict(row, references_table="")), [])
        with self.assertRaisesRegex(ValueError, "differ in length"):
            read_references(dict(row, references_column="clone_id"))

    def test_references_resolve_within_the_release_on_rows_meeting_their_condition(self):
        bundle = Path(self.directory.name) / "bundle"
        (bundle / "metrics").mkdir(parents=True)
        (bundle / "metrics/clones.tsv").write_text("clone_id\nc1\n")
        (bundle / "metrics/shm.tsv").write_text(
            "clone_id\tunit\nc1\tclone\nc2\tother\n\tclone\n")
        conditional = Reference("clone_id", "clones", "clone_id", ("unit", "=", "clone"))
        contracts = {"clones": TableContract(), "shm": TableContract(references=[conditional])}
        self.assertEqual(check_references(bundle, contracts), [])
        contracts["shm"].references = [Reference("clone_id", "clones", "clone_id")]
        self.assertIn("1 value(s) absent from clones.clone_id",
                      "\n".join(check_references(bundle, contracts)))
        contracts["shm"].references = [Reference("clone_id", "Repertoire.json", "clone_id")]
        self.assertEqual(check_references(bundle, contracts), [])

    def test_light_chain_abundance_and_curve_rows_meet_their_contracts(self):
        # Light chains have support vectors at every library and curves when a
        # target names them.
        dictionary = load_dictionary(Path(__file__).resolve().parents[4] / "docs/data_dictionary.csv")
        for table, measurements in (
                ("group_abundance", dict(group_key="junction_v_j", group_rank=1, n_members=1,
                                               consensus_count=3, support_fraction=1)),
                ("rarefaction_curves", dict(support_depth=3, richness=1, hill_shannon=1,
                                                inverse_simpson=1, is_target="TRUE", is_endpoint="TRUE"))):
            contract = dictionary[table]
            for locus in ("IGK", "IGL", "OTHER"):
                with self.subTest(table=table, locus=locus):
                    row = dict(scope_id="scope", library_id="0040", donor_id="donor",
                               tissue_id="tissue", locus=locus, group_unit="clonotype", **measurements)
                    with gzip.open(self.path, "wt", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=contract.columns, delimiter="\t", lineterminator="\n")
                        writer.writeheader()
                        writer.writerow(row)
                    errors = check_values(self.path, contract, contract.columns)
                    if locus == "OTHER":
                        self.assertIn("vocabulary", "\n".join(errors))
                    else:
                        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
