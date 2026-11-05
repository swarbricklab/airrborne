"""Scope collection preserves library rows and never requires another scope.

Run: python test/run_tests.py --component scopes --level integration
"""

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path
import tempfile
import unittest

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
from workflow.lib import scope_tables as MODULE


class ScopeCollectionTest(unittest.TestCase):
    def test_collects_only_supplied_members_with_string_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.tsv"
            manifest.write_text("scope_id\tlibrary_id\tdonor_id\ttissue_id\n"
                                "project\t0460\t001\t001-1\nproject\tNA\t002\t002-1\n")
            first = root / "0460.tsv"
            first.write_text("library_id\tvalue\n0460\t0\n0460\t\n")
            empty = root / "NA.tsv"
            empty.write_text("library_id\tvalue\n")
            output = root / "collected.tsv"
            MODULE.collect_scope_metric(manifest, [first, empty], output,
                                        ("scope_id", "tissue_id", "donor_id"))
            self.assertEqual(output.read_text(),
                             "scope_id\tlibrary_id\ttissue_id\tdonor_id\tvalue\n"
                             "project\t0460\t001-1\t001\t0\nproject\t0460\t001-1\t001\t\n")

    def test_sorts_a_scope_arm_into_the_library_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.tsv"
            manifest.write_text("scope_id\tlibrary_id\tdonor_id\ttissue_id\n"
                                "cohort\tL1\tD1\tT1\ncohort\tL2\tD2\tT2\n")
            first = root / "L1.tsv"
            first.write_text("library_id\tlocus\tsupport_depth\nL1\tTRB\t0\nL1\tTRB\t5\n")
            second = root / "L2.tsv"
            second.write_text("library_id\tlocus\tsupport_depth\nL2\tTRB\t0\n")
            scope = root / "IGH.tsv"
            scope.write_text("library_id\tlocus\tsupport_depth\nL2\tIGH\t0\nL1\tIGH\t5\nL1\tIGH\t0\n")
            output = root / "collected.tsv"
            MODULE.collect_scope_metric(manifest, [first, second, scope], output,
                                        ("scope_id", "tissue_id", "donor_id"),
                                        sort_columns=("library_id", "locus"))
            with output.open(newline="") as handle:
                rows = [(row["library_id"], row["locus"], row["support_depth"])
                        for row in csv.DictReader(handle, delimiter="\t")]
            # Each library's IGH rows come first, and a curve keeps its support_depth order.
            self.assertEqual(rows, [("L1", "IGH", "5"), ("L1", "IGH", "0"), ("L1", "TRB", "0"),
                                    ("L1", "TRB", "5"), ("L2", "IGH", "0"), ("L2", "TRB", "0")])

    def test_labels_each_entity_beside_its_column(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.tsv"
            manifest.write_text("scope_id\tlibrary_id\tdonor_id\ttissue_id\ncohort\tL1\tD1\tT1\n")
            reference = root / "metadata.tsv"
            reference.write_text("allele\tidentity_id\tidentity_label\n"
                                 "IGKV1-39*01\tVSEQ_a\tIGKV1-39,IGKV1D-39\n"
                                 "IGKV1D-39*01\tVSEQ_a\tIGKV1-39,IGKV1D-39\n"
                                 "IGKJ1*01\t\t\n")
            source = root / "L1.tsv"
            source.write_text("library_id\tv_gene\tj_gene\tn\nL1\tVSEQ_a\tIGKJ1\t1\nL1\tIGKV3-20\tIGKJ1\t2\n")
            output = root / "collected.tsv"
            MODULE.collect_scope_metric(manifest, [source], output, ("scope_id", "tissue_id", "donor_id"),
                                        label_columns={"v_gene": "v_gene_names", "j_gene": "j_gene_names"},
                                        entity_labels=MODULE.read_entity_labels(reference))
            header, *rows = [line.split("\t") for line in output.read_text().splitlines()]
            self.assertEqual(header, ["scope_id", "library_id", "tissue_id", "donor_id", "v_gene",
                                      "v_gene_names", "j_gene", "j_gene_names", "n"])
            # A composite reads as its member genes; a single gene is its own label.
            self.assertEqual([row[4:8] for row in rows],
                             [["VSEQ_a", "IGKV1-39,IGKV1D-39", "IGKJ1", "IGKJ1"],
                              ["IGKV3-20", "IGKV3-20", "IGKJ1", "IGKJ1"]])

    def test_rejects_an_out_of_scope_observation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.tsv"
            manifest.write_text("library_id\tdonor_id\ttissue_id\nA\tD\tT\n")
            source = root / "source.tsv"
            source.write_text("library_id\tvalue\nB\t1\n")
            output = root / "output.tsv"
            with self.assertRaisesRegex(ValueError, "outside the scope"):
                MODULE.collect_scope_metric(manifest, [source], output, ("tissue_id", "donor_id"))


    @unittest.skipUnless(shutil.which(os.environ.get("RSCRIPT", "Rscript")),
                         "Executed overlap headers require base R.")
    def test_empty_overlap_headers_match_producer(self):
        module = Path(__file__).resolve().parents[4] / "workflow/lib/overlap_metrics.R"
        command = ("source(" + json.dumps(str(module)) + "); "
                   'x <- make_empty_overlap_tables(); '
                   'for(n in names(x)) cat(n,paste(names(x[[n]]),collapse="\\t"),"\\n",sep="\\t")')
        output = subprocess.check_output([os.environ.get("RSCRIPT", "Rscript"), "-e", command], text=True)
        names = []
        for line in output.splitlines():
            name, *columns = line.rstrip("\t").split("\t")
            names.append(name)
            self.assertEqual(("scope_id", *columns), MODULE.OVERLAP_HEADERS[name])
        self.assertEqual(names, list(MODULE.OVERLAP_HEADERS))

    def test_overlap_keeps_supplied_rows_and_scope_first(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first, second, empty, output = [root / p for p in ("a", "b", "empty", "out")]
            header = "donor_id\tlocus\tn_groups_shared\n"
            first.write_text(header + "001\tTRB\t4\n001\tIGH\t0\n")
            second.write_text(header + "NA\tTRB\t9\n")
            empty.write_text(header)
            MODULE.collect_scope_overlap([first, second, empty], output, "S", "donor_overlap_pairs")
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({(row["scope_id"], row["donor_id"], row["locus"], row["n_groups_shared"])
                              for row in rows},
                             {("S", "001", "TRB", "4"), ("S", "001", "IGH", "0"),
                              ("S", "NA", "TRB", "9")})

    def test_overlap_leads_with_its_identifiers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / "shard", root / "out"
            source.write_text("library_id\ttissue_id\tdonor_id\tlocus\nL1\tT1\tD1\tTRB\n")
            MODULE.collect_scope_overlap([source], output, "S", "shared_groups")
            self.assertEqual(output.read_text(), "scope_id\tlibrary_id\ttissue_id\tdonor_id\tlocus\n"
                                                 "S\tL1\tT1\tD1\tTRB\n")

    def test_selection_preserves_source_counts_and_missing_values(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first, second, output = [root / p for p in ("a", "b", "out")]
            header = "library_id\tisotype_group\tvalue\n"
            first.write_text(header + "NA\tall_IGH\t0\nNA\tIgM\t\n")
            second.write_text(header + "0460\tall_IGH\tNA\n0460\tIgM\t1\n")
            MODULE.collect_scope_selection([first, second], output)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({(row["library_id"], row["isotype_group"], row["value"]) for row in rows},
                             {("0460", "all_IGH", "NA"), ("0460", "IgM", "1"),
                              ("NA", "all_IGH", "0"), ("NA", "IgM", "")})

    def test_group_abundance_merges_by_library_then_locus(self):
        header = "library_id\tlocus\tgroup_unit\tgroup_key\tn_members\tconsensus_count\tsupport_fraction\n"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.tsv"
            manifest.write_text("scope_id\tlibrary_id\ttissue_id\tdonor_id\n"
                                "S\t0460\tT1\tD1\nS\tL2\tT2\tD1\nS\tL3\tT3\tD2\n")
            first, second, third = [root / f"{name}.tsv" for name in ("0460", "L2", "L3")]
            first.write_text(header + "0460\tTRA\tclonotype\ta\t1\t2\t1\n"
                                      "0460\tTRB\tclonotype\tz\t1\t9\t0.45\n"
                                      "0460\tTRB\tclonotype\tb\t1\t10\t0.5\n"
                                      "0460\tTRB\tclonotype\tc\t1\t1\t0.05\n")
            second.write_text(header)
            third.write_text(header + "L3\tTRB\tclonotype\td\t1\t4\t1\n")
            scope = root / "IGH.tsv"
            scope.write_text(header + "0460\tIGH\tclone\tc2\t1\t10\t0.5\n"
                                      "0460\tIGH\tclone\tc10\t2\t10\t0.5\n"
                                      "L2\tIGH\tclone\tc1\t2\t3\t1\n")
            output = root / "out.tsv"
            MODULE.collect_scope_group_abundance(manifest, [first, second, third], scope, output, "S")
            header_row, *rows = [line.split("\t") for line in output.read_text().splitlines()]
            self.assertEqual(header_row[:5], ["scope_id", "library_id", "tissue_id", "donor_id", "locus"])
            # Each library leads with its IGH clones, and every locus keeps its
            # producer's row order, which neither the keys nor the counts as text
            # would give.
            self.assertEqual([(row[1], row[4], row[6], row[8]) for row in rows],
                             [("0460", "IGH", "c2", "10"), ("0460", "IGH", "c10", "10"),
                              ("0460", "TRA", "a", "2"), ("0460", "TRB", "z", "9"),
                              ("0460", "TRB", "b", "10"), ("0460", "TRB", "c", "1"),
                              ("L2", "IGH", "c1", "3"), ("L3", "TRB", "d", "4")])
            self.assertEqual({tuple(row[:4]) for row in rows if row[1] == "L3"}, {("S", "L3", "T3", "D2")})
            collected = output.read_text()
            unsorted = root / "unsorted.tsv"
            unsorted.write_text(header + "L2\tIGH\tclone\tc1\t2\t3\t1\n"
                                         "0460\tIGH\tclone\tc2\t1\t10\t0.5\n")
            with self.assertRaisesRegex(ValueError, "must be sorted by library_id, then locus"):
                MODULE.collect_scope_group_abundance(manifest, [first, second, third], unsorted, output, "S")
            repeated = root / "repeated.tsv"
            repeated.write_text(header + "L3\tTRB\tclonotype\te\t1\t1\t1\n")
            with self.assertRaisesRegex(ValueError, "must come from one table"):
                MODULE.collect_scope_group_abundance(manifest, [first, second, third], repeated, output, "S")
            self.assertEqual(output.read_text(), collected)


if __name__ == "__main__":
    unittest.main()
