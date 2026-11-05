"""Pair a library's FASTQs by Illumina and SRA names before assembly.

Run: python test/run_tests.py --component inputs --level unit
"""

import os
from pathlib import Path
import re
import tempfile
from types import SimpleNamespace
import unittest

import pandas as pd


ASSEMBLY_RULES = Path(__file__).resolve().parents[4] / "workflow/rules/assembly.smk"


def load_find_all_fastq_files(data_dir: Path):
    """Return find_all_fastq_files() from assembly.smk with deps.data set.

    The function sits among Snakemake rules, so its top-level block is read
    from the rule file and run with the names the Snakefile supplies.
    """
    lines = ASSEMBLY_RULES.read_text().splitlines()
    start = lines.index("def find_all_fastq_files(wildcards, samplesheet):")
    end = next(index for index in range(start + 1, len(lines))
               if lines[index] and not lines[index][0].isspace())
    namespace = {"config": {"deps": {"data": str(data_dir)}}, "os": os, "re": re, "Path": Path}
    exec("\n".join(lines[start:end]), namespace)
    return namespace["find_all_fastq_files"]


class FastqDiscoveryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.data = Path(self.temporary.name) / "data"
        self.library = self.data / "library_a"
        self.library.mkdir(parents=True)
        self.find = load_find_all_fastq_files(self.data)
        self.samplesheet = pd.DataFrame({"library_id": ["library_a"], "path": ["library_a"]})

    def write(self, *names):
        for name in names:
            (self.library / name).write_bytes(b"")

    def pairs(self):
        found = self.find(SimpleNamespace(library_id="library_a"), self.samplesheet)
        return [(Path(r1).name, Path(r2).name) for r1, r2 in zip(found["r1"], found["r2"])]

    def test_illumina_lanes_pair_in_name_order(self):
        self.write("library_a_L002_R2_001.fastq.gz", "library_a_L001_R1_001.fastq.gz",
                   "library_a_L002_R1_001.fastq.gz", "library_a_L001_R2_001.fastq.gz")
        self.assertEqual(self.pairs(), [
            ("library_a_L001_R1_001.fastq.gz", "library_a_L001_R2_001.fastq.gz"),
            ("library_a_L002_R1_001.fastq.gz", "library_a_L002_R2_001.fastq.gz")])

    def test_sra_runs_pair_by_trailing_mate_number(self):
        self.write("run_a_1.fastq.gz", "run_a_2.fastq.gz", "run_b_1.fq.gz", "run_b_2.fq.gz")
        self.assertEqual(self.pairs(), [("run_a_1.fastq.gz", "run_a_2.fastq.gz"),
                                        ("run_b_1.fq.gz", "run_b_2.fq.gz")])

    def test_illumina_form_takes_precedence_over_a_trailing_number(self):
        self.write("library_a_R1_1.fastq.gz", "library_a_R2_1.fastq.gz")
        self.assertEqual(self.pairs(), [("library_a_R1_1.fastq.gz", "library_a_R2_1.fastq.gz")])

    def test_one_mate_named_both_ways_is_a_duplicate(self):
        self.write("library_a_R1.fastq.gz", "library_a_1.fastq.gz", "library_a_R2.fastq.gz")
        with self.assertRaisesRegex(ValueError, "Duplicate R1 FASTQ member"):
            self.pairs()

    def test_a_lone_sra_mate_stops_discovery(self):
        self.write("run_a_1.fastq.gz")
        with self.assertRaisesRegex(ValueError, "Mismatched R1/R2 FASTQ pairs"):
            self.pairs()

    def test_other_names_are_ignored(self):
        self.write("library_a_R1_001.fastq.gz", "library_a_R2_001.fastq.gz",
                   "library_a_I1_001.fastq.gz", "library_a_001.fastq.gz",
                   "library_a_3.fastq.gz", "library_a_1.fastq", "notes.txt")
        self.assertEqual(self.pairs(), [("library_a_R1_001.fastq.gz", "library_a_R2_001.fastq.gz")])

    def test_link_names_decide_the_match(self):
        raw = Path(self.temporary.name) / "raw"
        raw.mkdir()
        for mate in (1, 2):
            (raw / f"delivery_mate{mate}.fastq.gz").write_bytes(b"")
            (self.library / f"run_a_{mate}.fastq.gz").symlink_to(raw / f"delivery_mate{mate}.fastq.gz")
        self.assertEqual(self.pairs(), [("run_a_1.fastq.gz", "run_a_2.fastq.gz")])


if __name__ == "__main__":
    unittest.main()
