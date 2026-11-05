"""The public fixture validator detects specific corrupted evidence.

Run: python test/run_tests.py --component inputs --level unit
"""

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from test.system.check_inputs import ROOT, backbones, check_pair


class SyntheticFixtureTest(unittest.TestCase):
    def check_fixture(self, root: Path) -> subprocess.CompletedProcess:
        return subprocess.run([sys.executable, str(Path(__file__).resolve().parents[3] / "system/check_inputs.py"), "--root", str(root)],
                              text=True, capture_output=True)

    def test_frozen_fixture(self):
        result = self.check_fixture(ROOT)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_wrong_mate_name(self):
        fragment = dict(fragment_id="fragment", start=0, end=250)
        with self.assertRaisesRegex(AssertionError, "R2 name differs"):
            check_pair(("fragment/1", "A" * 150, "I" * 150),
                       ("wrong/2", "T" * 150, "I" * 150), fragment, "A" * 250)

    def test_short_vote_array(self):
        source = next((ROOT / "post_trust4").glob("libraries/*/trust4/*_final.out"))
        lines = source.read_text().splitlines()
        lines[2] = " ".join(lines[2].split()[:-1])
        with self.assertRaisesRegex(AssertionError, "Vote array length differs"):
            backbones(lines)

    def test_wrong_mutation_coordinate(self):
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "fixture"
            shutil.copytree(ROOT, copied)
            path = copied / "construction/transcripts.json"
            records = json.loads(path.read_text())
            altered = next(row for row in records if row["edits"])
            altered["edits"][0]["position"] += 1
            path.write_text(json.dumps(records))
            result = self.check_fixture(copied)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Mutation coordinates", result.stderr)

    def test_broken_native_backbone_link(self):
        with tempfile.TemporaryDirectory() as temporary:
            copied = Path(temporary) / "fixture"
            shutil.copytree(ROOT, copied)
            assembly = copied / "post_trust4/libraries/SYN_D1_T1/trust4"
            for name in ("SYN_D1_T1_airr.tsv", "SYN_D1_T1_cdr3.out"):
                path = assembly / name
                path.write_text(path.read_text().replace("assemble0", "assemble_missing"))
            result = self.check_fixture(copied)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Native backbone link", result.stderr)


if __name__ == "__main__":
    unittest.main()
