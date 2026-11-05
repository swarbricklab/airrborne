"""Check native bundle completeness and isolation at both fixture boundaries.

Run: python test/run_tests.py --component inputs --level unit
"""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import yaml

from test.system import stage_inputs as STAGING


class FixtureStagingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.manifest = self.root / "samplesheet.csv"
        self.manifest.write_text("library_id,tissue_id,donor_id,path\n001,T1,D1,reads/001\n002,T2,D1,reads/002\n")
        self.config = self.root / "config.yaml"
        self.config.write_text(yaml.safe_dump({
            "deps": {"samplesheet": str(self.manifest), "data": "test/data/reads",
                     "resources": "test/resources/trust4",
                     "immcantation_docker": "example@sha256:" + "a" * 64},
            "membership_csv": "test/data/membership.csv",
            "scope_metadata": {"project_a": {"tissue": {"file": "test/data/tissue_metadata.csv"}}},
            "outs": {"results_dir": "unused", "log_dir": "unused"}}))
        self.seeds = self.root / "seeds"
        self.seed_files = []
        for library in ("001", "002"):
            for directory, suffix in (("trust4", "airr.tsv"), ("trust4", "cdr3.out"),
                                      ("trust4", "final.out"), ("trust4", "assembly_provenance.json"),
                                      ("metrics", "seq_depth.tsv")):
                path = self.seeds / "libraries" / library / directory / f"{library}_{suffix}"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps({"producer": "constructed test producer", "library_id": library}) + "\n")
                self.seed_files.append(path)

    def prepare(self, name="run", seeds=True):
        work = self.root / name
        with contextlib.redirect_stdout(io.StringIO()):
            return STAGING.prepare(self.config, work, work / "config.yaml",
                                   self.seeds if seeds else None)

    def test_every_native_file_and_producer_record_survive_staging(self):
        unrelated = self.seeds / "libraries/unused/trust4/unused_airr.tsv"
        unrelated.parent.mkdir(parents=True)
        unrelated.write_text("Unlisted library.\n")
        config = self.prepare()
        actual = Path(config["outs"]["results_dir"])
        self.assertEqual(len(list(actual.rglob("*.*"))), 10)
        for source in self.seed_files:
            staged = actual / source.relative_to(self.seeds)
            self.assertEqual(staged.read_bytes(), source.read_bytes())
            self.assertNotEqual(staged.stat().st_ino, source.stat().st_ino)
        self.assertFalse((actual / "libraries/unused").exists())
        self.assertEqual(config["deps"]["data"], str((Path.cwd() / "test/data/reads").resolve()))
        self.assertTrue(Path(config["scope_metadata"]["project_a"]["tissue"]["file"]).is_absolute())

    def test_missing_native_input_rejects_bundle_before_any_copy(self):
        source = self.seed_files[1]
        source.unlink()
        with self.assertRaisesRegex(ValueError, source.name):
            self.prepare()
        self.assertFalse((self.root / "run" / "results").exists())

    def test_full_run_has_fresh_outputs_and_preserves_input_paths(self):
        post = self.prepare(name="post")
        full = self.prepare(name="full", seeds=False)
        self.assertNotEqual(post["outs"], full["outs"])
        self.assertEqual(post["deps"], full["deps"])
        self.assertEqual(list(Path(full["outs"]["results_dir"]).iterdir()), [])
        for source in self.seed_files:
            self.assertTrue(source.exists())

    def test_existing_outputs_and_source_config_cannot_be_overwritten(self):
        self.prepare()
        with self.assertRaisesRegex(ValueError, "already exist"):
            self.prepare()
        with self.assertRaisesRegex(ValueError, "must differ"):
            STAGING.prepare(self.config, self.root / "new", self.config)


if __name__ == "__main__":
    unittest.main()
