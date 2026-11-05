"""The Perl assembly recorder writes the same producer contract as Python.

Run: python test/run_tests.py --component trust4 --level integration
"""

import copy
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from test.support.assembly_record import build_assembly_provenance
from workflow.lib.provenance import load_assembly_provenance


SCRIPT = Path(__file__).resolve().parents[4] / "workflow/scripts/record_assembly_provenance.pl"


class AssemblyRecorderTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.launcher = self.bin / "run-trust4"
        self.launcher.write_text('#!/usr/bin/env perl\nmy $version = "v9.8.7-fixture" ;\n')
        self.launcher.chmod(0o755)
        self.environment = {**os.environ, "PATH": str(self.bin) + os.pathsep + os.environ["PATH"]}
        self.perl = shutil.which("perl")
        self.reference = self.root / "reference.fa"
        self.reference.write_bytes(b">ref\nACGT\n")
        self.outputs = {name: self.root / ("0460_" + name + ".txt") for name in ("airr", "cdr3", "final_out")}
        for name, path in self.outputs.items():
            path.write_bytes((name + "\nNA\n").encode())
        self.sidecar = self.root / "producer/0460_producer.json"
        self.context = {"library_id": "0460", "source": {"revision": "a" * 40, "dirty": True,
                                                        "state_sha256": "b" * 64},
                        "image": {"format": "oci", "runtime": "singularity", "digest": "sha256:" + "c" * 64,
                                  "configured_reference": "example/trust4@sha256:" + "c" * 64},
                        "settings": {"k_mer_size": 9, "contig_min_cov": 0, "cleanup": 0},
                        "read_pairs": [{"r1": "0460_lane1_R1.fastq.gz", "r2": "0460_lane1_R2.fastq.gz"},
                                       {"r1": "0460_lane2_R1.fastq.gz", "r2": "0460_lane2_R2.fastq.gz"}]}

    def run_recorder(self, context=None):
        return subprocess.run([self.perl, str(SCRIPT), "--context", json.dumps(self.context if context is None else context),
                               "--reference", str(self.reference), "--output", str(self.sidecar),
                               "--airr", str(self.outputs["airr"]), "--cdr3", str(self.outputs["cdr3"]),
                               "--final", str(self.outputs["final_out"])], env=self.environment,
                              capture_output=True, text=True)


    def test_record_matches_python_contract_and_preserves_string_ids(self):
        result = self.run_recorder()
        self.assertEqual(result.returncode, 0, result.stderr)
        actual = json.loads(self.sidecar.read_text())
        software = [{"name": "TRUST4", "version": "v9.8.7-fixture", "status": "measured",
                     "method": "installed run-trust4 version declaration"}]
        expected = build_assembly_provenance(**self.context, reference=self.reference,
                                             outputs=self.outputs, software=software)
        self.assertEqual(actual, expected)
        self.assertEqual(actual["library_id"], "0460")
        self.assertIs(actual["producer"]["source"]["dirty"], True)
        self.assertEqual(load_assembly_provenance("0460", self.sidecar, self.outputs)["status"], "recorded")
        self.assertEqual(list(self.sidecar.parent.glob(".assembly-provenance-*")), [])

    def test_literal_na_library_id_is_never_missing(self):
        self.context["library_id"] = "NA"
        result = self.run_recorder()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(self.sidecar.read_text())["library_id"], "NA")

    def test_altered_output_is_rejected_by_python_loader(self):
        self.assertEqual(self.run_recorder().returncode, 0)
        self.outputs["cdr3"].write_text("changed\n")
        with self.assertRaisesRegex(ValueError, "does not match"):
            load_assembly_provenance("0460", self.sidecar, self.outputs)

    def test_context_must_have_text_ids_and_explicit_mate_relationships(self):
        numeric = copy.deepcopy(self.context)
        numeric["library_id"] = 460
        bad_pair = copy.deepcopy(self.context)
        del bad_pair["read_pairs"][0]["r2"]
        extra = {**self.context, "another_scope": {"private": "label"}}
        for context in (numeric, bad_pair, extra):
            with self.subTest(context=context):
                result = self.run_recorder(context)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.sidecar.exists())


if __name__ == "__main__":
    unittest.main()
