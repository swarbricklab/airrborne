"""Reject publication when scope libraries used different scientific processing.

Run: python test/run_tests.py --component release --level integration
"""

import json
from pathlib import Path
import runpy
import tempfile
from types import SimpleNamespace
import unittest


class ScopeRepertoireAdapterTests(unittest.TestCase):
    def test_scope_libraries_require_complete_matching_processing_identities(self):
        repository = Path(__file__).resolve().parents[4]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "manifest.tsv"
            manifest.write_text("scope_id\tlibrary_id\tdonor_id\ttissue_id\n"
                                "scope\t0460\tD1\tT1\nscope\tNA\tD1\tT2\n")
            provenance = root / "provenance.json"
            provenance.write_text(json.dumps({"software": {"analysis": {"tools": [
                                                  {"name": "IgBLAST", "version": "1.22.0", "status": "measured"}]},
                                              "configured_assembly": {"tools": [
                                                  {"name": "TRUST4", "version": "packaging-version", "status": "measured"}]}},
                                              "assembly_producers": {"libraries": {"0460": "producer"}, "producers": {
                                                  "producer": {"settings": {"k_mer_size": 11},
                                                               "software": [{"name": "TRUST4", "version": "v1.1.9-r619", "status": "measured"},
                                                                           {"name": "Perl", "version": None, "status": "unknown"}]}}},
                                              "references": {"imgt_reference": {
                                                  "sha256": "a" * 64,
                                                  "imgt_snapshot": {"source": "https://www.imgt.org/genedb",
                                                                    "date": "2026.06.07"}}}}))
            config = root / "effective.yaml"
            config.write_text("params:\n  airr_qc:\n    max_alignment_n_fraction: 0\n")
            output = root / "repertoires.json"
            context = SimpleNamespace(input=SimpleNamespace(
                module=str(repository / "workflow/lib/repertoire_metadata.py"),
                manifest=str(manifest), provenance=str(provenance), effective_config=str(config)),
                output=SimpleNamespace(repertoire=str(output)),
                wildcards=SimpleNamespace(scope="scope"), log=[str(root / "metadata.log")])

            def execute(records):
                context.input.identities = []
                for index, record in enumerate(records):
                    path = root / f"identity-{index}.json"
                    path.write_text(json.dumps(record))
                    context.input.identities.append(str(path))
                runpy.run_path(str(repository / "workflow/scripts/build_scope_repertoires.py"),
                               init_globals={"snakemake": context})

            first = {"library_id": "0460", "processing_identity": "a" * 64}
            second = {"library_id": "NA", "processing_identity": "a" * 64}
            execute([second, first])
            result = json.loads(output.read_text())
            self.assertEqual([row["repertoire_name"] for row in result["Repertoire"]], ["0460", "NA"])
            self.assertEqual(result["scope_id"], "scope")
            with self.assertRaisesRegex(ValueError, "differ from scope"):
                execute([first])
            conflicting = {**first, "processing_identity": "b" * 64}
            with self.assertRaisesRegex(ValueError, "different processing settings"):
                execute([conflicting, second])


if __name__ == "__main__":
    unittest.main()
