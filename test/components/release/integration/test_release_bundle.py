"""Keep complete scope observations through publication and validation retry.

Run: python test/run_tests.py --component release --level integration
"""

import csv
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from test.components.release.integration.check_release import check_folder
from workflow.lib.release_bundle import _readme_context, build_release_bundle, release_path, table_for


ROOT = Path(__file__).resolve().parents[4]


from test.support.release_context import write_release_context


class ReleaseBundleTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name).resolve()
        self.scope = self.root / "scope"
        self.bundle = self.scope / "release"
        for name in ("manifest", "metrics", "report/assets"):
            (self.scope / name).mkdir(parents=True)
        (self.scope / "manifest/scope_manifest.tsv").write_text(
            "scope_id\tlibrary_id\tdonor_id\ttissue_id\nproject\t0460\t001\tT1\nproject\tNA\t001\tT1\n")
        self.headline = self.scope / "metrics/headline.tsv"
        self.headline.write_text("scope_id\tlibrary_id\tdonor_id\ttissue_id\tmetric\n"
                                 "project\t0460\t001\tT1\t0\nproject\tNA\t001\tT1\t\n")
        (self.scope / "report/index.html").write_text('<html><script src="assets/data.js"></script></html>\n')
        (self.scope / "report/assets/data.js").write_text('window.scope = "project";\n')
        (self.scope / "report/.snakemake_timestamp").touch()
        (self.scope / "private.tsv").write_text("library_id\tsecret\nunpublished\tprivate\n")
        self.dictionary = self.root / "dictionary.csv"
        self.dictionary.write_text("table,column,description,type,nullable,unit,source,mapping_note,implementation\n"
                                   "scope_manifest,library_id,Library,string,FALSE,,,,Snakefile\n"
                                   "headline,metric,Measurement,number,TRUE,,,,workflow/lib/headline.py\n"
                                   "reannotated_AIRR,library_id,Library,string,FALSE,,,,\n"
                                   "release_files,path,Delivered file,string,FALSE,,,,workflow/lib/release_bundle.py\n"
                                   "unreleased,private,Unpublished field,string,TRUE,,,,\n")
        self.airr = {}
        self.tables = {"manifest/scope_manifest.tsv": "scope_manifest",
                       "metrics/headline.tsv": "headline"}
        for library in ("0460", "NA"):
            source = self.root / (library + ".tsv.gz")
            with gzip.open(source, "wt") as handle:
                handle.write(f"library_id\tsequence_id\tconsensus_count\n{library}\tseq\t7\n")
            relative = "airr/" + library + ".tsv.gz"
            self.airr[library] = {source: relative}
            self.tables[relative] = "reannotated_AIRR"
        self.provenance, _, config_path, threshold = write_release_context(self.scope)
        self.arguments = {
            "scope_dir": self.scope, "bundle": self.bundle,
            "scope_files": ["manifest/scope_manifest.tsv", "metrics/headline.tsv"],
            "scope_directories": ["report"], "airr_files": self.airr, "dictionary": self.dictionary,
            "table_names": self.tables,
            "provenance": self.provenance,
            "effective_config": config_path, "calibration": threshold,
        }

    def test_readme_describes_scope_assignment_only_when_released(self):
        # The README names manifest/scope_assignment.tsv only in a release that ships it.
        staged = self.root / "staged"
        (staged / "manifest").mkdir(parents=True)
        manifest = [{"donor_id": "001"}]
        self.assertEqual(_readme_context(staged, manifest, None)["scope_assignment"], "")
        (staged / "manifest/scope_assignment.tsv").write_text("library_id\tproject\n")
        self.assertIn("`manifest/scope_assignment.tsv`",
                      _readme_context(staged, manifest, None)["scope_assignment"])
        template = ROOT / "workflow/release/README_template.md"
        build_release_bundle(**self.arguments, readme_template=template)
        readme = (self.bundle / "README.md").read_text()
        self.assertNotIn("scope_assignment", readme)
        self.assertIn("The scope manifest defines this release's population.", readme)

    def files(self, directory):
        return {str(path.relative_to(directory)): path.read_bytes() for path in directory.rglob("*") if path.is_file()}

    def test_complete_bundle_dictionary_metadata_and_standalone_folder(self):
        result = build_release_bundle(**self.arguments)
        expected = {"manifest/scope_manifest.tsv", "metrics/headline.tsv", "report/index.html",
                    "report/assets/data.js", "airr/0460.tsv.gz", "airr/NA.tsv.gz",
                    "manifest/release_files.tsv", "README.md", "data_dictionary.csv", "provenance.json",
                    "effective_config.yaml"}
        self.assertEqual(set(self.files(self.bundle)), expected)
        self.assertEqual(result["files"], len(expected))
        self.assertEqual(result["libraries"], 2)
        self.assertFalse((self.bundle / "metadata_descriptors.json").exists())
        with (self.bundle / "data_dictionary.csv").open() as handle:
            reader = csv.DictReader(handle)
            dictionary = list(reader)
        self.assertNotIn("implementation", reader.fieldnames)
        self.assertEqual({row["table"] for row in dictionary},
                         {"scope_manifest", "headline", "reannotated_AIRR", "release_files"})
        self.assertEqual(result["tables"], 4)
        self.assertIn("defines the 4 shipped scientific tables", (self.bundle / "README.md").read_text())
        self.assertFalse(any(row["column"].startswith("meta_") for row in dictionary))
        self.assertIn("Assembly producer history is recorded for 0 libraries and unavailable for 2.", (self.bundle / "README.md").read_text())
        # The record names the threshold files where the release places them.
        released = json.loads((self.bundle / "provenance.json").read_text())
        self.assertEqual(released["files"], {**self.provenance["files"],
            "threshold": "igh_clone_threshold/scope_threshold.tsv",
            "donor_thresholds": "igh_clone_threshold/donor_fits.tsv",
            "distance_diagnostics": "igh_clone_threshold/donor_distance_histograms.tsv"})
        self.assertNotIn("private", "\n".join(path for path in self.files(self.bundle)))
        # No entry is a symlink, so copying the folder somewhere unrelated
        # delivers the same bytes and reaches nothing outside itself.
        for path in self.bundle.rglob("*"):
            self.assertFalse(path.is_symlink(), path)
        published = self.files(self.bundle)
        moved = self.root / "copied-away"
        shutil.copytree(self.bundle, moved, symlinks=True)
        shutil.rmtree(self.scope)
        self.assertEqual(self.files(moved), published)
        self.assertIn('assets/data.js', (moved / "report/index.html").read_text())
        self.assertTrue((moved / "report/assets/data.js").is_file())


    def test_file_manifest_lists_every_other_file_and_the_receipt_holds_its_hash(self):
        def validate(staged):
            (staged / "airr/validation.json").write_text('{"actual_staged_files": true}\n')

        result = build_release_bundle(**self.arguments, validate_staged=validate)
        published = self.files(self.bundle)
        manifest = published.pop("manifest/release_files.tsv")
        self.assertEqual(result["release_files_sha256"], hashlib.sha256(manifest).hexdigest())
        self.assertEqual(result["files"], len(published) + 1)
        lines = manifest.decode().splitlines()
        self.assertEqual(lines[0], "path\tn_bytes\tsha256")
        rows = [line.split("\t") for line in lines[1:]]
        self.assertEqual([row[0] for row in rows], sorted(published))
        self.assertIn("airr/validation.json", [row[0] for row in rows])
        self.assertIn("report/assets/data.js", [row[0] for row in rows])
        for path, n_bytes, sha256 in rows:
            self.assertEqual(int(n_bytes), len(published[path]), path)
            self.assertEqual(sha256, hashlib.sha256(published[path]).hexdigest(), path)
        self.assertEqual(table_for(self.bundle / "manifest/release_files.tsv"), "release_files")
        with (ROOT / "docs/data_dictionary.csv").open(newline="") as handle:
            declared = [row["column"] for row in csv.DictReader(handle) if row["table"] == "release_files"]
        self.assertEqual(lines[0].split("\t"), declared)
        (self.bundle / ".complete.json").write_text(json.dumps(result))
        self.assertEqual(check_folder(self.bundle), result["files"])

    def test_release_check_rejects_a_changed_folder_or_receipt(self):
        receipt = build_release_bundle(**self.arguments)
        target = self.bundle / "report/index.html"
        original = target.read_bytes()

        def replace(content):
            # The staged file can be a hard link to its source, so replace it.
            target.unlink()
            target.write_bytes(content)

        (self.bundle / ".complete.json").write_text(json.dumps(receipt))
        check_folder(self.bundle)
        for label, change, restore, message in (
                ("hash", lambda: replace(original[::-1]), lambda: replace(original), "hash differs"),
                ("size", lambda: replace(original + b"\n"), lambda: replace(original), "size differs"),
                ("missing", target.unlink, lambda: target.write_bytes(original), "lists absent files"),
                ("extra", (self.bundle / "extra.txt").touch, (self.bundle / "extra.txt").unlink, "omits delivered files"),
                ("receipt", lambda: (self.bundle / ".complete.json").write_text(
                    json.dumps({**receipt, "release_files_sha256": "0" * 64})),
                 lambda: (self.bundle / ".complete.json").write_text(json.dumps(receipt)), "receipt hash"),
                ("count", lambda: (self.bundle / ".complete.json").write_text(json.dumps({**receipt, "files": 1})),
                 lambda: (self.bundle / ".complete.json").write_text(json.dumps(receipt)), "file count")):
            with self.subTest(label):
                change()
                with self.assertRaisesRegex(AssertionError, message):
                    check_folder(self.bundle)
                restore()
                check_folder(self.bundle)

    def test_file_manifest_leaves_source_tables_unchanged(self):
        before = (self.headline.read_bytes(), self.headline.stat().st_mtime_ns)
        build_release_bundle(**self.arguments)
        self.assertEqual((self.headline.read_bytes(), self.headline.stat().st_mtime_ns), before)
        self.assertEqual((self.bundle / "metrics/headline.tsv").read_bytes(), before[0])
        # A source cannot take the manifest's place in the staged folder.
        (self.scope / "manifest/release_files.tsv").write_text("path\tn_bytes\tsha256\n")
        arguments = {**self.arguments, "scope_files": [*self.arguments["scope_files"], "manifest/release_files.tsv"]}
        with self.assertRaisesRegex(ValueError, "reserved release destination"):
            build_release_bundle(**arguments)

    def test_publisher_generated_validation_is_published_and_failure_preserves_release(self):
        # This callback tests publication and rollback with a constructed marker.
        # AIRR schema and file validation are exercised in test_airr_exports.py.
        def validate(staged):
            # The publisher generates this file from the staged copies, and the
            # previous release must not be visible to the validation step.
            self.assertFalse((staged / "airr/validation.json").exists())
            self.assertTrue((staged / "airr/0460.tsv.gz").is_file())
            (staged / "airr/validation.json").write_text('{"actual_staged_files": true}\n')

        arguments = {**self.arguments, "validate_staged": validate}
        result = build_release_bundle(**arguments)
        self.assertEqual(result["files"], len(self.files(self.bundle)))
        self.assertEqual((self.bundle / "airr/validation.json").read_text(), '{"actual_staged_files": true}\n')
        before = self.files(self.bundle)

        def fail(staged):
            validate(staged)
            raise ValueError("actual staged validation failed")

        with self.assertRaisesRegex(ValueError, "actual staged"):
            build_release_bundle(**{**arguments, "validate_staged": fail})
        self.assertEqual(self.files(self.bundle), before)
        build_release_bundle(**arguments)
        self.assertEqual(self.files(self.bundle), before)


    def test_scope_privacy_rejects_other_libraries_and_wrong_airr_source(self):
        for defect in ("missing", "extra", "wrong_source", "wrong_metric", "wrong_provenance"):
            with self.subTest(defect=defect):
                arguments = dict(self.arguments)
                if defect == "missing":
                    arguments["airr_files"] = {"0460": self.airr["0460"]}
                elif defect == "extra":
                    arguments["airr_files"] = {**self.airr, "unpublished": self.airr["0460"]}
                elif defect == "wrong_source":
                    arguments["airr_files"] = {**self.airr, "0460": {next(iter(self.airr["NA"])): "airr/wrong.tsv.gz"}}
                elif defect == "wrong_metric":
                    arguments["scope_files"] = [*arguments["scope_files"], "private.tsv"]
                    arguments["table_names"] = {**self.tables, "private.tsv": "headline"}
                else:
                    arguments["provenance"] = {**self.provenance, "scope_id": "other"}
                with self.assertRaisesRegex(ValueError, "in-scope|another library|out-of-scope|another scope"):
                    build_release_bundle(**arguments)
        self.assertFalse(self.bundle.exists())

    def test_threshold_files_move_to_their_release_folder(self):
        for source, released, table in (
                ("clones/project_clone_threshold.tsv", "igh_clone_threshold/scope_threshold.tsv",
                 "scope_threshold"),
                ("clones/project_donor_clone_thresholds.tsv", "igh_clone_threshold/donor_fits.tsv",
                 "donor_fits"),
                ("clones/project_clone_distance_diagnostics.tsv",
                 "igh_clone_threshold/donor_distance_histograms.tsv", "donor_distance_histograms")):
            with self.subTest(source=source):
                self.assertEqual(release_path(source), released)
                self.assertEqual(table_for(Path("scopes/project") / source), table)
                self.assertEqual(table_for(Path("scopes/project/release") / released), table)
        # Donor clone membership stays in the working tree under its own name.
        self.assertEqual(release_path("clones/D1_IGH_clone_membership.tsv"),
                         "clones/D1_IGH_clone_membership.tsv")


if __name__ == "__main__":
    unittest.main()
