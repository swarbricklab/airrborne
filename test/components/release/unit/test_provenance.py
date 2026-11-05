"""Scientific settings and reference evidence describe the measurements as run.

Run: python test/run_tests.py --component release --level unit
"""

import copy
import json
from pathlib import Path
import runpy
from types import SimpleNamespace
import unittest

import yaml

from test.support.provenance_cases import ProvenanceCase, image
from workflow.lib.provenance import (
    capture_source, effective_scope_config, file_identity, library_processing_identity,
    write_scope_provenance,
)
from workflow.lib.release_methods import FIXED_METHODS, methods_paragraphs

ROOT = Path(__file__).resolve().parents[4]
# The source file holding each fixed method's call.
FIXED_METHOD_SOURCES = {
    "shazam::distToNearest": "workflow/lib/donor_distance_threshold.R",
    "shazam::findThreshold": "workflow/lib/donor_distance_threshold.R",
    "scoper::hierarchicalClones": "workflow/scripts/assign_scope_donor_clones.R",
    "shazam::calcObservedMutations": "workflow/lib/shm_regions.R",
    "shazam::collapseClones": "workflow/lib/baseline_selection.R",
    "shazam::calcBaseline": "workflow/lib/baseline_selection.R",
}


def r_call(source: str, function: str) -> str:
    """Return one multi-line R call, from its function name to its closing parenthesis."""
    start = source.index(function.partition("::")[2] + "(\n")
    depth = 0
    for position in range(start, len(source)):
        depth += {"(": 1, ")": -1}.get(source[position], 0)
        if source[position] == ")" and depth == 0:
            return source[start:position + 1]
    raise ValueError(f"{function} has no closing parenthesis")


def r_argument(name: str, value: str | bool | None) -> str:
    """Write a recorded argument as the R call writes it."""
    if value is None:
        return f"{name} = NULL"
    if isinstance(value, bool):
        return f"{name} = {str(value).upper()}"
    return f"{name} = shazam::{value}" if name in ("targetingModel", "regionDefinition") else f'{name} = "{value}"'


class ProvenanceTests(ProvenanceCase, unittest.TestCase):
    def test_reference_metadata_is_identified_without_copying_shared_audit_details(self):
        metadata = self.root / "reference.metadata.tsv"
        metadata.write_text("allele\tcompleteness\nIGLJ2*01\tcomplete\n")
        receipt = copy.deepcopy(self.receipt)
        receipt["reference_metadata"] = {**file_identity(metadata), "rebuild_matches": True}
        receipt["identity_classes"] = [{"identity_id": "JSEQ_constructed", "label": "IGLJ2/IGLJ3",
                                        "region": "J-REGION", "alleles": ["IGLJ2*01", "IGLJ3*01"]}]
        args = {"references": {"assembly_reference": self.reference, "reference_metadata": metadata},
                "reference_receipt": receipt}
        record = self.scope(**args)
        self.assertEqual(record["references"]["reference_metadata"]["sha256"], file_identity(metadata)["sha256"])
        audit = record["references"]["assembly_reference"]["audit"]
        self.assertNotIn("identity_classes", audit)
        self.assertNotIn("reference_metadata", audit)
        self.assertEqual(len(receipt["identity_classes"]), 1)
        metadata.write_text("altered\n")
        with self.assertRaisesRegex(ValueError, "Reference metadata differs"):
            self.scope(**args)


    def test_source_capture_and_packaging_do_not_touch_scientific_inputs(self):
        scientific = [self.repo / "workflow/lib" / "analysis.R", self.reference, self.output]
        before = [(path.read_bytes(), path.stat().st_mtime_ns) for path in scientific]
        first = self.scope()
        write_scope_provenance(first, effective_scope_config(self.config, "project_b"), self.root / "release/provenance.json", self.root / "release/effective_config.yaml")
        self.assertEqual([(path.read_bytes(), path.stat().st_mtime_ns) for path in scientific], before)
        (self.repo / "workflow/lib" / "new.R").write_text("new_function <- function() 1\n")
        second = self.scope()
        self.assertNotIn("state_sha256", first["source"])
        self.assertEqual(len(second["source"]["state_sha256"]), 64)
        self.assertEqual(first["files"]["effective_config"], second["files"]["effective_config"])
        self.assertEqual([(path.read_bytes(), path.stat().st_mtime_ns) for path in scientific], before)


    def test_library_processing_identity_ignores_report_metadata_and_repository_revision(self):
        code = {"annotation": self.repo / "workflow/lib" / "analysis.R"}
        settings = {"airr_qc": self.config["params"]["airr_qc"]}
        references = {"imgt": self.reference}
        before = library_processing_identity(code, settings, image(), references)
        (self.repo / "workflow").mkdir(exist_ok=True)
        (self.repo / "workflow" / "report.html").write_text("new layout\n")
        self.config["scope_metadata"]["project_b"]["tissue"]["file"] = "different-labels.csv"
        self.config["params"]["lineage"]["clonal_threshold"] = 0.2
        self.assertTrue(capture_source(self.repo)["dirty"])
        self.assertEqual(library_processing_identity(code, settings, image(), references), before)
        code["annotation"].write_text("value <- 2\n")
        self.assertNotEqual(library_processing_identity(code, settings, image(), references), before)
        code["annotation"].write_text("value <- 1\n")
        settings["airr_qc"]["max_alignment_n_fraction"] = 0.2
        self.assertNotEqual(library_processing_identity(code, settings, image(), references), before)


    def test_effective_config_excludes_report_metadata_and_preserves_scientific_values(self):
        before = copy.deepcopy(self.config)
        effective = effective_scope_config(self.config, "project_b")
        self.assertEqual(effective["params"], self.config["params"])
        self.assertNotIn("scope_metadata", effective)
        self.assertNotIn("library_metadata", effective)
        self.assertEqual(effective["seed"], 1337)
        text = json.dumps(effective)
        for value in ("private_project", "private_label", "do-not-release", "all-libraries", "global-results"):
            self.assertNotIn(value, text)
        effective["params"]["lineage"]["max_n"] = 1
        self.assertEqual(self.config, before)


    def test_fixed_methods_match_the_calls_that_make_them(self):
        entries = [*(FIXED_METHODS[key] for key in ("clone_distance", "clone_threshold_fit", "clone_assignment", "shm_counting")),
                   *FIXED_METHODS["selection"].values()]
        self.assertEqual({entry["function"] for entry in entries}, set(FIXED_METHOD_SOURCES))
        for entry in entries:
            call = r_call((ROOT / FIXED_METHOD_SOURCES[entry["function"]]).read_text(), entry["function"])
            for name, value in entry.items():
                if name != "function":
                    with self.subTest(function=entry["function"], argument=name):
                        self.assertIn(r_argument(name, value), call)

    def test_scope_script_records_fixed_methods_beside_the_configured_settings(self):
        results = self.root / "results"
        for library in ("0460", "NA"):
            directory = results / "libraries" / library / "trust4"
            directory.mkdir(parents=True)
            for suffix in ("_airr.tsv", "_cdr3.out", "_final.out"):
                (directory / (library + suffix)).write_text(library + suffix + "\n")
        manifest = self.root / "manifest.tsv"
        manifest.write_text("library_id\ttissue_id\tdonor_id\n0460\tNA\t001\nNA\tNA\t001\n")
        versions = {role: self.root / f"{role}_versions.json" for role in self.versions}
        for role, path in versions.items():
            path.write_text(json.dumps(self.versions[role]))
        metadata, auxiliary = self.root / "reference.metadata.tsv", self.root / "human_gl.aux"
        metadata.write_text("allele\n")
        auxiliary.write_text("allele\n#\n")
        receipt = self.root / "receipt.json"
        receipt.write_text(json.dumps({**self.receipt, "reference_metadata": file_identity(metadata)}))
        output = SimpleNamespace(provenance=str(self.root / "scope/provenance.json"),
                                 config=str(self.root / "scope/effective_config.yaml"))
        context = SimpleNamespace(
            input=SimpleNamespace(module=str(ROOT / "workflow/lib/provenance.py"), manifest=str(manifest),
                                  analysis_versions=str(versions["analysis"]), assembly_versions=str(versions["assembly"]),
                                  reference=str(self.reference), reference_metadata=str(metadata),
                                  auxiliary=str(auxiliary), reference_receipt=str(receipt)),
            output=output, wildcards=SimpleNamespace(scope="project_b"), log=[str(self.root / "scope.log")],
            params=SimpleNamespace(effective_config=self.config, source=capture_source(self.repo), repo=str(self.repo),
                                   source_evidence=str(self.root / "evidence"), results_dir=str(results)))
        runpy.run_path(str(ROOT / "workflow/scripts/build_scope_provenance.py"), init_globals={"snakemake": context})
        written = yaml.safe_load(Path(output.config).read_text())
        self.assertEqual(written["fixed_methods"], FIXED_METHODS)
        expected = effective_scope_config(self.config, "project_b")
        self.assertEqual({key: value for key, value in written.items() if key != "fixed_methods"}, expected)
        self.assertNotIn("fixed_methods", json.loads(Path(output.provenance).read_text()))

    def test_scope_output_is_compact_and_methods_use_authoritative_inputs(self):
        record = self.scope()
        self.assertEqual(set(record), {"schema_version", "scope_id", "source", "images",
                         "software", "references", "files", "assembly_producers"})
        self.assertEqual(record["references"]["assembly_reference"], {
            **{key: file_identity(self.reference)[key] for key in ("name", "sha256")},
            "imgt_snapshot": {"source": "https://www.imgt.org/genedb", "date": "2026.06.07"},
            "audit": {"agreement": {"reference_rebuild": True, "auxiliary_rebuild": True, "igblast_databases": True},
                      "auxiliary_data": {
                          "identical_pairs_checked": 3, "sibling_pairs_checked": 34,
                          "supplemented": [{"allele": "IGHJ5*04", "rule": "same_gene_equal_length",
                                            "copied_from": ["IGHJ5*01", "IGHJ5*02", "IGHJ5*03"]},
                                           {"allele": "IGLJ2A*01", "rule": "identical_sequence", "copied_from": ["IGLJ2*01"]}],
                          "j_alleles_without_entry": [{"allele": "TRAJ55*01", "functionality": "P"}]}}})
        config = effective_scope_config(self.config, "project_b")
        config["params"]["lineage"].update(fit_edge=0.9, fit_specificity=0.99, fit_subsample=1000)
        config["params"].update(rarefaction={"targets": {"IGH": 1000}},
                                diversity={"min_consensus_count": {"IGH": 0}})
        libraries = [dict(row, scope_id="project_b") for row in self.libraries]
        threshold = {"scope_id": "project_b", "clone_threshold": 0.1, "threshold_source": "fallback"}
        methods = " ".join(methods_paragraphs(record, config, libraries, threshold))
        for text in ("threshold 0.1 (fallback)", "40 informative", "IgBLAST 1.0-fixture", "unavailable for 2",
                     "human.fa reference was curated from the IMGT/GENE-DB download dated 2026.06.07",
                     "entries added for IGHJ5*04 (from IGHJ5*01, IGHJ5*02, IGHJ5*03), IGLJ2A*01 (identical to IGLJ2*01)",
                     "reproduced all 3 identical-sequence pairs and all 34 same-gene pairs of equal length in the file",
                     "No auxiliary entry exists for J alleles TRAJ55*01 (P). A record whose first-listed top J hit"):
            self.assertIn(text, methods)
        self.assertNotIn("differs from the annotation reference", methods)
        self.assertNotIn("TRUST4 assembly was configured", methods)
        config["params"]["airr_qc"]["min_informative_c_alignment_positions"] = 41
        self.assertIn("41 informative", " ".join(methods_paragraphs(record, config, libraries, threshold)))
        payload = json.dumps(record)
        for removed in ("command", "method", "analysis_details", "effective_config\\\": {", "size_bytes"):
            self.assertNotIn(removed, payload)


    def test_receipt_for_another_reference_or_image_fails(self):
        receipt = copy.deepcopy(self.receipt)
        receipt["reference"]["sha256"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "reference receipt"):
            self.scope(reference_receipt=receipt)
        receipt = copy.deepcopy(self.receipt)
        receipt["image"] = image("c")
        with self.assertRaisesRegex(ValueError, "reference receipt"):
            self.scope(reference_receipt=receipt)
        receipt = copy.deepcopy(self.receipt)
        receipt["agreement"]["igblast_databases"] = False
        with self.assertRaisesRegex(ValueError, "must record agreement"):
            self.scope(reference_receipt=receipt)


    def test_different_image_tool_record_and_duplicate_scope_libraries_fail(self):
        versions = copy.deepcopy(self.versions)
        versions["analysis"]["image"] = image("c")
        with self.assertRaisesRegex(ValueError, "do not belong"):
            self.scope(tool_versions=versions)
        self.libraries.append(self.libraries[0])
        with self.assertRaisesRegex(ValueError, "one row per string library ID"):
            self.scope()


if __name__ == "__main__":
    unittest.main()
