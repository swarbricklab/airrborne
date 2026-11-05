"""Native observations retain their actual producer and reference history.

Run: python test/run_tests.py --component trust4 --level unit
"""

import json
import unittest

from test.support.assembly_record import build_assembly_provenance
from test.support.provenance_cases import ProvenanceCase, image
from workflow.lib.provenance import (
    ASSEMBLY_SCHEMA_VERSION, SCOPE_SCHEMA_VERSION, capture_source, effective_scope_config,
    file_identity, load_assembly_provenance,
)
from workflow.lib.release_methods import methods_paragraphs


class AssemblyHistoryTests(ProvenanceCase, unittest.TestCase):
    def test_assembly_record_preserves_paired_filenames_and_exact_outputs_without_hashing_fastqs(self):
        record = self.assembly()
        self.assertEqual(record["library_id"], "0460")
        self.assertEqual(record["inputs"]["read_pairs"][1]["r2"], "lane2_R2.fastq.gz")
        self.assertEqual(record["outputs"]["airr"], file_identity(self.output))
        sidecar = self.root / "assembly.json"
        sidecar.write_text(json.dumps(record))
        loaded = load_assembly_provenance("0460", sidecar, self.outputs)
        self.assertEqual(loaded["status"], "recorded")
        self.assertEqual(loaded["record"], record)
        self.output.write_text("changed output\n")
        with self.assertRaisesRegex(ValueError, "does not match"):
            load_assembly_provenance("0460", sidecar, self.outputs)


    def test_unknown_assembly_history_never_inherits_packaging_source_or_image(self):
        unknown = load_assembly_provenance("0460", self.root / "missing.json", self.outputs)
        self.assertEqual(unknown["status"], "unknown")
        self.assertEqual(unknown["outputs"]["airr"], file_identity(self.output))
        self.assertNotIn("producer", unknown)
        self.assertNotIn("source", unknown)
        self.assertNotIn("image", unknown)
        scoped = self.scope(assembly_records={"0460": unknown})
        self.assertEqual(scoped["assembly_producers"]["unknown"]["libraries"], ["0460", "NA"])
        self.assertEqual(scoped["assembly_producers"]["libraries"], {})
        self.assertNotIn('"outputs":', json.dumps(scoped["assembly_producers"]))


    def test_known_assembly_producer_history_survives_different_packaging_code(self):
        record = self.assembly()
        record["producer"]["settings"] = {"k_mer_size": 11, "contig_min_cov": 2, "clean": 0}
        producer_revision = record["producer"]["source"]["revision"]
        self.git("-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid",
                 "commit", "--allow-empty", "-qm", "new packaging revision")
        entry = {"library_id": "0460", "status": "recorded", "record": record}
        scope = self.scope(assembly_records={"0460": entry, "private_library": {"label": "private_label"}})
        self.assertNotEqual(scope["source"]["revision"], producer_revision)
        identity = scope["assembly_producers"]["libraries"]["0460"]
        self.assertEqual(scope["assembly_producers"]["producers"][identity]["source"]["revision"], producer_revision)
        self.assertEqual(scope["assembly_producers"]["producers"][identity]["settings"],
                         {"k_mer_size": 11, "contig_min_cov": 2, "clean": 0})
        self.assertNotEqual(scope["assembly_producers"]["producers"][identity]["settings"], self.config["params"]["trust4"])
        scope["assembly_producers"]["producers"][identity]["settings"]["k_mer_size"] = 13
        self.assertEqual(record["producer"]["settings"]["k_mer_size"], 11)
        self.assertEqual(scope["assembly_producers"]["unknown"]["libraries"], ["NA"])
        self.assertEqual(record["schema_version"], ASSEMBLY_SCHEMA_VERSION)
        self.assertEqual(scope["schema_version"], SCOPE_SCHEMA_VERSION)
        self.assertNotIn("private_library", json.dumps(scope))
        self.assertNotIn("private_label", json.dumps(scope))


    def test_assembly_reference_is_compared_with_the_annotation_reference(self):
        record = self.assembly()
        same = {"library_id": "0460", "status": "recorded", "record": record}
        scoped = self.scope(assembly_records={"0460": same})
        producer = scoped["assembly_producers"]["producers"][scoped["assembly_producers"]["libraries"]["0460"]]
        self.assertTrue(producer["reference"]["same_as_annotation_reference"])
        other = self.root / "older.fa"
        other.write_text(">gene\nACGA\n")
        older = build_assembly_provenance("NA", source=capture_source(self.repo), image=image("b"),
                                          settings={"k_mer_size": 11, "contig_min_cov": 2}, reference=other,
                                          read_pairs=[{"r1": "lane1_R1.fastq.gz", "r2": "lane1_R2.fastq.gz"}],
                                          outputs=self.outputs)
        scoped = self.scope(assembly_records={"0460": same,
                                              "NA": {"library_id": "NA", "status": "recorded", "record": older}})
        history = scoped["assembly_producers"]
        self.assertFalse(history["producers"][history["libraries"]["NA"]]["reference"]["same_as_annotation_reference"])
        self.assertTrue(history["producers"][history["libraries"]["0460"]]["reference"]["same_as_annotation_reference"])
        config = effective_scope_config(self.config, "project_b")
        config["params"]["lineage"].update(fit_edge=0.9, fit_specificity=0.99, fit_subsample=1000)
        config["params"].update(rarefaction={"targets": {"IGH": 1000}},
                                diversity={"min_consensus_count": {"IGH": 0}})
        libraries = [dict(row, scope_id="project_b") for row in self.libraries]
        threshold = {"scope_id": "project_b", "clone_threshold": 0.1, "threshold_source": "fallback"}
        methods = " ".join(methods_paragraphs(scoped, config, libraries, threshold))
        self.assertIn("Assembly of 1 library read a reference file that differs from the annotation "
                      "reference: NA.", methods)
        self.assertIn("reference older.fa (SHA256 " + file_identity(other)["sha256"] +
                      ") and settings contig_min_cov=2, k_mer_size=11.", methods)
        self.assertIn("and settings k_mer_size=9.", methods)
        self.assertNotIn("an image other than the configured assembly image", methods)
        self.assertNotIn("could not be queried", methods)


    def test_methods_name_reused_assembly_images_and_unqueried_versions(self):
        same = {"library_id": "0460", "status": "recorded", "record": self.assembly()}
        other = build_assembly_provenance("NA", source=capture_source(self.repo), image=image("c"),
                                          settings={"k_mer_size": 11, "contig_min_cov": 2},
                                          reference=self.reference,
                                          read_pairs=[{"r1": "lane1_R1.fastq.gz", "r2": "lane1_R2.fastq.gz"}],
                                          outputs=self.outputs)
        versions = json.loads(json.dumps(self.versions))
        versions["assembly"]["software"][0].update(version=None, status="unknown")
        scoped = self.scope(tool_versions=versions, assembly_records={
            "0460": same, "NA": {"library_id": "NA", "status": "recorded", "record": other}})
        config = effective_scope_config(self.config, "project_b")
        config["params"]["lineage"].update(fit_edge=0.9, fit_specificity=0.99, fit_subsample=1000)
        config["params"].update(rarefaction={"targets": {"IGH": 1000}},
                                diversity={"min_consensus_count": {"IGH": 0}})
        libraries = [dict(row, scope_id="project_b") for row in self.libraries]
        threshold = {"scope_id": "project_b", "clone_threshold": 0.1, "threshold_source": "fallback"}
        methods = " ".join(methods_paragraphs(scoped, config, libraries, threshold))
        self.assertIn("The TRUST4 version of the configured assembly image could not be queried.", methods)
        self.assertIn("Assembly of 1 library used an image other than the configured assembly image.",
                      methods)
        self.assertNotIn("in the analysis image could not be queried", methods)
        self.assertNotIn("The queried assembly image contains", methods)


if __name__ == "__main__":
    unittest.main()
