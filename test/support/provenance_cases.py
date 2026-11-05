"""Construct recorded assembly and annotation evidence for scientific provenance tests."""

import subprocess
import tempfile
from pathlib import Path

from test.support.assembly_record import build_assembly_provenance
from workflow.lib.provenance import build_scope_provenance, capture_source, file_identity


def image(digest="a", **changes):
    value = {"runtime": "docker", "format": "docker", "digest": "sha256:" + digest * 64,
             "configured_reference": "example/image@sha256:" + digest * 64,
             "image_id": "sha256:" + digest * 64, "architecture": "amd64", "os": "linux"}
    value.update(changes)
    return value


class ProvenanceCase:
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.repo = self.root / "repo"
        self.repo.mkdir()
        self.git("init", "-q")
        (self.repo / "workflow/lib").mkdir(parents=True)
        (self.repo / "workflow/lib" / "analysis.R").write_text("value <- 1\n")
        (self.repo / ".gitignore").write_text(".snakemake/\n")
        self.git("add", ".")
        self.git("-c", "user.name=Fixture", "-c", "user.email=fixture@example.invalid", "commit", "-qm", "fixture")
        self.reference = self.root / "human.fa"
        self.reference.write_text(">gene\nACGT\n")
        self.output = self.root / "0460_airr.tsv"
        self.output.write_text("sequence_id\nNA\n")
        self.outputs = {"airr": self.output}
        self.config = {"params": {"trust4": {"k_mer_size": 9, "contig_min_cov": 0, "cleanup": 0},
                                  "airr_qc": {"max_alignment_n_fraction": 0.1,
                                              "min_informative_c_alignment_positions": 40},
                                  "lineage": {"clonal_threshold": "auto", "min_fit_distances": 500, "max_n": 0},
                                  "igblast": {"num_alignments_v": 20, "num_alignments_j": 20,
                                              "num_alignments_d": 3, "num_alignments_c": 125},
                                  "shm": {"min_position_votes": 2,
                                          "min_supported_subregion_fraction": 1.0,
                                          "min_germline_occupied_v_positions": 270},
                                  "selection": {"collapse_min_freq": 0.6},
                                  "overlap": {"loci": ["IGH", "TRB"], "min_consensus_count": {"IGH": 0, "TRB": 0}}},
                       "refs": {"imgt_reference_name": "human.fa"},
                       "deps": {"immcantation_docker": image()["configured_reference"],
                                "trust4_container": "docker://" + image("b")["configured_reference"],
                                "samplesheet": "all-libraries.csv"},
                       "scope_metadata": {"project_b": {"tissue": {"file": "b.csv", "columns": ["response"]}},
                                          "private_project": {"library": {"file": "private.csv", "columns": ["batch"]}}},
                       "membership_csv": "all-projects.csv", "outs": {"results_dir": "global-results"}}
        self.libraries = [{"library_id": "0460", "tissue_id": "NA", "donor_id": "001"},
                          {"library_id": "NA", "tissue_id": "NA", "donor_id": "001"}]
        self.images = {"analysis": image(), "assembly": image("b")}
        # The receipt audit_imgt_reference.py writes for self.reference.
        self.receipt = {"schema_version": "1.0", "image": image(),
                        "reference": file_identity(self.reference),
                        "imgt_snapshot": {"path": "/germlines/imgt/IMGT.yaml",
                                          "source": "https://www.imgt.org/genedb", "date": "2026.06.07",
                                          "species": ["human:Homo+sapiens"]},
                        "germline_files": [], "curation": {"alleles_written": 1},
                        "igblast_databases": {},
                        "agreement": {"reference_rebuild": True, "auxiliary_rebuild": True, "igblast_databases": True},
                        "auxiliary_data": {
                            "shipped": {"path": "/igblast/optional_file/human_gl.aux", "sha256": "c" * 64,
                                        "entries": 132, "names": 131},
                            "written": {"name": "human_gl.aux", "sha256": "d" * 64, "entries": 140},
                            "rebuild_matches_written": True, "identical_pairs_checked": 3, "sibling_pairs_checked": 34,
                            "supplemented": [{"allele": "IGHJ5*04", "functionality": "F", "rule": "same_gene_equal_length",
                                              "copied_from": ["IGHJ5*01", "IGHJ5*02", "IGHJ5*03"], "frame": 2,
                                              "chain_type": "JH", "cdr3_stop": 16, "extra_bps": 1},
                                             {"allele": "IGLJ2A*01", "functionality": "F", "rule": "identical_sequence",
                                              "copied_from": ["IGLJ2*01"], "frame": 1, "chain_type": "JL",
                                              "cdr3_stop": 6, "extra_bps": 1}],
                            "j_alleles_without_entry": [{"allele": "TRAJ55*01", "functionality": "P",
                                                         "reason": "no same-gene allele of equal length"}],
                            "entries_without_j_allele": []}}
        self.versions = {role: {"image": identity, "software": [
            {"name": "IgBLAST" if role == "analysis" else "TRUST4", "version": "1.0-fixture",
             "status": "measured", "command": ["version-command"], "method": "fixture query"}]}
                         for role, identity in self.images.items()}


    def git(self, *arguments):
        return subprocess.run(["git", "-C", str(self.repo), *arguments], check=True,
                              capture_output=True, text=True).stdout.strip()


    def scope(self, **changes):
        args = {"source": capture_source(self.repo), "images": self.images, "tool_versions": self.versions,
                "references": {"assembly_reference": self.reference}, "reference_receipt": self.receipt}
        args.update(changes)
        return build_scope_provenance("project_b", self.libraries, self.config, **args)


    def assembly(self):
        return build_assembly_provenance("0460", source=capture_source(self.repo), image=image("b"),
                                          settings={"k_mer_size": 9}, reference=self.reference,
                                          read_pairs=[{"r1": "lane1_R1.fastq.gz", "r2": "lane1_R2.fastq.gz"},
                                                      {"r1": "lane2_R1.fastq.gz", "r2": "lane2_R2.fastq.gz"}],
                                          outputs=self.outputs)
