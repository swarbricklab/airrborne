"""Reference curation keeps the first record per allele, and the audit names every disagreement.

Run: python test/run_tests.py --component annotation --level unit
"""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from workflow.lib.germline_reference import (
    audit_reference,
    build_auxiliary_data,
    build_reference,
    build_reference_metadata,
    render_reference_metadata,
    reference_identity_classes,
    compare_alleles,
    expected_database_alleles,
    keep_first,
    parse_auxiliary_data,
    read_imgt_entries,
    read_imgt_records,
    read_imgt_snapshot,
    read_j_alleles,
    reference_alleles,
    supplement_auxiliary_data,
    write_reference,
)

REPOSITORY = Path(__file__).resolve().parents[4]
IMAGE = {"runtime": "docker", "format": "docker", "digest": "sha256:" + "a" * 64,
         "configured_reference": "example/suite@sha256:" + "a" * 64}


def imgt_header(allele: str, region: str = "V-REGION") -> str:
    return f">X00001|{allele}|Homo_sapiens|F|{region}|1..10|10 nt|1| | | | |10+0=10| | |"


# The public IGHJ1*01 J-REGION, so that the auxiliary entry NCBI ships for it
# (frame 0, CDR3 end 17, one extra base) describes a real sequence. IGHJ1*02 is
# a constructed allele differing at position 30, outside the conserved codon.
IGHJ1 = "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG"
IGHJ1_VARIANT = "GCTGAATACTTCCAGCACTGGGGCCAGGGCCCCCTGGTCACCGTCTCCTCAG"

# A constructed tree with the cases the keep-first rule decides: a dual-locus
# gene listed in two files with one sequence, a constant allele listed twice
# with two spliced forms, and lower-case bases with IMGT gap dots.
TREE = {
    "vdj/imgt_human_IGHV.fasta": [(imgt_header("IGHV1-2*01"), "caa...ggt"), (imgt_header("IGHV1-2*02"), "CAA...GGA")],
    "vdj/imgt_human_IGHD.fasta": [(imgt_header("IGHD1-1*01", "D-REGION"), "GGTACAAC")],
    "vdj/imgt_human_IGHJ.fasta": [(imgt_header("IGHJ1*01", "J-REGION"), IGHJ1),
                                  (imgt_header("IGHJ1*02", "J-REGION"), IGHJ1_VARIANT),
                                  (imgt_header("IGHJ1A*01", "J-REGION"), IGHJ1)],
    "vdj/imgt_human_TRAV.fasta": [(imgt_header("TRAV14/DV4*01"), "ACGT")],
    "vdj/imgt_human_TRDV.fasta": [(imgt_header("TRAV14/DV4*01"), "ACGT"), (imgt_header("TRDV1*01"), "GGCC")],
    "vdj/imgt_human_TRAJ.fasta": [(imgt_header("TRAJ1*01", "J-REGION"), "TTCGG")],
    "vdj/imgt_human_TRBD.fasta": [(imgt_header("TRBD1*01", "D-REGION"), "GGGACAGG")],
    "constant/imgt_human_IGHC.fasta": [(imgt_header("IGHM*01", "CH1+CH2+M"), "AAAATTTT"),
                                       (imgt_header("IGHM*01", "CH1+CH2-CHS"), "AAAA")],
    "constant/imgt_human_TRAC.fasta": [(imgt_header("TRAC*01", "EX1"), "CCCC")],
}
AUXILIARY = "# name, frame, chain, cdr3 end, extra\n\nIGHJ1*01\t0\tJH\t17\t1\nIGHJ9*01\t0\tJH\t17\t1\n"


def write_tree(root: Path) -> tuple[Path, Path, Path]:
    for name, records in TREE.items():
        path = root / "germlines/imgt/human" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(f"{header}\n{sequence}\n" for header, sequence in records))
    snapshot = root / "germlines/imgt/IMGT.yaml"
    snapshot.write_text("source:  https://www.imgt.org/genedb\ndate:    2026.06.07\nspecies:\n"
                        "    - human:Homo+sapiens\n")
    return root / "germlines/imgt/human/vdj", root / "germlines/imgt/human/constant", snapshot


def database_dumps() -> dict[str, dict[str, str]]:
    """The eight databases imgt2igblast.sh would build from TREE."""
    dumps = {name: {} for name in ("ig_v", "ig_d", "ig_j", "ig_c", "tr_v", "tr_d", "tr_j", "tr_c")}
    dumps["ig_v"] = {"IGHV1-2*01": "CAAGGT", "IGHV1-2*02": "CAAGGA"}
    dumps["ig_d"] = {"IGHD1-1*01": "GGTACAAC"}
    dumps["ig_j"] = {"IGHJ1*01": IGHJ1, "IGHJ1*02": IGHJ1_VARIANT, "IGHJ1A*01": IGHJ1}
    dumps["ig_c"] = {"IGHM*01": "AAAATTTT"}
    dumps["tr_v"] = {"TRAV14/DV4*01": "ACGT", "TRDV1*01": "GGCC"}
    dumps["tr_d"] = {"TRBD1*01": "GGGACAGG"}
    dumps["tr_j"] = {"TRAJ1*01": "TTCGG"}
    dumps["tr_c"] = {"TRAC*01": "CCCC"}
    return dumps


class CurationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.vdj, self.constant, self.snapshot = write_tree(self.root)

    def test_records_take_the_allele_from_the_second_header_field_in_file_order(self):
        records = read_imgt_records(self.constant / "imgt_human_IGHC.fasta")
        self.assertEqual(records, [("IGHM*01", "AAAATTTT"), ("IGHM*01", "AAAA")])
        bad = self.root / "bad.fasta"
        bad.write_text(">IGHV1-2*01\nACGT\n")
        with self.assertRaisesRegex(ValueError, "pipe-delimited"):
            read_imgt_records(bad)

    def test_keep_first_normalises_case_keeps_dots_and_describes_skipped_records(self):
        kept, skipped = keep_first([("A*01", "acg.t", "one.fasta"), ("A*01", "ACGT", "two.fasta"),
                                    ("B*01", "ggg", "two.fasta"), ("B*01", "gga", "three.fasta")])
        self.assertEqual(kept, [("A*01", "ACG.T"), ("B*01", "GGG")])
        self.assertEqual(skipped, [
            {"allele": "A*01", "file": "two.fasta", "kept_file": "one.fasta", "identical_sequence": True},
            {"allele": "B*01", "file": "three.fasta", "kept_file": "two.fasta", "identical_sequence": False}])

    def test_reference_lists_vdj_then_constant_files_with_one_record_per_allele(self):
        text, curation = build_reference(self.vdj, self.constant)
        self.assertEqual(text, (
            ">IGHD1-1*01\nGGTACAAC\n"
            f">IGHJ1*01\n{IGHJ1}\n>IGHJ1*02\n{IGHJ1_VARIANT}\n>IGHJ1A*01\n{IGHJ1}\n"
            ">IGHV1-2*01\nCAA...GGT\n>IGHV1-2*02\nCAA...GGA\n"
            ">TRAJ1*01\nTTCGG\n>TRAV14/DV4*01\nACGT\n>TRBD1*01\nGGGACAGG\n"
            ">TRDV1*01\nGGCC\n>IGHM*01\nAAAATTTT\n>TRAC*01\nCCCC\n"))
        self.assertEqual((curation["records_read"], curation["alleles_written"]), (14, 12))
        self.assertEqual(curation["duplicate_records_skipped"], [
            {"allele": "TRAV14/DV4*01", "file": "imgt_human_TRDV.fasta", "kept_file": "imgt_human_TRAV.fasta",
             "identical_sequence": True},
            {"allele": "IGHM*01", "file": "imgt_human_IGHC.fasta", "kept_file": "imgt_human_IGHC.fasta",
             "identical_sequence": False}])
        self.assertEqual([Path(entry["path"]).name for entry in curation["germline_files"]],
                         ["imgt_human_IGHD.fasta", "imgt_human_IGHJ.fasta", "imgt_human_IGHV.fasta",
                          "imgt_human_TRAJ.fasta", "imgt_human_TRAV.fasta", "imgt_human_TRBD.fasta",
                          "imgt_human_TRDV.fasta", "imgt_human_IGHC.fasta", "imgt_human_TRAC.fasta"])
        for entry in curation["germline_files"]:
            self.assertEqual(entry["sha256"], hashlib.sha256(Path(entry["path"]).read_bytes()).hexdigest())
        # The written file is what reannotated_airr.R reads: bare headers,
        # dots removable, no allele twice.
        output = self.root / "resources" / "human_IMGT+C.fa"
        write_reference(text, output)
        self.assertEqual(output.read_text(), text)
        self.assertFalse(output.with_name(output.name + ".tmp").exists())
        self.assertEqual(reference_alleles(text)["IGHV1-2*01"], "CAAGGT")

    def test_empty_tree_and_missing_directories_fail(self):
        empty = self.root / "empty"
        (empty / "vdj").mkdir(parents=True)
        (empty / "constant").mkdir()
        with self.assertRaisesRegex(ValueError, "No IMGT germline FASTA"):
            build_reference(empty / "vdj", empty / "constant")
        (empty / "vdj" / "imgt_human_IGHV.fasta").write_text("")
        with self.assertRaisesRegex(ValueError, "No IMGT record"):
            build_reference(empty / "vdj", empty / "constant")

    def test_snapshot_and_auxiliary_data_are_read_as_written(self):
        snapshot = read_imgt_snapshot(self.snapshot)
        self.assertEqual((snapshot["source"], snapshot["date"], snapshot["species"]),
                         ("https://www.imgt.org/genedb", "2026.06.07", ["human:Homo+sapiens"]))
        self.snapshot.write_text("date: 2026.06.07\n")
        with self.assertRaisesRegex(ValueError, "source, date and species"):
            read_imgt_snapshot(self.snapshot)
        aux = self.root / "human_gl.aux"
        aux.write_text(AUXILIARY)
        self.assertEqual(list(parse_auxiliary_data(aux.read_text())), ["IGHJ1*01", "IGHJ9*01"])

    def test_expected_database_alleles_follow_the_image_build_script(self):
        files = sorted(self.vdj.glob("*.fasta")) + sorted(self.constant.glob("*.fasta"))
        self.assertEqual(expected_database_alleles(files, "tr", "v"), {"TRAV14/DV4*01": "ACGT", "TRDV1*01": "GGCC"})
        self.assertEqual(expected_database_alleles(files, "ig", "c"), {"IGHM*01": "AAAATTTT"})
        self.assertEqual(expected_database_alleles(files, "ig", "d"), {"IGHD1-1*01": "GGTACAAC"})
        self.assertEqual(compare_alleles({"A": "X", "B": "Y"}, {"B": "Z", "C": "W"}),
                         {"missing_alleles": ["A"], "unexpected_alleles": ["C"], "sequence_mismatches": ["B"]})


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.vdj, self.constant, self.snapshot = write_tree(self.root)
        self.igblast = self.root / "igblast"
        (self.igblast / "database").mkdir(parents=True)
        (self.igblast / "optional_file").mkdir()
        (self.igblast / "optional_file/human_gl.aux").write_text(AUXILIARY)
        self.reference = self.root / "resources/human_IMGT+C.fa"
        text, curation = build_reference(self.vdj, self.constant)
        write_reference(text, self.reference)
        self.metadata = self.root / "resources/reference.metadata.tsv"
        write_reference(render_reference_metadata(curation["allele_metadata"]), self.metadata)
        self.auxiliary = self.root / "resources/human_gl.aux"
        auxiliary, _ = build_auxiliary_data(self.vdj, self.igblast / "optional_file/human_gl.aux")
        write_reference(auxiliary, self.auxiliary)
        self.dumps = database_dumps()

    def read_database(self, path: Path) -> dict[str, str]:
        return dict(self.dumps[path.name.removeprefix("imgt_human_")])

    def audit(self):
        return audit_reference(reference=self.reference, auxiliary=self.auxiliary, metadata=self.metadata, germline_dir=self.vdj,
                               constant_dir=self.constant, snapshot=self.snapshot, igblast_db=self.igblast,
                               image=IMAGE, read_database=self.read_database)

    def test_receipt_records_agreement_sources_and_auxiliary_gaps(self):
        receipt = self.audit()
        self.assertEqual(receipt["agreement"], {"reference_rebuild": True, "metadata_rebuild": True, "auxiliary_rebuild": True,
                                                "igblast_databases": True})
        self.assertEqual(receipt["reference"], {
            "name": "human_IMGT+C.fa", "sha256": hashlib.sha256(self.reference.read_bytes()).hexdigest(),
            "size_bytes": self.reference.stat().st_size})
        self.assertEqual(receipt["image"], IMAGE)
        self.assertEqual(receipt["imgt_snapshot"]["date"], "2026.06.07")
        self.assertEqual(receipt["curation"]["alleles_written"], 12)
        self.assertEqual(receipt["curation"]["rebuild_matches_reference"], True)
        self.assertEqual(set(receipt["igblast_databases"]), set(self.dumps))
        self.assertEqual(receipt["igblast_databases"]["tr_v"]["alleles"], 2)
        self.assertEqual(receipt["igblast_databases"]["ig_c"]["against_reference"],
                         {"missing_from_reference": [], "sequence_mismatches": []})
        auxiliary = receipt["auxiliary_data"]
        self.assertEqual(auxiliary["shipped"], {
            "path": str(self.igblast / "optional_file/human_gl.aux"),
            "sha256": hashlib.sha256((self.igblast / "optional_file/human_gl.aux").read_bytes()).hexdigest(),
            "entries": 2, "names": 2})
        self.assertEqual(auxiliary["written"], {
            "name": "human_gl.aux", "sha256": hashlib.sha256(self.auxiliary.read_bytes()).hexdigest(), "entries": 4})
        self.assertEqual((auxiliary["rebuild_matches_written"], auxiliary["identical_pairs_checked"],
                          auxiliary["sibling_pairs_checked"]), (True, 0, 0))
        self.assertEqual(auxiliary["supplemented"], [
            {"allele": "IGHJ1*02", "functionality": "F", "rule": "same_gene_equal_length", "copied_from": ["IGHJ1*01"],
             "frame": 0, "chain_type": "JH", "cdr3_stop": 17, "extra_bps": 1},
            {"allele": "IGHJ1A*01", "functionality": "F", "rule": "identical_sequence", "copied_from": ["IGHJ1*01"],
             "frame": 0, "chain_type": "JH", "cdr3_stop": 17, "extra_bps": 1}])
        self.assertEqual([(entry["allele"], entry["functionality"]) for entry in auxiliary["j_alleles_without_entry"]],
                         [("TRAJ1*01", "F")])
        self.assertEqual(auxiliary["entries_without_j_allele"], ["IGHJ9*01"])
        json.dumps(receipt, allow_nan=False)

    def test_changed_reference_or_database_fails_with_the_disagreement_named(self):
        self.reference.write_text(self.reference.read_text().replace("CAA...GGT", "CAA...GGC"))
        with self.assertRaisesRegex(ValueError, r"(?s)differs from a fresh curation.*IGHV1-2\*01") as failure:
            self.audit()
        self.assertIn("sequence_mismatches", str(failure.exception))
        text, curation = build_reference(self.vdj, self.constant)
        write_reference(text, self.reference)
        self.metadata = self.root / "resources/reference.metadata.tsv"
        write_reference(render_reference_metadata(curation["allele_metadata"]), self.metadata)
        self.dumps["ig_c"]["IGHM*01"] = "AAAA"
        with self.assertRaisesRegex(ValueError, "imgt_human_ig_c disagrees with the reference"):
            self.audit()
        self.dumps["ig_c"] = {"IGHM*01": "AAAATTTT", "IGHG1*01": "GGGG"}
        with self.assertRaisesRegex(ValueError, "missing_from_reference.*IGHG1\\*01"):
            self.audit()
        self.dumps["ig_c"] = {}
        with self.assertRaisesRegex(ValueError, "missing_alleles.*IGHM\\*01"):
            self.audit()
        self.dumps = database_dumps()
        self.auxiliary.write_text(self.auxiliary.read_text().replace("IGHJ1*02\t0", "IGHJ1*02\t1"))
        with self.assertRaisesRegex(ValueError, "differs from a fresh supplement"):
            self.audit()


class IdentityMetadataTests(unittest.TestCase):
    def test_partial_unknown_and_ambiguous_sequences_never_establish_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference.fasta"
            fields = imgt_header("IGHV1-4*01").split("|")
            fields[13] = "partial in 5'"
            path.write_text(imgt_header("IGHV1-2*01") + "\nac..gt\n" +
                imgt_header("IGHV1-3*01") + "\nACGT\n" + "|".join(fields) + "\nACGT\n" +
                imgt_header("IGHV1-5*01") + "\nACNT\n" +
                ">X|IGHV1-6*01|human|F|V-REGION\nACGT\n" +
                imgt_header("IGHV1-2*02") + "\nCCCC\n")
            rows = build_reference_metadata([path])
            by_allele = {row["allele"]: row for row in rows}
            identity = "VSEQ_" + hashlib.sha256(b"ACGT").hexdigest()
            self.assertEqual(by_allele["IGHV1-2*01"]["identity_id"], identity)
            self.assertEqual(by_allele["IGHV1-3*01"]["identity_id"], identity)
            self.assertNotEqual(by_allele["IGHV1-2*02"]["identity_id"], identity)
            for allele in ("IGHV1-4*01", "IGHV1-5*01", "IGHV1-6*01"):
                self.assertEqual(by_allele[allele]["identity_id"], "")
            self.assertEqual(by_allele["IGHV1-4*01"]["completeness"], "partial in 5'")
            classes = reference_identity_classes(rows)
            shared = next(row for row in classes if row["identity_id"] == identity)
            self.assertEqual(shared["alleles"], ["IGHV1-2*01", "IGHV1-3*01"])
            before = shared["identity_id"]
            path.write_text(path.read_text() + imgt_header("IGHV1-7*01") + "\nACGT\n")
            changed = reference_identity_classes(build_reference_metadata([path]))
            shared = next(row for row in changed if row["identity_id"] == before)
            self.assertEqual(shared["alleles"], ["IGHV1-2*01", "IGHV1-3*01", "IGHV1-7*01"])


class SupplementTests(unittest.TestCase):
    """Supplement entries from identical sequences or compatible same-gene alleles."""

    SHIPPED = ("# name, frame, chain, cdr3 end, extra\n\nIGHJ1*01\t0\tJH\t17\t1\nIGHJ9*01\t0\tJH\t17\t1\n"
               "IGHJ1P*01\t2\tJH\n")

    def alleles(self, **extra):
        alleles = {"IGHJ1*01": {"sequence": IGHJ1, "functionality": "F"},
                   "IGHJ1*02": {"sequence": IGHJ1_VARIANT, "functionality": "F"},
                   # A substitution inside the conserved codon (positions 18 to 20).
                   "IGHJ1*03": {"sequence": IGHJ1[:19] + "A" + IGHJ1[20:], "functionality": "F"},
                   # One base shorter, so the shipped positions do not carry over.
                   "IGHJ1*04": {"sequence": IGHJ1[1:], "functionality": "F"},
                   "IGHJ1P*01": {"sequence": IGHJ1[2:] + "GG", "functionality": "P"},
                   # Identical to IGHJ1*01 under another gene name, as IMGT lists IGLJ2A beside IGLJ2.
                   "IGHJ1A*01": {"sequence": IGHJ1, "functionality": "F"},
                   "TRAJ7*01": {"sequence": "GGGACCAAGCTGACCGTCCTAG", "functionality": "P"}}
        alleles.update(extra)
        return alleles

    def test_entries_are_copied_only_when_length_and_conserved_codon_match(self):
        text, report = supplement_auxiliary_data(self.SHIPPED, self.alleles())
        self.assertTrue(text.startswith(self.SHIPPED))
        self.assertEqual(parse_auxiliary_data(text)["IGHJ1*02"], ("0", "JH", "17", "1"))
        self.assertEqual(report["supplemented"], [
            {"allele": "IGHJ1*02", "functionality": "F", "rule": "same_gene_equal_length", "copied_from": ["IGHJ1*01"],
             "frame": 0, "chain_type": "JH", "cdr3_stop": 17, "extra_bps": 1},
            {"allele": "IGHJ1A*01", "functionality": "F", "rule": "identical_sequence", "copied_from": ["IGHJ1*01"],
             "frame": 0, "chain_type": "JH", "cdr3_stop": 17, "extra_bps": 1}])
        self.assertEqual([(entry["allele"], entry["functionality"]) for entry in report["j_alleles_without_entry"]],
                         [("IGHJ1*03", "F"), ("IGHJ1*04", "F"), ("TRAJ7*01", "P")])
        self.assertEqual((report["shipped_entries"], report["shipped_names"], report["identical_pairs_checked"],
                          report["sibling_pairs_checked"]), (3, 3, 0, 0))
        # The frame-only pseudogene entry is kept as shipped and never copied from.
        self.assertEqual(parse_auxiliary_data(text)["IGHJ1P*01"], ("2", "JH"))
        self.assertEqual(text.count("IGHJ1P*01"), 1)

    def test_nothing_is_appended_when_every_allele_has_an_entry(self):
        alleles = {"IGHJ1*01": {"sequence": IGHJ1, "functionality": "F"}}
        text, report = supplement_auxiliary_data(self.SHIPPED, alleles)
        self.assertEqual(text, self.SHIPPED)
        self.assertEqual((report["supplemented"], report["j_alleles_without_entry"]), ([], []))

    def test_inconsistent_shipped_siblings_stop_the_build(self):
        shipped = self.SHIPPED + "IGHJ1*02\t1\tJH\t17\t1\n"
        with self.assertRaisesRegex(ValueError, r"IGHJ1\*01 and IGHJ1\*02, same-gene alleles.*different entries"):
            supplement_auxiliary_data(shipped, self.alleles())
        identical = self.SHIPPED + "IGHJ1A*01\t1\tJH\t17\t1\n"
        with self.assertRaisesRegex(ValueError, r"IGHJ1\*01 and IGHJ1A\*01, alleles with identical sequences"):
            supplement_auxiliary_data(identical, self.alleles())
        consistent = self.SHIPPED + "IGHJ1A*01\t0\tJH\t17\t1\n"
        _, report = supplement_auxiliary_data(consistent, self.alleles())
        self.assertEqual(report["identical_pairs_checked"], 1)
        # The same two entries are consistent once the alleles differ in length.
        alleles = self.alleles(**{"IGHJ1*02": {"sequence": IGHJ1_VARIANT[1:], "functionality": "F"}})
        _, report = supplement_auxiliary_data(shipped, alleles)
        self.assertEqual(report["sibling_pairs_checked"], 0)

    def test_repeated_and_malformed_shipped_lines(self):
        self.assertEqual(parse_auxiliary_data("A*01\t0\tJH\t17\t1\nA*01\t0\tJH\t17\t1\n"),
                         {"A*01": ("0", "JH", "17", "1")})
        with self.assertRaisesRegex(ValueError, "twice with different fields"):
            parse_auxiliary_data("A*01\t0\tJH\t17\t1\nA*01\t1\tJH\t17\t1\n")
        with self.assertRaisesRegex(ValueError, "must hold a frame"):
            parse_auxiliary_data("A*01\tx\tJH\n")

    def test_j_alleles_come_from_the_j_files_with_their_functionality(self):
        with tempfile.TemporaryDirectory() as temporary:
            vdj, _, _ = write_tree(Path(temporary))
            alleles = read_j_alleles(vdj)
            self.assertEqual(set(alleles), {"IGHJ1*01", "IGHJ1*02", "IGHJ1A*01", "TRAJ1*01"})
            self.assertEqual(alleles["IGHJ1*01"], {"sequence": IGHJ1, "functionality": "F", "file": "imgt_human_IGHJ.fasta"})
            self.assertEqual(read_imgt_entries(vdj / "imgt_human_IGHJ.fasta")[1]["allele"], "IGHJ1*02")


if __name__ == "__main__":
    unittest.main()
