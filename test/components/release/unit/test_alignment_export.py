"""Check reference intervals against constructed alignments and reverse-D evidence."""
import copy
import gzip
from pathlib import Path
import tempfile
import unittest

from workflow.lib.alignment_export import (
    ALIGNMENT_FIELDS, add_alignment_evidence, cigar_interval, normalise_alignment, validate_alignment)


class AlignmentExportTests(unittest.TestCase):
    def test_partial_v_indels_and_trailing_offsets_keep_aligned_endpoints(self):
        row = dict.fromkeys(ALIGNMENT_FIELDS, "")
        row.update(sequence_id="synthetic", locus="IGH", sequence="TTACTGTTGGGG",
                   v_call="V*01", igblast_v_subject_id="V*01", v_cigar="2S3N2=1I1=1D2=4S5N",
                   v_sequence_start="3", v_sequence_end="8", v_germline_start="1", v_germline_end="40")
        native = copy.deepcopy(row)
        normalise_alignment(row)
        self.assertEqual((row["v_germline_start"], row["v_germline_end"]), ("4", "9"))
        self.assertEqual((row["changeo_v_germline_start"], row["changeo_v_germline_end"]), ("1", "40"))
        self.assertEqual(native["v_germline_end"], "40")
        validate_alignment(row, {"V*01": "GGGACGATT"})
        with self.assertRaisesRegex(ValueError, "reference bases"):
            validate_alignment(row, {"V*01": "GGGACGATA"})
        with self.assertRaisesRegex(ValueError, "exceeds"):
            validate_alignment(row, {"V*01": "GGGACGAT"})

    def test_reverse_d_keeps_native_hit_and_masks_unsupported_standard_fields(self):
        row = dict.fromkeys(ALIGNMENT_FIELDS, "")
        row.update(sequence_id="reverse", locus="IGH", sequence="TTCAT", d_call="D*01",
                   d_cigar="2S2N3=", d_germline_start="3", d_germline_end="2",
                   d_sequence_start="3", d_sequence_end="5", changeo_d_cigar="2S2N3=",
                   igblast_d_subject_id="D*01", igblast_d_subject_start="3", igblast_d_subject_end="1",
                   igblast_d_query_start="3", igblast_d_query_end="5",
                   igblast_d_query_alignment="CAT", igblast_d_germline_alignment="CAT", igblast_d_btop="3")
        normalise_alignment(row)
        self.assertEqual(row["d_alignment_status"], "unavailable_reverse_orientation")
        self.assertEqual(row["germline_reconstruction_status"], "unavailable_reverse_d")
        self.assertEqual(row["changeo_d_germline_end"], "2")
        self.assertEqual(row["changeo_d_cigar"], "2S2N3=")
        self.assertTrue(all(row[field] == "" for field in ("d_cigar", "d_germline_start", "d_germline_end")))
        validate_alignment(row, {"D*01": "ATG"})
        row["d_cigar"] = "2S2N3="
        with self.assertRaises(ValueError):
            validate_alignment(row, {"D*01": "ATG"})

    def test_forward_d_selected_allele_and_native_query_are_checked(self):
        row = dict.fromkeys(ALIGNMENT_FIELDS, "")
        row.update(sequence_id="forward", locus="TRB", sequence="TTATG", d_call="D*01,D*02",
                   d_cigar="2S3=", d_germline_start="1", d_germline_end="3",
                   d_sequence_start="3", d_sequence_end="5", igblast_d_subject_id="D*01",
                   igblast_d_subject_start="1", igblast_d_subject_end="3",
                   igblast_d_query_start="3", igblast_d_query_end="5",
                   igblast_d_query_alignment="ATG", igblast_d_germline_alignment="ATG")
        normalise_alignment(row)
        self.assertEqual(row["d_alignment_status"], "available")
        self.assertEqual(row["germline_reconstruction_status"], "not_applicable")
        validate_alignment(row, {"D*01": "ATG", "D*02": "CGA"})
        row["igblast_d_query_start"] = "2"
        with self.assertRaisesRegex(ValueError, "query bases"):
            validate_alignment(row, {"D*01": "ATG"})

    def test_native_first_hit_survives_equal_score_alternative(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "hits.gz"
            first = ["D", "reversed|s", "D*01", "100", "3", "0", "0", "0", "2", "4", "3", "1", "0", "10", "CAT", "CAT", "3"]
            second = first.copy(); second[2] = "D*02"
            with gzip.open(path, "wt") as handle:
                handle.write("# native hits\n" + "\t".join(first) + "\n" + "\t".join(second) + "\n")
            annotations = {"s": {}}
            self.assertEqual(add_alignment_evidence(annotations, [], path), list(ALIGNMENT_FIELDS))
            self.assertEqual(annotations["s"]["igblast_d_subject_id"], "D*01")
            self.assertEqual(annotations["s"]["igblast_d_subject_end"], "1")

    def test_empty_or_malformed_cigar_is_not_an_alignment(self):
        for cigar in ("", "3N", "0M", "2M2N2M", "1I0D", "2Mgarbage"):
            with self.subTest(cigar=cigar), self.assertRaises(ValueError):
                cigar_interval(cigar)
