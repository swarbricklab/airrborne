"""Export linkage and real AIRR 2.0 validation on small, literal TSV frames.

Run in an environment with airr==2.0.0. The analysis image bundles AIRR 1.6,
which is intentionally insufficient evidence for the declared export schema.

Run: python test/run_tests.py --component release --level integration
"""

import copy
import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import airr

from test.components.release.integration.check_release import check_semantic_counts
from workflow.lib.repertoire_metadata import (
    SEMANTIC_COUNTS, build_scope_repertoires, validate_scope_repertoires,
)

from workflow.lib.airr_exports import (
    EXPORT_DROPPED,
    EXPORT_NAMES,
    EXPORT_ORDER,
    GERMLINE_FIELDS,
    REGION_FIELDS,
    add_makedb_fields,
    annotate_member_assignments,
    export_library_airr,
    read_library_annotations,
    QUALIFIED_CALL_COLUMNS,
    validate_export_content,
    validate_library_airr,
    validate_cigar_query_spans,
)


ROOT = Path(__file__).resolve().parents[4]
FIELDS = ["sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call", "j_call",
          "sequence_alignment", "germline_alignment", "junction", "junction_aa", "cdr3", "cdr3_aa",
          "v_cigar", "d_cigar", "j_cigar", "library_id", "consensus_count",
          "v_gene", "j_gene", "v_gene_status", "j_gene_status", "c_call",
          "v_sequence_start", "v_sequence_end", "d_sequence_start", "d_sequence_end",
          "j_sequence_start", "j_sequence_end"]


# The isotype-call sidecar qualify_isotype_calls writes for every IGH record.
SIDECAR_COLUMNS = ["library_id", "sequence_id", "c_call", "qualified_c_call", "c_call_present",
                   "c_call_locus_match", "c_call_enough_evidence", "c_call_candidates_complete",
                   "c_call_callable", "isotype_class", "isotype_subclass"]
# An export of FIELDS without annotations, in EXPORT_ORDER.
UNANNOTATED_ORDER = ["sequence_id", "library_id", "source_sequence_id", "repertoire_id",
                     "sample_processing_id", "data_processing_id", "productive", "rev_comp", "changeo_rev_comp",
                     "v_call", "d_call", "j_call", "c_call", "v_gene_key", "j_gene_key",
                     "junction", "junction_aa", "cdr3", "cdr3_aa", "consensus_count",
                     "sequence", "sequence_alignment", "germline_alignment",
                     "v_cigar", "d_cigar", "j_cigar", "v_sequence_start", "v_sequence_end",
                     "d_sequence_start", "d_sequence_end", "j_sequence_start", "j_sequence_end",
                     "v_gene_status", "j_gene_status"]


def source_row(library_id, sequence_id="0460", **changes):
    row = dict.fromkeys(FIELDS, "")
    row.update(sequence_id=sequence_id, sequence="ACTG", rev_comp="F", productive="T",
               v_call="IGHV1-2*01", j_call="IGHJ1*01", library_id=library_id,
               consensus_count="10", junction="TGTGCC", junction_aa="CA",
               v_gene="IGHV1-2", j_gene="IGHJ1", v_gene_status="resolved",
               j_gene_status="resolved", c_call="NA")
    row.update(changes)
    return row


class AirrExportTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.processing = "recorded-post-qc-settings-and-reference-sha256"
        self.libraries = [
            {"library_id": "0460", "tissue_id": "001", "donor_id": "NA"},
            {"library_id": "NA", "tissue_id": "001", "donor_id": "NA"},
            {"library_id": "outside", "tissue_id": "002", "donor_id": "007"},
        ]
        self.files = {row["library_id"]: "airr/" + row["library_id"] + ".tsv.gz"
                      for row in self.libraries}
        for row in self.libraries:
            library_id = row["library_id"]
            self.write_source(library_id, [source_row(library_id),
                                           source_row(library_id, "NA", consensus_count="0")])
            export_library_airr(self.source(library_id), self.root / self.files[library_id],
                                library_id, self.processing)

    def test_reverse_source_exports_a_coherent_stored_query_orientation(self):
        original_query = "CCGTTAGG"
        oriented_query = "CCTAACGG"
        self.write_source("0460", [source_row(
            "0460", sequence=oriented_query, rev_comp="T", v_cigar="2S4=",
            v_sequence_start="3", v_sequence_end="6", sequence_alignment="TAAC")])
        export_library_airr(self.source("0460"), self.root / self.files["0460"],
                            "0460", self.processing)
        _, rows = self.read_export()
        row = rows[0]
        self.assertEqual(row["sequence"], oriented_query)
        self.assertEqual(row["changeo_rev_comp"], "T")
        self.assertEqual(row["rev_comp"], "F")
        reverse = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        self.assertEqual(reverse(row["sequence"]), original_query)
        query = reverse(row["sequence"]) if row["rev_comp"] == "T" else row["sequence"]
        self.assertEqual(query[int(row["v_sequence_start"]) - 1:int(row["v_sequence_end"])],
                         row["sequence_alignment"])
        self.assertEqual(query[2:6], "TAAC")
        self.assertEqual(row["v_cigar"], "2S4=")
        self.assertEqual(validate_library_airr(self.root / self.files["0460"],
                                              "0460", self.processing), 1)
        row["rev_comp"] = "T"
        with gzip.open(self.root / self.files["0460"], "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=row, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        with self.assertRaisesRegex(ValueError, "orientation"):
            validate_library_airr(self.root / self.files["0460"], "0460", self.processing)

    def source(self, library_id):
        return self.root / (library_id + ".tsv")

    def write_source(self, library_id, rows, fields=FIELDS):
        with self.source(library_id).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def metadata(self, libraries=None, **kwargs):
        return build_scope_repertoires("project_b", self.libraries[:2] if libraries is None else libraries,
                                       self.files, self.processing, **kwargs)

    def validate(self, data):
        return validate_scope_repertoires(data, self.root, self.libraries[:2], self.files, self.processing)

    def read_export(self, library_id="0460"):
        with gzip.open(self.root / self.files[library_id], "rt", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return reader.fieldnames, list(reader)

    def rewrite_export(self, fields, rows, library_id="0460"):
        with gzip.open(self.root / self.files[library_id], "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def test_reference_context_leaves_the_export_and_gene_names_are_renamed(self):
        fields = ["v_gene_label", "j_gene_label"] + [
            f"{segment}_{suffix}" for segment in ("v", "d", "j", "c")
            for suffix in ("reference_functionality", "reference_orphon", "call_candidates_complete")]
        context = {field: "" for field in fields}
        context.update(v_gene_label="IGHV1-2", v_reference_functionality="IGHV1-2*01:F;IGHV1-2*02:ORF",
                       v_reference_orphon="F", v_call_candidates_complete="TRUE",
                       c_reference_functionality="IGHM*01:F", c_call_candidates_complete="FALSE")
        source = source_row("0460", **context)
        self.write_source("0460", [source], FIELDS + fields)
        export_library_airr(self.source("0460"), self.root / self.files["0460"], "0460", self.processing)
        exported, rows = self.read_export()
        # The reference metadata holds functionality and orphon status per allele.
        for field in fields:
            if field in EXPORT_DROPPED:
                self.assertNotIn(field, exported)
            else:
                self.assertEqual(rows[0][EXPORT_NAMES.get(field, field)], source[field])
        self.assertNotIn("v_gene_label", exported)
        self.assertEqual(rows[0]["v_gene_names"], "IGHV1-2")
        self.assertEqual(set(EXPORT_DROPPED) - set(fields), set())
        self.assertEqual(validate_export_content(self.source("0460"), self.root / self.files["0460"],
                                                "0460", self.processing), 1)


    def test_shared_export_preserves_cells_and_qualifies_identifiers(self):
        fields, rows = self.read_export()
        self.assertEqual(fields, UNANNOTATED_ORDER)
        self.assertNotIn("v_gene", fields)
        self.assertEqual([rows[0]["v_gene_key"], rows[0]["j_gene_key"]], ["IGHV1-2", "IGHJ1"])
        self.assertEqual([row["sequence_id"] for row in rows], ["0460~0460", "0460~NA"])
        self.assertEqual([row["source_sequence_id"] for row in rows], ["0460", "NA"])
        self.assertEqual([row["library_id"] for row in rows], ["0460", "0460"])
        self.assertEqual([row["consensus_count"] for row in rows], ["10", "0"])
        self.assertEqual(rows[0]["c_call"], "NA")
        self.assertEqual(rows[0]["d_call"], "")
        self.assertNotIn("clone_id", fields)
        count = validate_export_content(self.source("0460"), self.root / self.files["0460"],
                                        "0460", self.processing)
        self.assertEqual(count, 2)
        # The official public file API also validates the decompressed export.
        plain = self.root / "decompressed.tsv"
        plain.write_bytes(gzip.decompress((self.root / self.files["0460"]).read_bytes()))
        self.assertTrue(airr.validate_rearrangement(str(plain), debug=True))


    def test_content_comparison_catches_scientific_and_missing_value_changes(self):
        fields, rows = self.read_export()
        for column, value in [("consensus_count", "11"), ("c_call", ""), ("library_id", "460")]:
            with self.subTest(column=column):
                altered = copy.deepcopy(rows)
                altered[0][column] = value
                self.rewrite_export(fields, altered)
                with self.assertRaisesRegex(ValueError, "differs from its source"):
                    validate_export_content(self.source("0460"), self.root / self.files["0460"],
                                            "0460", self.processing)

    def test_official_rearrangement_validator_rejects_missing_fields_and_bad_types(self):
        fields, rows = self.read_export()
        without = [field for field in fields if field != "junction"]
        self.rewrite_export(without, [{key: value for key, value in row.items() if key in without}
                                     for row in rows])
        with self.assertRaises(airr.ValidationError):
            validate_library_airr(self.root / self.files["0460"], "0460", self.processing)
        rows[0]["consensus_count"] = "not-a-count"
        self.rewrite_export(fields, rows)
        with self.assertRaises(airr.ValidationError):
            validate_library_airr(self.root / self.files["0460"], "0460", self.processing)

    def test_duplicate_native_rows_and_wrong_library_fail(self):
        for rows in [[source_row("0460"), source_row("0460")], [source_row("460")]]:
            self.write_source("0460", rows)
            with self.assertRaises(ValueError):
                export_library_airr(self.source("0460"), self.root / "bad.tsv.gz", "0460", self.processing)

    def test_cigar_semantics_check_all_segments_with_offsets_and_optional_tails(self):
        for segment in ("v", "d", "j"):
            for cigar, start, end in [("3=1I2=", "1", "6"), ("3=1D2=", "1", "5"),
                                      ("4S3N3=1I2=4S2N", "5", "10"),
                                      ("2=1I2=1D2=", "1", "7")]:
                with self.subTest(segment=segment, cigar=cigar):
                    row = source_row("0460", sequence="A" * 14, **{
                        f"{segment}_cigar": cigar, f"{segment}_sequence_start": start,
                        f"{segment}_sequence_end": end})
                    validate_cigar_query_spans(row)
                    self.write_source("0460", [row])
                    export_library_airr(self.source("0460"), self.root / self.files["0460"],
                                        "0460", self.processing)
                    self.assertTrue(self.validate(self.metadata())["cigar_query_spans_valid"])

    def test_staged_cigar_corruption_is_rejected_after_schema_validation(self):
        fields, rows = self.read_export()
        for segment in ("v", "d", "j"):
            for cigar in ("2=1D1=", "1S4=", "4Q", "2=1S2="):
                with self.subTest(segment=segment, cigar=cigar):
                    altered = copy.deepcopy(rows)
                    altered[0].update({f"{segment}_cigar": cigar,
                                       f"{segment}_sequence_start": "1",
                                       f"{segment}_sequence_end": "4"})
                    airr.schema.RearrangementSchema.validate_row(altered[0])
                    self.rewrite_export(fields, altered)
                    with self.assertRaisesRegex(ValueError, f"0460.*{segment}_cigar"):
                        self.validate(self.metadata())

    def test_present_cigar_requires_query_coordinates_and_sequence_bounds(self):
        for start, end in (("", ""), ("1", "5")):
            row = source_row("0460", v_cigar="5=", v_sequence_start=start, v_sequence_end=end)
            with self.assertRaisesRegex(ValueError, "0460 v_cigar"):
                validate_cigar_query_spans(row)

    def test_empty_library_still_has_valid_headers_and_unambiguous_file_links(self):
        self.write_source("0460", [])
        export_library_airr(self.source("0460"), self.root / self.files["0460"], "0460", self.processing)
        self.assertEqual(validate_export_content(self.source("0460"), self.root / self.files["0460"],
                                                "0460", self.processing), 0)
        self.assertEqual(validate_library_airr(self.root / self.files["0460"], "0460", self.processing), 0)
        self.assertTrue(self.validate(self.metadata())["structural_valid"])
        self.write_source("NA", [])
        export_library_airr(self.source("NA"), self.root / self.files["NA"], "NA", self.processing)
        data = self.metadata()
        data["Repertoire"][0]["data_processing"][0]["data_processing_files"] = [self.files["NA"]]
        with self.assertRaisesRegex(ValueError, "file link differs"):
            self.validate(data)

    def test_clone_assignments_are_not_accepted_in_shared_annotation(self):
        self.write_source("0460", [source_row("0460", clone_id="project_a|donor|1")], FIELDS + ["clone_id"])
        with self.assertRaisesRegex(ValueError, "must not contain clone assignments"):
            export_library_airr(self.source("0460"), self.root / "bad.tsv.gz", "0460", self.processing)

    def test_authority_annotations_cover_source_records_and_cannot_replace_native_values(self):
        annotations = {"0460": {"backbone_id": "contig1", "inference_consensus_count": "10"},
                       "NA": {"backbone_id": "contig1", "inference_consensus_count": ""}}
        fields = ["backbone_id", "inference_consensus_count"]
        output = self.root / "authority.tsv.gz"
        export_library_airr(self.source("0460"), output, "0460", self.processing,
                            annotations=annotations, annotation_fields=fields)
        self.assertEqual(validate_export_content(self.source("0460"), output, "0460", self.processing,
                                                annotations=annotations, annotation_fields=fields), 2)
        for change in [{"0460": annotations["0460"]}, {**annotations, "extra": annotations["NA"]}]:
            with self.assertRaisesRegex(ValueError, "authority annotation|absent from post-QC"):
                export_library_airr(self.source("0460"), output, "0460", self.processing,
                                    annotations=change, annotation_fields=fields)
        with self.assertRaisesRegex(ValueError, "overwrite source"):
            export_library_airr(self.source("0460"), output, "0460", self.processing,
                                annotations=annotations, annotation_fields=["consensus_count"])
        with self.assertRaisesRegex(ValueError, "must not contain clone assignments"):
            export_library_airr(self.source("0460"), output, "0460", self.processing,
                                annotations=annotations, annotation_fields=["clone_id"])

    def test_qualified_call_export_preserves_source_and_marks_nonigh_inapplicable(self):
        self.write_source("0460", [
            source_row("0460", c_call="IGHM*01", locus="IGH", c_call_candidates_complete="T"),
            source_row("0460", "NA", c_call="TRBC1", locus="TRB",
                       v_call="TRBV1*01", j_call="TRBJ1-1*01",
                       v_gene="TRBV1", j_gene="TRBJ1-1", c_call_candidates_complete="T"),
        ], FIELDS + ["locus", "c_call_candidates_complete"])
        columns = ["library_id", "sequence_id", "locus", "assembly_backbone_id",
                   "grouping_eligible", "grouping_status",
                   "backbone_representative_sequence_id", "is_backbone_representative",
                   "inference_representative_sequence_id", "is_inference_representative",
                   "inference_consensus_count"]
        paths = [self.root / "igh.tsv", self.root / "nonigh.tsv"]
        for path, key, locus in zip(paths, ("0460", "NA"), ("IGH", "TRB")):
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(columns)
                writer.writerow(["0460", key, locus, key, "TRUE", "eligible", key, "TRUE", key, "TRUE", "10"])
        calls = self.root / "calls.tsv"
        with calls.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(SIDECAR_COLUMNS)
            writer.writerow(["0460", "0460", "IGHM*01", "", "TRUE", "TRUE", "FALSE", "TRUE", "FALSE", "", ""])
        annotations, fields = read_library_annotations(paths, calls, "0460")
        output = self.root / self.files["0460"]
        export_library_airr(self.source("0460"), output, "0460", self.processing,
                            annotations=annotations, annotation_fields=fields)
        self.assertEqual(validate_export_content(self.source("0460"), output, "0460",
                         self.processing, annotations=annotations, annotation_fields=fields), 2)
        exported, rows = self.read_export()
        self.assertEqual([row["c_call"] for row in rows], ["IGHM*01", "TRBC1"])
        # The gated call is c_call where c_call_callable is T, so the sidecar's
        # qualified_c_call stays out of the export.
        self.assertNotIn("qualified_c_call", exported)
        start = exported.index("v_gene_key")
        self.assertEqual(exported[start:start + 4], ["v_gene_key", "j_gene_key", "isotype_class", "isotype_subclass"])
        start = exported.index("v_gene_status")
        self.assertEqual(exported[start:], ["v_gene_status", "j_gene_status", "c_call_present", "c_call_locus_match",
                                            "c_call_enough_evidence", "c_call_candidates_complete", "c_call_callable"])
        self.assertEqual(rows[0]["c_call_callable"], "F")
        self.assertEqual([rows[0]["isotype_class"], rows[0]["isotype_subclass"]], ["", ""])
        self.assertEqual([row["contig_id"] for row in rows], ["0460~0460", "0460~NA"])
        self.assertEqual([row["contig_representative_sequence_id"] for row in rows],
                         ["0460~0460", "0460~NA"])
        self.assertEqual([row["is_contig_representative"] for row in rows], ["T", "T"])
        for internal in ("assembly_backbone_id", "backbone_representative_sequence_id",
                         "is_backbone_representative", "v_gene", "j_gene"):
            self.assertNotIn(internal, exported)
        self.assertTrue(all(rows[1][column] == "" for column in QUALIFIED_CALL_COLUMNS))
        calls.write_text("\t".join(SIDECAR_COLUMNS) + "\n")
        with self.assertRaisesRegex(ValueError, "cover every common-QC IGH"):
            read_library_annotations(paths, calls, "0460")

    def member_export(self):
        """A two-member export built through the grouping annotations."""
        self.write_source("0460", [source_row("0460", locus="IGH"), source_row("0460", "NA", locus="TRB")],
                          FIELDS + ["locus"])
        columns = ["library_id", "sequence_id", "locus", "assembly_backbone_id", "grouping_eligible",
                   "grouping_status", "backbone_representative_sequence_id", "is_backbone_representative",
                   "inference_representative_sequence_id", "is_inference_representative",
                   "inference_consensus_count"]
        paths = [self.root / "igh.tsv", self.root / "nonigh.tsv"]
        for path, key, locus in zip(paths, ("0460", "NA"), ("IGH", "TRB")):
            path.write_text("\t".join(columns) + "\n" + "\t".join(
                ["0460", key, locus, key, "TRUE", "eligible", key, "TRUE", key, "TRUE", "10"]) + "\n")
        calls = self.root / "calls.tsv"
        calls.write_text("\t".join(SIDECAR_COLUMNS) + "\n"
                         + "\t".join(["0460", "0460", "IGHM*01", "IGHM*01", "TRUE", "TRUE", "TRUE", "TRUE",
                                      "TRUE", "IgM", "IgM"]) + "\n")
        return read_library_annotations(paths, calls, "0460")

    def test_blocked_sibling_status_takes_the_contig_name(self):
        annotations, _ = self.member_export()
        igh = self.root / "igh.tsv"
        header, row = igh.read_text().splitlines()
        cells = row.split("\t")
        cells[4:11] = ["FALSE", "ineligible_backbone_representative", "", "FALSE", "", "FALSE", "0"]
        igh.write_text(header + "\n" + "\t".join(cells) + "\n")
        annotations, _ = read_library_annotations([igh, self.root / "nonigh.tsv"],
                                                  self.root / "calls.tsv", "0460")
        self.assertEqual(annotations["0460"]["grouping_status"], "ineligible_contig_representative")
        self.assertEqual(annotations["NA"]["grouping_status"], "eligible")

    def test_member_fields_take_the_member_names_and_makedb_fields(self):
        annotations, fields = self.member_export()
        self.assertIn("member_consensus_count", fields)
        self.assertIn("member_representative_sequence_id", fields)
        self.assertNotIn("inference_consensus_count", fields)
        self.assertEqual(annotations["0460"]["is_member_representative"], "T")
        self.assertEqual([annotations["0460"]["isotype_class"], annotations["NA"]["isotype_class"]], ["IgM", ""])
        # MakeDb writes the regions IMGT-gapped, with deletions and N kept. The
        # TRB record's FWR1 is empty and stays empty. Its CDR1 and CDR2 hold
        # only '.' and '-', no sequenced base, so they are written empty, while
        # its FWR2 keeps a lone N.
        regions = {"0460": ["CAG...GTN", "GGA.TTC", "ATG-AGC", "ATT.AGT", "TACTGT", "TGGGGC"],
                   "NA": ["", "......", "..N", "-.", "TGT", "TTCGGC"]}
        written = {"0460": regions["0460"], "NA": ["", "", "..N", "", "TGT", "TTCGGC"]}
        makedb = self.root / "makedb.tsv.gz"
        with gzip.open(makedb, "wt", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sequence_id", *GERMLINE_FIELDS, *REGION_FIELDS, "cdr3"])
            writer.writerow(["0460", "3", "", "1", "12", "4", "48", *regions["0460"], "GCC"])
            writer.writerow(["NA", "5", "", "", "", "2", "40", *regions["NA"], "GCC"])
        fields = add_makedb_fields(annotations, fields, makedb)
        self.assertEqual(fields[-len(GERMLINE_FIELDS) - len(REGION_FIELDS):],
                         [*GERMLINE_FIELDS, *REGION_FIELDS])
        self.assertEqual(annotations["0460"]["j_germline_end"], "48")
        output = self.root / self.files["0460"]
        export_library_airr(self.source("0460"), output, "0460", self.processing,
                            annotations=annotations, annotation_fields=fields)
        self.assertEqual(validate_export_content(self.source("0460"), output, "0460", self.processing,
                                                annotations=annotations, annotation_fields=fields), 2)
        exported, rows = self.read_export()
        # The regions are written in V-domain order after the alignments. MakeDb's
        # germline coordinates follow each segment's query coordinates, and
        # np1_length and np2_length follow the J coordinates.
        start = exported.index("germline_alignment")
        self.assertEqual(exported[start:start + 22],
                         ["germline_alignment", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "fwr4",
                          "v_cigar", "d_cigar", "j_cigar", "v_sequence_start", "v_sequence_end",
                          "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end",
                          "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end",
                          "np1_length", "np2_length"])
        for row, key in zip(rows, ("0460", "NA")):
            self.assertEqual([row[column] for column in REGION_FIELDS], written[key])
        with gzip.open(makedb, "wt", newline="") as handle:
            handle.write("\t".join(["sequence_id", *GERMLINE_FIELDS, *REGION_FIELDS]) + "\n")
        with self.assertRaisesRegex(ValueError, "does not cover every post-QC record"):
            add_makedb_fields(self.member_export()[0], [], makedb)
        with gzip.open(makedb, "wt", newline="") as handle:
            handle.write("\t".join(["sequence_id", *GERMLINE_FIELDS]) + "\n0460" + "\t" * len(GERMLINE_FIELDS) + "\n")
        with self.assertRaisesRegex(ValueError, "lacks germline coordinates or V-domain regions"):
            add_makedb_fields(self.member_export()[0], [], makedb)
        with gzip.open(makedb, "wt", newline="") as handle:
            handle.write("sequence_id\tv_cigar\td_cigar\tj_cigar\n")
        self.assertEqual(add_makedb_fields({}, ["sequence_id"], makedb),
                         ["sequence_id", *GERMLINE_FIELDS, *REGION_FIELDS])

    def test_export_order_is_the_dictionary_order(self):
        with (ROOT / "docs/data_dictionary.csv").open(newline="") as handle:
            declared = [row["column"] for row in csv.DictReader(handle) if row["table"] == "exported_AIRR"]
        self.assertEqual(list(EXPORT_ORDER), declared)

    def test_export_order_skips_absent_columns_and_keeps_others_in_source_order_after_it(self):
        fields = [field for field in FIELDS if field != "cdr3"] + ["zeta_note", "alpha_note"]
        row = {key: value for key, value in source_row("0460").items() if key in fields}
        self.write_source("0460", [{**row, "zeta_note": "z", "alpha_note": "a"}], fields)
        output = self.root / self.files["0460"]
        annotations = {"0460": {"fwr1": ""}}
        export_library_airr(self.source("0460"), output, "0460", self.processing,
                            annotations=annotations, annotation_fields=["fwr1"])
        self.assertEqual(validate_export_content(self.source("0460"), output, "0460", self.processing,
                                                annotations=annotations, annotation_fields=["fwr1"]), 1)
        exported, rows = self.read_export()
        listed = [column for column in UNANNOTATED_ORDER if column != "cdr3"]
        listed.insert(listed.index("v_cigar"), "fwr1")
        self.assertEqual(exported, [*listed, "zeta_note", "alpha_note"])
        self.assertEqual((rows[0]["zeta_note"], rows[0]["alpha_note"]), ("z", "a"))

    def test_scope_copy_names_each_member_clone_and_nothing_else(self):
        annotations, fields = self.member_export()
        shared = self.root / "shared.tsv.gz"
        export_library_airr(self.source("0460"), shared, "0460", self.processing,
                            annotations=annotations, annotation_fields=fields)
        scoped = self.root / self.files["0460"]
        assignments = {"0460~0460": "project_b~NA~1", "0460~NA": "0460_TRB_ct00001"}
        self.assertEqual(annotate_member_assignments(shared, scoped, assignments), 2)
        exported, rows = self.read_export()
        # clone_id takes its listed place in the scope copy and only the copy.
        with gzip.open(shared, "rt", newline="") as handle:
            without = next(csv.reader(handle, delimiter="\t"))
        self.assertNotIn("clone_id", without)
        position = without.index("member_consensus_count") + 1
        self.assertEqual(exported, [*without[:position], "clone_id", *without[position:]])
        self.assertEqual([row["clone_id"] for row in rows], ["project_b~NA~1", "0460_TRB_ct00001"])
        self.assertEqual(validate_library_airr(scoped, "0460", self.processing, member_assignments=True), 2)
        with self.assertRaisesRegex(ValueError, "must not contain clone assignments"):
            validate_library_airr(scoped, "0460", self.processing)
        with self.assertRaisesRegex(ValueError, "differs from the member flag"):
            annotate_member_assignments(shared, self.root / "partial.tsv.gz", {"0460~0460": "project_b~NA~1"})

    def scope_copy(self, library_id, records):
        """Write one library's scope copy from (ID, locus, member link, source cells, flags) records.

        An empty member link marks a record outside grouping. Flags name
        contig_representative, member and member_consensus_count.
        """
        fields = FIELDS + ["locus", "germline_alignment_d_mask"]
        rows, annotations = [], {}
        for native_id, locus, link, source, flags in records:
            qualified = library_id + "~" + link if link else ""
            rows.append(source_row(library_id, native_id, **{
                "locus": locus, "germline_alignment": "ACTG", "germline_alignment_d_mask": "ACTG",
                "junction_aa": "CARW", **source}))
            annotations[native_id] = {
                "contig_id": library_id + "~" + native_id, "contig_representative_sequence_id": qualified,
                "is_contig_representative": "T" if flags.get("contig_representative") else "F",
                "member_representative_sequence_id": qualified,
                "is_member_representative": "T" if flags.get("member") else "F",
                "grouping_eligible": "T" if link else "F", "grouping_status": "eligible" if link else "ineligible_junction",
                "member_consensus_count": flags.get("member_consensus_count", "0")}
        self.write_source(library_id, rows, fields)
        shared = self.root / (library_id + "_shared.tsv.gz")
        export_library_airr(self.source(library_id), shared, library_id, self.processing,
                            annotations=annotations, annotation_fields=list(next(iter(annotations.values()))))
        assignments = {library_id + "~" + native_id: "clone_" + native_id
                       for native_id, _, _, _, flags in records if flags.get("member")}
        annotate_member_assignments(shared, self.root / self.files[library_id], assignments)

    def test_receipt_counts_each_locus_semantics_from_the_scope_copies(self):
        self.scope_copy("0460", [
            # A member with two contig sources, whose support exceeds its own record's.
            ("A", "IGH", "A", {"consensus_count": "10"},
             {"contig_representative": True, "member": True, "member_consensus_count": "15"}),
            ("B", "IGH", "A", {"consensus_count": "5"}, {"contig_representative": True}),
            # A contig sibling with blank germlines and a noncanonical first anchor.
            ("C", "IGH", "A", {"junction_aa": "GARW", "germline_alignment": "", "germline_alignment_d_mask": ""}, {}),
            # A member with a noncanonical last anchor and no transferred support.
            ("D", "IGH", "D", {"consensus_count": "3", "junction_aa": "CARR"},
             {"contig_representative": True, "member": True, "member_consensus_count": "3"}),
            # A record outside grouping, at another locus, with a canonical F anchor.
            ("E", "TRB", "", {"junction_aa": "CASSF", "germline_alignment_d_mask": ""}, {}),
        ])
        self.scope_copy("NA", [("A", "TRB", "A", {"consensus_count": "4"},
                                {"contig_representative": True, "member": True, "member_consensus_count": "4"})])
        report = validate_scope_repertoires(self.metadata(), self.root, self.libraries[:2], self.files,
                                            self.processing, member_assignments=True)
        self.assertEqual(report["linked_records"], {"0460": 5, "NA": 1})
        checks = report["semantic_checks"]
        zero = dict.fromkeys(SEMANTIC_COUNTS, 0)
        self.assertEqual(checks["by_locus"], {
            **{locus: zero for locus in ("IGK", "IGL", "TRA", "TRD", "TRG")},
            "IGH": {"n_records": 4, "n_member_representatives": 2, "n_members_with_several_contig_sources": 1,
                    "n_members_with_support_above_source_record": 1, "n_blank_germline_alignment": 1,
                    "n_blank_germline_alignment_d_mask": 1, "n_noncanonical_junction_aa_records": 2,
                    "n_noncanonical_junction_aa_members": 1},
            "TRB": {**zero, "n_records": 2, "n_member_representatives": 1, "n_blank_germline_alignment_d_mask": 1}})
        self.assertIn("ungapped selected reference", checks["v_germline_coordinates"])
        self.assertIn("separately", checks["cigar_query_spans_scope"])
        self.assertIn("^C.*[FW]$", checks["noncanonical_junction_aa"])
        self.assertIn("no record is filtered", checks["noncanonical_junction_aa"])
        # The receipt is JSON, and its counts stay integers there.
        written = json.loads(json.dumps(report, allow_nan=False))
        self.assertEqual(written, report)
        check_semantic_counts(written)
        for altered, message in (({"linked_records": {"0460": 5, "NA": 2}}, "does not sum"),
                                 ({"semantic_checks": {"by_locus": {"IGH": zero}}}, "lacks semantic counts"),
                                 ({"semantic_checks": {"by_locus": {**checks["by_locus"], "IGK": {"n_records": 0}}}},
                                  "incomplete at IGK")):
            with self.subTest(message=message), self.assertRaisesRegex(AssertionError, message):
                check_semantic_counts({**written, **altered})

    def test_shared_exports_receive_no_semantic_counts(self):
        self.assertNotIn("semantic_checks", self.validate(self.metadata()))

    def test_scope_metadata_contains_only_scope_members_with_literal_ids(self):
        # The manifest has already excluded this library. Its shared file need
        # not exist, even if a global export map still contains the path.
        (self.root / self.files["outside"]).unlink()
        data = self.metadata()
        self.assertEqual(data["scope_id"], "project_b")
        self.assertEqual([record["repertoire_name"] for record in data["Repertoire"]], ["0460", "NA"])
        self.assertNotIn("outside", json.dumps(data))
        for record in data["Repertoire"]:
            self.assertEqual(record["subject"]["subject_id"], "NA")
            sample = record["sample"][0]
            self.assertEqual(sample["sample_id"], "001")
            self.assertEqual(sample["template_class"], "RNA")
            self.assertFalse(sample["single_cell"])
            self.assertEqual(sample["physical_linkage"], "none")
            self.assertIsNone(sample["sequencing_run_id"])
            self.assertIsNone(sample["total_reads_passing_qc_filter"])
        self.assertEqual(self.validate(data)["linked_libraries"], 2)


    def test_missing_preparation_is_explicit_and_not_claimed_strict_valid(self):
        data = self.metadata()
        sample = data["Repertoire"][0]["sample"][0]
        self.assertIs(data["Repertoire"][0]["subject"]["synthetic"], False)
        self.assertIsNone(sample["library_generation_method"])
        self.assertIsNone(sample["library_generation_protocol"])
        self.assertIsNone(sample["complete_sequences"])
        self.assertNotIn("pcr_target", sample)
        report = self.validate(data)
        self.assertTrue(report["structural_valid"])
        self.assertFalse(report["strict_valid"])
        self.assertFalse(report["miairr_complete"])
        gaps = {gap["field"] for gap in report["strict_metadata_gaps"]}
        self.assertIn("sample[0].library_generation_method", gaps)
        self.assertIn("sample[0].complete_sequences", gaps)
        self.assertIn("sample[0].pcr_target", gaps)
        self.assertIn("study.contributors", gaps)
        self.assertIn("not supplied to this export", data["Info"]["description"])
        self.assertIn("do not denote negative findings", data["Info"]["description"])
        # Exercise the official public APIs as well as the version-matched Schema API.
        self.assertTrue(airr.validate_airr(data, check_nullable=False))
        with self.assertRaises(airr.ValidationError):
            airr.validate_airr(data, check_nullable=True)

    def test_unavailable_preparation_is_an_enum_constraint_not_free_text(self):
        data = self.metadata()
        data["Repertoire"][0]["sample"][0]["library_generation_method"] = "not recorded"
        with self.assertRaisesRegex(airr.ValidationError, "enumeration values"):
            self.validate(data)

    def test_official_optional_target_nullability_conflicts_with_non_targeted_omission(self):
        data = self.metadata()
        record = data["Repertoire"][0]
        record["study"]["contributors"] = []
        record["subject"]["diagnosis"] = []
        record["subject"]["synthetic"] = False
        sample = record["sample"][0]
        sample["library_generation_method"] = "RT(oligo-dT)+PCR"
        sample["complete_sequences"] = "partial"
        with self.assertRaisesRegex(airr.ValidationError, "pcr_target"):
            airr.schema.Schema("Repertoire").validate_object(record, nonairr=False)

    def test_project_annotations_do_not_enter_airr_repertoires(self):
        libraries = [{**row, "meta_group": "Private report annotation", "timepoint": "private_visit",
                      "timepoint_order": 27, "sex": "female", "synthetic": True,
                      "library_generation_method": "PCR"} for row in self.libraries[:2]]
        data = self.metadata(libraries=libraries)
        self.assertNotIn("Private report annotation", json.dumps(data))
        self.assertNotIn("private_visit", json.dumps(data))
        self.assertTrue(all("supplied_metadata" not in row for row in data["Repertoire"]))
        for record in data["Repertoire"]:
            self.assertIsNone(record["subject"]["sex"])
            self.assertIs(record["subject"]["synthetic"], False)
            self.assertIsNone(record["sample"][0]["library_generation_method"])
        self.assertTrue(self.validate(data)["structural_valid"])

    def test_unavailable_optional_metadata_is_omitted(self):
        data = self.metadata()
        for record in data["Repertoire"]:
            self.assertEqual(record["repertoire_type"], "observed")
            self.assertNotIn("diagnosis", record["subject"])
            self.assertNotIn("sequencing_files", record["sample"][0])
            self.assertIsNone(record["study"]["contributors"])
        report = self.validate(data)
        self.assertNotIn("required_metadata_gaps", report)
        gaps = {gap["field"] for gap in report["strict_metadata_gaps"]}
        self.assertIn("sample[0].sequencing_files", gaps)
        self.assertIn("subject.diagnosis", gaps)

    def test_processing_details_remain_recorded_values_and_ids_are_not_overridable(self):
        data = self.metadata(processing_fields={"software_versions": "IgBLAST 1.22.0",
                                                "germline_database": "recorded snapshot sha256:example"})
        processing = data["Repertoire"][0]["data_processing"][0]
        self.assertEqual(processing["software_versions"], "IgBLAST 1.22.0")
        self.assertIsNone(processing["paired_reads_assembly"])
        with self.assertRaisesRegex(ValueError, "declared AIRR DataProcessing"):
            self.metadata(processing_fields={"data_processing_id": "fake"})

    def test_cross_file_links_reject_wrong_donor_specimen_processing_and_repertoire(self):
        for section, field in [("subject", "subject_id"), ("sample", "sample_id"),
                               ("sample", "sample_processing_id"),
                               ("data_processing", "data_processing_id"), (None, "repertoire_id")]:
            with self.subTest(field=field):
                data = self.metadata()
                record = data["Repertoire"][0]
                target = record if section is None else record[section]
                if isinstance(target, list):
                    target = target[0]
                target[field] = "wrong"
                with self.assertRaisesRegex(ValueError, "linkage differs"):
                    self.validate(data)


    def test_wrong_export_row_identity_fails_even_when_official_schema_accepts_it(self):
        fields, rows = self.read_export()
        rows[0]["data_processing_id"] = "wrong"
        self.rewrite_export(fields, rows)
        with self.assertRaisesRegex(ValueError, "AIRR row linkage differs"):
            self.validate(self.metadata())


if __name__ == "__main__":
    unittest.main()
