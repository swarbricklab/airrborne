"""Protect identifier-preserving imports and the documented join cardinalities."""
import tempfile
from pathlib import Path
import unittest

from examples.audit_vote_mapping import audit
from examples.read_release import add_metadata, annotate_vj_pairs, join_member_airr, number, read_rows


class ReadExamplesTests(unittest.TestCase):
    def test_import_and_metadata_preserve_literal_ids_and_unmatched_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "headline.tsv"
            path.write_text("library_id\ttissue_id\tinput_read_pairs\n0460\tNA\t\nNA\t001\t0\n")
            rows = list(read_rows(path))
        self.assertEqual([row["library_id"] for row in rows], ["0460", "NA"])
        self.assertIsNone(number(rows[0]["input_read_pairs"]))
        self.assertEqual(number(rows[1]["input_read_pairs"]), 0)
        metadata = [{"tissue_id": "NA", "group": "observed"}]
        joined = add_metadata(rows, metadata)
        self.assertEqual([row["group"] for row in joined], ["observed", ""])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            add_metadata(rows, metadata * 2)
        with self.assertRaisesRegex(ValueError, "conflicts"):
            add_metadata(rows, [{"tissue_id": "NA", "library_id": "wrong"}])

    def test_pair_lookup_uses_library_locus_and_segment_context(self):
        pairs = [dict(scope_id="cohort", library_id=library, locus="IGH", v_gene_key="V",
                      j_gene_key="J", consensus_count="10") for library in ("0460", "NA")]
        genes = [dict(scope_id="cohort", library_id=library, locus="IGH", gene_type=segment,
                      gene_key=segment, gene_equivalence_group=library + segment)
                 for library in ("0460", "NA") for segment in ("V", "J")]
        joined = annotate_vj_pairs(pairs, genes)
        self.assertEqual([row["v_gene_equivalence_group"] for row in joined], ["0460V", "NAV"])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            annotate_vj_pairs(pairs, genes + genes[:1])

    def test_member_join_never_counts_transferred_support_on_sibling_rows(self):
        members = [dict(scope_id="cohort", library_id="0460", sequence_id="0460~assemble0_0",
                        member_consensus_count="12")]
        airr = [dict(library_id="0460", sequence_id="0460~assemble0_0", sequence="ACTG",
                     is_member_representative="T", member_consensus_count="12", consensus_count="10"),
                dict(library_id="0460", sequence_id="0460~assemble0_1", sequence="ACTA",
                     is_member_representative="F", member_consensus_count="0", consensus_count="2")]
        joined = join_member_airr(members, airr)
        self.assertEqual(len(joined), 1)
        self.assertEqual(joined[0]["member_consensus_count"], "12")
        self.assertEqual(joined[0]["source_consensus_count"], "10")


class VoteAuditTests(unittest.TestCase):
    def test_counts_sources_and_rejects_inconsistent_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "libraries/0460/airr/0460_igh_consensus_votes.tsv"
            path.parent.mkdir(parents=True)
            header = "sequence_id\timgt_v_positions_mapped\tjunction_vote_status\n"
            path.write_text(header + "a\t0\topposite_contig_orientation\nb\t23\tmapped\n")
            result = audit(directory)
            self.assertEqual(result["counts"], {"igh_source_records": 2, "opposite_contig_orientation": 1})
            self.assertEqual(result["libraries"][0]["library_id"], "0460")
            path.write_text(header + "a\t23\topposite_contig_orientation\n")
            with self.assertRaisesRegex(ValueError, "retained mapped"):
                audit(directory)
