"""Check independent SHM evidence selection without pooling partial sources.
Constructed cases anchor the independent system-result calculations.

Run: python test/run_tests.py --component shm --level unit
"""

import unittest
import csv
import json
from pathlib import Path
from test.components.shm.integration.reference_calculations import REGIONS, sequence_shm_evidence


class SequenceEvidenceReferenceTests(unittest.TestCase):
    def setUp(self):
        self.settings = {"min_supported_subregion_fraction": 1,
                         "min_germline_occupied_v_positions": 50}
        self.records = {(identity,): dict(sequence="ACGT" if identity != "B_1" else "TGCA",
                                         consensus_count=str(count))
                        for identity, count in (("A_0", 100), ("B_0", 50), ("B_1", 1000))}
        self.mapping = {key: dict(is_backbone_representative="FALSE" if key == ("B_1",) else "TRUE",
                                 inference_representative_sequence_id="A_0")
                        for key in self.records}
        feature = dict(qualified_c_call="IGHM*01", shm_frequency="0.1", informative_v="50")
        for region in REGIONS:
            feature.update({f"shm_frequency_{region}": "0.1", f"mutations_{region}": "1",
                            f"informative_{region}": "10", f"germline_occupied_{region}": "10",
                            f"informative_supported_{region}": "10", f"median_majority_votes_{region}": "3"})
        self.features = {key: dict(feature) for key in self.records}

    def evidence(self):
        return sequence_shm_evidence("A_0", self.records, self.mapping, self.features, self.settings)

    def test_complementary_sources_cannot_make_whole_v_eligible(self):
        self.features[("A_0",)]["informative_supported_fwr1"] = "0"
        self.features[("B_0",)]["informative_supported_cdr1"] = "0"
        selected = self.evidence()
        # Neither source passes, and their supported totals tie at 40, so the supported
        # count is that of A_0, the source with the larger original record support.
        self.assertEqual(selected["whole_v"], ("", (.1, False, 5, None), 40))
        self.assertEqual(selected["fwr1"][0], "B_0")
        self.assertEqual(selected["cdr1"][0], "A_0")
        # The discarded B_1 sibling is fully covered but has a different DNA
        # sequence. Its larger support supplies neither a vote nor eligibility.
        self.assertNotIn("B_1", [source for source, *_ in selected.values()])

    def test_supported_count_follows_the_source_ranked_first(self):
        self.features[("A_0",)]["informative_supported_fwr1"] = "0"
        self.features[("B_0",)]["informative_supported_fwr1"] = "10"
        self.features[("B_0",)]["informative_supported_cdr1"] = "4"
        selected = self.evidence()
        # FWR1 selects B_0 and its 10 positions over A_0's zero, CDR1 selects A_0,
        # and whole V ranks B_0's 44 above A_0's 40 with neither source passing.
        self.assertEqual((selected["fwr1"][0], selected["fwr1"][2]), ("B_0", 10))
        self.assertEqual((selected["cdr1"][0], selected["cdr1"][2]), ("A_0", 10))
        self.assertEqual((selected["whole_v"][0], selected["whole_v"][2]), ("", 44))
        # A source with no count for a subregion has no total and ranks last.
        self.features[("B_0",)]["informative_supported_cdr1"] = ""
        self.assertEqual(self.evidence()["whole_v"][2], 40)

    def test_coverage_then_original_support_then_literal_id_choose_evidence(self):
        self.settings["min_supported_subregion_fraction"] = .8
        self.features[("A_0",)]["informative_supported_cdr2"] = "8"
        self.assertEqual(self.evidence()["cdr2"][0], "B_0")
        self.features[("A_0",)]["informative_supported_cdr2"] = "10"
        self.assertEqual(self.evidence()["cdr2"][0], "A_0")
        self.records[("B_0",)]["consensus_count"] = "100"
        self.mapping = dict(reversed(list(self.mapping.items())))
        self.assertEqual(self.evidence()["cdr2"][0], "A_0")

    def test_candidate_dna_and_sequence_derived_features_must_agree(self):
        self.records[("B_0",)]["sequence"] = "TGCA"
        self.assertRaisesRegex(AssertionError, "evidence DNA", self.evidence)
        self.records[("B_0",)]["sequence"] = "ACGT"
        self.features[("B_0",)]["qualified_c_call"] = "IGHG1*01"
        self.assertRaisesRegex(AssertionError, "sequence-derived features", self.evidence)
