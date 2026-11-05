#!/usr/bin/env python3
"""Check editable scientific thresholds against the configuration schema.

Run: python test/run_tests.py --component inputs --level unit
"""

import copy
from pathlib import Path
import unittest

import jsonschema
import yaml


ROOT = Path(__file__).resolve().parents[4]
SCHEMA = yaml.safe_load((ROOT / "workflow/schemas/config.schema.yaml").read_text())
EXAMPLE = yaml.safe_load((ROOT / "config/airrborne/config.yaml").read_text())


class ConfigurationThresholds(unittest.TestCase):
    def configuration(self, family, key, value):
        config = copy.deepcopy(EXAMPLE)
        config["params"][family][key] = value
        return config

    def test_examples_use_valid_matching_scientific_defaults(self):
        fixture = yaml.safe_load((ROOT / "test/system/config.yaml").read_text())
        for config in (EXAMPLE, fixture):
            jsonschema.validate(config, SCHEMA)
        # Rarefaction targets differ for the fixture; these families match.
        for family in ("trust4", "lineage", "igblast", "airr_qc", "shm", "selection", "overlap"):
            self.assertEqual(EXAMPLE["params"][family], fixture["params"][family])

    def test_empirical_floors_can_be_lowered(self):
        cases = (
            ("lineage", "min_fit_distances", 2),
            ("lineage", "fit_subsample", None),
            ("lineage", "max_n", 2),
            ("airr_qc", "min_informative_c_alignment_positions", 0),
            ("shm", "min_position_votes", 0),
            ("shm", "min_supported_subregion_fraction", 0),
            ("shm", "min_germline_occupied_v_positions", 0),
            ("shm", "min_germline_occupied_v_positions", 312),
        )
        for family, key, value in cases:
            with self.subTest(family=family, key=key, value=value):
                jsonschema.validate(self.configuration(family, key, value), SCHEMA)

    def test_selection_frequency_requires_a_positive_fraction(self):
        for value in (1e-12, 0.5, 0.6, 1):
            with self.subTest(value=value):
                jsonschema.validate(
                    self.configuration("selection", "collapse_min_freq", value), SCHEMA)
        for value in (0, -0.1, 1.01, None, True, "0.6"):
            with self.subTest(value=value), self.assertRaises(jsonschema.ValidationError):
                jsonschema.validate(
                    self.configuration("selection", "collapse_min_freq", value), SCHEMA)

    def test_pair_support_minima_allow_disabled_and_per_locus_values(self):
        for value in ({}, {"IGH": 0, "TRB": 100}, {"TRA": 1, "TRG": 0}):
            jsonschema.validate(self.configuration("overlap", "min_consensus_count", value), SCHEMA)
        for value in ({"TRB": -1}, {"IGH": 1.5}, {"IGH": None}, {"IGK": 10}, 100):
            with self.subTest(value=value), self.assertRaises(jsonschema.ValidationError):
                jsonschema.validate(self.configuration("overlap", "min_consensus_count", value), SCHEMA)
        properties = SCHEMA["properties"]["params"]["properties"]["overlap"]["properties"]
        self.assertTrue(all(value == 0 for value in properties["min_consensus_count"]["default"].values()))

    def test_algorithm_limits_are_rejected(self):
        cases = (
            ("lineage", "min_fit_distances", 1),
            ("lineage", "fit_subsample", 1),
            ("lineage", "fit_edge", 0),
            ("lineage", "fit_edge", 1),
            ("lineage", "fit_specificity", 1.01),
            ("lineage", "max_n", -1),
            ("igblast", "num_alignments_c", 0),
            ("shm", "min_position_votes", -1),
            ("shm", "min_supported_subregion_fraction", 1.01),
            ("shm", "min_germline_occupied_v_positions", 313),
        )
        for family, key, value in cases:
            with self.subTest(family=family, key=key, value=value):
                with self.assertRaises(jsonschema.ValidationError):
                    jsonschema.validate(self.configuration(family, key, value), SCHEMA)


if __name__ == "__main__":
    unittest.main()
