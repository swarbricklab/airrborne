"""Check independent diversity arithmetic against exact small-population cases.
Constructed cases anchor the independent system-result calculations.

Run: python test/run_tests.py --component diversity --level unit
"""

import unittest

from collections import Counter
from itertools import combinations
import statistics
import math
from test.components.diversity.integration import reference_calculations as DIVERSITY


class AnalyticalReferenceTests(unittest.TestCase):
    def test_exact_enumeration(self):
            for labels in ((0, 0, 1), (0, 0, 1, 1, 2)):
                counts = tuple(Counter(labels).values())
                for depth in range(1, len(labels) + 1):
                    samples = [Counter(labels[i] for i in positions) for positions in combinations(range(len(labels)), depth)]
                    richness = statistics.mean(len(sample) for sample in samples)
                    entropy = statistics.mean(-math.fsum(n / depth * math.log(n / depth) for n in sample.values()) for sample in samples)
                    concentration = statistics.mean(math.fsum((n / depth) ** 2 for n in sample.values()) for sample in samples)
                    actual = DIVERSITY.interpolate(counts, depth)
                    self.assertAlmostEqual(actual["richness"], richness, places=12)
                    self.assertAlmostEqual(actual["hill_shannon"], math.exp(entropy), places=12)
                    self.assertAlmostEqual(actual["inverse_simpson"], 1 / concentration, places=12)

    def test_empty_singleton_minimum_and_target_statuses(self):
            empty = DIVERSITY.expected_detail((), 0, 100)
            self.assertEqual(empty["rarefaction_status"], "empty_repertoire")
            self.assertIsNone(empty["inverse_simpson"])
            single = DIVERSITY.expected_detail((3,), 0, 2)
            self.assertEqual(single["top_10_support_fraction"], 1)
            self.assertIsNone(single["clonality"])
            below = DIVERSITY.expected_detail((2, 1), 4, 2)
            self.assertEqual(below["grouped_consensus_count"], 3)
            self.assertEqual(below["richness"], 2)
            self.assertEqual(below["rarefaction_status"], "below_minimum_support")
            self.assertIsNone(below["rarefied_inverse_simpson"])
            self.assertEqual(DIVERSITY.expected_detail((2, 1), 0, 4)["rarefaction_status"], "below_target_support")
            self.assertEqual(DIVERSITY.expected_detail((), 0, None)["rarefaction_status"], "not_configured")

    def test_tied_support_shares_the_smaller_rank(self):
            def vector(*rows):
                return [(count, {"group_key": key, "group_rank": str(rank)}) for key, count, rank in rows]
            DIVERSITY.verify_ranks(vector(("b", 10, 1), ("a", 9, 2), ("c", 9, 2), ("d", 1, 4)), {}, "tied")
            with self.assertRaises(AssertionError):
                DIVERSITY.verify_ranks(vector(("b", 10, 1), ("a", 9, 2), ("c", 9, 3)), {}, "dense")
            with self.assertRaises(AssertionError):
                DIVERSITY.verify_ranks(vector(("b", 10, 1), ("c", 9, 2), ("a", 9, 2)), {}, "key order")
            # Tied IGH clones follow the smallest member junction before clone_id.
            junctions = {"c": "TGTAAA", "a": "TGTCCC"}
            DIVERSITY.verify_ranks(vector(("c", 9, 1), ("a", 9, 1)), junctions, "junction order")
