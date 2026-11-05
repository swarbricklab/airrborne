"""V and J pair usage reproduces gene usage for each segment.

Run: python test/run_tests.py --component gene_usage --level unit
"""

import unittest

from test.components.release.integration.check_release import check_vj_pair_marginals


def pair(v_key, j_key, n_groups, support):
    return dict(library_id="L", locus="TRB", group_unit="clonotype", v_gene_key=v_key, v_gene_names=v_key,
                j_gene_key=j_key, j_gene_names=j_key, n_groups=str(n_groups), group_fraction=str(n_groups / 4),
                consensus_count=str(support), support_fraction=str(support / 40))


def usage(gene_type, key, n_groups, support):
    return dict(library_id="L", locus="TRB", group_unit="clonotype", gene_type=gene_type, gene_key=key, gene_names=key,
                n_groups=str(n_groups), group_fraction=str(n_groups / 4),
                consensus_count=str(support), support_fraction=str(support / 40))


class PairMarginalTests(unittest.TestCase):
    def setUp(self):
        self.pairs = [pair("V1", "J1", 2, 30), pair("V1", "J2", 1, 5), pair("V2", "J1", 1, 5)]
        self.usage = [usage("V", "V1", 3, 35), usage("V", "V2", 1, 5),
                      usage("J", "J1", 3, 35), usage("J", "J2", 1, 5)]

    def test_marginals_reproduce_both_segments(self):
        check_vj_pair_marginals(self.pairs, self.usage)

    def test_a_differing_marginal_fails(self):
        for column, value in (("n_groups", "3"), ("consensus_count", "31"), ("group_fraction", "0.6"),
                              ("support_fraction", "0.8"), ("group_unit", "clone")):
            with self.subTest(column=column), self.assertRaises(AssertionError):
                check_vj_pair_marginals([{**self.pairs[0], column: value}, *self.pairs[1:]], self.usage)

    def test_an_entity_in_one_table_only_fails(self):
        with self.assertRaises(AssertionError):
            check_vj_pair_marginals(self.pairs[:-1], self.usage)
        with self.assertRaises(AssertionError):
            check_vj_pair_marginals(self.pairs, [*self.usage, usage("V", "V3", 1, 1)])

    def test_one_entity_with_two_labels_fails(self):
        with self.assertRaises(AssertionError):
            check_vj_pair_marginals([{**self.pairs[1], "v_gene_names": "other"}, self.pairs[0], self.pairs[2]], self.usage)
