"""Check that available selection results carry estimates and intervals.
Constructed cases anchor the independent system-result calculations.

Run: python test/run_tests.py --component selection --level unit
"""

import unittest
import csv
import json
from pathlib import Path
from test.components.selection.integration.reference_calculations import verify_selection_estimates


class SelectionAvailabilityTests(unittest.TestCase):
    def test_success_requires_measured_estimates_and_intervals(self):
        row = dict(library_id="L1", isotype_group="all_IGH", n_clone_consensus_sequences="2",
                   n_clone_consensus_sequences_used_cdr="1", sigma_cdr="1", sigma_cdr_ci_low="-1",
                   sigma_cdr_ci_high="2", sigma_cdr_status="estimated", n_clone_consensus_sequences_used_fwr="0",
                   sigma_fwr="", sigma_fwr_ci_low="", sigma_fwr_ci_high="",
                   sigma_fwr_status="no_usable_mutation_pdfs")
        self.assertEqual(verify_selection_estimates([row]), 1)
        for column, value in (("sigma_cdr", ""), ("sigma_cdr", "inf"),
                              ("sigma_cdr_ci_high", "0"), ("n_clone_consensus_sequences_used_cdr", "0"),
                              ("sigma_fwr", "0")):
            with self.subTest(column=column, value=value):
                self.assertRaises(AssertionError, verify_selection_estimates,
                                  [{**row, column: value}])
