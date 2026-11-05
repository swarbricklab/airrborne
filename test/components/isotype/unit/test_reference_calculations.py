"""Check independent isotype reference arithmetic on callable populations.
Constructed cases anchor the independent system-result calculations.

Run: python test/run_tests.py --component isotype --level unit
"""

import unittest

from test.components.isotype.integration.reference_calculations import resolve_call


class AnalyticalReferenceTests(unittest.TestCase):
    def test_callable_populations(self):
            self.assertEqual(resolve_call("IGHG1*01,IGHG2*02"), ("IgG", None))
            self.assertEqual(resolve_call("IGHG1*01,IGHA1*01"), (None, None))
            self.assertEqual(resolve_call("IGHG4A*01"), ("IgG", "IgG4"))
            self.assertEqual(resolve_call("IGHG1*01,IGHGP*01"), (None, None))
