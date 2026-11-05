"""A fitted donor's released mixture reproduces its histogram's component densities.

The fit and bin cells are those the pinned analysis image writes for the separated
two-mode fixture of test/support/helpers.R, with distances from 1 to 18 mismatches
over a 60 nt junction. They span the fitted peaks and both far tails.

Run: python test/run_tests.py --component release --level unit
"""

import unittest

from test.components.release.integration.check_release import check_threshold_components

FIT = dict(scope_id="fit_scope", donor_id="fit_donor", n_libraries="1", n_members="600", n_distances="600",
           n_distances_fitted="600", fit_status="fitted", fitted_threshold="0.205151255137562",
           near_weight="0.502490377207745", near_shape="4.18334126708814", near_scale="0.0100277439765795",
           far_weight="0.497509622792255", far_mean="0.266621028216678", far_sd="0.0236055525204918",
           absolute_log_likelihood="1054.2", fit_sensitivity="0.999996767503661")
PARAMETERS = ("near_weight", "near_shape", "near_scale", "far_weight", "far_mean", "far_sd",
              "absolute_log_likelihood", "fit_sensitivity")


def histogram_bin(midpoint, near, far, donor_id="fit_donor"):
    return dict(donor_id=donor_id, distance_bin_midpoint=midpoint, near_component_density=near,
                far_component_density=far)


BINS = [histogram_bin("0.005", "0.437583744664219", "1.78539295591612e-26"),
        histogram_bin("0.125", "0.0783549395003468", "0.000000128455579697292"),
        histogram_bin("0.245", "0.00000423937927884471", "5.52747292744109"),
        histogram_bin("0.485", "0.00000000000000150374415733437", "2.18946198443868e-18"),
        histogram_bin("0.965", "2.18700978182958e-35", "7.19470367181696e-190")]
UNFITTED = {**FIT, "donor_id": "sparse", "n_distances_fitted": "", "fit_status": "too_few_distances",
            "fitted_threshold": "", **dict.fromkeys(PARAMETERS, "")}
UNFITTED_BINS = [histogram_bin("0.005", "", "", "sparse")]


class ThresholdComponentTests(unittest.TestCase):
    def test_the_written_mixture_reproduces_the_written_densities(self):
        check_threshold_components([FIT, UNFITTED], [*BINS, *UNFITTED_BINS])
        check_threshold_components([], [])

    def test_a_changed_parameter_fails(self):
        for column in PARAMETERS[:6]:
            changed = {**FIT, column: str(float(FIT[column]) * (1 + 1e-6))}
            with self.subTest(column=column), self.assertRaisesRegex(AssertionError, "does not reproduce"):
                check_threshold_components([changed], BINS)

    def test_a_changed_or_empty_density_fails(self):
        for column in ("near_component_density", "far_component_density"):
            for value in (str(float(BINS[0][column]) * (1 + 1e-6)), ""):
                with self.subTest(column=column, value=value), self.assertRaisesRegex(AssertionError, "does not reproduce"):
                    check_threshold_components([FIT], [{**BINS[0], column: value}, *BINS[1:]])

    def test_a_tolerance_of_one_part_in_a_billion_applies(self):
        near = float(BINS[0]["near_component_density"])
        check_threshold_components([FIT], [{**BINS[0], "near_component_density": repr(near * (1 + 1e-10))}, *BINS[1:]])
        with self.assertRaisesRegex(AssertionError, "does not reproduce"):
            check_threshold_components([FIT], [{**BINS[0], "near_component_density": repr(near * (1 + 1e-8))}, *BINS[1:]])

    def test_subnormal_values_compare_by_their_difference(self):
        # At this standard deviation the far component underflows to zero at the midpoint.
        narrow = {**FIT, "far_sd": "0.018"}
        midpoint = histogram_bin("0.965", BINS[4]["near_component_density"], "1e-310")
        check_threshold_components([narrow], [midpoint])
        with self.assertRaisesRegex(AssertionError, "does not reproduce"):
            check_threshold_components([narrow], [{**midpoint, "far_component_density": "1e-300"}])

    def test_a_donor_without_a_fit_holds_no_mixture(self):
        for column in PARAMETERS:
            with self.subTest(column=column), self.assertRaisesRegex(AssertionError, "holds mixture parameters"):
                check_threshold_components([{**UNFITTED, column: "1"}], UNFITTED_BINS)
        with self.assertRaisesRegex(AssertionError, "holds component densities"):
            check_threshold_components([UNFITTED], [{**UNFITTED_BINS[0], "far_component_density": "1"}])

    def test_a_fitted_donor_lacking_a_parameter_fails(self):
        for column in PARAMETERS:
            with self.subTest(column=column), self.assertRaisesRegex(AssertionError, "lacks mixture parameters"):
                check_threshold_components([{**FIT, column: ""}], BINS)

    def test_tables_naming_other_donors_fail(self):
        with self.assertRaisesRegex(AssertionError, "name other donors"):
            check_threshold_components([FIT], UNFITTED_BINS)


if __name__ == "__main__":
    unittest.main()
