/* Display transforms preserve supplied measurements and biological identities. */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const selection = require("../../../../workflow/report/assets/selection_controller.js");
const {SelectionController, derivePairMeasures} = selection;
global.window = {ReportSelection: selection};
require("../../../../workflow/report/assets/report_view.js");
const {axis, format, orderedComparisons, violinSummary} = window.ReportView;

function fixture() {
    return {
        libraries: [
            {library_id: "0460", donor_id: "d1", tissue_id: "t1", group: "a", score: 0, diversity: 5, shm: 0, timepoint_order: 1},
            {library_id: "NA", donor_id: "d1", tissue_id: "t1", group: "a", score: null, diversity: null, diversity_status: "below_target", shm: 0.1, timepoint_order: 1},
            {library_id: "third", donor_id: "d1", tissue_id: "t3", group: "b", score: 5, diversity: 10, shm: null, timepoint_order: 3},
            {library_id: "other", donor_id: "d2", tissue_id: "t4", group: null, score: 9, diversity: 1, shm: 0.2, timepoint_order: null},
        ],
        pairs: [
            {pair_id: "p1", library_id_a: "0460", library_id_b: "NA"},
            {pair_id: "p2", library_id_a: "NA", library_id_b: "third"},
            {pair_id: "p3", library_id_a: "third", library_id_b: "0460"},
        ],
        default_context: {locus: "IGH", group_field: "group", y_metric: "diversity"},
    };
}

// Released pair columns only; the report derives the rest before any pair is
// oriented. Simpson concentration is 1 / inverse_simpson: 0.5 for 0460, 0.4
// for NA and 0.25 for third.
function longitudinalFixture() {
    const data = fixture();
    for (const pair of data.pairs) Object.assign(pair, {
        locus: "IGH", donor_id: "d1", library_id_1: pair.library_id_a, library_id_2: pair.library_id_b,
        tissue_id_1: "later_tissue", tissue_id_2: "earlier_tissue",
        pair_status: "below_minimum_support", min_consensus_count_per_library: 120,
        n_groups_1: 10, n_groups_2: 20,
        grouped_consensus_count_1: 200, grouped_consensus_count_2: 100,
        n_groups_shared: 5, shared_consensus_count_1: 80, shared_consensus_count_2: 60,
        shared_support_fraction_1: .4, shared_support_fraction_2: .6,
        bray_curtis_proportion: .45, n_shared_higher_in_1: 4, n_shared_higher_in_2: 1, n_shared_equal: 0,
    });
    Object.assign(data.pairs[1], {grouped_consensus_count_2: 150, pair_status: "available"});
    data.diversity = [["0460", 2], ["NA", 2.5], ["third", 4]].map(([library_id, inverse_simpson]) =>
        ({library_id, locus: "IGH", inverse_simpson, diversity_status: "available"}));
    data.pairs = derivePairMeasures(data.pairs, data.diversity);
    return data;
}

const close = (actual, expected, context) => assert.ok(Math.abs(actual - expected) < 1e-12, `${context}: ${actual} != ${expected}`);


test("reference comparisons retain repeated libraries and orient reversed source pairs", () => {
    const data = fixture();
    data.libraries.forEach(row => { row.timepoint = row.timepoint_order === 1 ? "Before" : row.timepoint_order === 3 ? "After" : null; });
    const controller = new SelectionController(data);
    assert.deepEqual(controller.orderedTimepoints(), [{label: "Before", order: 1}, {label: "After", order: 3}]);
    const pairs = controller.referencePairs("Before");
    assert.deepEqual(pairs.map(row => [row.pair_id, row.reference_id, row.followup_id]), [["p2", "NA", "third"], ["p3", "0460", "third"]]);
    assert.equal(pairs[1].library_id_a, "third", "Source A/B abundance denominators retain their order.");
    assert.deepEqual(controller.referencePairs("After"), []);
    controller.setCategoricalFilter("group", ["a"]);
    assert.deepEqual(controller.referencePairs("Before"), [], "Both libraries must satisfy the filters.");
});

test("reference visits require consistent labels and do not impute missing orders", () => {
    const data = fixture();
    data.libraries.forEach(row => { row.timepoint = row.timepoint_order === 1 ? "Before" : row.timepoint_order === 3 ? "After" : null; });
    data.libraries[0].timepoint_order = null;
    let controller = new SelectionController(data);
    assert.deepEqual(controller.referencePairs("Before").map(row => row.pair_id), ["p2"]);
    data.libraries[0].timepoint_order = 2;
    controller = new SelectionController(data);
    assert.deepEqual(controller.orderedTimepoints(), [{label: "After", order: 3}]);
    assert.deepEqual(controller.referencePairs("Before"), []);
    data.libraries[0].timepoint_order = 1;
    data.libraries[2].timepoint_order = 1;
    assert.deepEqual(new SelectionController(data).orderedTimepoints(), [], "Distinct labels with the same order are ambiguous.");
});

test("missing metrics remain selected and zero is measured", () => {
    const controller = new SelectionController(fixture());
    controller.selectLibraries(["0460", "NA", "third"]);
    const shm = controller.availability(["shm"], controller.selectedRows());
    assert.equal(shm.eligible.libraries, 2);
    assert.equal(shm.zeros, 1);
    assert.deepEqual(shm.missing_library_ids, ["third"]);
    const joint = controller.availability(["diversity", "shm"], controller.selectedRows());
    assert.equal(joint.eligible.libraries, 1);
    assert.deepEqual(joint.missing_library_ids, ["NA", "third"]);
});

test("display orientation swaps directional measurements without changing released identities", () => {
    const data = longitudinalFixture(), original = JSON.stringify(data.pairs);
    const controller = new SelectionController(data), row = controller.displayPair("p3");
    assert.equal(row.library_id_a, "0460");
    assert.equal(row.library_id_b, "third");
    assert.equal(row.library_id_1, "0460");
    assert.equal(row.grouped_consensus_count_1, 100);
    assert.equal(row.grouped_consensus_count_2, 200);
    assert.equal(row.n_groups_1, 20);
    assert.equal(row.shared_consensus_count_1, 60);
    assert.equal(row.shared_support_fraction_1, .6);
    // Canonical third to 0460 is 0.5 - 0.25; 0460 to third reverses it.
    close(row.simpson_concentration_difference, -.25, "Simpson difference");
    close(row.shared_difference_balance, .6, "balance");
    assert.equal(row.n_shared_higher_in_2, 4);
    assert.equal(row.n_shared_higher_in_1, 1);
    close(row.shared_higher_in_2_fraction, .8, "higher in 2");
    close(row.shared_higher_in_1_fraction, .2, "higher in 1");
    close(row.support_fraction_only_in_2, .6, "only in 2");
    close(row.support_fraction_only_in_1, .4, "only in 1");
    assert.equal(row.bray_curtis_proportion, .45);
    assert.equal(row.shared_equal_fraction, 0);
    assert.equal(row.pair_id, "p3");
    assert.equal(row.has_visit_order, true);
    assert.equal(JSON.stringify(data.pairs), original);
    assert.equal(controller.displayPair("p1").has_visit_order, false);
    assert.equal(controller.displayPair("p1").library_id_a, "0460");
});

test("reversing empty-endpoint status preserves unavailable values and genuine zeros", () => {
    const data = longitudinalFixture();
    Object.assign(data.pairs[2], {grouped_consensus_count_1: 0, pair_status: "empty_first",
        shared_support_fraction_1: null, bray_curtis_proportion: null, n_shared_higher_in_1: null,
        n_shared_higher_in_2: null, n_shared_equal: null, n_groups_shared: 0});
    data.pairs = derivePairMeasures(data.pairs, data.diversity);
    const controller = new SelectionController(data);
    const pair = controller.pairs.get("p3");
    for (const field of ["overlap_status", "simpson_difference_status", "shared_difference_status"]) assert.equal(pair[field], "empty_first");
    const row = controller.displayPair("p3");
    for (const field of ["overlap_status", "pair_status", "simpson_difference_status", "shared_difference_status"]) assert.equal(row[field], "empty_second");
    for (const field of ["simpson_concentration_difference", "shared_difference_balance", "support_fraction_only_in_1"]) assert.equal(row[field], null);
    pair.simpson_concentration_difference = 0;
    assert.equal(controller.displayPair("p3").simpson_concentration_difference, 0);
});

// The reviewed fixture vectors: IGH support [1392, 480] in SYN_D1_T1 and
// [1255, 446, 60] in SYN_D1_T2, with only the first group shared and higher in
// T1. The pipeline released these four values before the report derived them.
test("derived pair measurements reproduce the reviewed released values", () => {
    const [p, q] = [[1392, 480], [1255, 446, 60]];
    const inverse = counts => { const total = counts.reduce((a, b) => a + b, 0); return 1 / counts.reduce((sum, n) => sum + (n / total) ** 2, 0); };
    const diversity = [["SYN_D1_T1", p], ["SYN_D1_T2", q]].map(([library_id, counts]) =>
        ({library_id, locus: "IGH", inverse_simpson: inverse(counts), diversity_status: "available"}));
    const [row] = derivePairMeasures([{library_id_1: "SYN_D1_T1", library_id_2: "SYN_D1_T2", locus: "IGH",
        grouped_consensus_count_1: 1872, grouped_consensus_count_2: 1761,
        n_groups_shared: 1, shared_support_fraction_1: 1392 / 1872, shared_support_fraction_2: 1255 / 1761,
        n_shared_higher_in_1: 1, n_shared_higher_in_2: 0, n_shared_equal: 0}], diversity);
    close(row.simpson_concentration_difference, -0.045478877785994125, "Simpson difference");
    assert.equal(row.shared_difference_balance, -1);
    close(row.support_fraction_only_in_2, 0.2873367404883589, "only in 2");
    close(row.support_fraction_only_in_1, 0.2564102564102564, "only in 1");
    assert.deepEqual([row.shared_higher_in_1_fraction, row.shared_higher_in_2_fraction, row.shared_equal_fraction], [1, 0, 0]);
    for (const field of ["overlap_status", "simpson_difference_status", "shared_difference_status"]) assert.equal(row[field], "available");
});

test("derived pair statuses follow empty endpoints, shared groups and the diversity minimum", () => {
    const base = {library_id_1: "a", library_id_2: "b", locus: "TRB", grouped_consensus_count_1: 10,
        grouped_consensus_count_2: 5, n_groups_shared: 0, shared_support_fraction_1: 0,
        shared_support_fraction_2: 0, n_shared_higher_in_1: 0, n_shared_higher_in_2: 0, n_shared_equal: 0};
    const diversity = [{library_id: "a", locus: "TRB", inverse_simpson: 2, diversity_status: "available"},
        {library_id: "b", locus: "TRB", inverse_simpson: null, diversity_status: "below_minimum_support"}];
    const [disjoint, empty, both, missing] = derivePairMeasures([
        base, {...base, grouped_consensus_count_2: 0}, {...base, grouped_consensus_count_1: 0, grouped_consensus_count_2: 0},
        {...base, locus: "TRA"}], diversity);
    // Nonempty endpoints without shared groups: zero counts, unavailable fractions.
    assert.equal(disjoint.overlap_status, "available");
    assert.equal(disjoint.shared_difference_status, "no_shared_groups");
    assert.equal(disjoint.shared_difference_balance, null);
    assert.equal(disjoint.support_fraction_only_in_2, 1);
    assert.equal(disjoint.simpson_difference_status, "below_minimum_support");
    assert.equal(disjoint.simpson_concentration_difference, null);
    assert.equal(empty.overlap_status, "empty_second");
    assert.equal(empty.support_fraction_only_in_1, null);
    assert.equal(both.overlap_status, "empty_both");
    assert.equal(both.simpson_difference_status, "empty_both");
    assert.equal(missing.simpson_difference_status, "unavailable");
});

test("pair support overrides are per locus, reversible and leave library visibility and flags intact", () => {
    const controller = new SelectionController(longitudinalFixture());
    controller.selectPairs(["p3"]);
    assert.equal(controller.pairMinimum(), 120);
    assert.equal(controller.pairMeetsSupport(controller.pairs.get("p3")), false);
    assert.equal(controller.pairMeetsSupport(controller.pairs.get("p2")), true);
    const visible = controller.visibleRows();
    controller.setPairMinimum(100);
    assert.equal(controller.pairMeetsSupport(controller.pairs.get("p3")), true);
    assert.equal(controller.pairs.get("p3").pair_status, "below_minimum_support");
    controller.setContext({locus: "TRB"});
    assert.equal(controller.pairMinimum(), 0);
    controller.setPairMinimum(40);
    controller.setContext({locus: "IGH"});
    assert.equal(controller.pairMinimum(), 100);
    assert.deepEqual(controller.visibleRows(), visible);
    assert.deepEqual([...controller.selected_pairs], ["p3"]);
    for (const value of [-1, 1.5, NaN, Infinity]) assert.throws(() => controller.setPairMinimum(value));
    controller.reset();
    assert.equal(controller.pairMinimum(), 120);
});

test("plot copies separate support exclusion, missing visit order and measurement availability", () => {
    const controller = new SelectionController(longitudinalFixture());
    const descriptor = {column: "shared_difference_balance", directional: true, measurement_reason_column: "shared_difference_status"};
    const rows = [...controller.pairs.keys()].map(id => controller.displayPair(id));
    let plotted = controller.pairMetricRows(rows, descriptor);
    assert.deepEqual(plotted.map(row => row.shared_difference_balance), [null, -.6, null]);
    assert.equal(plotted[2].pair_plot_status, "below_active_support_minimum");
    controller.setPairMinimum(0);
    plotted = controller.pairMetricRows(rows, descriptor);
    assert.equal(plotted[0].pair_plot_status, "visit_order_unavailable");
    assert.equal(plotted[2].shared_difference_balance, .6);
    assert.equal(controller.pairs.get("p3").shared_difference_balance, -.6);
    assert.equal(controller.pairMetricRows(rows, {column: "bray_curtis_proportion"})[0].bray_curtis_proportion, .45);
    rows[1].shared_difference_status = "no_shared_groups";
    rows[1].shared_difference_balance = null;
    assert.equal(controller.pairMetricRows(rows, descriptor)[1].pair_plot_status, "no_shared_groups");
    rows[1].overlap_status = "empty_first";
    controller.setPairMinimum(1000);
    assert.equal(controller.pairMetricRows(rows, descriptor)[1].pair_plot_status, "empty_first");
});

test("log violin median uses the central original values of the positive population", () => {
    const result = violinSummary([null, 0, 1, 2, 3, 9, 10, 100], true);
    assert.equal(result.median, Math.log10(6));
    assert.notEqual(result.median, (Math.log10(3) + Math.log10(9)) / 2);
    assert.deepEqual(result.values, [1, 2, 3, 9, 10, 100].map(Math.log10));
    assert.equal(violinSummary([0, 0, 1, 2, 3, 4], true), null);
});

test("formatting distinguishes missing, zero and tiny nonzero proportions", () => {
    assert.equal(format(null), "Unavailable");
    assert.equal(format(0, {format: "percent"}), "0%");
    assert.equal(format(.075, {format: "percent"}), "7.5%");
    assert.notEqual(format(1e-8, {format: "percent"}), "0%");
    assert.equal(format(-.125), "-0.125");
    const small = axis("SHM", [0, .00004], {format: "percent"}, "linear");
    assert.ok(Number(small.tickformat.match(/\d+/)[0]) >= 3);
});

test("comparison choices follow supplied chronology before lexical labels", () => {
    const rows = [{timepoint: "Z baseline", timepoint_order: 1},
        {timepoint: "M middle", timepoint_order: 2}, {timepoint: "A end", timepoint_order: 3}];
    const values = [["Z baseline", "A end"], ["M middle", "A end"], ["Z baseline", "M middle"]].map(JSON.stringify);
    assert.deepEqual(orderedComparisons(values, rows), [values[2], values[0], values[1]]);
});

