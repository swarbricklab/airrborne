/* Shared report state. Filters control visibility; selections retain identities.
 * All measurements arrive from pipeline tables. This controller selects rows,
 * counts their availability for display and derives the pair measurements that
 * follow exactly from released columns.
 */
(function (root) {
    "use strict";

    const present = value => typeof value === "number" && Number.isFinite(value);
    const counts = rows => ({
        libraries: rows.length,
        tissues: new Set(rows.map(row => row.tissue_id)).size,
        donors: new Set(rows.map(row => row.donor_id)).size,
    });
    const category = value => value === null || value === undefined || value === "" ? null : String(value);

    // Pair measurements the released pair table leaves out because they follow
    // exactly from its columns or from diversity.tsv. Each is derived for the
    // canonical endpoints, so displayPair() reorients it with the released ones.
    // Endpoint 1 is p and endpoint 2 is q. An empty endpoint has zero
    // grouping-eligible support, and its status replaces every comparison.
    function derivePairMeasures(pairs, diversity = []) {
        const diversityRows = new Map(diversity.map(row => [JSON.stringify([row.library_id, row.locus]), row]));
        return pairs.map(pair => {
            const first = pair.grouped_consensus_count_1, second = pair.grouped_consensus_count_2;
            const overlap = !first && !second ? "empty_both" : !first ? "empty_first" : !second ? "empty_second" : "available";
            const shared = overlap !== "available" ? overlap : pair.n_groups_shared > 0 ? "available" : "no_shared_groups";
            const endpoints = [pair.library_id_1, pair.library_id_2].map(id => diversityRows.get(JSON.stringify([id, pair.locus])));
            // diversity.tsv withholds inverse Simpson below the diversity minimum,
            // and 1 / inverse_simpson is that endpoint's Simpson concentration.
            const simpson = overlap !== "available" ? overlap
                : endpoints.some(row => !row) ? "unavailable"
                : endpoints.some(row => row.diversity_status !== "available") ? "below_minimum_support" : "available";
            const fraction = count => shared === "available" ? count / pair.n_groups_shared : null;
            return {
                ...pair,
                overlap_status: overlap,
                simpson_difference_status: simpson,
                simpson_concentration_difference: simpson === "available"
                    ? 1 / endpoints[1].inverse_simpson - 1 / endpoints[0].inverse_simpson : null,
                shared_difference_status: shared,
                shared_higher_in_2_fraction: fraction(pair.n_shared_higher_in_2),
                shared_higher_in_1_fraction: fraction(pair.n_shared_higher_in_1),
                shared_equal_fraction: fraction(pair.n_shared_equal),
                shared_difference_balance: shared === "available"
                    ? (pair.n_shared_higher_in_2 - pair.n_shared_higher_in_1) / pair.n_groups_shared : null,
                support_fraction_only_in_2: overlap === "available" ? 1 - pair.shared_support_fraction_2 : null,
                support_fraction_only_in_1: overlap === "available" ? 1 - pair.shared_support_fraction_1 : null,
            };
        });
    }

    function keyed(rows, key) {
        const result = new Map();
        for (const row of rows) {
            if (typeof row[key] !== "string" || !row[key] || result.has(row[key])) {
                throw new Error(`Expected a unique string ${key}.`);
            }
            result.set(row[key], row);
        }
        return result;
    }

    class SelectionController {
        constructor(data) {
            this.libraries = keyed(data.libraries, "library_id");
            this.pairs = keyed(data.pairs || [], "pair_id");
            for (const pair of this.pairs.values()) {
                const a = this.libraries.get(pair.library_id_a);
                const b = this.libraries.get(pair.library_id_b);
                if (!a || !b || a.library_id === b.library_id || a.donor_id !== b.donor_id) {
                    throw new Error(`Invalid in-scope donor pair: ${pair.pair_id}.`);
                }
            }
            this.defaults = {...data.default_context};
            this.listeners = new Set();
            this.reset();
        }

        subscribe(listener) {
            this.listeners.add(listener);
            return () => this.listeners.delete(listener);
        }

        notify() {
            for (const listener of this.listeners) listener(this);
        }

        reset() {
            this.filters = new Map();
            this.selected = new Set();
            this.selected_pairs = new Set();
            this.selected_donors = new Set();
            this.selection_only = false;
            this.detail_library_id = null;
            this.context = {...this.defaults};
            this.scales = new Map();
            this.notify();
        }

        clearSelection() {
            this.selected.clear();
            this.selected_pairs.clear();
            this.selected_donors.clear();
            this.selection_only = false;
            this.context.pair_detail_id = "";
            this.notify();
        }

        toggleLibrary(id, additive = false) {
            this.requireIds([id], this.libraries);
            if (!additive && this.selected.size === 1 && this.selected.has(id)) {
                this.clearSelection();
                return;
            }
            const ids = additive ? new Set(this.selected) : new Set();
            if (ids.has(id)) ids.delete(id);
            else ids.add(id);
            if (ids.size) this.selectLibraries([...ids]);
            else this.clearSelection();
        }

        togglePair(id, additive = false) {
            this.requireIds([id], this.pairs);
            if (!additive && this.selected_pairs.size === 1 && this.selected_pairs.has(id)) {
                this.clearSelection();
                return;
            }
            const ids = additive ? new Set(this.selected_pairs) : new Set();
            if (ids.has(id)) ids.delete(id);
            else ids.add(id);
            if (ids.size) this.selectPairs([...ids]);
            else this.clearSelection();
        }

        setScale(metric, scale) {
            this.scales.set(metric, scale);
            this.notify();
        }

        setContext(changes) {
            this.context = {...this.context, ...changes};
            this.notify();
        }

        setCategoricalFilter(field, values) {
            if (values === null) this.filters.delete(field);
            else this.filters.set(field, {type: "category", values: new Set(values.map(category))});
            this.notify();
        }

        setNumericFilter(field, minimum, maximum, include_missing = false) {
            if (minimum === null && maximum === null && include_missing) this.filters.delete(field);
            else this.filters.set(field, {type: "numeric", minimum, maximum, include_missing});
            this.notify();
        }

        matches(row) {
            for (const [field, filter] of this.filters) {
                const value = row[field];
                if (filter.type === "category") {
                    if (!filter.values.has(category(value))) return false;
                } else if (!present(value)) {
                    if (!filter.include_missing) return false;
                } else if ((filter.minimum !== null && value < filter.minimum) ||
                           (filter.maximum !== null && value > filter.maximum)) return false;
            }
            return true;
        }

        filteredRows() {
            return [...this.libraries.values()].filter(row => this.matches(row));
        }

        visibleRows() {
            return this.filteredRows().filter(row => !this.selection_only || this.selected.has(row.library_id));
        }

        selectedRows() {
            return [...this.selected].map(id => this.libraries.get(id));
        }

        highlightedPairs() {
            if (this.selected_pairs.size) return this.selected_pairs;
            return new Set([...this.pairs.values()].filter(pair =>
                this.selected.has(pair.library_id_a) && this.selected.has(pair.library_id_b)
            ).map(pair => pair.pair_id));
        }

        visiblePairs() {
            const visible = new Set(this.visibleRows().map(row => row.library_id));
            return [...this.pairs.values()].filter(pair => visible.has(pair.library_id_a) && visible.has(pair.library_id_b)
                && (!pair.locus || !this.context.locus || pair.locus === this.context.locus)
                && (!this.context.timepoint_pair || pair.timepoint_pair === this.context.timepoint_pair));
        }

        requireIds(ids, source) {
            for (const id of ids) {
                if (!source.has(id)) throw new Error(`Unknown in-scope identity: ${id}.`);
            }
        }

        selectLibraries(ids) {
            this.requireIds(ids, this.libraries);
            this.selected = new Set(ids);
            this.selected_pairs.clear();
            this.selected_donors.clear();
            this.notify();
        }

        selectPairs(ids) {
            this.requireIds(ids, this.pairs);
            this.selected_pairs = new Set(ids);
            this.selected = new Set();
            this.selected_donors = new Set();
            for (const id of ids) {
                const pair = this.pairs.get(id);
                this.selected.add(pair.library_id_a);
                this.selected.add(pair.library_id_b);
                this.selected_donors.add(this.libraries.get(pair.library_id_a).donor_id);
            }
            this.notify();
        }

        selectDonors(ids) {
            const donors = new Set(ids);
            const known = new Map([...this.libraries.values()].map(row => [row.donor_id, true]));
            this.requireIds(ids, known);
            this.selected = new Set([...this.libraries.values()].filter(row => donors.has(row.donor_id)).map(row => row.library_id));
            this.selected_donors = donors;
            this.selected_pairs.clear();
            this.notify();
        }

        selectTissues(ids) {
            const tissues = new Set(ids);
            this.selectLibraries([...this.libraries.values()].filter(row => tissues.has(row.tissue_id)).map(row => row.library_id));
        }

        expandDonors() {
            this.selectDonors([...new Set(this.selectedRows().map(row => row.donor_id))]);
        }

        setSelectionOnly(value) {
            this.selection_only = Boolean(value);
            this.notify();
        }

        openDetail(id) {
            if (id !== null) this.requireIds([id], this.libraries);
            this.detail_library_id = id;
            this.notify();
        }

        availability(fields, rows = this.visibleRows(), reasons = {}) {
            const eligible = rows.filter(row => fields.every(field => present(row[field])));
            const missing = rows.filter(row => !fields.every(field => present(row[field])));
            const unavailable = {};
            for (const row of missing) {
                for (const field of fields.filter(field => !present(row[field]))) {
                    const reason = row[reasons[field]] || "unavailable";
                    unavailable[reason] = (unavailable[reason] || 0) + 1;
                }
            }
            return {
                total: counts(rows), eligible: counts(eligible),
                zeros: rows.filter(row => fields.some(field => row[field] === 0)).length,
                missing_library_ids: missing.map(row => row.library_id), reasons: unavailable,
            };
        }

        orderedPair(id) {
            const pair = this.pairs.get(id);
            const a = this.libraries.get(pair.library_id_a);
            const b = this.libraries.get(pair.library_id_b);
            if (!present(a.timepoint_order) || !present(b.timepoint_order) || a.timepoint_order === b.timepoint_order) return null;
            return a.timepoint_order < b.timepoint_order ? [a, b] : [b, a];
        }

        // Return a view of supplied measurements. Canonical endpoints and pair
        // identity stay unchanged in this.pairs and in the released payload.
        displayPair(id) {
            const source = this.pairs.get(id), ordered = this.orderedPair(id);
            const row = {...source, has_visit_order: Boolean(ordered)};
            if (!ordered || ordered[0].library_id === source.library_id_a) return row;
            for (const [first, second] of [
                ["library_id_a", "library_id_b"], ["library_id_1", "library_id_2"],
                ["tissue_id_1", "tissue_id_2"], ["n_groups_1", "n_groups_2"],
                ["grouped_consensus_count_1", "grouped_consensus_count_2"],
                ["shared_consensus_count_1", "shared_consensus_count_2"],
                ["shared_support_fraction_1", "shared_support_fraction_2"],
                ["n_shared_higher_in_2", "n_shared_higher_in_1"],
                ["shared_higher_in_2_fraction", "shared_higher_in_1_fraction"],
                ["support_fraction_only_in_2", "support_fraction_only_in_1"],
            ]) [row[first], row[second]] = [source[second], source[first]];
            for (const field of ["simpson_concentration_difference", "shared_difference_balance"]) {
                if (present(row[field])) row[field] = row[field] === 0 ? 0 : -row[field];
            }
            for (const field of ["overlap_status", "pair_status", "simpson_difference_status", "shared_difference_status"]) {
                if (row[field] === "empty_first") row[field] = "empty_second";
                else if (row[field] === "empty_second") row[field] = "empty_first";
            }
            return row;
        }

        configuredPairMinimum(locus = this.context.locus) {
            return [...this.pairs.values()].find(pair => pair.locus === locus)?.min_consensus_count_per_library ?? 0;
        }

        pairMinimum(locus = this.context.locus) {
            return this.context.pair_support_minima?.[locus] ?? this.configuredPairMinimum(locus);
        }

        setPairMinimum(value) {
            if (!Number.isSafeInteger(value) || value < 0) throw new Error("Pair support minimum must be a non-negative integer.");
            this.setContext({pair_support_minima: {...this.context.pair_support_minima, [this.context.locus]: value}});
        }

        pairMeetsSupport(pair) {
            return pair.grouped_consensus_count_1 >= this.pairMinimum(pair.locus) &&
                pair.grouped_consensus_count_2 >= this.pairMinimum(pair.locus);
        }

        // Mask only the plotted copy. Excluded values remain available for
        // inspection and exports, with their original measurement statuses.
        pairMetricRows(rows, descriptor) {
            return rows.map(row => {
                let reason = row[descriptor.measurement_reason_column] || "available";
                if (row.overlap_status && row.overlap_status !== "available") reason = row.overlap_status;
                else if (!this.pairMeetsSupport(row)) reason = "below_active_support_minimum";
                else if (descriptor.directional && !row.has_visit_order) reason = "visit_order_unavailable";
                return {...row, [descriptor.column]: reason === "available" ? row[descriptor.column] : null,
                    pair_plot_status: reason};
            });
        }

        orderedTimepoints() {
            const orders = new Map();
            for (const row of this.libraries.values()) {
                if (category(row.timepoint) === null || !present(row.timepoint_order)) continue;
                if (!orders.has(row.timepoint)) orders.set(row.timepoint, new Set());
                orders.get(row.timepoint).add(row.timepoint_order);
            }
            return [...orders].filter(([, values]) => values.size === 1).map(([label, values]) => ({label, order: [...values][0]}))
                .filter(visit => [...orders].every(([label, values]) => label === visit.label || !values.has(visit.order)))
                .sort((a, b) => a.order - b.order);
        }

        referencePairs(reference) {
            const visits = new Map(this.orderedTimepoints().map(visit => [visit.label, visit.order]));
            return this.visiblePairs().flatMap(pair => {
                const endpoints = this.orderedPair(pair.pair_id);
                if (!endpoints || endpoints[0].timepoint !== reference || endpoints.some(row => visits.get(row.timepoint) !== row.timepoint_order)) return [];
                return [{...pair, reference_id: endpoints[0].library_id, followup_id: endpoints[1].library_id,
                    followup_timepoint: endpoints[1].timepoint, followup_order: endpoints[1].timepoint_order}];
            });
        }
    }

    root.ReportSelection = {SelectionController, derivePairMeasures, present, counts, category};
    if (typeof module !== "undefined" && module.exports) module.exports = root.ReportSelection;
})(typeof globalThis !== "undefined" ? globalThis : this);
