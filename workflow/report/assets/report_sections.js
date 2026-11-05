/* Chapter views read released measurements. Library and gene-use scripts
 * load independently, so cohort gene comparisons do not load sequence details.
 */
(function (root) {
    "use strict";
    const {present, category, derivePairMeasures} = root.ReportSelection;
    const {byId, escape, label, options, number, format, axis, palette, compare} = root.ReportView;
    const unpack = table => table.rows.map(values => Object.fromEntries(table.columns.map((column, index) => [column, values[index]])));
    const loci = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
    const open = id => byId(id).open;

    function mount(data) {
        data.libraries = unpack(data.libraries);
        const tables = Object.fromEntries(Object.entries(data.tables).map(([key, table]) => [key, unpack(table)]));
        // Derived before the controller orients any pair, from the canonical endpoints.
        data.pairs = derivePairMeasures(unpack(data.pairs), tables.diversity);
        const calibration = Object.fromEntries(Object.entries(data.calibration).map(([key, table]) => [key, unpack(table)]));
        const densityRows = Object.values(calibration).flat().filter(row => "distance_bin_midpoint" in row);
        const donors = [...new Set(densityRows.map(row => row.donor_id))].sort(compare);
        const selectionDonors = [...new Set(data.libraries.map(row => row.donor_id))].sort(compare);
        Object.assign(data.default_context, {antibody_tab: "composition", isotype_resolution: "main", dominance: "top_5_support_fraction", shm_isotype: "all_IGH", shm_region: "whole_v", curve_library_id: "", curve_order: "inverse_simpson", cohort_curve_order: "inverse_simpson", overlap_metric: "morisita_horn", reference_timepoint: "", pair_abundance_scale: "support_fraction", pair_group_search: "", pair_detection: "all", pair_group_sort: "difference", qc_view: "common", qc_stage: "", timepoint_pair: "", pair_detail_id: "", secondary_diversity: "richness", depth_metric: "input_read_pairs", gene_type: "V", gene_group: "", gene_search: "", gene_share: "group_fraction", calibration_donor: donors[0] || ""});
        // Selection estimates come from their detail table, not the headline.
        for (const region of ["cdr", "fwr"]) data.metrics.push({column: `sigma_${region}`, label: `${region.toUpperCase()} BASELINe`, unit: "BASELINe sigma", source: "selection_pressure.tsv", format: "number", scale: "linear", reason_column: `sigma_${region}_status`, interval_low: `sigma_${region}_ci_low`, interval_high: `sigma_${region}_ci_high`, count_column: `n_clone_consensus_sequences_used_${region}`});
        const report = root.ReportView.mount(data), controller = report.controller;
        const m = report.metric, distribution = report.distribution, scatter = report.scatter;
        const details = new Map(), genes = new Map(), curves = new Map(), loading = new Map(), signatures = new Map();
        root.REPORT_DETAIL = {}; root.REPORT_GENES = {}; root.REPORT_CURVES = {};
        function loadScript(kind, id, filename, target, unpackPayload, cache) {
            if (cache.has(id)) return Promise.resolve(cache.get(id));
            const key = `${kind}:${id}`;
            if (!loading.has(key)) loading.set(key, new Promise((resolve, reject) => {
                const script = document.createElement("script"); script.src = filename;
                script.onload = () => { const value = unpackPayload(target[id]); cache.set(id, value); delete target[id]; script.remove(); resolve(value); };
                script.onerror = () => reject(new Error(`Could not load ${kind} detail for ${id}. Keep the report folder's files together.`));
                document.head.append(script);
            }));
            return loading.get(key);
        }
        const loadDetail = id => loadScript("library", id, data.detail_files[id], root.REPORT_DETAIL, payload => Object.fromEntries(Object.entries(payload).map(([key, table]) => [key, unpack(table)])), details);
        const loadGenes = locus => loadScript("gene", locus, data.gene_usage_files[locus], root.REPORT_GENES, unpack, genes);
        const loadCurves = () => loadScript("curves", "cohort", data.curve_file, root.REPORT_CURVES, unpack, curves);
        Object.assign(report, {load_detail: loadDetail, load_genes: loadGenes, loaded_detail_ids: () => [...details.keys()], loaded_gene_loci: () => [...genes.keys()]});
        const joinRows = (source, predicate, libraries = controller.visibleRows()) => { const rows = new Map(source.filter(predicate).map(row => [row.library_id, row])); return libraries.map(row => ({...row, ...rows.get(row.library_id)})); };
        const supportExtra = row => [["Retained read support", number(row.grouped_consensus_count)], ["Observed groups", number(row.richness)], ["Target read support", number(row.rarefaction_target)]];
        const shmExtra = row => [["Eligible clones", number(row.n_clones_eligible)], ["Eligible members", number(row.n_members_eligible)], ["Members before coverage gate", number(row.n_members_measured)]];
        const qcView = row => ["assembly", "annotation", "common_qc"].includes(row.analysis) ? (row.kind === "criterion" ? "criteria" : "common") : row.analysis;
        const qcKey = row => JSON.stringify([row.stage, row.locus_source, row.isotype_group, row.v_region]);
        const qcLabel = row => [row.stage.replaceAll("_", " "), row.locus === "all" ? "all loci" : row.locus, row.isotype_group, row.v_region].filter(Boolean).join("; ");
        // qc.tsv fractions, each with the support count in its row that divides it.
        const qcFractions = [["grouping_eligible_fraction_of_common_qc", "Grouping-eligible share of common-QC support", "common_qc_consensus_count"], ["common_qc_fraction_of_annotated", "Common-QC share of annotated support", "igblast_consensus_count"], ["nonproductive_fraction", "Non-productive share of annotated support", "igblast_consensus_count"], ["v_call_unresolved_fraction", "Unresolved V call share of common-QC support", "common_qc_consensus_count"], ["j_call_unresolved_fraction", "Unresolved J call share of common-QC support", "common_qc_consensus_count"], ["junction_ineligible_fraction", "Ineligible junction share of common-QC support", "common_qc_consensus_count"], ["c_call_ineligible_fraction", "Ineligible constant-region call share of common-QC support", "common_qc_consensus_count"], ["contig_representative_ineligible_fraction", "Share of common-QC support left out with its contig's representative", "common_qc_consensus_count"], ["class_callable_fraction", "IGH class-callable share of common-QC support", "common_qc_consensus_count"], ["subclass_callable_fraction", "IGH subclass-callable share of common-QC support", "common_qc_consensus_count"]];
        function sourceRows(locus = controller.context.locus) { return joinRows(tables.diversity, row => row.locus === locus); }
        function metricForDiversity(column, overrides = {}) { return m(column, {source: "diversity.tsv", comparison_column: `${controller.context.locus.toLowerCase()}_${column}`, reason_column: column.startsWith("rarefied_") ? "rarefaction_status" : "diversity_status", ...overrides}); }

        async function recovery() {
            const rows = controller.visibleRows();
            for (const [id, chains, title] of [["signal-bcr", loci.slice(0, 3), "Antibody receptor chains"], ["signal-tcr", loci.slice(3), "T-cell receptor chains"]]) {
                const values = rows.flatMap(row => chains.map(locus => ({...row, series: locus, value: row[`${locus.toLowerCase()}_yield_per_million_read_pairs`], measurement_label: `${locus} receptor yield`})));
                await distribution(id, values, m("receptor_yield", {column: "value", label: "Recovered receptor yield", unit: "Read support per million input pairs", scale: "log"}), {title, series: chains, source: "yield.tsv", description: "Each point is one library at one receptor chain. A violin describes the distribution of libraries; it is not a pooled repertoire.", extra: row => [["Input read pairs", number(row.input_read_pairs)]]});
            }
            if (open("recovery-depth")) await scatter("input-depth", rows, m("input_read_pairs"), m(`${controller.context.locus.toLowerCase()}_yield_per_million_read_pairs`), {title: "Receptor recovery and input depth"});
            if (open("qc-section")) await qc();
        }
        async function qc() {
            const locus = controller.context.locus;
            const [fraction, fractionLabel, denominator] = qcFractions.find(([column]) => column === controller.context.qc_fraction) || qcFractions[0];
            options("qc-fraction", qcFractions.map(([column, label]) => [column, label]), fraction);
            await distribution("qc-fractions", joinRows(tables.qc, row => row.locus === locus), m(fraction, {label: `${locus} ${fractionLabel.charAt(0).toLowerCase()}${fractionLabel.slice(1)}`, unit: "Share of read support", bounded: true}), {source: "qc.tsv", description: "Each point is one library at the selected receptor chain. A fraction is unavailable when its read-support denominator is zero. Grouping reasons can overlap.", extra: row => [["Annotated read support", number(row.igblast_consensus_count)], ["Common-QC read support", number(row.common_qc_consensus_count)], ["Denominator", denominator === "igblast_consensus_count" ? "Annotated read support" : "Common-QC read support"]]});
            const trail = tables.qc_trail;
            const locusStages = new Set(trail.filter(row => row.locus === locus).map(qcKey));
            const candidates = trail.filter(row => qcView(row) === controller.context.qc_view && (row.locus === locus || row.locus === "all" && !locusStages.has(qcKey(row))));
            const choices = [...new Map([...candidates].sort((a, b) => a.stage_order - b.stage_order).map(row => [qcKey(row), row])).entries()];
            if (!choices.some(([key]) => key === controller.context.qc_stage)) controller.context.qc_stage = (choices.find(([, row]) => row.stage === "common_qc") || choices[0] || [""])[0];
            options("qc-stage", choices.map(([key, row]) => [key, qcLabel(row)]), controller.context.qc_stage);
            const selected = candidates.filter(row => qcKey(row) === controller.context.qc_stage), first = selected[0];
            const rows = joinRows(selected, () => true);
            const countLabel = key => first ? first[key].replaceAll("_", " ") : "Units";
            await scatter("qc-records", rows, m("n_units_in", {label: `${countLabel("input_unit")} entering`, scale: "log"}), m("n_units_out", {label: `${countLabel("output_unit")} leaving`, scale: "log"}), {title: "Counted units", source: "qc_trail.tsv", extra: row => [["Decision", row.decision_status || "Unavailable"]]});
            await scatter("qc-support", rows, m("consensus_count_in", {label: "Read support entering"}), m("consensus_count_out", {label: "Read support leaving"}), {title: "Read support", source: "qc_trail.tsv"});
            byId("qc-note").textContent = first ? `${qcLabel(first)}. Counted units: ${countLabel("input_unit")} entering; ${countLabel("output_unit")} leaving. ${first.kind === "criterion" ? "Failures can overlap with other criteria." : ""}` : "No applicable population is supplied for this chain and view.";
        }
        async function diversity() {
            const rows = sourceRows(), locus = controller.context.locus;
            await distribution("diversity-main", rows, metricForDiversity("rarefied_inverse_simpson", {label: `${locus} target-depth inverse Simpson`, unit: "Effective groups", scale: "log"}), {description: "The effective number of equally supported groups at the supplied target depth. It depends on both group count and the distribution of support.", extra: supportExtra});
            const dominance = controller.context.dominance;
            await distribution("dominance-main", rows, metricForDiversity(dominance, {label: `${locus} observed top-${dominance.includes("10") ? "10" : "5"} support share`, unit: "Share of retained read support"}), {extra: supportExtra});
            if (open("diversity-diagnostics")) {
                const paired = rows.filter(row => present(row.inverse_simpson) && present(row.rarefied_inverse_simpson)).map(row => ({...row, difference: row.rarefied_inverse_simpson - row.inverse_simpson}));
                await scatter("observed-target", paired, m("inverse_simpson", {label: `${locus} observed inverse Simpson`}), m("rarefied_inverse_simpson", {label: "Target-depth inverse Simpson"}), {title: "Observed and target-depth diversity", source: metricForDiversity("inverse_simpson").source, extra: supportExtra});
                await scatter("target-change", paired, m("inverse_simpson", {label: "Observed inverse Simpson"}), m("difference", {label: "Target minus observed", scale: "linear"}), {title: "Change at the target depth", source: metricForDiversity("inverse_simpson").source});
                byId("diversity-availability").textContent = `${paired.length}/${rows.length} libraries have both values. ` + [...new Set(rows.map(row => `Target: ${number(row.rarefaction_target)} read support; configured minimum: ${number(row.min_consensus_count_for_diversity)}.`))].join(" ");
                const column = controller.context.depth_metric === "yield" ? `${locus.toLowerCase()}_yield_per_million_read_pairs` : controller.context.depth_metric;
                await scatter("depth-comparison", rows, m(column, {label: column === "grouped_consensus_count" ? "Retained read support" : m(column).label}), m("rarefied_inverse_simpson", {label: `${locus} target-depth inverse Simpson`}), {title: "Diversity and depth", source: metricForDiversity("inverse_simpson").source, extra: supportExtra});
            }
            if (open("curve-section")) {
                if (!controller.context.curve_library_id) controller.context.curve_library_id = controller.selected.size === 1 ? [...controller.selected][0] : controller.detail_library_id || controller.visibleRows()[0]?.library_id || data.libraries[0].library_id;
                byId("curve-library").value = controller.context.curve_library_id;
                const id = controller.context.curve_library_id, payload = await loadDetail(id), curves = payload.curves.filter(row => row.locus === locus);
                const names = {richness: "Richness (q = 0)", hill_shannon: "Effective Shannon (q = 1)", inverse_simpson: "Inverse Simpson (q = 2)"};
                const columns = controller.context.curve_order === "all" ? Object.keys(names) : [controller.context.curve_order];
                await report.simple("rarefaction-curves", m("curve", {label: `${locus} recovered-support rarefaction`}), columns.map((column, index) => ({type: "scatter", mode: "lines+markers", name: names[column], x: curves.map(row => row.support_depth), y: curves.map(row => row[column]), customdata: curves.map(row => [id, `Library ${escape(id)}<br>${names[column]}: ${number(row[column])}<br>Subsample read support: ${number(row.support_depth)}${row.is_target ? "<br>Target depth" : ""}${row.is_endpoint ? "<br>Observed endpoint" : ""}`, row[column]]), hovertemplate: "%{customdata[1]}<extra></extra>", meta: {kind: "library", metric: column}, marker: {size: curves.map(row => row.is_target ? 10 : row.is_endpoint ? 8 : 3), color: palette[index]}, line: {color: palette[index]}})), {xaxis: {title: {text: "Subsample read support"}}, yaxis: {title: {text: "Effective groups"}}, showlegend: true, legend: {orientation: "h", y: -.25}}, {linked: true, source: "rarefaction_curves.tsv", description: "Curves interpolate the observed support vector and stop at its full depth."});
                byId("curve-note").textContent = `${id}. ${curves.length ? "Target points: " + (curves.filter(row => row.is_target).map(row => number(row.support_depth)).join(", ") || "none available") : "No curve is available for this library and chain"}. ${controller.matches(controller.libraries.get(id)) ? "" : "This library is outside the active filters."}`;
            }
            if (open("cohort-curve-section")) await cohortCurves();
            if (open("gene-section")) await geneComparison();
            if (open("secondary-section")) {
                const column = controller.context.secondary_diversity;
                const titles = {richness: "Observed group count", hill_shannon: "Observed effective Shannon diversity", shannon_entropy: "Shannon entropy", clonality: "Normalised Shannon clonality", top_1_support_fraction: "Top-1 support share", rarefied_richness: "Target-depth expected group count", rarefied_hill_shannon: "Target-depth effective Shannon diversity", inverse_simpson_asymptote: "Estimated q = 2 asymptote"};
                await distribution("secondary-diversity-plot", rows, metricForDiversity(column, {label: titles[column], reason_column: column === "inverse_simpson_asymptote" ? "asymptote_status" : column.startsWith("rarefied_") ? "rarefaction_status" : "diversity_status"}), {extra: supportExtra});
            }
        }
        async function cohortCurves() {
            const libraries = controller.visibleRows(), locus = controller.context.locus;
            byId("cohort-curve-note").textContent = `Loading curves for ${libraries.length} visible libraries.`;
            const visible = new Map(libraries.map(row => [row.library_id, row]));
            const rows = (await loadCurves()).filter(row => row.locus === locus && visible.has(row.library_id)).map(row => ({...visible.get(row.library_id), ...row}));
            const column = controller.context.cohort_curve_order;
            const names = {richness: "Richness (q = 0)", hill_shannon: "Effective Shannon (q = 1)", inverse_simpson: "Inverse Simpson (q = 2)"};
            await report.cohortCurves(libraries, rows, m(column, {label: `${locus} ${names[column]}`, scale: "linear"}));
            const available = new Set(rows.filter(row => present(row[column])).map(row => row.library_id));
            byId("cohort-curve-note").textContent = `${available.size}/${libraries.length} visible libraries have ${locus} curves. ${libraries.length - available.size} have no available curve for this measure.`;
        }
        async function geneComparison() {
            const locus = controller.context.locus, type = controller.context.gene_type, share = controller.context.gene_share;
            const rows = (await loadGenes(locus)).filter(row => row.gene_type.toUpperCase() === type);
            const search = controller.context.gene_search.toLowerCase();
            const matches = [...new Set(rows.map(row => row.gene_equivalence_group))].sort(compare).filter(group => group.toLowerCase().includes(search));
            if (!matches.includes(controller.context.gene_group)) controller.context.gene_group = matches[0] || "";
            options("gene-group", matches.map(group => [group, group]), controller.context.gene_group);
            const values = controller.visibleRows().map(library => {
                const calls = rows.filter(row => row.library_id === library.library_id), chosen = calls.filter(row => row.gene_equivalence_group === controller.context.gene_group);
                const assigned = calls.reduce((sum, row) => sum + row.n_groups, 0), selected = chosen.reduce((sum, row) => sum + row.n_groups, 0);
                return {...library, gene_group_fraction: !controller.context.gene_group || !assigned ? null : chosen.reduce((sum, row) => sum + row[share], 0),
                    assigned_units: assigned, selected_units: assigned ? selected : null,
                    availability_reason: controller.context.gene_group ? "No assigned receptor groups for this gene type" : "No matching gene equivalence group"};
            });
            await distribution("gene-distribution", values, m("gene_group_fraction", {label: `${locus} ${type}: ${controller.context.gene_group || "no matching gene equivalence group"}`, unit: geneShareUnit(share), bounded: true}), {source: "gene_usage.tsv", description: `${share === "support_fraction" ? "Each clone or clonotype contributes its read support." : "Each clone or clonotype contributes once."} An equivalence group sums its single-gene and composite rows. An absent equivalence group is zero only when the library has assigned receptor groups of this gene type.`, extra: row => [["Receptor groups in this equivalence group", number(row.selected_units)], ["All assigned receptor groups", number(row.assigned_units)]]});
            byId("gene-status").textContent = `${values.filter(row => row.assigned_units > 0).length}/${values.length} libraries have assigned ${type} genes.`;
        }
        function geneShareUnit(share) { return share === "support_fraction" ? "Share of assigned read support" : "Share of assigned receptor groups"; }
        async function antibody() {
            const visible = controller.visibleRows();
            if (controller.context.antibody_tab === "composition") {
                const resolution = controller.context.isotype_resolution;
                const names = resolution === "main" ? ["IgM", "IgG", "IgA"] : resolution === "class" ? ["IgM", "IgD", "IgG", "IgA", "IgE"] : ["IgM", "IgD", "IgG1", "IgG2", "IgG3", "IgG4", "IgA1", "IgA2", "IgE"];
                const rows = names.flatMap(name => joinRows(tables.isotype, row => row.isotype_level === (resolution === "subclass" ? "subclass" : "class") && row.isotype === name).map(row => ({...row, series: name, measurement_label: `${name} support share`})));
                await distribution("class-fractions", rows, m("support_fraction", {label: "Isotype proportions across libraries", unit: `Share of ${resolution === "subclass" ? "subclass" : "class"}-callable support`}), {series: names, source: "isotype_usage.tsv", description: "Each point is one library's proportion. Libraries stay separate when grouped by clinical labels or timepoints.", extra: row => [["Isotype read support", number(row.consensus_count)], ["Callable read support", number(row.callable_consensus_count)]]});
                await distribution("ratio-plot", visible, m("log2_igg_iga", {label: "IgG:IgA balance", unit: "Log2 IgG:IgA support ratio"}), {description: "Zero means equal IgG and IgA support. Positive values favour IgG; negative values favour IgA. Both classes require positive support.", extra: row => [["IgG read support", number(row.igg_consensus_count)], ["IgA read support", number(row.iga_consensus_count)]]});
            } else {
                const rows = joinRows(tables.shm, row => row.isotype_group === controller.context.shm_isotype && row.v_region === controller.context.shm_region);
                const comparisonPrefix = controller.context.shm_region === "whole_v" ? controller.context.shm_isotype.toLowerCase() : null;
                await distribution("shm-plot", rows, m("mean_clone_median_shm", {label: "Mean clone-median SHM", unit: "Mutation frequency", bounded: false, comparison_column: comparisonPrefix ? `${comparisonPrefix}_mean_clone_median_shm_whole_v` : null}), {source: "shm.tsv", description: "Mean of the within-library clone medians. Every eligible clone receives equal weight, independent of read support.", extra: shmExtra});
                await distribution("shm-zero", rows, m("fraction_clones_zero_median_shm", {label: "Clones with zero median SHM", unit: "Share of eligible clones", bounded: true, comparison_column: comparisonPrefix ? `${comparisonPrefix}_fraction_clones_zero_median_shm_whole_v` : null}), {source: "shm.tsv", extra: row => [["Clones with zero median SHM", number(row.n_clones_zero_median_shm)], ...shmExtra(row)]});
                byId("shm-counts").textContent = `${rows.reduce((sum, row) => sum + (row.n_members_eligible || 0), 0)} eligible members from ${rows.reduce((sum, row) => sum + (row.n_members_measured || 0), 0)} measured members; ${rows.reduce((sum, row) => sum + (row.n_clones_eligible || 0), 0)} library-local clone observations. Changing region changes eligibility.`;
                if (open("shm-coverage-section")) await scatter("shm-coverage", rows, m("n_members_measured", {label: "Members before coverage gate", scale: "log"}), m("n_members_eligible", {label: "Eligible members", scale: "log"}), {title: "SHM coverage eligibility", source: "shm.tsv"});
                if (open("selection-section")) {
                    const selection = joinRows(tables.selection, row => row.isotype_group === controller.context.shm_isotype);
                    for (const region of ["cdr", "fwr"]) await distribution(`selection-${region}`, selection, m(`sigma_${region}`, {label: region === "cdr" ? "CDR1 and CDR2 selection" : "FWR1 through FWR3 selection", unit: "BASELINe sigma", reason_column: `sigma_${region}_status`}), {low: `sigma_${region}_ci_low`, high: `sigma_${region}_ci_high`, source: "selection_pressure.tsv", extra: row => [["95% posterior interval", `${number(row[`sigma_${region}_ci_low`])} to ${number(row[`sigma_${region}_ci_high`])}`], ["Usable clone densities", number(row[`n_clone_consensus_sequences_used_${region}`])], ["Status", (row[`sigma_${region}_status`] || "Unavailable").replaceAll("_", " ")]]});
                    byId("selection-counts").textContent = "Each region reports its usable clone densities and supplied 95% interval. Sparse estimates remain available. Of a library's effective clone sequences, the table counts those with a callable position in the region and those that gave a usable density.";
                    const table = byId("selection-evidence"); table.replaceChildren();
                    const head = table.createTHead().insertRow();
                    for (const title of ["Library", "Region", "Status", "Effective clones", "Clones with callable positions", "Usable densities", "Sigma (95% posterior interval)"]) {
                        const cell = document.createElement("th"); cell.textContent = title; head.append(cell);
                    }
                    const body = table.createTBody();
                    for (const row of selection) for (const region of ["cdr", "fwr"]) {
                        const line = body.insertRow(); line.dataset.library = row.library_id; line.dataset.region = region;
                        const values = [row.library_id, region.toUpperCase(), (row[`sigma_${region}_status`] || "Unavailable").replaceAll("_", " "), number(row.n_clone_consensus_sequences), number(row[`n_clone_consensus_sequences_with_positions_${region}`]), number(row[`n_clone_consensus_sequences_used_${region}`]), `${number(row[`sigma_${region}`])} (${number(row[`sigma_${region}_ci_low`])} to ${number(row[`sigma_${region}_ci_high`])})`];
                        for (const value of values) line.insertCell().textContent = value;
                    }
                }
            }
        }
        function pairRows() {
            return controller.visiblePairs().map(source => {
                const pair = controller.displayPair(source.pair_id);
                const a = controller.libraries.get(pair.library_id_a), b = controller.libraries.get(pair.library_id_b);
                const group = controller.context.group_field, colour = controller.context.colour_field;
                const values = [a, b].map(row => category(row[group]));
                const separator = pair.has_visit_order ? " to " : " / ";
                const groupLabel = values[0] === values[1] ? label(values[0]) : values.map(label).join(separator);
                return {...pair, [colour]: category(a[colour]) === category(b[colour]) ? a[colour] : null,
                    mixed_colour: category(a[colour]) !== category(b[colour]), display_group: groupLabel,
                    display_visit: a.timepoint && b.timepoint ? [a.timepoint, b.timepoint].join(separator) + (!pair.has_visit_order && a.timepoint !== b.timepoint ? " (order unavailable)" : "") : null};
            });
        }
        const pairExtra = (row, shownMetric) => [
            ["Direction", `${row.library_id_a} to ${row.library_id_b}${row.has_visit_order ? " (supplied visit order)" : " (endpoint order; chronology unavailable)"}`],
            [`${row.library_id_a} read support`, number(row.grouped_consensus_count_1)],
            [`${row.library_id_b} read support`, number(row.grouped_consensus_count_2)],
            ["Groups detected in A / B / both", `${number(row.n_groups_1)} / ${number(row.n_groups_2)} / ${number(row.n_groups_shared)}`],
            ["Shared groups higher in B / equal / higher in A", `${number(row.n_shared_higher_in_2)} / ${number(row.n_shared_equal)} / ${number(row.n_shared_higher_in_1)}`],
            ["Configured pair minimum", number(row.min_consensus_count_per_library)],
            ["Released support status", (row.pair_status || "unavailable").replaceAll("_", " ")],
            ...(shownMetric === "morisita_horn" ? [] : [["Morisita-Horn", number(row.morisita_horn)]]),
            ...(shownMetric === "jaccard" ? [] : [["Jaccard", number(row.jaccard)]]),
            [`${row.library_id_a} support in shared groups`, format(row.shared_support_fraction_1, {format: "percent"})],
            [`${row.library_id_b} support in shared groups`, format(row.shared_support_fraction_2, {format: "percent"})]
        ];
        function overlapMetric(column = controller.context.overlap_metric) {
            const definitions = {
                morisita_horn: {label: "Morisita-Horn overlap", bounded: true},
                jaccard: {label: "Jaccard overlap", bounded: true},
                bray_curtis_proportion: {label: "Bray-Curtis redistribution", bounded: true, description: "Half the total absolute difference in support proportions over all receptor groups detected in either library. Zero means identical proportions; one means no detected groups in common."},
                simpson_concentration_difference: {label: "Difference in Simpson concentration", directional: true, measurement_reason_column: "simpson_difference_status", description: "Sum of squared support proportions in B minus that in A. Positive means B is more concentrated than A. Both endpoints must meet the diversity minimum."},
                shared_difference_balance: {label: "Balance of shared groups higher in B and in A", directional: true, measurement_reason_column: "shared_difference_status", range: [-1.04, 1.04], description: "(Number higher in B - number higher in A) / number shared. Every shared group has equal weight. Positive means more shared groups have a higher observed proportion in B; this is not absolute population expansion."},
                shared_higher_in_2_fraction: {label: "Shared groups with a higher share in B", directional: true, bounded: true, format: "percent", measurement_reason_column: "shared_difference_status"},
                shared_higher_in_1_fraction: {label: "Shared groups with a higher share in A", directional: true, bounded: true, format: "percent", measurement_reason_column: "shared_difference_status"},
            };
            return m(column, {format: "number", scale: "linear", bounded: false, measurement_reason_column: "overlap_status", reason_column: "pair_plot_status", ...definitions[column]});
        }
        function pairInclusion(rows, descriptor) {
            const minimum = controller.pairMinimum(), configured = controller.configuredPairMinimum();
            byId("pair-minimum").value = minimum;
            byId("pair-minimum-note").textContent = `Both endpoints need at least ${number(minimum)} retained support units at ${controller.context.locus}. Configured minimum: ${number(configured)}.${minimum === 0 ? " No positive support minimum is active." : ""}${minimum !== configured ? " Exploratory override; released flags are unchanged." : ""} Open circles identify pairs below the configured minimum. Filters apply only to pair views.`;
            const plotted = rows.filter(row => present(row[descriptor.column]));
            const donors = items => new Set(items.map(row => row.donor_id)).size;
            byId("pair-counts").textContent = `${descriptor.label}: ${plotted.length}/${rows.length} pairs plotted; ${donors(plotted)}/${donors(rows)} donors represented. Repeated pairs from a donor remain separate observations.`;
            const table = byId("pair-inclusion-table"); table.replaceChildren();
            const head = table.createTHead().insertRow();
            for (const text of ["Group", "Candidate pairs / donors", "Plotted pairs / donors", "Excluded pairs / donors", "Below active support minimum", "Visit order unavailable", "Measurement unavailable"]) {
                const cell = document.createElement("th"); cell.scope = "col"; cell.textContent = text; head.append(cell);
            }
            const body = table.createTBody();
            for (const group of [...new Set(rows.map(row => row.display_group))].sort(compare)) {
                const members = rows.filter(row => row.display_group === group), included = members.filter(row => present(row[descriptor.column]));
                const excluded = members.filter(row => !present(row[descriptor.column]));
                const low = excluded.filter(row => row.pair_plot_status === "below_active_support_minimum");
                const unordered = excluded.filter(row => row.pair_plot_status === "visit_order_unavailable");
                const tr = body.insertRow();
                for (const value of [group, `${members.length} / ${donors(members)}`, `${included.length} / ${donors(included)}`,
                    `${excluded.length} / ${donors(excluded)}`, low.length, unordered.length, excluded.length - low.length - unordered.length]) tr.insertCell().textContent = value;
            }
        }
        async function sharing() {
            const candidates = pairRows(), descriptor = overlapMetric();
            const rows = controller.pairMetricRows(candidates, descriptor);
            const overlapExtra = row => pairExtra(row, descriptor.column);
            // The Simpson difference is derived from each endpoint's diversity row.
            const sourceFor = column => column === "simpson_concentration_difference" ? "donor_overlap_pairs.tsv; diversity.tsv" : "donor_overlap_pairs.tsv";
            const source = sourceFor(descriptor.column);
            for (const [id, column] of [["pair-redistribution", "bray_curtis_proportion"], ["pair-dominance", "simpson_concentration_difference"]]) {
                const measure = overlapMetric(column);
                await distribution(id, controller.pairMetricRows(candidates, measure), measure, {kind: "pair", source: sourceFor(column), extra: row => pairExtra(row, column)});
            }
            await distribution("sharing", rows, descriptor, {title: `${controller.context.locus} ${descriptor.label}`, kind: "pair", source,
                description: `${descriptor.description || "Each point is one within-donor library pair."} Group labels retain both endpoints when metadata differ. Pairs from one donor are not independent biological replicates.`, extra: overlapExtra});
            pairInclusion(rows, descriptor);
            if (open("donor-trend-section")) {
                const visits = controller.orderedTimepoints();
                if (!visits.some(visit => visit.label === controller.context.reference_timepoint)) controller.context.reference_timepoint = visits[0]?.label || "";
                const reference = controller.context.reference_timepoint;
                options("reference-timepoint", visits.length ? visits.map(visit => [visit.label, visit.label]) : [["", "No ordered timepoints supplied"]], reference);
                byId("reference-timepoint").disabled = !visits.length;
                const comparisons = controller.referencePairs(reference).map(pair => ({
                    ...controller.libraries.get(pair.reference_id), ...controller.displayPair(pair.pair_id),
                    reference_id: pair.reference_id, followup_id: pair.followup_id,
                    followup_timepoint: pair.followup_timepoint, followup_order: pair.followup_order,
                }));
                await report.donorTrends(controller.pairMetricRows(comparisons, descriptor), visits.filter(visit => visit.order > visits.find(item => item.label === reference)?.order), descriptor, overlapExtra);
                byId("donor-trend-note").textContent = visits.length ? `${comparisons.length} candidate comparisons use a ${reference} reference and a later library with supplied visit order. Active pair-support and library filters apply. ${rows.length - comparisons.length} other visible pairs are outside this reference comparison.` : "This view needs labelled timepoints with consistent, distinct numeric order. Unordered pairs remain inspectable with explicit endpoint direction.";
            }
            if (open("sharing-depth-section")) {
                await scatter("sharing-depth", rows, m("smaller_endpoint_consensus_count", {label: "Smaller endpoint read support"}), descriptor, {kind: "pair", title: "Pair measurements and recovered support", source, extra: overlapExtra});
                const shareMeasure = m("support_fraction", {label: "Read support in shared groups", unit: "Share of library read support", measurement_reason_column: "overlap_status", reason_column: "pair_plot_status"});
                const shares = candidates.flatMap(row => [1, 2].map(endpoint => ({...row, series: `Library ${endpoint === 1 ? "A" : "B"}`, support_fraction: row[`shared_support_fraction_${endpoint}`]})));
                await distribution("shared-support", controller.pairMetricRows(shares, shareMeasure), shareMeasure, {kind: "pair", series: ["Library A", "Library B"], source: sourceFor("shared_support_fraction"), extra: pairExtra, description: "A and B identify the named endpoints. Supplied visit order sets before-to-after direction when available; otherwise source endpoint order is retained."});
            }
            if (open("pair-detail-section")) await pairDetail();
        }
        async function pairDetail() {
            const pair = controller.pairs.has(controller.context.pair_detail_id) ? controller.displayPair(controller.context.pair_detail_id) : null;
            const available = pair && controller.selected_pairs.has(pair.pair_id);
            byId("pair-detail-content").hidden = !available;
            if (!available) { byId("pair-support").textContent = "Select a point in the overlap or donor comparison plot to inspect that pair."; return; }
            const measurements = byId("pair-measurements"); measurements.replaceChildren();
            const metricHead = measurements.createTHead().insertRow();
            for (const text of ["Pair measurement", "Value", "Measurement status"]) {
                const cell = document.createElement("th"); cell.scope = "col"; cell.textContent = text; metricHead.append(cell);
            }
            const metricBody = measurements.createTBody();
            for (const column of ["bray_curtis_proportion", "simpson_concentration_difference", "shared_difference_balance", "shared_higher_in_2_fraction", "shared_higher_in_1_fraction", "support_fraction_only_in_2", "support_fraction_only_in_1"]) {
                const measure = column === "support_fraction_only_in_2" || column === "support_fraction_only_in_1" ? m(column, {
                    label: column === "support_fraction_only_in_2" ? "Support share in groups only in B" : "Support share in groups only in A",
                    measurement_reason_column: "overlap_status"}) : overlapMetric(column);
                const tr = metricBody.insertRow();
                for (const value of [measure.label, format(pair[column], measure), pair[measure.measurement_reason_column].replaceAll("_", " ")]) tr.insertCell().textContent = value;
            }
            const bar = byId("pair-change-bar"); bar.replaceChildren();
            const directions = [["higher_in_2", "Higher in B", "#177e89"], ["equal", "Equal", "#a4afb5"], ["higher_in_1", "Higher in A", "#bc6435"]];
            const changeText = pair.shared_difference_status === "available" ? directions.map(([key, name]) => `${name}: ${number(pair[`n_shared_${key}`])} of ${number(pair.n_groups_shared)} (${format(pair[`shared_${key}_fraction`], {format: "percent"})})`).join("; ") : `Shared-group fractions unavailable: ${pair.shared_difference_status.replaceAll("_", " ")}.`;
            for (const [key, name, colour] of directions) {
                const value = pair[`shared_${key}_fraction`];
                if (!present(value) || value === 0) continue;
                const segment = document.createElement("span"); segment.style.width = `${100 * value}%`; segment.style.backgroundColor = colour;
                segment.title = `${name}: ${number(pair[`n_shared_${key}`])} of ${number(pair.n_groups_shared)}`; bar.append(segment);
            }
            bar.hidden = pair.shared_difference_status !== "available";
            bar.setAttribute("aria-label", changeText);
            byId("pair-change-counts").textContent = changeText;
            const payloads = await Promise.all([loadDetail(pair.library_id_a), loadDetail(pair.library_id_b)]);
            const vectors = payloads.map(payload => new Map(payload.abundance.filter(row => row.locus === pair.locus).map(row => [row.group_key, row])));
            const keys = [...new Set([...vectors[0].keys(), ...vectors[1].keys()])].sort(compare);
            const rows = keys.map(key => {
                const values = vectors.map(vector => vector.get(key));
                // Endpoint absence is a display zero only for a measured population.
                const fractions = values.map((row, index) => row ? row.support_fraction : pair[`grouped_consensus_count_${index + 1}`] > 0 ? 0 : null);
                return {key, support: values.map(row => row?.consensus_count || 0), fractions,
                    detection: values.every(Boolean) ? "shared" : values[0] ? "a" : "b",
                    difference: fractions.every(present) ? Math.abs(fractions[1] - fractions[0]) : null};
            });
            const detectionNames = {shared: "Both libraries", a: `${pair.library_id_a} only`, b: `${pair.library_id_b} only`};
            const fraction = controller.context.pair_abundance_scale === "support_fraction";
            const shown = rows.filter(row => !fraction || row.fractions.every(present));
            const values = row => fraction ? row.fractions : row.support;
            const traces = ["shared", "a", "b"].map((detection, index) => {
                const subset = shown.filter(row => row.detection === detection);
                return {type: "scatter", mode: "markers", name: detectionNames[detection], x: subset.map(row => values(row)[0]), y: subset.map(row => values(row)[1]),
                    customdata: subset.map(row => [row.key, `${escape(row.key)}<br>${escape(detectionNames[detection])}<br>${escape(pair.library_id_a)}: ${number(row.support[0])} read support (${format(row.fractions[0], {format: "percent"})})<br>${escape(pair.library_id_b)}: ${number(row.support[1])} read support (${format(row.fractions[1], {format: "percent"})})`]),
                    hovertemplate: "%{customdata[1]}<extra></extra>", marker: {size: 7, color: palette[index]}};
            });
            const unit = fraction ? "share of read support" : "read support";
            const descriptor = {format: fraction ? "percent" : "number", bounded: fraction};
            await report.simple("pair-abundance", m("pair_abundance", {label: `${pair.locus} groups in the selected pair`}), traces, {xaxis: axis(`${pair.library_id_a} ${unit}`, shown.map(row => values(row)[0]), descriptor, "linear"), yaxis: axis(`${pair.library_id_b} ${unit}`, shown.map(row => values(row)[1]), descriptor, "linear"), shapes: fraction && shown.length ? [{type: "line", x0: 0, y0: 0, x1: 1, y1: 1, line: {color: "#89969d", dash: "dot", width: 1}, layer: "below"}] : [], showlegend: true, legend: {orientation: "h", y: -.3}}, {source: "group_abundance.tsv", description: "The complete union of the two released abundance vectors. Relative support uses each library's own retained support denominator. The dotted line marks equal proportions."});
            const a = controller.libraries.get(pair.library_id_a), b = controller.libraries.get(pair.library_id_b);
            const inView = controller.visiblePairs().some(row => row.pair_id === pair.pair_id) && controller.pairMeetsSupport(pair);
            byId("pair-support").textContent = `${pair.locus}, donor ${a.donor_id}. Library A: ${a.library_id} (${label(a.timepoint)}); library B: ${b.library_id} (${label(b.timepoint)}). ${keys.length} groups detected in either library; ${number(pair.n_groups_shared)} in both. Morisita-Horn: ${number(pair.morisita_horn)}; Jaccard: ${number(pair.jaccard)}. Shared groups account for ${format(pair.shared_support_fraction_1, {format: "percent"})} of A support and ${format(pair.shared_support_fraction_2, {format: "percent"})} of B support.${inView ? "" : " This selected pair is outside the active pair-view filters."}${shown.length < rows.length ? " Relative support is unavailable for an empty library; use read support to see all groups." : ""}`;
            byId("pair-support").textContent += ` Direction: ${pair.has_visit_order ? "earlier A to later B, using supplied visit order" : "named endpoint A to B; chronology is unavailable"}. Retained support A / B: ${number(pair.grouped_consensus_count_1)} / ${number(pair.grouped_consensus_count_2)}. Groups A / B: ${number(pair.n_groups_1)} / ${number(pair.n_groups_2)}. Configured minimum: ${number(pair.min_consensus_count_per_library)}; ${pair.pair_status.replaceAll("_", " ")}. Inspector values retain measurements excluded by the active pair filter.`;
            options("pair-detection", [["all", "All groups"], ...Object.entries(detectionNames)], controller.context.pair_detection);
            const filtered = rows.filter(row => (controller.context.pair_detection === "all" || row.detection === controller.context.pair_detection) && row.key.toLowerCase().includes(controller.context.pair_group_search.toLowerCase()));
            const sort = controller.context.pair_group_sort;
            filtered.sort((a, b) => sort === "id" ? compare(a.key, b.key) : (sort === "support" ? Math.max(...b.support) - Math.max(...a.support) : (b.difference ?? -1) - (a.difference ?? -1)) || compare(a.key, b.key));
            const table = byId("pair-group-table"); table.replaceChildren();
            const head = table.createTHead().insertRow();
            for (const text of [pair.locus === "IGH" ? "Clone" : "Clonotype", "Detection", "A read support", "B read support", "A support share", "B support share"]) { const th = document.createElement("th"); th.scope = "col"; th.textContent = text; head.append(th); }
            const body = table.createTBody();
            for (const row of filtered) {
                const tr = body.insertRow(); tr.dataset.group = row.key;
                for (const text of [row.key, detectionNames[row.detection], ...row.support.map(number), ...row.fractions.map(value => format(value, {format: "percent"}))]) tr.insertCell().textContent = text;
            }
            byId("pair-group-counts").textContent = `${filtered.length}/${rows.length} groups listed. A and B refer to the named libraries above. Differences in support share are descriptive and can reflect recovery or expression.`;
        }
        async function inspection() {
            const id = controller.detail_library_id;
            byId("library-detail").value = id || ""; byId("detail-panel").hidden = id === null;
            if (id === null) return;
            const library = controller.libraries.get(id), payload = await loadDetail(id);
            byId("detail-title").textContent = `${id} | ${library.tissue_id} | ${library.donor_id}${controller.matches(library) ? "" : " (outside filters)"}`;
            // clone_shm.tsv has a whole-V all-IGH row only for a clone with an
            // SHM-eligible member, so the IGH clones of the abundance table
            // without one are those whose whole-V SHM is unavailable.
            const clones = payload.clone_shm.filter(row => row.v_region === "whole_v" && row.isotype_group === "all_IGH");
            const measured = clones.filter(row => present(row.median_shm_frequency)), scored = new Set(measured.map(row => row.clone_id));
            const missing = payload.abundance.filter(row => row.locus === "IGH" && !scored.has(row.group_key));
            const compositions = new Map();
            for (const row of payload.clone_shm.filter(row => row.v_region === "whole_v" && ["IgM", "IgG", "IgA"].includes(row.isotype_group))) { if (!compositions.has(row.clone_id)) compositions.set(row.clone_id, []); compositions.get(row.clone_id).push(`${row.isotype_group}: ${number(row.eligible_consensus_count)} eligible read support`); }
            const descriptor = m("median_shm_frequency", {label: "Clone support and whole-V SHM", format: "percent", bounded: false});
            const support = m("clone_consensus_count", {label: "Clone read support"});
            await report.simple("clone-detail", descriptor, [{type: "scatter", mode: "markers", x: measured.map(row => row.consensus_count), y: measured.map(row => row.median_shm_frequency), customdata: measured.map(row => [row.clone_id, `${escape(row.clone_id)}<br>Read support: ${number(row.consensus_count)}<br>Median SHM: ${format(row.median_shm_frequency, descriptor)}<br>${number(row.n_members_eligible)} eligible members<br>${number(row.n_members_zero_shm)} zero; ${number(row.n_members_positive_shm)} positive<br>${escape(row.mutation_state.replaceAll("_", " "))}<br>Classes with eligible members: ${escape((compositions.get(row.clone_id) || []).join("; "))}`, row.median_shm_frequency]), hovertemplate: "%{customdata[1]}<extra></extra>", marker: {size: 6, color: palette[0]}}], {xaxis: axis("Clone read support (log scale)", measured.map(row => row.consensus_count), support, "log"), yaxis: axis("Median whole-V SHM", measured.map(row => row.median_shm_frequency), descriptor, "linear")}, {source: "clone_shm.tsv", description: "Each point is this library's observation of one inferred IGH clone."});
            byId("detail-unavailable").textContent = `${missing.length} ${missing.length === 1 ? "clone has" : "clones have"} unavailable whole-V SHM. ${missing.slice().sort((a, b) => b.consensus_count - a.consensus_count).slice(0, 10).map(row => `${row.group_key} (support ${number(row.consensus_count)})`).join("; ")}${missing.length > 10 ? ". The remaining clones are the IGH rows of group_abundance.tsv without a whole-V row in clone_shm.tsv." : ""}`;
            const isotypes = tables.isotype.filter(row => row.library_id === id);
            await report.simple("library-isotypes", m("library_composition", {label: "Isotype composition in this library"}), ["class", "subclass"].map((resolution, index) => { const rows = isotypes.filter(row => row.isotype_level === resolution); return {type: "bar", name: resolution, x: rows.map(row => row.isotype), y: rows.map(row => row.support_fraction), customdata: rows.map(row => `${escape(row.isotype)} (${resolution})<br>Share: ${format(row.support_fraction, {format: "percent"})}<br>Callable support: ${number(row.callable_consensus_count)}`), hovertemplate: "%{customdata}<extra></extra>", marker: {color: palette[index]}}; }), {barmode: "group", showlegend: true, legend: {orientation: "h"}, yaxis: {title: {text: "Share of callable support at this resolution"}, tickformat: ".0%", range: [0, 1.04]}}, {source: "isotype_usage.tsv"});
            if (open("library-gene-section")) {
                const rows = (await loadGenes(controller.context.locus)).filter(row => row.library_id === id);
                const share = controller.context.gene_share;
                for (const type of ["V", "J"]) {
                    const members = new Map();
                    for (const row of rows.filter(row => row.gene_type.toUpperCase() === type && present(row[share]))) { if (!members.has(row.gene_equivalence_group)) members.set(row.gene_equivalence_group, []); members.get(row.gene_equivalence_group).push(row); }
                    const calls = [...members].map(([group, entries]) => ({group, entries, value: entries.reduce((sum, row) => sum + row[share], 0)})).sort((a, b) => b.value - a.value || compare(a.group, b.group));
                    const plot = await report.simple(`library-${type.toLowerCase()}`, m(`library_${type}`, {label: `${controller.context.locus} ${type} gene equivalence group usage`}), [{type: "bar", orientation: "h", x: calls.map(row => row.value), y: calls.map(row => row.group), customdata: calls.map(row => `${escape(row.group)}<br>Share: ${format(row.value, {format: "percent"})}${row.entries.map(entry => `<br>${escape(entry.gene_names)}${entry.identical_alleles ? " (" + escape(entry.identical_alleles) + ")" : ""}: ${number(entry.n_groups)} receptor groups, ${number(entry.consensus_count)} read support`).join("")}`), hovertemplate: "%{customdata}<extra></extra>", marker: {color: palette[0]}}], {height: Math.max(320, calls.length * 23), margin: {l: 150, r: 22, t: 12, b: 60}, xaxis: {title: {text: geneShareUnit(share)}, tickformat: ".0%"}, yaxis: {autorange: "reversed", automargin: true}}, {source: "gene_usage.tsv"});
                    plot.querySelector(".plot").style.height = `${Math.max(320, calls.length * 23)}px`;
                }
            }
        }
        async function methods() {
            if (!open("calibration-section")) return;
            const density = densityRows.filter(row => row.donor_id === controller.context.calibration_donor);
            await report.simple("calibration-density", m("calibration", {label: "IGH threshold calibration"}), [["observed_density", "Observed distances"], ["near_component_density", "Near fitted component"], ["far_component_density", "Far fitted component"]].map(([column, name], index) => ({type: "scatter", mode: "lines", name, x: density.map(row => row.distance_bin_midpoint), y: density.map(row => row[column]), customdata: density.map(row => `${name}<br>Distance: ${number(row.distance_bin_midpoint)}<br>Density: ${number(row[column])}`), hovertemplate: "%{customdata}<extra></extra>", line: {color: palette[index]}})), {xaxis: {title: {text: "Distance to nearest"}}, yaxis: {title: {text: "Density"}}, showlegend: true, legend: {orientation: "h", y: -.25}}, {description: "Nearest-neighbour junction distances are fitted separately for each donor. The near and far components model related and unrelated sequences. They are not known biological labels. Display filters do not refit them."});
            byId("calibration-note").textContent = Object.entries(calibration).filter(([, rows]) => rows.some(row => present(row.clone_threshold))).map(([source, rows]) => `${source}: ${rows.filter(row => present(row.clone_threshold)).map(row => number(row.clone_threshold)).join(", ")}`).join(". ") || "Full calibration results are included with the release.";
        }

        options("shm-isotype", [...new Set(tables.shm.map(row => row.isotype_group))].map(value => [value, value === "all_IGH" ? "All IGH" : value]), "all_IGH");
        options("comparison-isotype", [...new Set(tables.selection.map(row => row.isotype_group))].map(value => [value, value === "all_IGH" ? "All IGH" : value]), controller.context.shm_isotype);
        byId("comparison-isotype").onchange = event => controller.setContext({shm_isotype: event.target.value});
        options("shm-region", [...new Set(tables.shm.map(row => row.v_region))].map(value => [value, value === "whole_v" ? "Whole V" : value.toUpperCase()]), "whole_v");
        options("curve-library", data.libraries.map(row => [row.library_id, row.library_id]), controller.context.curve_library_id);
        options("library-detail", [["", "Choose a library"], ...data.libraries.map(row => [row.library_id, `${row.library_id} | ${row.tissue_id}`])], "");
        options("calibration-donor", donors.map(id => [id, id]), controller.context.calibration_donor);
        for (const select of document.querySelectorAll(".locus-control")) { options(select, loci.map(locus => [locus, locus]), "IGH"); select.onchange = event => controller.setContext({locus: event.target.value, gene_group: "", gene_search: ""}); }
        for (const select of document.querySelectorAll(".gene-share-control")) { options(select, [["group_fraction", "Receptor groups"], ["support_fraction", "Read support"]], "group_fraction"); select.onchange = event => controller.setContext({gene_share: event.target.value}); }
        options("qc-view", [["common", "Assembly, annotation and common QC"], ["criteria", "Common QC criteria"], ["grouping", "Grouping"], ["isotype", "Isotype calls"], ["shm", "SHM eligibility"], ["selection", "Selection inputs"]], "common");
        byId("qc-fraction").onchange = event => controller.setContext({qc_fraction: event.target.value});
        const comparisons = root.ReportView.orderedComparisons(data.pairs.map(row => row.timepoint_pair).filter(Boolean), data.libraries);
        options("timepoint-pair", [["", "All supplied timepoints"], ...comparisons.map(value => [value, JSON.parse(value).join(" / ")])], "");
        byId("timepoint-pair").parentElement.hidden = comparisons.length === 0;
        for (const field of ["donor", "tissue"]) { options(`${field}-selection`, [["", `Choose a ${field}`], ...[...new Set(data.libraries.map(row => row[`${field}_id`]))].sort(compare).map(id => [id, id])], ""); byId(`${field}-selection`).onchange = event => { if (event.target.value) controller[field === "donor" ? "selectDonors" : "selectTissues"]([event.target.value]); }; }
        byId("donor-search").oninput = () => {
            const query = byId("donor-search").value.trim().toLowerCase();
            const selected = byId("donor-selection").value;
            const matches = selectionDonors.filter(id => id.toLowerCase().includes(query));
            options("donor-selection", [["", matches.length ? "Choose a donor" : "No matching donors"], ...matches.map(id => [id, id])], selected);
        };
        const fields = [["shm-isotype", "shm_isotype"], ["shm-region", "shm_region"], ["curve-library", "curve_library_id"], ["curve-order", "curve_order"], ["qc-view", "qc_view"], ["qc-stage", "qc_stage"], ["timepoint-pair", "timepoint_pair"], ["pair-detail", "pair_detail_id"], ["secondary-diversity", "secondary_diversity"], ["calibration-donor", "calibration_donor"], ["dominance", "dominance"], ["depth-metric", "depth_metric"], ["gene-type", "gene_type"], ["gene-group", "gene_group"], ["isotype-resolution", "isotype_resolution"], ["cohort-curve-order", "cohort_curve_order"], ["overlap-metric", "overlap_metric"], ["reference-timepoint", "reference_timepoint"], ["pair-abundance-scale", "pair_abundance_scale"], ["pair-detection", "pair_detection"], ["pair-group-sort", "pair_group_sort"]];
        for (const [id, field] of fields) byId(id).onchange = event => controller.setContext({[field]: event.target.value});
        byId("pair-minimum").onchange = event => {
            if (!event.target.reportValidity() || event.target.value === "") return;
            controller.setPairMinimum(Number(event.target.value));
        };
        byId("pair-minimum-reset").onclick = () => controller.setPairMinimum(controller.configuredPairMinimum());
        byId("pair-group-search").oninput = event => controller.setContext({pair_group_search: event.target.value});
        byId("gene-search").oninput = event => controller.setContext({gene_search: event.target.value});
        byId("library-detail").onchange = event => controller.openDetail(event.target.value || null);
        byId("inspect-pair").onclick = async () => {
            byId("pair-detail-section").open = true;
            report.schedule();
            await report.ready;
            byId("pair-detail").focus();
            byId("pair-detail-section").scrollIntoView({block: "start"});
        };
        byId("curve-use-selection").onclick = () => { if (controller.selected.size === 1) controller.setContext({curve_library_id: [...controller.selected][0]}); };
        for (const tab of document.querySelectorAll("[data-antibody-tab]")) tab.onclick = () => controller.setContext({antibody_tab: tab.dataset.antibodyTab});
        for (const element of document.querySelectorAll("details.supporting")) element.addEventListener("toggle", report.schedule);
        byId("methods-narrative").replaceChildren(...data.methods.map(text => { const paragraph = document.createElement("p"); paragraph.textContent = text; return paragraph; }));
        byId("methods-provenance").textContent = JSON.stringify(data.provenance, null, 2);
        byId("metadata-definitions").replaceChildren(...data.metadata.filter(field => field.column !== "report_group").flatMap(field => { const term = document.createElement("dt"), definition = document.createElement("dd"); term.textContent = field.label; definition.textContent = `${field.source_level === "tissue" ? "Tissue" : "Library"} column; ${field.type === "number" ? "numeric" : "categorical"}.`; return [term, definition]; }));
        const renderers = {signal: recovery, diversity, antibody, sharing, methods};
        report.renderViews = async () => {
            for (const [id, field] of fields) byId(id).value = controller.context[field];
            byId("gene-search").value = controller.context.gene_search;
            byId("pair-group-search").value = controller.context.pair_group_search;
            byId("curve-use-selection").disabled = controller.selected.size !== 1;
            byId("curve-use-selection").title = controller.selected.size === 1 ? `Show ${[...controller.selected][0]}` : "Select one library to use it here.";
            for (const element of document.querySelectorAll(".locus-control")) element.value = controller.context.locus;
            for (const element of document.querySelectorAll(".gene-share-control")) element.value = controller.context.gene_share;
            for (const tab of document.querySelectorAll("[data-antibody-tab]")) tab.setAttribute("aria-pressed", String(tab.dataset.antibodyTab === controller.context.antibody_tab));
            byId("composition-view").hidden = controller.context.antibody_tab !== "composition";
            byId("mutation-view").hidden = controller.context.antibody_tab !== "mutation";
            const selectedPairs = [...controller.selected_pairs].map(id => controller.pairs.get(id));
            byId("inspect-pair").disabled = !selectedPairs.length;
            for (const field of ["donor", "tissue"]) {
                const selected = controller.selectedRows(), id = selected[0]?.[`${field}_id`];
                byId(`${field}-selection`).value = id && selected.every(row => row[`${field}_id`] === id) && data.libraries.filter(row => row[`${field}_id`] === id).length === selected.length ? id : "";
            }
            if (!controller.selected_pairs.has(controller.context.pair_detail_id)) controller.context.pair_detail_id = selectedPairs[0]?.pair_id || "";
            options("pair-detail", [["", "Choose a selected pair"], ...selectedPairs.map(pair => [pair.pair_id, pair.pair_label])], controller.context.pair_detail_id);
            const chapter = controller.context.chapter;
            const opened = [...byId(`${chapter}-section`).querySelectorAll("details[id]")].filter(node => node.open).map(node => node.id);
            const signature = report.visibleSignature() + JSON.stringify(opened);
            if (signatures.get(chapter) !== signature) { await renderers[chapter](); signatures.set(chapter, report.visibleSignature() + JSON.stringify(opened)); }
            if (byId("comparison-dialog").open) {
                const baseline = [controller.context.x_metric, controller.context.y_metric].some(column => column.startsWith("sigma_"));
                byId("comparison-isotype-control").hidden = !baseline;
                byId("comparison-isotype").value = controller.context.shm_isotype;
                const needsPairs = controller.context.trajectories && !controller.selected_pairs.size;
                byId("comparison-choose-pairs").hidden = !needsPairs;
                byId("comparison-note").textContent = [baseline ? "BASELINe uses the antibody subset chosen here. Bars show supplied 95% posterior intervals; hover shows the clone counts each estimate used. These exploratory estimates do not measure affinity." : "", needsPairs ? "Choose a pair in Sharing to connect its libraries across visits." : ""].filter(Boolean).join(" ");
                const key = report.visibleSignature() + (controller.context.trajectories ? JSON.stringify([...controller.selected_pairs]) : "");
                if (signatures.get("comparison") !== key) {
                    const rows = baseline ? joinRows(tables.selection, row => row.isotype_group === controller.context.shm_isotype) : controller.visibleRows();
                    const descriptor = column => m(column, column.startsWith("sigma_") ? {label: `${controller.context.shm_isotype === "all_IGH" ? "All IGH" : controller.context.shm_isotype} ${m(column).label}`} : {});
                    await scatter("explorer", rows, descriptor(controller.context.x_metric), descriptor(controller.context.y_metric), {title: "Library measurements", trajectories: true}); signatures.set("comparison", key);
                }
            }
            if (byId("inspection-dialog").open) {
                const key = JSON.stringify([controller.detail_library_id, controller.context.locus, open("library-gene-section"), report.visibleSignature()]);
                if (signatures.get("inspection") !== key) { await inspection(); signatures.set("inspection", key); }
            }
        };
        report.schedule();
        return report;
    }
    root.FinalReport = {mount};
})(window);
