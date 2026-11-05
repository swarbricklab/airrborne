/* Shared presentation for supplied measurements. Density and log coordinates
 * describe the display; raw values and string identities stay on every point.
 */
(function (root) {
    "use strict";
    const {SelectionController, present, category, counts} = root.ReportSelection;
    const byId = id => document.getElementById(id);
    const escape = value => String(value).replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
    const label = value => category(value) === null ? "Unlabelled" : String(value);
    const palette = ["#177e89", "#c36b37", "#7768a6", "#668c45", "#c65b7c", "#4b81b4", "#9b8036", "#586974"];
    // The same Viridis stops colour points, lines and the visible legend.
    const viridis = ["#440154", "#482878", "#3e4989", "#31688e", "#26828e", "#1f9e89", "#35b779", "#6ece58", "#b5de2b", "#fde725"];
    const numericScale = viridis.map((colour, index) => [index / (viridis.length - 1), colour]);
    const compare = (a, b) => a === null ? 1 : b === null ? -1 : String(a).localeCompare(String(b), "en", {numeric: true});
    const unique = values => [...new Set(values)];
    function colourRange(values) {
        const finite = values.filter(present);
        return finite.length ? finite.reduce(([low, high], value) => [Math.min(low, value), Math.max(high, value)], [Infinity, -Infinity]) : null;
    }
    function continuousColour(value, range) {
        if (!present(value) || !range) return "#89969d";
        const fraction = range[0] === range[1] ? .5 : Math.max(0, Math.min(1, (value - range[0]) / (range[1] - range[0])));
        const position = fraction * (viridis.length - 1), index = Math.min(Math.floor(position), viridis.length - 2);
        const blend = position - index;
        const channels = [1, 3, 5].map(offset => {
            const a = parseInt(viridis[index].slice(offset, offset + 2), 16), b = parseInt(viridis[index + 1].slice(offset, offset + 2), 16);
            return Math.round(a + (b - a) * blend).toString(16).padStart(2, "0");
        });
        return "#" + channels.join("");
    }
    function violinSummary(values, logarithmic = false) {
        const eligible = values.filter(value => present(value) && (!logarithmic || value > 0)).sort((a, b) => a - b);
        // Five observations and three distinct values govern smoothing only.
        if (eligible.length < 5 || unique(eligible).length < 3) return null;
        const middle = Math.floor(eligible.length / 2);
        const median = eligible.length % 2 ? eligible[middle] : eligible[middle - 1] / 2 + eligible[middle] / 2;
        return {values: logarithmic ? eligible.map(Math.log10) : eligible, median: logarithmic ? Math.log10(median) : median};
    }
    const population = rows => counts([...new Map(rows.map(row => [row.library_id, row])).values()]);
    const countText = value => `${value.libraries} ${value.libraries === 1 ? "library" : "libraries"}, ${value.tissues} ${value.tissues === 1 ? "tissue" : "tissues"}, ${value.donors} ${value.donors === 1 ? "donor" : "donors"}`;
    const number = value => !present(value) ? "Unavailable" : (value !== 0 && (Math.abs(value) < .001 || Math.abs(value) >= 1e7) ? value.toExponential(2).replace("e+", "e") : value.toLocaleString("en-US", {maximumSignificantDigits: 4}));
    const format = (value, metric = {}) => !present(value) ? "Unavailable" : metric.format === "percent" ? `${number(value * 100)}%` : number(value);
    const jitter = id => { let hash = 0; for (const c of id) hash = (hash * 31 + c.charCodeAt(0)) % 997; return (hash / 997 - .5) * .48; };
    function options(element, choices, selected) {
        if (typeof element === "string") element = byId(element);
        element.replaceChildren(...choices.map(([value, text]) => { const option = document.createElement("option"); option.value = value; option.textContent = text; option.selected = value === selected; return option; }));
    }
    function button(text, action, className = "") { const node = document.createElement("button"); node.textContent = text; node.className = className; node.onclick = action; return node; }
    function orderedVisits(rows) {
        const values = unique(rows.map(row => category(row.timepoint)));
        const orders = new Map(values.map(value => [value, unique(rows.filter(row => category(row.timepoint) === value && present(row.timepoint_order)).map(row => row.timepoint_order))]));
        const known = values.filter(value => value !== null);
        const chronological = known.every(value => orders.get(value).length === 1) && unique(known.map(value => orders.get(value)[0])).length === known.length;
        return values.sort((a, b) => a === null ? 1 : b === null ? -1 : chronological ? orders.get(a)[0] - orders.get(b)[0] : compare(a, b));
    }
    function orderedComparisons(values, rows) {
        const ranks = new Map(orderedVisits(rows).map((value, index) => [value, index]));
        return unique(values).sort((a, b) => {
            const left = JSON.parse(a), right = JSON.parse(b);
            return ranks.get(left[0]) - ranks.get(right[0]) || ranks.get(left[1]) - ranks.get(right[1]) || compare(a, b);
        });
    }
    function metric(column, overrides = {}) {
        const percent = /fraction|_shm|incidence|support_fraction/.test(column);
        const log = /consensus_count|input_read_pairs|yield|inverse_simpson|richness|hill_shannon/.test(column);
        return {column, label: column.replaceAll("_", " "), unit: percent ? "Proportion" : "", format: percent ? "percent" : "number", scale: log ? "log" : "linear", bounded: /fraction/.test(column), ...overrides};
    }
    function axis(title, values, descriptor, scale, transformed = false) {
        const result = {title: {text: title, standoff: 12}, automargin: true, gridcolor: "#edf1f3", zerolinecolor: "#bdccd2", tickfont: {size: 12}, rangemode: "tozero"};
        if (scale === "log") {
            const positive = values.filter(value => present(value) && value > 0);
            let low = positive.length ? Math.min(...positive) : 1, high = positive.length ? Math.max(...positive) : 10;
            if (low === high) { low /= 2; high *= 2; }
            const lo = Math.log10(low), hi = Math.log10(high), pad = Math.max(.07, (hi - lo) * .07);
            let ticks = [];
            for (let power = Math.floor(lo) - 1; power <= Math.ceil(hi) + 1; power++) for (const base of [1, 2, 5]) { const value = base * 10 ** power; if (value >= 10 ** (lo - pad) && value <= 10 ** (hi + pad)) ticks.push(value); }
            if (ticks.length > 9) ticks = ticks.filter(value => Math.abs(Math.log10(value) - Math.round(Math.log10(value))) < 1e-8);
            if (ticks.length < 2) ticks = [low, high];
            Object.assign(result, {type: transformed ? "linear" : "log", range: [lo - pad, hi + pad], tickmode: "array", tickvals: ticks.map(value => transformed ? Math.log10(value) : value), ticktext: ticks.map(value => format(value, descriptor)), minorloglabels: "none"});
        } else {
            if (descriptor.format === "percent") {
                const maximum = Math.max(0, ...values.filter(present).map(value => Math.abs(value)));
                const digits = descriptor.bounded ? 0 : Math.min(6, Math.max(1, Math.ceil(-Math.log10(Math.max(maximum * 100, 1e-6))) + 1));
                result.tickformat = `.${digits}%`;
            }
            if (descriptor.bounded) result.range = [0, 1.04];
            if (descriptor.range) result.range = descriptor.range;
        }
        return result;
    }

    function mount(data) {
        const controller = new SelectionController(data);
        const metrics = new Map(data.metrics.map(item => [item.column, item]));
        const categorical = data.metadata.filter(field => field.type === "category");
        const numeric = data.metadata.filter(field => field.type === "number");
        const colourMaps = new Map(categorical.map(field => [field.column, new Map(unique(data.libraries.map(row => category(row[field.column]))).sort(compare).map((value, index) => [value, value === null ? "#89969d" : palette[index % palette.length]]))]));
        const colourRanges = new Map(numeric.map(field => [field.column, colourRange(data.libraries.map(row => row[field.column]))]));
        const report = {controller, data, metrics, ready: Promise.resolve(), timings: [], render_count: 0, renderViews: async () => {}};
        root.trust4_report = report;
        let queued = false;
        let hadSelection = false;
        function schedule() {
            if (queued) return;
            queued = true;
            report.ready = report.ready.then(async () => {
                queued = false;
                const started = performance.now();
                renderControls();
                await report.renderViews();
                await restyleSelection();
                report.render_count++;
                report.timings.push(performance.now() - started);
            });
            report.ready.catch(error => { byId("report-error").textContent = error.message; });
        }
        report.schedule = schedule;
        controller.subscribe(schedule);
        report.metric = (column, overrides = {}) => metric(column, {...metrics.get(column), ...overrides});
        report.scale = descriptor => controller.scales.get(descriptor.column) || descriptor.scale;
        report.visibleSignature = () => JSON.stringify({context: {...controller.context, pair_detail_id: controller.context.chapter === "sharing" ? controller.context.pair_detail_id : null}, filters: [...controller.filters].map(([key, filter]) => [key, {...filter, values: filter.values ? [...filter.values] : undefined}]), scales: [...controller.scales], selection: controller.selection_only ? [...controller.selected].sort() : null});

        function pointRows(rows) { return rows.map(row => controller.libraries.get(row.library_id)).filter(Boolean); }
        function choose(ids, kind) { if (!ids.length) controller.clearSelection(); else if (kind === "pair") controller.selectPairs(unique(ids)); else controller.selectLibraries(unique(ids)); }
        function toggle(id, kind, additive) { if (kind === "pair") controller.togglePair(id, additive); else controller.toggleLibrary(id, additive); }
        function colour(row) {
            const field = controller.context.colour_field;
            if (row.mixed_colour) return "#536971";
            if (colourRanges.has(field)) return continuousColour(row[field], colourRanges.get(field));
            return colourMaps.get(field)?.get(category(row[field])) || "#89969d";
        }
        function marker(rows, kind) {
            const field = controller.context.colour_field;
            const range = colourRanges.get(field);
            const numericColour = kind !== "pair" && range && range[0] !== range[1] && rows.every(row => present(row[field]));
            const result = {size: 6, color: rows.map(colour), opacity: .8, line: {width: .4, color: "white"}};
            if (numericColour) Object.assign(result, {color: rows.map(row => row[field]), colorscale: numericScale, cmin: range[0], cmax: range[1], showscale: false});
            // Pair rows carry endpoint labels, not an invented averaged numeric score.
            if (kind === "pair" && numeric.some(item => item.column === field)) result.color = "#177e89";
            if (kind === "pair") result.symbol = rows.map(row => row.pair_status === "below_minimum_support" ? "circle-open" : "circle");
            return result;
        }
        function hover(row, descriptor, value, extra = [], kind = "library") {
            const lines = [];
            const group = data.metadata.find(item => item.column === controller.context.group_field);
            if (kind === "pair") {
                const a = controller.libraries.get(row.library_id_a), b = controller.libraries.get(row.library_id_b);
                lines.push(`<b>${escape(row.locus)} | Donor ${escape(a.donor_id)}</b>`);
                for (const [index, item] of [a, b].entries()) {
                    lines.push(`Library ${index ? "B" : "A"}: ${escape(item.library_id)} | Tissue: ${escape(item.tissue_id)}`);
                    if (data.has_metadata) lines.push(`${escape(group.label)}: ${escape(label(item[group.column]))} | Timepoint: ${escape(label(item.timepoint))}`);
                }
                if (controller.context.colour_field !== controller.context.group_field) {
                    const field = data.metadata.find(item => item.column === controller.context.colour_field);
                    for (const item of [a, b]) lines.push(`${escape(item.library_id)} ${escape(field.label)}: ${escape(label(item[field.column]))}${field.unit ? " " + escape(field.unit) : ""}`);
                }
            } else {
                const library = controller.libraries.get(row.library_id) || row;
                lines.push(`<b>Library ${escape(library.library_id)}</b>`, `Tissue: ${escape(library.tissue_id)} | Donor: ${escape(library.donor_id)}`);
                if (data.has_metadata) lines.push(`${escape(group.label)}: ${escape(label(library[group.column]))}`, `Timepoint: ${escape(label(library.timepoint))}`);
                if (controller.context.colour_field !== controller.context.group_field) { const field = data.metadata.find(item => item.column === controller.context.colour_field); lines.push(`${escape(field.label)}: ${escape(label(library[field.column]))}${field.unit ? " " + escape(field.unit) : ""}`); }
            }
            lines.push(`${escape(row.measurement_label || descriptor.label)}: ${escape(format(value, descriptor))}${descriptor.unit && descriptor.format !== "percent" ? " " + escape(descriptor.unit) : ""}`);
            for (const [name, value] of extra) lines.push(`${escape(name)}: ${escape(value)}`);
            return lines.join("<br>");
        }
        report.hover = hover;
        function pointTrace(rows, xs, ys, descriptor, settings = {}) {
            const kind = settings.kind || "library", key = kind === "pair" ? "pair_id" : "library_id";
            const selected = kind === "pair" ? controller.highlightedPairs() : controller.selected;
            const trace = {type: "scatter", mode: "markers", name: settings.name || "Libraries", ids: rows.map(row => row[key]), x: xs, y: ys,
                customdata: rows.map(row => [row[key], hover(row, descriptor, row[descriptor.column], settings.extra ? settings.extra(row) : [], kind), row[descriptor.column], settings.x ? row[settings.x.column] : null, row.series || null]),
                hovertemplate: "%{customdata[1]}<extra></extra>", meta: {kind, metric: descriptor.column}, marker: marker(rows, kind),
                selectedpoints: controller.selected.size ? rows.flatMap((row, index) => selected.has(row[key]) ? [index] : []) : null,
                selected: {marker: {opacity: 1, size: 9}}, unselected: {marker: {opacity: .14}}};
            for (const [direction, item] of [["x", settings.x], ["y", descriptor]]) if (item?.interval_low) {
                trace[`error_${direction}`] = {type: "data", symmetric: false,
                    array: rows.map(row => present(row[item.interval_high]) ? row[item.interval_high] - row[item.column] : null),
                    arrayminus: rows.map(row => present(row[item.interval_low]) ? row[item.column] - row[item.interval_low] : null), width: 2, thickness: .8, color: "#82969e"};
            }
            return trace;
        }
        function bind(div) {
            if (div.trust4_bound) return;
            div.trust4_bound = true;
            div.on("plotly_click", event => {
                const point = event.points?.find(point => point.data.meta?.kind && point.customdata);
                if (point) toggle(point.customdata[0], point.data.meta.kind, Boolean(event.event?.shiftKey));
            });
            div.on("plotly_selected", event => {
                if (!event) return;
                const points = (event.points || []).filter(point => point.data?.meta?.kind && point.customdata);
                choose(points.map(point => point.customdata[0]), points[0]?.data.meta.kind || "library");
            });
            div.on("plotly_deselect", () => controller.clearSelection());
            // A point click need not create Plotly's native lasso selection.
            div.on("plotly_doubleclick", () => controller.clearSelection());
            // Plotly's lasso drag cover can consume the DOM double-click event.
            // The second native click still reaches its plotting surface.
            div.addEventListener("click", event => {
                if (event.detail === 2 && event.target.closest(".nsewdrag")) controller.clearSelection();
            }, true);
            div.on("plotly_legendclick", () => false);
            div.on("plotly_legenddoubleclick", () => false);
        }
        async function plot(div, traces, layout = {}, linked = true) {
            if (typeof div === "string") div = byId(div);
            await Plotly.react(div, traces, {margin: {l: 68, r: 22, t: 14, b: 65}, font: {family: "Arial, sans-serif", size: 12, color: "#36515b"}, paper_bgcolor: "white", plot_bgcolor: "white", hovermode: "closest", hoverlabel: {font: {size: 12}, align: "left"}, dragmode: linked ? "lasso" : "zoom", showlegend: false, uirevision: div.id, ...layout}, {responsive: true, displaylogo: false, modeBarButtonsToRemove: ["toImage", "autoScale2d"]});
            if (linked) bind(div);
        }
        report.plot = plot;

        function frame(id, descriptor, settings = {}) {
            const container = byId(id);
            if (!container.querySelector(".figure-header")) {
                container.innerHTML = '<div class="figure-header"><h3></h3><div class="figure-tools"></div></div><div class="figure-facets"></div><div class="figure-foot"></div>';
            }
            container.querySelector("h3").textContent = settings.title || descriptor.label;
            const tools = container.querySelector(".figure-tools"); tools.replaceChildren();
            for (const [item, name] of [[settings.x, "X scale"], [descriptor, settings.x ? "Y scale" : "Scale"]]) {
                if (!item || item.scale !== "log" && !item.allow_log) continue;
                const select = document.createElement("select"); select.setAttribute("aria-label", `${settings.title || descriptor.label}: ${name}`); select.className = "scale-control";
                options(select, [["log", "Log scale"], ["linear", "Linear scale"]], report.scale(item));
                select.onchange = () => controller.setScale(item.column, select.value); tools.append(select);
            }
            const comparison = descriptor.comparison_column || descriptor.column;
            if (id !== "explorer" && metrics.has(comparison)) tools.append(button("Compare", () => {
                controller.setContext({y_metric: comparison, ...(settings.x && metrics.has(settings.x.column) ? {x_metric: settings.x.column} : {})});
                byId("comparison-dialog").showModal(); schedule();
            }));
            const foot = container.querySelector(".figure-foot"), helpOpen = foot.querySelector("details")?.open; foot.replaceChildren();
            const about = document.createElement("details"); about.className = "metric-help";
            about.open = Boolean(helpOpen);
            const summary = document.createElement("summary"); summary.textContent = "About this measurement"; about.append(summary);
            const paragraph = document.createElement("p"); paragraph.textContent = [settings.description || descriptor.description || `${descriptor.label}${descriptor.unit ? " is shown in " + descriptor.unit : ""}.`, settings.source || descriptor.detail_source || descriptor.source ? `Source: ${settings.source || descriptor.detail_source || descriptor.source}.` : ""].filter(Boolean).join(" "); about.append(paragraph); foot.append(about);
            return {container, facets: container.querySelector(".figure-facets"), foot};
        }
        function facets(frame, labels) {
            const previous = [...frame.facets.children];
            if (previous.length !== labels.length || previous.some((element, index) => element.dataset.label !== labels[index])) {
                for (const plot of frame.facets.querySelectorAll(".js-plotly-plot")) Plotly.purge(plot);
                frame.facets.replaceChildren(...labels.map((name, index) => {
                    const facet = document.createElement("div"); facet.className = "facet"; facet.dataset.label = name;
                    const heading = document.createElement("p"); heading.className = "facet-label"; heading.textContent = name;
                    const graph = document.createElement("div"); graph.className = "plot"; graph.id = `${frame.container.id}-plot-${index}`;
                    const caption = document.createElement("p"); caption.className = "plot-caption";
                    const zeros = document.createElement("div"); zeros.className = "zero-strip"; zeros.hidden = true;
                    const viewport = document.createElement("div"); viewport.className = "plot-viewport"; viewport.tabIndex = 0; viewport.setAttribute("aria-label", "Plot. Wide figures can be scrolled horizontally."); viewport.append(graph);
                    facet.append(heading, viewport, zeros, caption); return facet;
                }));
            }
            return [...frame.facets.children];
        }
        function zeroStrip(node, rows, descriptor, kind, settings = {}) {
            node.replaceChildren(); node.hidden = rows.length === 0;
            if (!rows.length) return;
            const key = kind === "pair" ? "pair_id" : "library_id";
            node.append(button(`Select ${rows.length} ${settings.x ? "points with zero coordinates" : "zero values"}`, () => choose(rows.map(row => row[key]), kind)));
            const points = document.createElement("div"); points.className = "zero-points";
            for (const row of rows) {
                const point = button("", event => toggle(row[key], kind, event.shiftKey), "zero-point");
                point.dataset.identity = row[key]; point.dataset.kind = kind;
                const text = document.createElement("span");
                const extra = [...(settings.x ? [[settings.x.label, format(row[settings.x.column], settings.x)]] : []), ...(settings.extra ? settings.extra(row) : [])];
                text.innerHTML = hover(row, descriptor, row[descriptor.column], extra, kind).replaceAll("<br>", "\n");
                point.title = text.textContent; point.setAttribute("aria-label", point.title); point.setAttribute("aria-pressed", "false"); points.append(point);
            }
            node.append(points);
        }
        function unavailable(frame, rows, descriptor, settings) {
            const key = settings.kind === "pair" ? "pair_id" : "library_id";
            const missing = rows.filter(row => !present(row[descriptor.column]));
            if (missing.length) {
                const ids = unique(missing.map(row => row[key]));
                const select = button(`Select ${ids.length} unavailable ${settings.kind === "pair" ? "pairs" : "libraries"}`, () => choose(ids, settings.kind));
                const reasons = unique(missing.map(row => row[descriptor.reason_column] || row.availability_reason || "No available measurement"));
                select.title = reasons.join("; "); frame.foot.append(select);
                const reason = document.createElement("span"); reason.className = "note"; reason.textContent = reasons.map(value => value.replaceAll("_", " ")).join("; "); frame.foot.append(reason);
            }
        }
        async function distribution(id, rows, descriptor, settings = {}) {
            const kind = settings.kind || "library", key = kind === "pair" ? "pair_id" : "library_id";
            const scale = report.scale(descriptor), logarithmic = scale === "log";
            const groupField = controller.context.group_field;
            const groupOf = row => kind === "pair" ? row.display_group : category(row[groupField]);
            const groupValues = kind === "pair" ? unique(rows.map(groupOf)).sort(compare) : unique(data.libraries.map(row => category(row[groupField]))).sort(compare);
            const series = settings.series || [null];
            const combinations = series.flatMap(name => groupValues.map(group => ({name, group})));
            const split = controller.context.facet && (kind === "pair" ? rows.some(row => row.display_visit !== null) : data.libraries.some(row => row.timepoint !== null));
            const visitOf = row => kind === "pair" ? row.display_visit : category(row.timepoint);
            const visits = split ? (kind === "pair" ? unique(rows.map(visitOf)).sort(compare) : orderedVisits(data.libraries)) : [null];
            if (!visits.length) visits.push(null);
            const prepared = frame(id, descriptor, settings);
            const panels = facets(prepared, visits.map(value => split ? (value === null ? "Timepoint unlabelled" : value) : ""));
            const allValues = rows.map(row => row[descriptor.column]).filter(value => present(value) && (!logarithmic || value > 0));
            for (let visitIndex = 0; visitIndex < visits.length; visitIndex++) {
                const subset = split ? rows.filter(row => visitOf(row) === visits[visitIndex]) : rows;
                const traces = [];
                combinations.forEach(({name, group}, index) => {
                    const eligible = subset.filter(row => groupOf(row) === group && (name === null || row.series === name) && present(row[descriptor.column]) && (!logarithmic || row[descriptor.column] > 0));
                    const summary = !settings.low && violinSummary(eligible.map(row => row[descriptor.column]), logarithmic);
                    if (summary) {
                        traces.push({type: "violin", x: summary.values.map(() => index), y: summary.values, width: .8, scalemode: "width", scalegroup: id, spanmode: "hard", points: false, hoverinfo: "skip", fillcolor: "#dce8eb", line: {color: "#a6bfc5", width: 1}, name: label(group)});
                        traces.push({type: "scatter", mode: "lines", x: [index - .25, index + .25], y: [summary.median, summary.median], hoverinfo: "skip", line: {color: "#304f5b", width: 2.5}, showlegend: false, name: "Median"});
                    }
                    const numericColour = colourRanges.has(controller.context.colour_field) && kind !== "pair";
                    const groups = numericColour ? [eligible.filter(row => present(row[controller.context.colour_field])), eligible.filter(row => !present(row[controller.context.colour_field]))] : [eligible];
                    for (const points of groups) {
                        if (!points.length) continue;
                        const trace = pointTrace(points, points.map(row => index + jitter(row[key])), points.map(row => logarithmic ? Math.log10(row[descriptor.column]) : row[descriptor.column]), descriptor, {...settings, name: label(group)});
                        if (settings.low) trace.error_y = {type: "data", symmetric: false, array: points.map(row => present(row[settings.high]) ? row[settings.high] - row[descriptor.column] : 0), arrayminus: points.map(row => present(row[settings.low]) ? row[descriptor.column] - row[settings.low] : 0), width: 2, thickness: .8};
                        traces.push(trace);
                    }
                });
                const title = `${descriptor.unit || descriptor.label}${logarithmic ? " (log scale)" : ""}`;
                const minWidth = Math.max(320, combinations.length * 90 + 90);
                prepared.facets.style.gridTemplateColumns = `repeat(auto-fit, minmax(min(100%, ${minWidth}px), 1fr))`;
                const graph = panels[visitIndex].querySelector(".plot"); graph.style.minWidth = `${minWidth}px`;
                const ticktext = combinations.map(({name, group}) => {
                    const members = subset.filter(row => groupOf(row) === group && (name === null || row.series === name));
                    const available = members.filter(row => present(row[descriptor.column]));
                    const plotted = available.filter(row => !logarithmic || row[descriptor.column] > 0);
                    const donors = new Set(plotted.map(row => kind === "pair" ? controller.libraries.get(row.library_id_a).donor_id : row.donor_id)).size;
                    const text = name === null ? escape(label(group)) : `${escape(name)}${groupValues.length > 1 ? "<br>" + escape(label(group)) : ""}`;
                    return `${text}<br>${available.length}/${members.length} available<br>${plotted.length} plotted<br>${donors} donors`;
                });
                await plot(graph, traces, {xaxis: {tickmode: "array", tickvals: combinations.map((_, index) => index), ticktext, range: [-.7, Math.max(.7, combinations.length - .3)], automargin: true, tickangle: 0, zeroline: false, tickfont: {size: 11}}, yaxis: axis(title, allValues, descriptor, scale, logarithmic), annotations: traces.length ? [] : [{text: "No available values", showarrow: false, xref: "paper", yref: "paper", x: .5, y: .5}]});
                zeroStrip(panels[visitIndex].querySelector(".zero-strip"), logarithmic ? subset.filter(row => row[descriptor.column] === 0) : [], descriptor, kind, settings);
                const available = subset.filter(row => present(row[descriptor.column]));
                panels[visitIndex].querySelector(".plot-caption").textContent = kind === "pair" ? `${available.length}/${subset.length} pairs available.` : `${unique(available.map(row => row.library_id)).length}/${unique(subset.map(row => row.library_id)).length} libraries available. ${countText(population(pointRows(subset)))} represented.`;
                panels[visitIndex].querySelector(".plot-caption").textContent += ` Group labels show available/total ${kind === "pair" ? "pairs" : "libraries"}.`;
                if (traces.some(trace => trace.type === "violin")) panels[visitIndex].querySelector(".plot-caption").textContent += ` Violins describe plotted observations; dark ticks mark their median${logarithmic ? " in original units, using positive values only" : ""}.`;
            }
            unavailable(prepared, rows, descriptor, settings);
            legend(prepared, rows, kind);
        }
        function legend(frame, rows, kind) {
            const field = data.metadata.find(item => item.column === controller.context.colour_field);
            if (!data.has_metadata) return;
            const node = document.createElement("span"); node.className = "plot-caption";
            if (kind === "pair" && field.type === "number") node.textContent = `${field.label} is shown for each endpoint on hover. Pair colours do not average their values.`;
            else if (field.type === "number") {
                node.classList.add("numeric-legend");
                const range = colourRanges.get(field.column);
                node.append(`${field.label}${field.unit ? " (" + field.unit + ")" : ""}: `);
                if (range) {
                    const scale = document.createElement("span"); scale.className = "colour-scale";
                    const gradient = document.createElement("span"); gradient.className = "colour-gradient"; gradient.setAttribute("aria-hidden", "true");
                    gradient.style.background = range[0] === range[1] ? continuousColour(range[0], range) : `linear-gradient(to right, ${viridis.join(", ")})`;
                    const ticks = document.createElement("span"); ticks.className = "colour-ticks";
                    for (const value of range[0] === range[1] ? [range[0]] : [range[0], range[0] / 2 + range[1] / 2, range[1]]) {
                        const tick = document.createElement("span"); tick.textContent = number(value); ticks.append(tick);
                    }
                    scale.append(gradient, ticks); node.append(scale);
                } else node.append("Unavailable. ");
                const missing = document.createElement("span"); missing.className = "missing-colour"; missing.textContent = "Grey: unavailable"; node.append(missing);
            }
            else {
                const entries = unique(rows.filter(row => !row.mixed_colour).map(row => category(row[field.column]))).sort(compare).map(value => [label(value), colourMaps.get(field.column)?.get(value) || "#89969d"]);
                if (rows.some(row => row.mixed_colour)) entries.push(["Different endpoint values", "#536971"]);
                for (const [text, colour] of entries) { const item = document.createElement("span"); item.textContent = text; item.style.borderLeft = `8px solid ${colour}`; item.style.paddingLeft = "5px"; item.style.marginRight = "12px"; node.append(item); }
            }
            frame.foot.append(node);
        }
        report.distribution = distribution;
        async function scatter(id, rows, x, y, settings = {}) {
            const xScale = report.scale(x), yScale = report.scale(y), kind = settings.kind || "library";
            const eligible = rows.filter(row => present(row[x.column]) && present(row[y.column]));
            const points = eligible.filter(row => (xScale !== "log" || row[x.column] > 0) && (yScale !== "log" || row[y.column] > 0));
            const visitOf = row => kind === "pair" ? row.display_visit : category(row.timepoint);
            const split = controller.context.facet && (kind === "pair" ? rows.some(row => row.display_visit !== null) : data.libraries.some(row => row.timepoint !== null));
            const visits = split ? (kind === "pair" ? unique(rows.map(visitOf)).sort(compare) : orderedVisits(data.libraries)) : [null];
            if (!visits.length) visits.push(null);
            const source = settings.source || unique([x, y].map(item => item.detail_source || item.source).filter(Boolean)).join("; ");
            const prepared = frame(id, y, {...settings, x, source});
            const panels = facets(prepared, visits.map(value => split ? value === null ? "Timepoint unlabelled" : value : ""));
            for (let index = 0; index < visits.length; index++) {
                const shown = split ? points.filter(row => visitOf(row) === visits[index]) : points;
                const subset = split ? rows.filter(row => visitOf(row) === visits[index]) : rows;
                const excluded = subset.filter(row => present(row[x.column]) && present(row[y.column]) && ((xScale === "log" && row[x.column] <= 0) || (yScale === "log" && row[y.column] <= 0)));
                const groups = colourRanges.has(controller.context.colour_field) && kind !== "pair" ? [shown.filter(row => present(row[controller.context.colour_field])), shown.filter(row => !present(row[controller.context.colour_field]))] : [shown];
                const extra = row => [...[x, y].flatMap(item => item.interval_low ? [[`${item.label} 95% interval`, `${number(row[item.interval_low])} to ${number(row[item.interval_high])}`], [`${item.label} evaluable clones`, number(row[item.count_column])]] : []), ...(settings.extra ? settings.extra(row) : [])];
                const traces = groups.filter(group => group.length).map(group => pointTrace(group, group.map(row => row[x.column]), group.map(row => row[y.column]), y, {...settings, x, extra: row => [[x.label, `${format(row[x.column], x)}${x.unit && x.format !== "percent" ? " " + x.unit : ""}`], ...extra(row)]}));
                let ordered = 0, connected = 0;
                if (settings.trajectories && controller.context.trajectories && !split) for (const id of controller.selected_pairs) {
                    const endpoints = controller.orderedPair(id);
                    if (!endpoints) continue;
                    ordered++;
                    const pair = endpoints.map(row => shown.find(point => point.library_id === row.library_id));
                    if (pair.every(Boolean)) { connected++; traces.unshift({type: "scatter", mode: "lines", x: pair.map(row => row[x.column]), y: pair.map(row => row[y.column]), hoverinfo: "skip", line: {color: "#9caeb5", width: 1}, showlegend: false}); }
                }
                await plot(panels[index].querySelector(".plot"), traces, {xaxis: axis(x.label + (xScale === "log" ? " (log scale)" : ""), points.map(row => row[x.column]), x, xScale), yaxis: axis(y.label + (yScale === "log" ? " (log scale)" : ""), points.map(row => row[y.column]), y, yScale), ...settings.layout});
                panels[index].querySelector(".plot-caption").textContent = `${shown.length}/${subset.length} ${kind === "pair" ? "pairs" : "libraries"} plotted. ${excluded.length ? `${excluded.length} zero coordinates are outside the logarithmic axes.` : ""}${settings.trajectories && controller.context.trajectories && split ? " Turn off timepoint panels to connect selected pairs across visits." : ""}`;
                if (settings.trajectories && controller.context.trajectories && !split && controller.selected_pairs.size) panels[index].querySelector(".plot-caption").textContent += ordered ? ` ${connected}/${ordered} selected ordered pairs connected.${connected < ordered ? " Other endpoints are filtered, unavailable or zero on a log axis." : ""}` : " Selected pairs have no distinct supplied visit order.";
                zeroStrip(panels[index].querySelector(".zero-strip"), excluded, y, kind, {...settings, x, extra});
            }
            unavailable(prepared, rows.map(row => ({...row, [y.column]: present(row[x.column]) ? row[y.column] : null})), y, settings);
            legend(prepared, rows, settings.kind);
        }
        report.scatter = scatter;
        report.simple = async (id, descriptor, traces, layout = {}, settings = {}) => { const prepared = frame(id, descriptor, settings); const panel = facets(prepared, [""])[0]; await plot(panel.querySelector(".plot"), traces, layout, Boolean(settings.linked)); return panel; };

        report.cohortCurves = async (libraries, rows, descriptor) => {
            const group = controller.context.group_field;
            const groups = unique(libraries.map(row => category(row[group]))).sort(compare);
            const visits = controller.context.facet ? orderedVisits(libraries) : [null];
            const combinations = groups.flatMap(value => visits.map(visit => ({value, visit})));
            const prepared = frame("cohort-curves", descriptor, {source: "rarefaction_curves.tsv", description: "Each line connects supplied depths for one library. Curve points remain available below the configured diversity minimum. No extrapolation or pooled estimate is shown."});
            const panels = facets(prepared, combinations.length ? combinations.map(({value, visit}) => `${label(value)}${controller.context.facet ? " | " + label(visit) : ""}`) : [""]);
            const maxDepth = rows.reduce((maximum, row) => Math.max(maximum, row.support_depth), 1);
            const maxValue = rows.reduce((maximum, row) => present(row[descriptor.column]) ? Math.max(maximum, row[descriptor.column]) : maximum, 1);
            for (const [index, panel] of panels.entries()) {
                const combination = combinations[index];
                const members = combination ? libraries.filter(row => category(row[group]) === combination.value && (!controller.context.facet || category(row.timepoint) === combination.visit)) : [];
                const traces = members.map(library => {
                    const curve = rows.filter(row => row.library_id === library.library_id && present(row[descriptor.column])).sort((a, b) => a.support_depth - b.support_depth);
                    const trace = pointTrace(curve, curve.map(row => row.support_depth), curve.map(row => row[descriptor.column]), descriptor, {extra: row => [["Subsample read support", number(row.support_depth)], ["Curve point", row.is_endpoint ? "Observed endpoint" : row.is_target ? "Configured target" : "Interpolation"]]});
                    return {...trace, mode: "lines+markers", name: library.library_id, line: {width: 1, color: colour(library)}, marker: {...trace.marker, size: curve.map(row => row.is_endpoint || row.is_target ? 6 : 3)}};
                }).filter(trace => trace.x.length);
                await plot(panel.querySelector(".plot"), traces, {xaxis: {title: {text: "Subsample read support"}, range: [0, maxDepth * 1.04], automargin: true}, yaxis: {title: {text: "Effective groups"}, range: [0, maxValue * 1.04], automargin: true}, annotations: traces.length ? [] : [{text: "No available curves", showarrow: false, xref: "paper", yref: "paper", x: .5, y: .5}]});
                panel.querySelector(".plot-caption").textContent = `${traces.length}/${members.length} libraries have curves. ${countText(population(members))} represented.`;
            }
            legend(prepared, libraries, "library");
        };

        report.donorTrends = async (rows, visits, descriptor, extra) => {
            const group = controller.context.group_field;
            const groups = unique(rows.map(row => category(row[group]))).sort(compare);
            const prepared = frame("donor-trends", descriptor, {title: `${controller.context.locus} ${descriptor.label} from ${controller.context.reference_timepoint || "a reference timepoint"}`, source: "donor_overlap_pairs.tsv", description: "Each measured pair retains its reference library and donor. Repeated libraries are separate observations, and pairs from one donor are not independent biological replicates."});
            const panels = facets(prepared, groups.length ? groups.map(value => `Reference: ${label(value)}`) : [""]);
            for (const [index, panel] of panels.entries()) {
                const subset = rows.filter(row => category(row[group]) === groups[index]);
                const traces = [];
                for (const reference of unique(subset.map(row => row.reference_id))) {
                    const pairs = subset.filter(row => row.reference_id === reference);
                    // A line needs one observation at each successive supplied visit.
                    // Replicates and missing measurements remain points or gaps.
                    const path = visits.map(visit => {
                        const matches = pairs.filter(row => row.followup_timepoint === visit.label);
                        return matches.length === 1 && present(matches[0][descriptor.column]) ? matches[0] : null;
                    });
                    const available = pairs.filter(row => present(row[descriptor.column]));
                    const lineColour = typeof colour(pairs[0]) === "string" ? colour(pairs[0]) : "#a6bfc5";
                    traces.push({type: "scatter", mode: "lines", name: reference, x: visits.map(visit => visit.order), y: path.map(row => row ? row[descriptor.column] : null), connectgaps: false, hoverinfo: "skip", line: {width: 1, color: lineColour}});
                    traces.push(pointTrace(available, available.map(row => row.followup_order), available.map(row => row[descriptor.column]), descriptor, {kind: "pair", extra: row => [["Reference library", row.reference_id], ...extra(row)]}));
                }
                await plot(panel.querySelector(".plot"), traces, {xaxis: {title: {text: "Follow-up timepoint (supplied order)"}, tickmode: "array", tickvals: visits.map(visit => visit.order), ticktext: visits.map(visit => escape(visit.label)), automargin: true}, yaxis: axis(descriptor.label, [], descriptor, "linear"), annotations: subset.some(row => present(row[descriptor.column])) ? [] : [{text: "No available ordered comparisons", showarrow: false, xref: "paper", yref: "paper", x: .5, y: .5}]});
                panel.querySelector(".plot-caption").textContent = `${subset.filter(row => present(row[descriptor.column])).length}/${subset.length} comparisons available. ${countText(population(subset.map(row => controller.libraries.get(row.reference_id))))} represented by the reference libraries.`;
            }
            unavailable(prepared, rows, descriptor, {kind: "pair"});
            legend(prepared, rows, "pair");
        };
        async function restyleSelection() {
            const updates = [];
            const highlightedPairs = controller.highlightedPairs();
            for (const div of document.querySelectorAll(".plot.js-plotly-plot")) {
                if (!div.checkVisibility()) continue;
                if (div._fullLayout.width !== div.clientWidth) await Plotly.Plots.resize(div);
                // Clear the native outline as well as shared point highlights.
                if (!controller.selected.size && div.layout.selections?.length) await Plotly.relayout(div, {selections: []});
                const indices = div.data.flatMap((trace, index) => trace.meta?.kind ? [index] : []);
                if (indices.length) updates.push(Plotly.restyle(div, {selectedpoints: indices.map(index => { const trace = div.data[index], selected = trace.meta.kind === "pair" ? highlightedPairs : controller.selected; return controller.selected.size ? trace.customdata.flatMap((value, point) => selected.has(value[0]) ? [point] : []) : null; })}, indices));
            }
            for (const point of document.querySelectorAll(".zero-point")) point.setAttribute("aria-pressed", (point.dataset.kind === "pair" ? highlightedPairs : controller.selected).has(point.dataset.identity) ? "true" : "false");
            await Promise.all(updates);
        }

        function navigate(chapter, fromHistory = false) {
            if (!chapters.includes(chapter)) chapter = "signal";
            if (!fromHistory && location.hash !== `#${chapter}`) history.pushState(null, "", `#${chapter}`);
            controller.setContext({chapter});
            window.scrollTo({top: 0, behavior: "instant"});
        }
        const chapters = ["signal", "diversity", "antibody", "sharing", "methods"];
        report.navigate = navigate;
        for (const link of document.querySelectorAll("[data-chapter][href]")) link.onclick = event => { event.preventDefault(); navigate(link.dataset.chapter); };
        root.addEventListener("popstate", () => navigate(location.hash.slice(1), true));
        if (chapters.includes(location.hash.slice(1))) controller.context.chapter = location.hash.slice(1);
        byId("previous-chapter").onclick = () => navigate(chapters[chapters.indexOf(controller.context.chapter) - 1]);
        byId("next-chapter").onclick = () => navigate(chapters[chapters.indexOf(controller.context.chapter) + 1]);
        byId("scope-title").textContent = data.title;
        document.title = data.title;
        byId("scope-counts").textContent = countText(counts(data.libraries));
        options("group-field", categorical.map(field => [field.column, field.label]), controller.context.group_field);
        options("colour-field", data.metadata.map(field => [field.column, field.label]), controller.context.colour_field);
        options("filter-field", categorical.map(field => [field.column, field.label]), controller.context.group_field);
        options("numeric-field", [["", "Choose a numeric field"], ...numeric.map(field => [field.column, `${field.label}${field.unit ? " (" + field.unit + ")" : ""}`])], "");
        const metricChoices = data.metrics.map(item => [item.column, item.label]);
        for (const [id, field] of [["x-metric", "x_metric"], ["y-metric", "y_metric"], ["group-field", "group_field"], ["colour-field", "colour_field"]]) {
            if (id.endsWith("metric")) options(id, metricChoices, controller.context[field]);
            byId(id).onchange = event => controller.setContext({[field]: event.target.value});
        }
        function filterOptions() {
            const field = byId("filter-field").value;
            const filter = controller.filters.get(field);
            byId("group-filter").replaceChildren(...unique(data.libraries.map(row => category(row[field]))).sort(compare).map(value => {
                const item = document.createElement("label"); item.className = "check";
                const input = document.createElement("input"); input.type = "checkbox"; input.value = JSON.stringify(value);
                input.checked = !filter || filter.values.has(value); item.append(input, document.createTextNode(label(value))); return item;
            }));
        }
        filterOptions(); byId("filter-field").onchange = filterOptions;
        for (const [id, value] of [["category-all", true], ["category-none", false]]) byId(id).onclick = () => { for (const input of byId("group-filter").querySelectorAll("input")) input.checked = value; };
        byId("apply-category-filter").onclick = () => {
            const inputs = [...byId("group-filter").querySelectorAll("input")], chosen = inputs.filter(input => input.checked);
            controller.setCategoricalFilter(byId("filter-field").value, chosen.length === inputs.length ? null : chosen.map(input => JSON.parse(input.value)));
        };
        function numericFilterOptions() {
            const filter = controller.filters.get(byId("numeric-field").value);
            byId("numeric-min").value = filter?.minimum ?? ""; byId("numeric-max").value = filter?.maximum ?? "";
            byId("numeric-missing").checked = filter?.include_missing || false; byId("numeric-max").setCustomValidity("");
        }
        byId("numeric-field").onchange = numericFilterOptions;
        byId("apply-numeric-filter").onclick = () => {
            const field = byId("numeric-field").value;
            if (!field) return;
            const low = byId("numeric-min").value === "" ? null : Number(byId("numeric-min").value), high = byId("numeric-max").value === "" ? null : Number(byId("numeric-max").value);
            byId("numeric-max").setCustomValidity(low !== null && high !== null && low > high ? "Maximum must be at least the minimum." : "");
            if (!byId("numeric-max").reportValidity()) return;
            controller.setNumericFilter(field, low, high, byId("numeric-missing").checked);
        };
        for (const [id, target] of [["filter-toggle", "filter-panel"], ["selection-toggle", "selection-panel"]]) byId(id).onclick = () => { byId(target).hidden = !byId(target).hidden; byId(id).setAttribute("aria-expanded", String(!byId(target).hidden)); if (target === "selection-panel") renderSelected(); };
        byId("facet").onchange = event => controller.setContext({facet: event.target.checked});
        byId("comparison-facet").onchange = event => controller.setContext({facet: event.target.checked});
        byId("trajectories").onchange = event => controller.setContext({trajectories: event.target.checked});
        byId("selection-only").onchange = event => controller.setSelectionOnly(event.target.checked);
        byId("clear-selection").onclick = () => controller.clearSelection();
        byId("expand-donors").onclick = () => controller.expandDonors();
        byId("selection-search").oninput = renderSelected;
        byId("reset").onclick = () => {
            controller.reset(); history.replaceState(null, "", "#signal");
            for (const id of ["numeric-field", "numeric-min", "numeric-max", "selection-search", "donor-search", "donor-selection", "tissue-selection"]) byId(id).value = "";
            byId("donor-search").dispatchEvent(new Event("input"));
            byId("numeric-missing").checked = false; byId("numeric-max").setCustomValidity("");
            byId("filter-field").value = controller.context.group_field; filterOptions();
            for (const dialog of document.querySelectorAll("dialog[open]")) dialog.close();
            for (const section of document.querySelectorAll("details[open]")) section.open = false;
            for (const [id, target] of [["filter-toggle", "filter-panel"], ["selection-toggle", "selection-panel"]]) { byId(target).hidden = true; byId(id).setAttribute("aria-expanded", "false"); }
            window.scrollTo({top: 0, behavior: "instant"});
        };
        root.addEventListener("keydown", event => { if (event.key === "Escape" && !document.querySelector("dialog[open]") && !["INPUT", "SELECT", "TEXTAREA"].includes(event.target.tagName)) controller.clearSelection(); });
        byId("open-comparison").onclick = () => { byId("comparison-dialog").showModal(); schedule(); };
        byId("open-inspection").onclick = () => {
            if (controller.selected.size === 1) controller.openDetail([...controller.selected][0]);
            byId("inspection-dialog").showModal(); schedule();
        };
        byId("comparison-choose-pairs").onclick = () => { byId("comparison-dialog").close(); navigate("sharing"); };
        for (const close of document.querySelectorAll("[data-close]")) close.onclick = () => byId(close.dataset.close).close();
        for (const dialog of document.querySelectorAll("dialog")) dialog.addEventListener("close", schedule);
        for (const node of document.querySelectorAll(".metadata-control")) node.hidden = !data.has_metadata;
        for (const node of document.querySelectorAll(".numeric-control")) node.hidden = !numeric.length;
        for (const node of document.querySelectorAll(".timepoint-control")) node.hidden = !data.libraries.some(row => row.timepoint !== null);
        byId("filter-toggle").hidden = !data.has_metadata;
        function renderSelected() {
            const search = byId("selection-search").value.toLowerCase();
            byId("selected-libraries").replaceChildren(...controller.selectedRows().filter(row => [row.library_id, row.tissue_id, row.donor_id].some(value => value.toLowerCase().includes(search))).map(row => {
                const item = document.createElement("li"), text = document.createElement("span"); text.textContent = `${row.library_id} | ${row.tissue_id} | ${row.donor_id}`;
                item.append(text, button("Inspect", () => { controller.openDetail(row.library_id); byId("inspection-dialog").showModal(); }), button("Remove", () => controller.toggleLibrary(row.library_id, true))); return item;
            }));
        }
        function renderControls() {
            const chapter = controller.context.chapter;
            for (const node of document.querySelectorAll(".chapter")) node.hidden = node.dataset.chapter !== chapter;
            for (const node of document.querySelectorAll(".chapter-nav > a")) { if (node.dataset.chapter === chapter) node.setAttribute("aria-current", "page"); else node.removeAttribute("aria-current"); }
            byId("previous-chapter").disabled = chapters.indexOf(chapter) === 0;
            byId("next-chapter").disabled = chapters.indexOf(chapter) === chapters.length - 1;
            const visible = controller.visibleRows(), selected = controller.selectedRows();
            if (hadSelection && !selected.length) {
                for (const id of ["donor-search", "donor-selection", "tissue-selection"]) byId(id).value = "";
                byId("donor-search").dispatchEvent(new Event("input"));
            }
            hadSelection = selected.length > 0;
            byId("view-counts").textContent = `${countText(counts(visible))} visible. ${selected.length} selected${selected.filter(row => !controller.matches(row)).length ? ` (${selected.filter(row => !controller.matches(row)).length} outside filters)` : ""}.`;
            byId("selected-pairs").textContent = controller.selected_pairs.size ? `${controller.selected_pairs.size} selected pairs. Overlapping endpoints count once.` : "No pairs selected.";
            byId("clear-selection").disabled = !selected.length && !controller.selection_only;
            byId("expand-donors").disabled = !selected.length;
            byId("selection-only").checked = controller.selection_only;
            byId("facet").checked = controller.context.facet;
            byId("comparison-facet").checked = controller.context.facet;
            byId("trajectories").checked = controller.context.trajectories;
            for (const [id, field] of [["x-metric", "x_metric"], ["y-metric", "y_metric"], ["group-field", "group_field"], ["colour-field", "colour_field"]]) byId(id).value = controller.context[field];
            byId("active-filters").replaceChildren(...[...controller.filters].map(([field, filter]) => {
                const name = data.metadata.find(item => item.column === field)?.label || field;
                const text = filter.type === "category" ? [...filter.values].map(label).join(", ") || "No categories" : filter.minimum === null && filter.maximum === null ? "Available values only" : `${filter.minimum ?? "Any"} to ${filter.maximum ?? "any"}; ${filter.include_missing ? "include" : "exclude"} missing`;
                return button(`${name}: ${text} | Clear`, () => { controller.setCategoricalFilter(field, null); filterOptions(); numericFilterOptions(); });
            }));
            if (!byId("selection-panel").hidden) renderSelected();
        }
        return report;
    }
    root.ReportView = {mount, byId, escape, label, options, metric, number, format, axis, countText, orderedVisits, orderedComparisons, compare, palette, colourRange, continuousColour, violinSummary};
})(window);
