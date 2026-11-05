"""Render Snakemake DOT from stdin as a simplified SVG on stdout.

Graphviz removes transitive edges for display only. Every workflow rule and
dependency path remains; the aggregate target is omitted. No jobs are run.
"""

import json
from html import escape
import subprocess
import sys
import textwrap
from xml.etree import ElementTree as ET


# Stage names and colours match the Mermaid overview in README.md. Assign by
# scientific purpose rather than source file: one rule file can span stages.
STAGES = (
    ("assembly", "Receptor assembly", "#e0f2fe", (
        "run_trust4",
    )),
    ("annotation", "Annotation and quality checks", "#fef3c7", (
        "reannotate_igblast", "filter_airr",
    )),
    ("igh", "IGH analysis", "#fce7f3", (
        "prepare_igh_grouping", "pool_scope_donor_igh", "fit_donor_threshold",
        "resolve_scope_clone_threshold", "assign_scope_donor_clones",
        "qualify_isotype_calls", "calculate_isotype", "build_igh_consensus_votes",
        "calculate_igh_features", "calculate_donor_igh_summaries",
        "calculate_donor_selection",
    )),
    ("other", "Other receptor chains", "#ffedd5", (
        "build_clonotypes",
    )),
    ("airr", "AIRR annotations", "#ccfbf1", (
        "export_library_airr", "annotate_scope_airr", "build_scope_repertoires",
    )),
    ("summaries", "Repertoire summaries", "#ede9fe", (
        "summarise_library_yield", "calculate_library_clonotype_diversity",
        "calculate_scope_clone_diversity", "calculate_donor_overlap",
        "collect_scope_library_metric", "collect_scope_donor_metric",
        "collect_scope_rarefaction_curves", "collect_scope_diversity",
        "collect_scope_gene_usage", "collect_scope_vj_pair_usage",
        "collect_scope_group_abundance", "collect_scope_qc_audit",
        "collect_scope_overlap", "collect_scope_selection", "build_scope_headline",
    )),
    ("release", "Release and report", "#dcfce7", (
        "prepare_scope_report", "render_scope_report", "assemble_scope_release",
        "write_scope_intermediate_readme",
    )),
    ("support", "Workflow support", "#eef1f5", (
        "count_sequencing_depth", "normalise_input_depth", "build_scope_manifest",
        "build_scope_assignment", "obtain_analysis_image", "obtain_assembly_image",
        "obtain_imgt_reference", "audit_imgt_reference", "capture_tool_versions",
        "build_scope_provenance",
    )),
)


def graphviz(command: list[str], source: str) -> str:
    return subprocess.run(
        command, input=source, text=True, stdout=subprocess.PIPE, check=True
    ).stdout


def render(source: str) -> str:
    graph = json.loads(graphviz(["dot", "-Tdot_json"], source))
    nodes = {
        node["_gvid"]: node["label"] for node in graph["objects"]
        if node["label"] not in {"all", "post_trust4_all"}
    }
    stages = {
        rule: (stage, label, colour)
        for stage, label, colour, rules in STAGES for rule in rules
    }
    missing = set(nodes.values()) - stages.keys()
    if missing:
        raise ValueError(f"Assign these rules to README overview stages: {sorted(missing)}")
    legend = []
    for start in (0, 4):
        legend.append("<TR>" + "".join(
            f'<TD BGCOLOR="{colour}"> {escape(label)} </TD>'
            for _, label, colour, _ in STAGES[start:start + 4]
        ) + "</TR>")

    lines = [
        "digraph airrborne {",
        'graph [rankdir=TB, bgcolor="white", pad=0.35, nodesep=0.22, ranksep=0.42,',
        'fontname="sans-serif", fontcolor="#243746", labelloc=t, label=<',
        '<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="8">',
        '<TR><TD COLSPAN="4"><FONT POINT-SIZE="24"><B>airrborne</B></FONT></TD></TR>',
        '<TR><TD COLSPAN="4">Workflow rules grouped by the README analysis stages</TD></TR>',
        *legend,
        '<TR><TD COLSPAN="4"><FONT POINT-SIZE="11">Redundant arrows omitted; all dependency paths retained.',
        '<BR/>Workflow support covers inputs, software, references and provenance.</FONT></TD></TR>',
        '</TABLE>>];',
        'node [shape=box, style="rounded,filled", color="#9aabb8",',
        'fontname="sans-serif", fontcolor="#243746", fontsize=11, margin="0.14,0.10"];',
        'edge [color="#8e9eac", arrowsize=0.6, penwidth=1];',
    ]
    for node_id, rule in sorted(nodes.items(), key=lambda item: item[1]):
        label = textwrap.fill(rule.replace("_", " "), width=25, break_long_words=False)
        lines.append(
            f'{node_id} [label={json.dumps(label)}, tooltip={json.dumps(rule)}, '
            f'fillcolor="{stages[rule][2]}"];'
        )
    for edge in graph.get("edges", []):
        if edge["tail"] in nodes and edge["head"] in nodes:
            lines.append(f'{edge["tail"]} -> {edge["head"]};')
    lines.append("}")
    reduced = graphviz(["tred"], "\n".join(lines))
    svg = ET.fromstring(graphviz(["dot", "-Tsvg"], reduced))
    namespace = "http://www.w3.org/2000/svg"
    ET.register_namespace("", namespace)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    for group in svg.findall(f".//{{{namespace}}}g"):
        if group.get("class") == "node":
            title = group.find(f"{{{namespace}}}title")
            rule = nodes[int(title.text)]
            # Graphviz's default <title> is a numeric DOT id. Supply the rule
            # name as the native SVG tooltip and accessible name instead.
            title.text = rule
            group.set("data-stage", stages[rule][0])
    return ET.tostring(svg, encoding="unicode") + "\n"


if __name__ == "__main__":
    sys.stdout.write(render(sys.stdin.read()))
