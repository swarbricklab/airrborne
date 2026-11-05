# Build the offline report from supplied scientific measurements and calibration.
# report_data.py prepares display data; report_metadata.py attaches optional
# project annotations only here. Rendering owns browser assets and methods text,
# so layout changes do not re-read measurement tables.
import json

REPORT_METRICS = (
    "headline", "shm", "selection_pressure", "isotype_usage", "diversity",
    "donor_overlap_pairs", "rarefaction_curves", "group_abundance", "qc", "qc_trail", "clone_shm",
    "gene_usage", "yield", "isotype_summary",
)

# prepare() checks table schemas and keys while retaining unavailable estimates.
# Metadata file and column-selection changes refresh preparation, rendering and
# publication without becoming inputs to scientific analyses or AIRR exports.
rule prepare_scope_report:
    input:
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        metrics = lambda wc: [str(results_dir / "scopes" / wc.scope / "metrics" /
                                 (metric + ".tsv")) for metric in REPORT_METRICS],
        metadata = lambda wc: [source["file"] for source in config["scope_metadata"].get(wc.scope, {}).values()],
        calibration = results_dir / "scopes/{scope}/clones/{scope}_clone_threshold.tsv",
        donor_calibration = results_dir / "scopes/{scope}/clones/{scope}_donor_clone_thresholds.tsv",
        diagnostics = results_dir / "scopes/{scope}/clones/{scope}_clone_distance_diagnostics.tsv",
        script = str(SCRIPTS_DIR / "prepare_scope_report.py"),
        module = str(Path(WORKFLOW_DIR) / "workflow/lib/report_data.py"),
        headline_module = str(Path(WORKFLOW_DIR) / "workflow/lib/headline.py"),
        metadata_module = str(Path(WORKFLOW_DIR) / "workflow/lib/report_metadata.py"),
        bundle_module = str(Path(WORKFLOW_DIR) / "workflow/lib/release_bundle.py"),
    output:
        data = directory(results_dir / "scopes/{scope}/intermediate/report_data"),
    params:
        scope_dir = lambda wc: str(results_dir / "scopes" / wc.scope),
        metadata_sources = lambda wc: json.dumps(config["scope_metadata"].get(wc.scope, {})),
    benchmark:
        log_dir / "benchmarks/prepare_scope_report/{scope}.tsv"
    log:
        log_dir / "prepare_scope_report/{scope}.log"
    shell:
        "python {input.script:q} --scope-dir {params.scope_dir:q} --output-dir {output.data:q} "
        "--metadata-sources {params.metadata_sources:q} --calibration {input.calibration:q} "
        "{input.donor_calibration:q} {input.diagnostics:q} > {log:q} 2>&1"


# render_scope_report.py copies prepared data and local assets into the report
# directory and generates methods from effective settings and calibration.
# The release publisher waits for this directory alongside the scientific tables.
rule render_scope_report:
    input:
        data = rules.prepare_scope_report.output.data,
        provenance = results_dir / "scopes/{scope}/provenance.json",
        script = str(SCRIPTS_DIR / "render_scope_report.py"),
        methods = str(Path(WORKFLOW_DIR) / "workflow/lib/release_methods.py"),
        bundle_module = str(Path(WORKFLOW_DIR) / "workflow/lib/release_bundle.py"),
        effective_config = results_dir / "scopes/{scope}/effective_config.yaml",
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        threshold = results_dir / "scopes/{scope}/clones/{scope}_clone_threshold.tsv",
        template = str(Path(WORKFLOW_DIR) / "workflow/report/index.html"),
        assets = [str(Path(WORKFLOW_DIR) / "workflow/report/assets" / name) for name in
                  ("plotly.min.js", "LICENSE.plotly.txt", "selection_controller.js",
                   "report_view.js", "report_sections.js", "report.css")],
    output:
        report = directory(results_dir / "scopes/{scope}/report"),
    conda:
        "../../envs/release.yaml"
    benchmark:
        log_dir / "benchmarks/render_scope_report/{scope}.tsv"
    log:
        log_dir / "render_scope_report/{scope}.log"
    shell:
        "python {input.script:q} --data-dir {input.data:q} --output-dir {output.report:q} "
        "--provenance {input.provenance:q} --effective-config {input.effective_config:q} "
        "--manifest {input.manifest:q} --threshold {input.threshold:q} > {log:q} 2>&1"
