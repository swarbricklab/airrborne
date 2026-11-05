# Record descriptive source state and the evidence for one scope's handover.
# provenance.py captures source changes, image and reference identities and
# assembly history. Library processing identity is defined separately by
# export_library_airr from consumed code, settings, image and resources.
# SOURCE_RECORD is a provenance parameter, outside scientific rule inputs.
from workflow.lib.provenance import capture_source

SOURCE_RECORD = capture_source(Path(WORKFLOW_DIR))


def _scope_provenance_config(wildcards):
    return {
        "params": config["params"], "refs": config["refs"],
        "deps": {key: config["deps"][key] for key in ("immcantation_docker", "trust4_container")},
    }


# Versions are measured once per actual image identity and cached under
# params.cache, which Snakemake does not declare. Changing version_queries()
# makes the cached record fail its query check: rerun capture_tool_versions.py
# with --refresh, or delete that cache directory, to measure the new queries.
rule capture_tool_versions:
    input:
        script = str(SCRIPTS_DIR / "capture_tool_versions.py"),
        module = str(Path(WORKFLOW_DIR) / "workflow/lib/provenance.py"),
        container_exec = CONTAINER_EXEC,
        # Each obtained image's receipt already holds its identity. A run that
        # assembles obtains the assembly image. A --post-trust4 run does not,
        # so the configured assembly reference is identified here instead and
        # may be unavailable.
        image_record = lambda wc: ([ANALYSIS_IMAGE_RECEIPT] if wc.role == "analysis"
                                   else [] if ANALYSIS_ONLY else [ASSEMBLY_IMAGE_RECEIPT]),
    output:
        versions = results_dir / "shared/provenance/{role}_versions.json",
    wildcard_constraints:
        role = "analysis|assembly",
    params:
        image = lambda wc: (str(ANALYSIS_SIF) if wc.role == "analysis"
                            else config["deps"]["trust4_container"] if ANALYSIS_ONLY
                            else str(ASSEMBLY_SIF)),
        docker_image = lambda wc: ANALYSIS_IMAGE_REFERENCE if wc.role == "analysis" else ASSEMBLY_IMAGE_REFERENCE,
        record_option = lambda wc, input: ("--image-record " + shlex.quote(input.image_record[0])
                                           if input.image_record else ""),
        cache = str(results_dir / "shared/provenance/tool_cache"),
        runtime = os.environ.get("AIRRBORNE_CONTAINER_RUNTIME", ""),
    benchmark:
        log_dir / "benchmarks/capture_tool_versions/{role}.tsv"
    log:
        log_dir / "capture_tool_versions/{role}.log"
    shell:
        "python {input.script:q} --image {params.image:q} --docker-image {params.docker_image:q} "
        "{params.record_option} --role {wildcards.role:q} --cache-dir {params.cache:q} --output {output.versions:q} "
        "--container-exec {input.container_exec:q} > {log:q} 2>&1"


# Record the scope's handover provenance: source state, image and reference
# identities, the reference receipt, measured tool versions and each library's
# assembly history in provenance.json, and the params, refs, image settings and
# fixed methods in effective_config.yaml.
rule build_scope_provenance:
    input:
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        analysis_versions = results_dir / "shared/provenance/analysis_versions.json",
        assembly_versions = results_dir / "shared/provenance/assembly_versions.json",
        reference = rules.obtain_imgt_reference.output.imgt_ref,
        reference_metadata = rules.obtain_imgt_reference.output.metadata,
        auxiliary = rules.obtain_imgt_reference.output.auxiliary,
        reference_receipt = IMGT_REFERENCE_RECEIPT,
        module = str(Path(WORKFLOW_DIR) / "workflow/lib/provenance.py"),
        methods = str(Path(WORKFLOW_DIR) / "workflow/lib/release_methods.py"),
        assembly = lambda wc: [str(results_dir / "libraries" / library / "trust4" / (library + suffix))
                               for library in sorted(SCOPE_MEMBERS[wc.scope])
                               for suffix in ("_airr.tsv", "_cdr3.out", "_final.out")],
        # A reused assembly has no sidecar, so this list holds only the records
        # that already exist when the DAG is built and never requests one. The
        # script looks the same paths up again at run time, which is how a
        # sidecar written by run_trust4 in this same run is still recorded.
        producer_records = lambda wc: [str(path) for library in sorted(SCOPE_MEMBERS[wc.scope])
                                       for path in [results_dir / "libraries" / library / "trust4" /
                                                    (library + "_assembly_provenance.json")]
                                       if path.exists()],
    output:
        provenance = results_dir / "scopes/{scope}/provenance.json",
        config = results_dir / "scopes/{scope}/effective_config.yaml",
    params:
        effective_config = _scope_provenance_config,
        source = SOURCE_RECORD,
        repo = str(Path(WORKFLOW_DIR)),
        source_evidence = str(results_dir / "shared/provenance/source_evidence"),
        results_dir = str(results_dir),
    benchmark:
        log_dir / "benchmarks/build_scope_provenance/{scope}.tsv"
    log:
        log_dir / "build_scope_provenance/{scope}.log"
    script:
        "../scripts/build_scope_provenance.py"
