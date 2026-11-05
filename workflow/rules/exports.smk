# Export shared library annotations, link scope metadata and publish release/.
# airr_exports.py owns export validation; repertoire_metadata.py owns linkage
# checks. release_bundle.py stages and replaces the complete handover.

# Write each library's common-QC AIRR table once as the compressed export every
# scope's release copies, with the grouping links, qualified C calls and their
# isotypes, and the MakeDb germline coordinates and V-domain regions added,
# and record its processing identity: the consumed code, the settings below, the
# analysis image, the reference files and the AIRR schema.
# export_library_airr() retains every common-QC record, including records without
# a grouping link or qualified C call. Scope clone assignments and project
# annotations do not enter this shared library export or its processing identity.
rule export_library_airr:
    input:
        airr = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        igblast_hits = results_dir / "libraries/{library_id}/airr/{library_id}_igblast_hits.fmt7.gz",
        makedb = results_dir / "libraries/{library_id}/airr/{library_id}_makedb_pass.tsv.gz",
        igh_grouping = results_dir / "libraries/{library_id}/metrics/{library_id}_igh_grouping.tsv",
        nonigh_grouping = results_dir / "libraries/{library_id}/metrics/{library_id}_nonigh_grouping.tsv",
        isotype_calls = results_dir / "libraries/{library_id}/metrics/{library_id}_isotype_calls.tsv",
        reference = rules.obtain_imgt_reference.output.imgt_ref,
        reference_metadata = rules.obtain_imgt_reference.output.metadata,
        auxiliary = rules.obtain_imgt_reference.output.auxiliary,
        export_module = str(Path(WORKFLOW_DIR) / "workflow/lib/airr_exports.py"),
        provenance_module = str(Path(WORKFLOW_DIR) / "workflow/lib/provenance.py"),
        schema = str(Path(WORKFLOW_DIR) / "workflow/schemas/airr-2.0.0.yaml"),
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        code = LIBRARY_PROCESSING_CODE,
    output:
        airr = results_dir / "libraries/{library_id}/export/{library_id}_rearrangements.tsv.gz",
        identity = results_dir / "libraries/{library_id}/export/{library_id}_processing.json",
    params:
        settings = {"igblast": config["params"]["igblast"],
                    "max_alignment_n_fraction": config["params"]["airr_qc"]["max_alignment_n_fraction"],
                    "min_informative_c_alignment_positions": config["params"]["airr_qc"]["min_informative_c_alignment_positions"],
                    "max_n": config["params"]["lineage"]["max_n"]},
    conda:
        "../../envs/release.yaml"
    benchmark:
        log_dir / "benchmarks/export_library_airr/{library_id}.tsv"
    log:
        log_dir / "export_library_airr/{library_id}.log"
    script:
        "../scripts/export_library_airr.py"


# build_scope_repertoires() writes one Repertoire entry per manifest library
# using scope-relative AIRR file names and that library's processing identity.
# Missing clinical and preparation fields remain missing and appear as strict
# metadata gaps. This job checks metadata; the publisher validates staged files.
rule build_scope_repertoires:
    input:
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        provenance = results_dir / "scopes/{scope}/provenance.json",
        identities = lambda wc: [str(results_dir / "libraries" / library / "export" /
                                    (library + "_processing.json"))
                                for library in sorted(SCOPE_MEMBERS[wc.scope])],
        module = str(Path(WORKFLOW_DIR) / "workflow/lib/repertoire_metadata.py"),
        alignment_module = str(Path(WORKFLOW_DIR) / "workflow/lib/alignment_export.py"),
        airr_module = str(Path(WORKFLOW_DIR) / "workflow/lib/airr_exports.py"),
        effective_config = results_dir / "scopes/{scope}/effective_config.yaml",
        schema = str(Path(WORKFLOW_DIR) / "workflow/schemas/airr-2.0.0.yaml"),
    output:
        repertoire = results_dir / "scopes/{scope}/airr/Repertoire.json",
    params:
        libraries = lambda wc: sorted(SCOPE_MEMBERS[wc.scope]),
    conda:
        "../../envs/release.yaml"
    benchmark:
        log_dir / "benchmarks/build_scope_repertoires/{scope}.tsv"
    log:
        log_dir / "build_scope_repertoires/{scope}.log"
    script:
        "../scripts/build_scope_repertoires.py"


# Write each scope's copy of a library export with the scope's member
# assignments: the IGH clone_id of every member and the clonotype_id of every
# non-IGH member, in clone_id. Contig siblings and records outside grouping
# keep an empty clone_id. The scope's release links these copies.
def _scope_library_donor(scope, library_id):
    return {library: donor for library, donor, _ in scope_manifest_by_scope[scope]}[library_id]


rule annotate_scope_airr:
    input:
        airr = results_dir / "libraries/{library_id}/export/{library_id}_rearrangements.tsv.gz",
        clone_members = lambda wc: str(results_dir / "scopes" / wc.scope / "intermediate" /
                                       _scope_library_donor(wc.scope, wc.library_id) /
                                       (_scope_library_donor(wc.scope, wc.library_id)
                                        + "_clone_members.tsv")),
        clonotype_members = results_dir / "libraries/{library_id}/metrics/{library_id}_clonotype_members.tsv",
        alignment_module = str(LIB_DIR / "alignment_export.py"),
        export_module = str(Path(WORKFLOW_DIR) / "workflow/lib/airr_exports.py"),
    output:
        airr = results_dir / "scopes/{scope}/airr/{library_id}_rearrangements.tsv.gz",
    benchmark:
        log_dir / "benchmarks/annotate_scope_airr/{scope}/{library_id}.tsv"
    log:
        log_dir / "annotate_scope_airr/{scope}/{library_id}.log"
    script:
        "../scripts/annotate_scope_airr.py"


# Snakemake removes declared outputs before a job starts. Declare only the
# completion receipt so the prior release folder survives a failed replacement.
# This rule alone owns that folder and publishes it as one recoverable
# transaction. Removing the receipt requests the folder's rebuild.
# assemble_scope_release.py requires a complete allow-list and matching library
# processing identities. Its staged validation writes airr/validation.json from
# the actual staged AIRR files before release_bundle.py replaces the folder.
# Raw inputs, pre-QC records and engineering evidence stay outside release/.
rule assemble_scope_release:
    input:
        products = lambda wc: [str(results_dir / "scopes" / wc.scope / relative)
                              for relative in scope_release_source_files(wc.scope, bool(membership_csv))],
        airr = lambda wc: [str(results_dir / "scopes" / wc.scope / "airr" /
                              (library + "_rearrangements.tsv.gz"))
                          for library in sorted(SCOPE_MEMBERS[wc.scope])],
        reference = rules.obtain_imgt_reference.output.imgt_ref,
        dictionary = DATA_DICTIONARY,
        identities = lambda wc: [str(results_dir / "libraries" / library / "export" /
                                     (library + "_processing.json"))
                                for library in sorted(SCOPE_MEMBERS[wc.scope])],
        alignment_module = str(Path(WORKFLOW_DIR) / "workflow/lib/alignment_export.py"),
        airr_module = str(Path(WORKFLOW_DIR) / "workflow/lib/airr_exports.py"),
        airr_schema = str(Path(WORKFLOW_DIR) / "workflow/schemas/airr-2.0.0.yaml"),
        readme = RELEASE_README,
        module = RELEASE_BUNDLE,
        methods = str(Path(WORKFLOW_DIR) / "workflow/lib/release_methods.py"),
        metadata_module = str(Path(WORKFLOW_DIR) / "workflow/lib/repertoire_metadata.py"),
    output:
        receipt = results_dir / "scopes/{scope}/release/.complete.json",
    params:
        scope_dir = lambda wc: str(results_dir / "scopes" / wc.scope),
        bundle = lambda wc: str(results_dir / "scopes" / wc.scope / "release"),
        libraries = lambda wc: sorted(SCOPE_MEMBERS[wc.scope]),
        membership = bool(membership_csv),
    conda:
        "../../envs/release.yaml"
    benchmark:
        log_dir / "benchmarks/assemble_scope_release/{scope}.tsv"
    log:
        log_dir / "assemble_scope_release/{scope}.log"
    script:
        "../scripts/assemble_scope_release.py"
