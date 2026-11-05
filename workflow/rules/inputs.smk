# Write scope membership and prepare the shared analysis image and references.
# Snakefile validates membership and builds scope_manifest_by_scope and
# SCOPE_ASSIGNMENT_ROWS; these rules write those decisions without reapplying
# membership policy. Image acquisition runs locally, while reference preparation
# and its independent audit use the analysis image through container_exec.sh.

# Write one row per library in the selected scope, including donor and tissue IDs.
# Only this scope's validated rows enter params.rows, so editing another scope's
# membership does not schedule this manifest or the jobs that read it.
localrules: build_scope_manifest
rule build_scope_manifest:
    output:
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv"
    params:
        rows = lambda wc: sorted(scope_manifest_by_scope[wc.scope]),
    run:
        # Written field by field rather than through csv.writer so this file uses
        # the same dialect as every table R writes through write_pipeline_table():
        # tab-separated, unquoted, newline-terminated. csv.writer would end each
        # record with a carriage return and line feed, which leaves a stray
        # carriage return on the last column for anything reading the file by
        # bytes rather than through a CSV parser.
        # The identifiers follow IDENTITY_ORDER in scope_tables.py, as every
        # collected table does.
        rows = params.rows
        lines = ["\t".join(["scope_id", "library_id", "tissue_id", "donor_id"])]
        lines += ["\t".join([wildcards.scope, library_id, tissue_id, donor_id])
                  for library_id, donor_id, tissue_id in sorted(rows)]
        rendered = "\n".join(lines) + "\n"
        # The validated rows are tracked as this scope's parameters. Changes to
        # another scope's input rows do not schedule this manifest or its consumers.
        manifest = Path(output.manifest)
        if manifest.exists() and manifest.read_text() == rendered:
            logger.info("scope manifest unchanged, leaving as is: %s", manifest)
        else:
            manifest.parent.mkdir(parents=True, exist_ok=True)
            manifest.write_text(rendered)


# Write the run-wide membership audit beside the reserved cohort manifest.
# SCOPE_ASSIGNMENT_ROWS contains each membership row and a row for each
# samplesheet library no project claims. Its analysed flag records project
# membership after usage filtering; it does not report downstream record QC.
# A samplesheet library remains in the cohort even when a project rejects it.
localrules: build_scope_assignment
rule build_scope_assignment:
    output:
        # Composed rather than formatted: Snakemake re-tokenises this file, and an
        # f-string comes back out with spaces around what it interpolated.
        assignment = results_dir / "scopes" / RESERVED_SCOPE / "manifest/scope_assignment.tsv"
    params:
        rows = sorted(SCOPE_ASSIGNMENT_ROWS),
    run:
        # Written field by field for the same reason as the scope manifest above:
        # tab-separated, unquoted, newline-terminated, like every table R writes.
        header = ["library_id", "project", "usage", "in_samplesheet", "analysed"]
        lines = ["\t".join(header)]
        lines += ["\t".join([library_id, project, usage,
                             str(in_samplesheet).upper(), str(analysed).upper()])
                  for library_id, project, usage, in_samplesheet, analysed
                  in params.rows]
        rendered = "\n".join(lines) + "\n"
        # An unchanged file is left exactly as it is, so that editing one project's
        # rows does not advance this file's timestamp and rebuild the cohort bundle
        # for nothing.
        assignment = Path(output.assignment)
        if assignment.exists() and assignment.read_text() == rendered:
            logger.info("scope assignment unchanged, leaving as is: %s", assignment)
        else:
            assignment.parent.mkdir(parents=True, exist_ok=True)
            assignment.write_text(rendered)


# Obtain the pinned analysis image and record its runtime and image identity.
# obtain_image.py calls provenance.obtain_image(), which follows the runtime
# selection in container_exec.sh unless AIRRBORNE_CONTAINER_RUNTIME overrides it.
# This local job pulls a SIF for Apptainer or Singularity, or uses Docker's daemon.
# Analysis jobs declare the receipt as their image input. setup.sh writes the
# same receipt beside the SIF before the first launch, so this job runs only
# when the image or its identification code changes. The helper preserves an
# identical receipt; Snakemake can still update output timestamps after a job.
localrules: obtain_analysis_image
rule obtain_analysis_image:
    input:
        script = str(SCRIPTS_DIR / "obtain_image.py"),
        module = str(Path(WORKFLOW_DIR) / "workflow/lib/provenance.py"),
    output:
        receipt = ANALYSIS_IMAGE_RECEIPT,
    params:
        reference = ANALYSIS_IMAGE_REFERENCE,
        sif = str(ANALYSIS_SIF),
        runtime_option = lambda wc: ("--runtime " + shlex.quote(os.environ["AIRRBORNE_CONTAINER_RUNTIME"])
                                     if os.environ.get("AIRRBORNE_CONTAINER_RUNTIME") else ""),
    log:
        log_dir / "obtain_analysis_image.log"
    shell:
        "python {input.script:q} --reference {params.reference:q} --sif {params.sif:q} "
        "--receipt {output.receipt:q} {params.runtime_option} > {log:q} 2>&1"


# Build the shared reference FASTA, allele metadata and IgBLAST auxiliary data.
# obtain_imgt_reference.py calls build_reference() and build_auxiliary_data() in
# workflow/lib/germline_reference.py using the analysis image's IMGT files.
# build_reference() keeps the first record per allele and writes bare allele
# headers. supplement_auxiliary_data() first tries an identical-sequence allele,
# then a same-gene allele of equal length with matching conserved-codon evidence.
# An allele with no qualifying source stays without an auxiliary entry.
# The FASTA serves TRUST4. Annotation reads the image's databases and germline
# files alongside this rule's auxiliary data and gene-call metadata.
rule obtain_imgt_reference:
    input:
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "obtain_imgt_reference.py"),
        module = GERMLINE_REFERENCE,
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
    output:
        imgt_ref = resources_dir / config['refs']['imgt_reference_name'],
        auxiliary = AUXILIARY_DATA_FILE,
        metadata = resources_dir / (config['refs']['imgt_reference_name'] + ".metadata.tsv"),
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        germline_dir = GERMLINE_DIR,
        constant_dir = CONSTANT_DIR,
        shipped_auxiliary = SHIPPED_AUXILIARY_DATA,
    log:
        log_dir / "obtain_imgt_reference.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            python3 {input.script:q} \
            --germline-dir {params.germline_dir:q} --constant-dir {params.constant_dir:q} \
            --auxiliary-data {params.shipped_auxiliary:q} \
            --output {output.imgt_ref:q} --auxiliary-output {output.auxiliary:q} \
            --metadata-output {output.metadata:q} \
            > {log:q} 2>&1
        """


# Audit the curated FASTA, allele metadata and auxiliary data against fresh builds.
# audit_reference() in workflow/lib/germline_reference.py reads each IG/TR V/D/J
# IgBLAST database with blastdbcmd. Each database must match its source allele set
# and sequences; every database allele must occur unchanged in the combined FASTA.
# The receipt lists supplementary auxiliary entries and J alleles still without
# one. Scope provenance requires this receipt, so disagreement stops publication.
# Annotation does not declare it: Snakemake touches outputs after a job, so an
# audit code change would otherwise rerun IgBLAST for every library.
rule audit_imgt_reference:
    input:
        container_exec = CONTAINER_EXEC,
        reference = rules.obtain_imgt_reference.output.imgt_ref,
        auxiliary = rules.obtain_imgt_reference.output.auxiliary,
        metadata = rules.obtain_imgt_reference.output.metadata,
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        script = str(SCRIPTS_DIR / "audit_imgt_reference.py"),
        module = GERMLINE_REFERENCE,
    output:
        receipt = IMGT_REFERENCE_RECEIPT,
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        germline_dir = GERMLINE_DIR,
        constant_dir = CONSTANT_DIR,
        snapshot = IMGT_SNAPSHOT,
        igblast_db = config['params']['igblast']['igblast_db'],
    log:
        log_dir / "audit_imgt_reference.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            python3 {input.script:q} \
            --reference {input.reference:q} --auxiliary {input.auxiliary:q} \
            --metadata {input.metadata:q} \
            --germline-dir {params.germline_dir:q} --constant-dir {params.constant_dir:q} \
            --snapshot {params.snapshot:q} --igblast-db {params.igblast_db:q} \
            --image-receipt {input.analysis_image:q} --receipt {output.receipt:q} \
            > {log:q} 2>&1
        """
