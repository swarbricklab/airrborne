# Run TRUST4 on each library's read pairs and retain its native sequence evidence.
# The run_trust4 rule below owns the tool command and its fixed arguments;
# params.trust4 in config/airrborne/config.yaml supplies the configurable assembly settings.
# parse_trust4_airr.R removes imputed CDR3 rows using flag_imputed_cdr3_records()
# in workflow/lib/trust4_parsing.R before IgBLAST annotation and common record QC.
# library.smk runs those later stages; assembly_backbone.R prepares grouping.
# Each library runs independently. Read support describes recovered reads, not cells.

def find_all_fastq_files(wildcards, samplesheet):
    """Return sorted, matched R1 and R2 FASTQ paths for one samplesheet library.

    Search only its directory beneath deps.data. Match names by prefix, suffix
    after R1/R2, or after 1/2 in SRA names, and .fastq.gz or .fq.gz extension.
    Missing directories, no pairs and incomplete pairs fail before either
    assembly or depth counting starts.
    TRUST4's launcher constructs an internal shell command, so paths are
    deliberately restricted to a small shell-safe character set even though
    Snakemake quotes them at the outer command boundary.
    """
    data_dir = config["deps"]["data"]

    row = samplesheet.loc[samplesheet['library_id'] == wildcards.library_id]

    if row.empty:
        raise ValueError("Could not find entry for library_id {} in samplesheet".format(wildcards.library_id))

    sample_dir = Path(data_dir) / row['path'].values[0]

    # The joined directory has to stay inside deps.data. The samplesheet schema
    # already rejects a dot segment, so this is the same rule applied to the
    # value that reaches the filesystem. It is lexical on purpose: the FASTQ
    # files here are symlinks into a content-addressed cache outside the data
    # directory, so resolving them would reject every library.
    relative_dir = os.path.relpath(os.path.normpath(sample_dir),
                                   os.path.normpath(data_dir))
    if relative_dir.split(os.sep)[0] == os.pardir:
        raise ValueError(
            "FASTQ directory for library_id {} leaves the configured data directory "
            "{}: {}".format(wildcards.library_id, data_dir, sample_dir))

    safe_path = re.compile(r"^[A-Za-z0-9_./:+-]+$")
    if not safe_path.fullmatch(str(sample_dir)):
        raise ValueError(
            "FASTQ directory for library_id {} contains whitespace or shell-special "
            "characters unsupported by run-trust4: {}".format(
                wildcards.library_id, sample_dir))
    if not sample_dir.is_dir():
        raise ValueError("Could not find FASTQ files for library_id {} in {}".format(wildcards.library_id, sample_dir))

    # Illumina names carry _R1 or _R2 and an optional suffix such as _001. SRA
    # and ENA downloads end in _1 or _2 directly before the extension. The
    # Illumina form is tried first, so x_R1_1.fastq.gz keeps its _1 suffix. Both
    # forms share one key, so x_R1.fastq.gz beside x_1.fastq.gz is a duplicate.
    fastq_names = (
        re.compile(r"^(?P<prefix>.+)_R(?P<mate>[12])(?P<suffix>(?:_[^/]*)?)"
                   r"(?P<extension>\.(?:fastq|fq)\.gz)$"),
        re.compile(r"^(?P<prefix>.+)_(?P<mate>[12])(?P<suffix>)"
                   r"(?P<extension>\.(?:fastq|fq)\.gz)$"),
    )
    pairs = {}
    for path in sorted(sample_dir.iterdir()):
        if not path.is_file():
            continue
        match = next((found for found in (name.match(path.name) for name in fastq_names)
                      if found is not None), None)
        if match is None:
            continue
        path_string = str(path)
        if not safe_path.fullmatch(path_string):
            raise ValueError(
                "FASTQ path contains whitespace or shell-special characters "
                "unsupported by run-trust4: {}".format(path_string))
        key = (match.group("prefix"), match.group("suffix"), match.group("extension"))
        mate = match.group("mate")
        if mate in pairs.setdefault(key, {}):
            raise ValueError("Duplicate R{} FASTQ member for {}: {}".format(
                mate, key, path_string))
        pairs[key][mate] = path_string

    if not pairs:
        raise ValueError("Could not find paired FASTQ files for library_id {} in {}".format(
            wildcards.library_id, sample_dir))

    incomplete = ["{} (found R{})".format(key, ",R".join(sorted(mates)))
                  for key, mates in pairs.items() if set(mates) != {"1", "2"}]
    if incomplete:
        raise ValueError(
            "Mismatched R1/R2 FASTQ pairs for library_id {} in {}: {}".format(
                wildcards.library_id, sample_dir, ", ".join(incomplete)))

    ordered = sorted(pairs)
    return {
        "r1": [pairs[key]["1"] for key in ordered],
        "r2": [pairs[key]["2"] for key in ordered],
    }

# Obtain the pinned TRUST4 image using the same helper as obtain_analysis_image.
# This local job pulls a SIF for Apptainer or Singularity, or uses Docker's daemon.
# obtain_image.py records the selected runtime and image identity in the receipt
# required by run_trust4. analysis_only omits this file and image acquisition.
localrules: obtain_assembly_image
rule obtain_assembly_image:
    input:
        script = str(SCRIPTS_DIR / "obtain_image.py"),
        module = str(Path(WORKFLOW_DIR) / "workflow/lib/provenance.py"),
    output:
        receipt = ASSEMBLY_IMAGE_RECEIPT,
    params:
        reference = ASSEMBLY_IMAGE_REFERENCE,
        sif = str(ASSEMBLY_SIF),
        runtime_option = lambda wc: ("--runtime " + shlex.quote(os.environ["AIRRBORNE_CONTAINER_RUNTIME"])
                                     if os.environ.get("AIRRBORNE_CONTAINER_RUNTIME") else ""),
    log:
        log_dir / "obtain_assembly_image.log"
    shell:
        "python {input.script:q} --reference {params.reference:q} --sif {params.sif:q} "
        "--receipt {output.receipt:q} {params.runtime_option} > {log:q} 2>&1"

# Run TRUST4 on one library's paired-end FASTQ files. Multi-lane, top-up and rerun
# inputs are passed as repeated -1/-2 flags; fastq-extractor combines them
# into a single _toassemble_{1,2}.fq for downstream assembly.
#
# Four retained outputs are declared: the three assembly files the pipeline
# consumes and the producer record written after the run. The run uses
# --clean 0 to retain consensus votes, so annot.fa, report.tsv and the AIRR
# alignments remain in the same directory as evidence for how each assembly
# was built.
#
# TRUST4's own cleanup cannot be used to discard the rest: --clean 1 would also
# delete final.out, whose per-position votes are the reason cleanup is 0. So the
# four bulky intermediates nothing downstream reads (the extracted candidate
# reads, the assembled reads and the raw consensus) are declared as temp()
# outputs instead and Snakemake removes them once the job finishes. All of them
# are reproducible from the source FASTQs by assembling again. Pass --notemp to
# keep them when debugging an assembly.
#
# params.trust4.k_mer_size supplies -k and contig_min_cov supplies --contigMinCov.
# The command fixes --clean 0 and gives TRUST4 the rule's allocated threads.
# run-trust4 options/defaults match the pinned TRUST4 1.1.10 image:
#   https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/README.md
# Its simplerep stage applies a within-assembly alternative-CDR3 threshold of
# 5% for TCR and 0 for BCR before run-trust4 converts the report to AIRR:
#   https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/trust-simplerep.pl#L193-L222
#   https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/trust-simplerep.pl#L346-L437
#   https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/run-trust4#L552-L565
# Method: Song et al., TRUST4, Nat Methods 2021 (PMC9328942; doi:10.1038/s41592-021-01142-2).
# Paired FASTQ input uses the curated IMGT reference for both -f and --ref.

# capture_assembly_context.py adds the image identity from the assembly image
# receipt, so the producer record names the image the job actually ran.
def _assembly_producer_context(wildcards, input):
    return json.dumps({
        "library_id": wildcards.library_id,
        # Record the flags the command actually used, including the literal
        # --clean 0 that keeps final.out.
        "settings": {**config["params"]["trust4"], "clean": 0},
        "read_pairs": [{"r1": first, "r2": second} for first, second in zip(input.r1, input.r2)],
    }, sort_keys=True)

rule run_trust4:
    input:
        container_exec = CONTAINER_EXEC,
        producer_script = str(SCRIPTS_DIR / "record_assembly_provenance.pl"),
        source_context_script = str(SCRIPTS_DIR / "capture_assembly_context.py"),
        provenance_module = str(Path(WORKFLOW_DIR) / "workflow/lib/provenance.py"),
        assembly_image = ASSEMBLY_IMAGE_RECEIPT,
        imgt_ref = rules.obtain_imgt_reference.output.imgt_ref,
        r1 = lambda wildcards: find_all_fastq_files(wildcards, samplesheet)["r1"],
        r2 = lambda wildcards: find_all_fastq_files(wildcards, samplesheet)["r2"]
    output:
        airr = results_dir / "libraries/{library_id}/trust4/{library_id}_airr.tsv",
        # Column 10 repeats the contig CDR3 annotation score on each row.
        # flag_imputed_cdr3_records() matches sequence_id to the CDR3 key and
        # flags the imputed score 0.01; parse_trust4_airr.R removes those rows.
        cdr3 = results_dir / "libraries/{library_id}/trust4/{library_id}_cdr3.out",
        # The assembly consensus with its per-position A/C/G/T read votes, which
        # are the only record of how much read evidence stands behind each
        # consensus base. build_igh_consensus_votes reads them for the IGH records.
        final_out = results_dir / "libraries/{library_id}/trust4/{library_id}_final.out",
        producer = results_dir / "libraries/{library_id}/trust4/{library_id}_assembly_provenance.json",
        # Assembly scratch, declared only so Snakemake can delete it. The candidate
        # reads fastq-extractor pulled out of the source FASTQs, the reads that went
        # into each contig, and the consensus before annotation and filtering.
        toassemble_r1 = temp(
            results_dir / "libraries/{library_id}/trust4/{library_id}_toassemble_1.fq"),
        toassemble_r2 = temp(
            results_dir / "libraries/{library_id}/trust4/{library_id}_toassemble_2.fq"),
        assembled_reads = temp(
            results_dir / "libraries/{library_id}/trust4/{library_id}_assembled_reads.fa"),
        raw_out = temp(
            results_dir / "libraries/{library_id}/trust4/{library_id}_raw.out")
    threads:
        8
    params:
        k_mer = config["params"]["trust4"]["k_mer_size"],
        contig_min_cov = config["params"]["trust4"]["contig_min_cov"],
        output_dir = lambda wildcards: str(
            results_dir / "libraries" / wildcards.library_id / "trust4"),
        output_prefix = lambda wildcards: wildcards.library_id,
        r1_args = lambda wildcards, input: " ".join("-1 " + shlex.quote(p) for p in input.r1),
        r2_args = lambda wildcards, input: " ".join("-2 " + shlex.quote(p) for p in input.r2),
        producer_context = _assembly_producer_context,
        workflow_root = str(Path(WORKFLOW_DIR)),
        assembly_sif = str(ASSEMBLY_SIF),
        docker_image = ASSEMBLY_IMAGE_REFERENCE,
    log:
        log_dir / "run_trust4/libraries/{library_id}.log"
    shell:
        """
        set -euo pipefail
        mkdir -p {params.output_dir:q}

        {input.container_exec:q} {params.assembly_sif:q} {params.docker_image:q} \
          run-trust4 \
          -f {input.imgt_ref:q} \
          --ref {input.imgt_ref:q} \
          {params.r1_args} {params.r2_args} \
          -t {threads} \
          -k {params.k_mer} \
          --contigMinCov {params.contig_min_cov} \
          -o {params.output_prefix:q} \
          --od {params.output_dir:q} \
          --clean 0 \
          > {log:q} 2>&1

        CONTEXT=$(python {input.source_context_script:q} \
          --repo {params.workflow_root:q} --context {params.producer_context:q} \
          --image-receipt {input.assembly_image:q})
        {input.container_exec:q} {params.assembly_sif:q} {params.docker_image:q} \
          perl {input.producer_script:q} --context "$CONTEXT" \
          --reference {input.imgt_ref:q} --output {output.producer:q} \
          --airr {output.airr:q} --cdr3 {output.cdr3:q} --final {output.final_out:q} \
          >> {log:q} 2>&1
        """

# Count raw reads and bases in the same R1/R2 FASTQs supplied to TRUST4. This
# host-side job runs independently of assembly. count_seq_depth.sh validates
# four-line FASTQ records and equal mate counts, then writes one library row with
# file-pair, read-pair, read and base counts. Empty paired files have zero depth.
rule count_sequencing_depth:
    input:
        r1 = lambda wildcards: find_all_fastq_files(wildcards, samplesheet)["r1"],
        r2 = lambda wildcards: find_all_fastq_files(wildcards, samplesheet)["r2"],
        helper = str(SCRIPTS_DIR / "count_seq_depth.sh"),
    output:
        seq_depth = results_dir / "libraries/{library_id}/metrics/{library_id}_seq_depth.tsv"
    params:
        # Quote each path before passing the list to the helper's --r1/--r2 loop.
        # Discovery also restricts paths for TRUST4's internal shell command.
        r1_args = lambda wildcards, input: " ".join(shlex.quote(p) for p in input.r1),
        r2_args = lambda wildcards, input: " ".join(shlex.quote(p) for p in input.r2),
    log:
        log_dir / "count_sequencing_depth/libraries/{library_id}.log"
    shell:
        """
        bash {input.helper:q} {output.seq_depth:q} {wildcards.library_id:q} \
          --r1 {params.r1_args} --r2 {params.r2_args} \
          > {log:q} 2>&1
        """
