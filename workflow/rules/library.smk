# Annotate and summarise each library independently of scope membership.
# AIRR rows are TRUST4-reported sequence variants with recovered read support.
# They do not identify cells. Common QC precedes separate grouping, isotype and
# SHM requirements; scope-dependent IGH clone inference lives in igh.smk.

# Reannotate native IG[HKL] and TR[ABDG] records with IgBLAST and MakeDb.
# run_igblast_all_loci.sh owns database routing, tool arguments and IG germline
# reconstruction. build_reannotated_airr.R combines the annotations with native
# IDs and support. Its candidate-list completeness and gene_calls.R resolution
# determine whether V/J calls name one accepted analytical gene entity.
rule reannotate_igblast:
    input:
        airr = results_dir / "libraries/{library_id}/trust4/{library_id}_airr.tsv",
        # TRUST4 CDR3 report; repeats the contig imputation score on each variant.
        cdr3 = results_dir / "libraries/{library_id}/trust4/{library_id}_cdr3.out",
        # Original allele metadata and complete-reference identity classes.
        # Later stages read the analytical entities attached by reannotation.
        reference_metadata = rules.obtain_imgt_reference.output.metadata,
        # IgBLAST auxiliary data written beside the reference: the image's NCBI
        # file plus supported supplementary J-allele entries from germline_reference.py.
        auxiliary = rules.obtain_imgt_reference.output.auxiliary,
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        igblast_script = str(SCRIPTS_DIR / "run_igblast_all_loci.sh"),
        parse_script = str(SCRIPTS_DIR / "parse_trust4_airr.R"),
        build_script = str(SCRIPTS_DIR / "build_reannotated_airr.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        airr_io = AIRR_IO,
        gene_calls = GENE_CALLS,
        reannotated_airr = REANNOTATED_AIRR,
        qc_audit = QC_AUDIT,
        trust4_parsing = TRUST4_PARSING,
    output:
        # Not temp(): this is the annotated repertoire before the common record-QC
        # gate runs, so it is the file that lets a post-annotation QC question be
        # re-asked, such as a different N threshold, a non-productive tier or the cost
        # of the D-call locus rule, without re-running IgBLAST over the whole cohort.
        # Earlier decisions are already applied here: the parser removed imputed CDR3
        # calls, and MakeDb dropped what it could not annotate. Identical sequences
        # keep their own records until grouping.
        igblast = results_dir / "libraries/{library_id}/airr/{library_id}_igblast.tsv",
        # What MakeDb rejected, and its own per-record reason for doing so. The
        # QC audit counts these; only these files say which records and why.
        makedb_fail = results_dir / "libraries/{library_id}/airr/{library_id}_makedb_fail.tsv",
        makedb_log = results_dir / "libraries/{library_id}/airr/{library_id}_makedb.log.gz",
        igblast_hits = results_dir / "libraries/{library_id}/airr/{library_id}_igblast_hits.fmt7.gz",
        makedb_pass = results_dir / "libraries/{library_id}/airr/{library_id}_makedb_pass.tsv.gz",
        qc_audit = results_dir / "libraries/{library_id}/qc/{library_id}_qc_audit_reannotate.tsv"
    threads:
        8
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        sample_label = lambda wc: wc.library_id,
        num_alignments_v = config['params']['igblast']['num_alignments_v'],
        num_alignments_j = config['params']['igblast']['num_alignments_j'],
        num_alignments_d = config['params']['igblast']['num_alignments_d'],
        num_alignments_c = config['params']['igblast']['num_alignments_c'],
        igblast_db = config['params']['igblast']['igblast_db'],
        germline_dir = config['params']['igblast']['germline_dir'],
        outdir = lambda wc: str(results_dir / "libraries" / wc.library_id / "airr")
    log:
        log_dir / "reannotate_igblast/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir:q}
        # Clear any scratch directory an interrupted earlier attempt left behind.
        rm -rf {params.outdir:q}/igblast.*
        IGB=$(mktemp -d {params.outdir:q}/igblast.XXXXXX)
        trap 'rm -rf "$IGB"' EXIT

        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            bash {input.igblast_script:q} \
            {input.airr:q} {output.igblast:q} "$IGB" \
            {params.sample_label:q} {threads} \
            {params.igblast_db:q} {params.germline_dir:q} \
            {input.cdr3:q} {output.qc_audit:q} \
            {output.makedb_fail:q} {output.makedb_log:q} \
            {params.num_alignments_v} {params.num_alignments_j} \
            {params.num_alignments_d} {params.num_alignments_c} \
            {output.makedb_pass:q} {input.auxiliary:q} {input.reference_metadata:q} {output.igblast_hits:q} \
            > {log:q} 2>&1

        if [ ! -s {output.igblast:q} ]; then
            echo "Expected all-locus IgBLAST output was not created: {output.igblast}" >> {log:q}
            exit 1
        fi

        # Validate the assembled AIRR Rearrangement table before common QC.
        # Publication separately validates the staged exports and their links.
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            airr-tools validate rearrangement -a {output.igblast:q} \
            >> {log:q} 2>&1
        """

# Keep all-locus records passing airr_qc_flags() and airr_qc_pass() in
# record_qc.R. params.airr_qc.max_alignment_n_fraction supplies max_n_fraction.
# filter_airr.R writes the retained rows and criterion counts; failures can
# overlap. A valid empty input retains the output schema. Grouping, isotype and
# SHM apply their further requirements in independent downstream jobs.
rule filter_airr:
    input:
        igblast = results_dir / "libraries/{library_id}/airr/{library_id}_igblast.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        filter_script = str(SCRIPTS_DIR / "filter_airr.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        airr_io = AIRR_IO,
        gene_calls = GENE_CALLS,
        qc_audit = QC_AUDIT,
        record_qc = RECORD_QC,
    output:
        qc_pass = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        qc_audit = results_dir / "libraries/{library_id}/qc/{library_id}_qc_audit_qc.tsv"
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        max_n_fraction = config['params']['airr_qc']['max_alignment_n_fraction'],
    log:
        log_dir / "filter_airr/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.filter_script:q} \
            --input_file={input.igblast:q} \
            --output_file={output.qc_pass:q} \
            --library_id={wildcards.library_id:q} \
            --max_n_fraction={params.max_n_fraction:q} \
            --qc_audit_file={output.qc_audit:q} \
            > {log:q} 2>&1
        """

# Measure one row per common-QC IGH record with
# build_igh_v_region_feature_table() in igh_features.R. min_position_votes marks
# supported IMGT V positions from the matching assembly votes. Missing regional
# measurements remain missing. Qualified C calls supply isotype context without
# changing clone eligibility. The retained features precede the region-specific
# coverage and evidence-source choices in calculate_donor_igh_summaries.
rule calculate_igh_features:
    input:
        reannotated = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        isotype_calls = results_dir / "libraries/{library_id}/metrics/{library_id}_isotype_calls.tsv",
        votes = results_dir / "libraries/{library_id}/airr/{library_id}_igh_consensus_votes.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        feature_script = str(SCRIPTS_DIR / "calculate_igh_features.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        consensus_votes = CONSENSUS_VOTES,
        gene_calls = GENE_CALLS,
        igh_features = IGH_FEATURES,
        isotype_calls_module = ISOTYPE_CALLS,
        shm_regions = SHM_REGIONS,
        qc_audit_module = QC_AUDIT,
    output:
        qc_audit = results_dir / "libraries/{library_id}/qc/{library_id}_qc_audit_shm_measurement.tsv",
        features = results_dir / "libraries/{library_id}/metrics/{library_id}_igh_features.tsv"
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        min_position_votes = config['params']['shm']['min_position_votes'],
    log:
        log_dir / "calculate_igh_features/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.feature_script:q} \
            --input_file={input.reannotated:q} \
            --votes_file={input.votes:q} \
            --min_position_votes={params.min_position_votes} \
            --output_file={output.features:q} \
            --qc_audit_file={output.qc_audit:q} \
            --library_id={wildcards.library_id:q} \
            --isotype_calls_file={input.isotype_calls:q} \
            > {log:q} 2>&1
        """

# Map TRUST4 _final.out read tallies to each common-QC IGH record's IMGT V
# positions through build_igh_consensus_votes.R and trust4_consensus_votes.R.
# The assembly and its reannotation are both required to match positions.
# The retained table records positional evidence; it does not pool assemblies
# or decide which source qualifies for a region's SHM summary.
rule build_igh_consensus_votes:
    input:
        reannotated = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        final_out = results_dir / "libraries/{library_id}/trust4/{library_id}_final.out",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        votes_script = str(SCRIPTS_DIR / "build_igh_consensus_votes.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        consensus_votes = CONSENSUS_VOTES,
        gene_calls = GENE_CALLS,
        assembly_backbone = ASSEMBLY_BACKBONE,
    output:
        votes = results_dir / "libraries/{library_id}/airr/{library_id}_igh_consensus_votes.tsv"
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
    log:
        log_dir / "build_igh_consensus_votes/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.votes_script:q} \
            --input_file={input.reannotated:q} \
            --final_out_file={input.final_out:q} \
            --output_file={output.votes:q} \
            --library_id={wildcards.library_id:q} \
            > {log:q} 2>&1
        """

# Build non-IGH groups with build_rearrangement_group_table() in
# rearrangement_groups.R: exact A/C/G/T nucleotide junction and resolved V/J
# entities define each library-local group. Constant-call compatibility also
# controls eligibility. IGK and IGL choose the dominant common-QC backbone
# variant before eligibility; an ineligible winner excludes its minor siblings.
# TCR loci retain eligible variants. prepare_assembly_mapping() then merges
# identical retained full DNA within library and locus, requiring agreeing
# annotations. The grouping and membership tables preserve representative and
# source links and add each eligible source's support once. Empty loci retain
# typed tables. These groups describe sequence observations, not cells.
rule build_clonotypes:
    input:
        reannotated = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "build_clonotypes.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        gene_calls = GENE_CALLS,
        assembly_backbone = ASSEMBLY_BACKBONE,
        rearrangement_keys = REARRANGEMENT_KEYS,
        rearrangement_groups = REARRANGEMENT_GROUPS,
        qc_audit_module = QC_AUDIT,
        pooled_records = POOLED_RECORDS,
    output:
        grouping = results_dir / "libraries/{library_id}/metrics/{library_id}_nonigh_grouping.tsv",
        membership = results_dir / "libraries/{library_id}/metrics/{library_id}_clonotype_members.tsv",
        qc_audit = results_dir / "libraries/{library_id}/qc/{library_id}_qc_audit_nonigh_grouping.tsv",
        groups_igk = results_dir / "libraries/{library_id}/clonotypes/{library_id}_IGK_clonotypes.tsv",
        groups_igl = results_dir / "libraries/{library_id}/clonotypes/{library_id}_IGL_clonotypes.tsv",
        groups_tra = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRA_clonotypes.tsv",
        groups_trb = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRB_clonotypes.tsv",
        groups_trd = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRD_clonotypes.tsv",
        groups_trg = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRG_clonotypes.tsv",
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
    log:
        log_dir / "build_clonotypes/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --input_file={input.reannotated:q} \
            --library_id={wildcards.library_id:q} \
            --membership_file={output.membership:q} \
            --grouping_file={output.grouping:q} \
            --qc_audit_file={output.qc_audit:q} \
            --igk_groups_file={output.groups_igk:q} \
            --igl_groups_file={output.groups_igl:q} \
            --tra_groups_file={output.groups_tra:q} \
            --trb_groups_file={output.groups_trb:q} \
            --trd_groups_file={output.groups_trd:q} \
            --trg_groups_file={output.groups_trg:q} \
            > {log:q} 2>&1
        """

# Calculate per-library non-IGH diversity and dominance, rarefaction curves,
# support vectors and V/J usage from the rearrangement-group tables.
# diversity.R applies min_consensus_count to entropy, Hill diversity, dominance
# and target-depth estimates. Observed richness, support and curves remain
# available. Every non-IGH locus receives a diversity row and a support vector;
# rarefaction.targets supplies target depths and curves. gene_usage.R counts V/J
# usage and reads gene names, gene groups and identical alleles from the
# reference metadata.
rule calculate_library_clonotype_diversity:
    input:
        igk = results_dir / "libraries/{library_id}/clonotypes/{library_id}_IGK_clonotypes.tsv",
        igl = results_dir / "libraries/{library_id}/clonotypes/{library_id}_IGL_clonotypes.tsv",
        tra = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRA_clonotypes.tsv",
        trb = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRB_clonotypes.tsv",
        trd = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRD_clonotypes.tsv",
        trg = results_dir / "libraries/{library_id}/clonotypes/{library_id}_TRG_clonotypes.tsv",
        reference_metadata = rules.obtain_imgt_reference.output.metadata,
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "calculate_library_clonotype_diversity.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        diversity = DIVERSITY,
        gene_calls = GENE_CALLS,
        gene_usage = GENE_USAGE,
    output:
        diversity = results_dir / "libraries/{library_id}/metrics/{library_id}_diversity.tsv",
        curves = results_dir / "libraries/{library_id}/metrics/{library_id}_rarefaction_curves.tsv",
        group_abundance = results_dir / "libraries/{library_id}/metrics/{library_id}_group_abundance.tsv",
        gene_usage = results_dir / "libraries/{library_id}/metrics/{library_id}_clonotype_gene_usage.tsv",
        vj_pair_usage = results_dir / "libraries/{library_id}/metrics/{library_id}_clonotype_vj_pair_usage.tsv",
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        rarefaction_targets = named_values({k: v for k, v in config["params"]["rarefaction"]["targets"].items() if k != "IGH"}),
        minimum_consensus_count = named_values({k: v for k, v in config["params"]["diversity"]["min_consensus_count"].items() if k != "IGH"}),
    log:
        log_dir / "calculate_library_clonotype_diversity/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --igk_file={input.igk:q} \
            --igl_file={input.igl:q} \
            --tra_file={input.tra:q} \
            --trb_file={input.trb:q} \
            --trd_file={input.trd:q} \
            --trg_file={input.trg:q} \
            --library_id={wildcards.library_id:q} \
            --reference_metadata_file={input.reference_metadata:q} \
            --rarefaction_targets={params.rarefaction_targets:q} \
            --minimum_consensus_count={params.minimum_consensus_count:q} \
            --diversity_file={output.diversity:q} \
            --curves_file={output.curves:q} \
            --group_abundance_file={output.group_abundance:q} \
            --gene_usage_file={output.gene_usage:q} \
            --vj_pair_usage_file={output.vj_pair_usage:q} \
            > {log:q} 2>&1
        """

# Calculate IGH class and subclass composition from the common QC-pass table,
# weighted by read support. A record contributes to a category only when its
# qualified constant-region call resolves there; support that does not resolve
# stays outside that level's callable support, and the summary reports it as
# class-uncallable support. build_isotype_metrics() in isotype_metrics.R owns
# the denominators. Empty IGH input gives zero support
# counts and missing fractions, while retaining the complete category rows.
rule calculate_isotype:
    input:
        reannotated = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        isotype_calls = results_dir / "libraries/{library_id}/metrics/{library_id}_isotype_calls.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "calculate_isotype.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        gene_calls = GENE_CALLS,
        isotype_calls_module = ISOTYPE_CALLS,
        isotype_resolution = ISOTYPE_RESOLUTION,
        isotype_metrics = ISOTYPE_METRICS,
        qc_audit_module = QC_AUDIT,
    output:
        qc_audit = results_dir / "libraries/{library_id}/qc/{library_id}_qc_audit_isotype.tsv",
        isotype = results_dir / "libraries/{library_id}/metrics/{library_id}_isotype_usage.tsv",
        summary = results_dir / "libraries/{library_id}/metrics/{library_id}_isotype_summary.tsv",
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
    log:
        log_dir / "calculate_isotype/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --input_file={input.reannotated:q} \
            --isotype_calls_file={input.isotype_calls:q} \
            --library_id={wildcards.library_id:q} \
            --qc_audit_file={output.qc_audit:q} \
            --isotype_file={output.isotype:q} \
            --summary_file={output.summary:q} \
            > {log:q} 2>&1
        """

# Summarise the support carried by TRUST4's reported full-CDR3 records. This is
# an output support quantity, not an assembly-recovery or unique-read count.
# repertoire_support.R writes one row per library before pipeline QC with the
# library's input depth, and divides locus support by input read pairs for yield
# per million. A header-only TRUST4 AIRR table gives zero support; absent input
# read pairs leave yield unavailable.
rule summarise_library_yield:
    input:
        depth = results_dir / "libraries/{library_id}/metrics/{library_id}_input_depth.tsv",
        airr = results_dir / "libraries/{library_id}/trust4/{library_id}_airr.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "summarise_library_yield.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        repertoire_support = REPERTOIRE_SUPPORT,
    output:
        library_yield = results_dir / "libraries/{library_id}/metrics/{library_id}_yield.tsv"
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
    log:
        log_dir / "summarise_library_yield/libraries/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --airr_file={input.airr:q} \
            --depth_file={input.depth:q} \
            --library_id={wildcards.library_id:q} \
            --output_file={output.library_yield:q} \
            > {log:q} 2>&1
        """

# Qualify IGH constant-region calls once for composition, mutation analyses,
# BASELINe selection and the released library AIRR export.
# build_qualified_isotype_calls() in isotype_calls.R uses min_c_positions for
# the informative C-alignment minimum. Failed calls are masked without deleting
# records. Original calls remain in the retained sidecar and post-QC AIRR table;
# the audit records overlapping qualification failures. resolve_isotype() in
# isotype_resolution.R names each qualified call's class and subclass for the
# AIRR export, as the isotype tables and clone members name them.
rule qualify_isotype_calls:
    input:
        reannotated = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "qualify_isotype_calls.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        gene_calls = GENE_CALLS,
        isotype_calls_module = ISOTYPE_CALLS,
        isotype_resolution = ISOTYPE_RESOLUTION,
        qc_audit_module = QC_AUDIT,
    output:
        qc_audit = results_dir / "libraries/{library_id}/qc/{library_id}_qc_audit_isotype_calls.tsv",
        isotype_calls = results_dir / "libraries/{library_id}/metrics/{library_id}_isotype_calls.tsv",
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        min_c_positions = config['params']['airr_qc']['min_informative_c_alignment_positions'],
    benchmark:
        log_dir / "benchmarks/qualify_isotype_calls/{library_id}.tsv"
    log:
        log_dir / "qualify_isotype_calls/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --input_file={input.reannotated:q} \
            --library_id={wildcards.library_id:q} \
            --min_c_positions={params.min_c_positions:q} \
            --qc_audit_file={output.qc_audit:q} \
            --output_file={output.isotype_calls:q} \
            > {log:q} 2>&1
        """

# Convert the native depth record count_sequencing_depth writes beside assembly
# to the analytical schema. This rule reads no FASTQ, so a --post-trust4 run
# reuses retained depth records and the assembly boundary stays reusable.
# Convert the retained sequencing-depth row with normalise_input_depth.py.
# It validates paired-read counts and names input_read_pairs for yield's
# denominator. It reads no FASTQ, so analysis-only runs can reuse that evidence.
rule normalise_input_depth:
    input:
        depth = results_dir / "libraries/{library_id}/metrics/{library_id}_seq_depth.tsv",
    output:
        depth = results_dir / "libraries/{library_id}/metrics/{library_id}_input_depth.tsv",
    benchmark:
        log_dir / "benchmarks/normalise_input_depth/{library_id}.tsv"
    log:
        log_dir / "normalise_input_depth/{library_id}.log"
    script:
        "../scripts/normalise_input_depth.py"

# Prepare IGH grouping beside the unchanged common-QC pool with
# prepare_igh_grouping() in igh_grouping.R. prepare_assembly_mapping() chooses
# the dominant backbone variant using original record support before testing
# eligibility, then merges identical retained full DNA with agreeing annotations.
# max_n supplies the junction N/? allowance. The retained mapping records every
# source, both representative links and support counted once for inference.
# Empty IGH input yields the typed mapping and audit without clone assignment.
rule prepare_igh_grouping:
    input:
        reannotated = results_dir / "libraries/{library_id}/airr/{library_id}_airr_qc_pass.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "prepare_igh_grouping.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        gene_calls = GENE_CALLS,
        assembly_backbone = ASSEMBLY_BACKBONE,
        igh_clones = IGH_CLONES,
        igh_grouping = IGH_GROUPING,
        qc_audit_module = QC_AUDIT,
    output:
        grouping = results_dir / "libraries/{library_id}/metrics/{library_id}_igh_grouping.tsv",
        qc_audit = results_dir / "libraries/{library_id}/qc/{library_id}_qc_audit_igh_grouping.tsv",
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        max_n = config["params"]["lineage"]["max_n"],
    benchmark:
        log_dir / "benchmarks/prepare_igh_grouping/{library_id}.tsv"
    log:
        log_dir / "prepare_igh_grouping/{library_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --input_file={input.reannotated:q} \
            --library_id={wildcards.library_id:q} \
            --max_n={params.max_n:q} \
            --grouping_file={output.grouping:q} \
            --qc_audit_file={output.qc_audit:q} \
            > {log:q} 2>&1
        """
