# Infer IGH clones separately within each donor and scope, then summarise them.
# pool_scope_donor_igh.R keeps unchanged common-QC rows beside prepared grouping.
# igh_grouping.R and assembly_backbone.R choose backbone and exact-sequence
# representatives before fitting or assignment. A backbone denotes a TRUST4
# contig with reported CDR3 alternatives, not a cell or confirmed rearrangement.
# Shared library features remain independent of these scope-dependent jobs.

def _scope_donor_reann_files(scope, donor):
    return [str(results_dir / "libraries/{}/airr/{}_airr_qc_pass.tsv".format(g, g))
            for g in sorted(g for t, g in scope_donor_to_libraries[(scope, donor)])]

def _scope_donor_igh_feature_files(scope, donor):
    return [str(results_dir / "libraries/{}/metrics/{}_igh_features.tsv".format(g, g))
            for g in sorted(g for t, g in scope_donor_to_libraries[(scope, donor)])]

def _scope_donor_isotype_call_files(scope, donor):
    return [str(results_dir / "libraries" / library / "metrics" / (library + "_isotype_calls.tsv"))
            for library in sorted(g for t, g in scope_donor_to_libraries[(scope, donor)])]

def _scope_donor_grouping_files(scope, donor):
    return [str(results_dir / "libraries" / library / "metrics" / (library + "_igh_grouping.tsv"))
            for library in sorted(g for t, g in scope_donor_to_libraries[(scope, donor)])]

# Pool one donor's in-scope IGH rows and grouping in sorted library order.
# pool_scope_donor_igh.R includes library IDs in sequence, backbone and
# representative IDs and validates the links. Both temporary tables remain
# available until their consumers finish, including selection and SHM summaries.
rule pool_scope_donor_igh:
    input:
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        reannotated = lambda wc: _scope_donor_reann_files(wc.scope, wc.donor_id),
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        grouping = lambda wc: _scope_donor_grouping_files(wc.scope, wc.donor_id),
        table_io = TABLE_IO,
        assembly_backbone = ASSEMBLY_BACKBONE,
        pool_script = str(SCRIPTS_DIR / "pool_scope_donor_igh.R"),
        runtime = RUNTIME,
        pooled_records = POOLED_RECORDS,
    output:
        pooled = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_pooled.tsv"),
        grouping = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_grouping.tsv"),
    threads: 1
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        reannotated_files = lambda wc, input: ",".join(input.reannotated),
        grouping_files = lambda wc, input: ",".join(input.grouping),
    log:
        log_dir / "pool_scope_donor_igh/scopes/{scope}/{donor_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.pooled:q})"
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.pool_script:q} \
            --out_file={output.pooled:q} \
            --grouping_files={params.grouping_files:q} \
            --grouping_file={output.grouping:q} \
            --manifest_file={input.manifest:q} \
            --scope_id={wildcards.scope:q} \
            --donor_id={wildcards.donor_id:q} \
            --reannotated_files={params.reannotated_files:q} \
            > {log:q} 2>&1
        """

# Assign prepared inference representatives with SCOPer hierarchicalClones in
# assign_scope_donor_clones.R. The saved scope cutoff limits each single-linkage
# joining distance; linked members can be farther apart. Donors never join.
# igh_clones.R builds membership with each eligible source's read support counted
# once. Clone IDs include donor and scope. An empty prepared pool writes empty
# membership and records plus the population audit.
rule assign_scope_donor_clones:
    input:
        grouping = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_grouping.tsv",
        pooled = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_pooled.tsv",
        # One scope-wide cutoff keeps clone assignment comparable across donors.
        scope_threshold = results_dir / "scopes/{scope}/clones/{scope}_clone_threshold.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        assign_script = str(SCRIPTS_DIR / "assign_scope_donor_clones.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        gene_calls = GENE_CALLS,
        assembly_backbone = ASSEMBLY_BACKBONE,
        igh_clones = IGH_CLONES,
        igh_grouping = IGH_GROUPING,
        qc_audit_module = QC_AUDIT,
        rearrangement_keys = REARRANGEMENT_KEYS,
        scope_clone_threshold = SCOPE_CLONE_THRESHOLD,
    output:
        qc_audit = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_qc_audit_clones.tsv",
        membership = results_dir / "scopes/{scope}/clones/{donor_id}_IGH_clone_membership.tsv",
        # Retain the contigs and exact-nucleotide keys for independent overlap
        # and summary reruns. This working shard remains outside the release.
        records = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_lineage_records.tsv",
    threads: 1
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        # Fitting and assignment consume the same prepared representatives and
        # exact V/J entities. A complete-sequence composite joins only itself.
        # Both stages use the same junction allowance.
        max_n = config['params']['lineage']['max_n'],
    log:
        log_dir / "assign_scope_donor_clones/scopes/{scope}/{donor_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.assign_script:q} \
            --input_file={input.pooled:q} \
            --grouping_file={input.grouping:q} \
            --threshold_file={input.scope_threshold:q} \
            --scope_id={wildcards.scope:q} \
            --donor_id={wildcards.donor_id:q} \
            --max_n={params.max_n:q} \
            --membership_file={output.membership:q} \
            --records_file={output.records:q} \
            --qc_audit_file={output.qc_audit:q} \
            > {log:q} 2>&1
        """

# --- Scope input-gathering helpers ---
# These use the validated scope mappings constructed in the Snakefile.
def _scope_donors(scope):
    return sorted({d for (s, d) in SCOPE_DONOR_PAIRS if s == scope})

def _scope_membership_files(scope):
    # One donor-pooled clone membership table for each donor in the scope.
    return [str(results_dir / "scopes/{}/clones/{}_IGH_clone_membership.tsv".format(scope, d))
            for d in _scope_donors(scope)]

def _scope_donor_threshold_files(scope):
    # Temporary donor fits used to derive the operative threshold and its
    # retained donor calibration table.
    return [str(results_dir / "scopes/{}/intermediate/{}/{}_IGH_donor_threshold.tsv".format(scope, d, d))
            for d in _scope_donors(scope)]

def _scope_donor_diagnostic_files(scope):
    return [str(results_dir / "scopes/{}/intermediate/{}/{}_IGH_distance_diagnostics.tsv".format(scope, d, d))
            for d in _scope_donors(scope)]

# Join donor membership, original records, grouping links and library features
# with write_donor_igh_summaries() in igh_summaries.R. choose_shm_evidence() in
# shm_summaries.R tests each retained backbone representative separately for
# each region, then chooses the qualifying source with most supported positions,
# original read support and sequence ID as successive tie-breaks. Whole V needs
# one source passing all five V sub-regions and the germline occupancy minimum.
# CDR3 and FWR4 lie outside whole V. The two shm settings below control these
# coverage requirements. Sources never pool partial coverage or vote matrices.
# Enriched members retain regional evidence IDs; absent qualifying evidence
# leaves SHM unavailable. The member shard keeps every measured rate, evidence
# ID, sequence and alignment, and collect_scope_donor_metric releases a subset
# of its columns. The member SHM shard holds one row per member and V region with
# a measured rate, and collection empties its rates and counts outside
# eligibility. The six shards and audit retain support and membership even when
# SHM cannot be scored; the released clone_shm.tsv leaves out the rows with no
# evaluable member.
rule calculate_donor_igh_summaries:
    input:
        pooled = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_pooled.tsv",
        grouping = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_grouping.tsv",
        assembly_backbone = ASSEMBLY_BACKBONE,
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        features = lambda wc: _scope_donor_igh_feature_files(wc.scope, wc.donor_id),
        membership = results_dir / "scopes/{scope}/clones/{donor_id}_IGH_clone_membership.tsv",
        records = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_lineage_records.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        shm_script = str(SCRIPTS_DIR / "calculate_donor_igh_summaries.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        igh_summaries = IGH_SUMMARIES,
        observed_representatives = str(LIB_DIR / "observed_representatives.R"),
        shm_summaries = SHM_SUMMARIES,
        gene_calls = GENE_CALLS,
        igh_features = IGH_FEATURES,
        igh_clones = IGH_CLONES,
        isotype_resolution = ISOTYPE_RESOLUTION,
        isotype_metrics = ISOTYPE_METRICS,
        pooled_records = POOLED_RECORDS,
        shm_regions = SHM_REGIONS,
        qc_audit_module = QC_AUDIT,
        gene_usage = GENE_USAGE,
        reference_metadata = rules.obtain_imgt_reference.output.metadata,
    output:
        qc_audit = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_qc_audit_shm.tsv",
        clone_shm = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_clone_shm.tsv",
        shm = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_shm.tsv",
        clones = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_clones.tsv",
        clone_isotype = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_clone_isotype.tsv",
        members = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_clone_members.tsv",
        clone_member_shm_coverage = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_clone_member_shm_coverage.tsv",
        clone_member_shm = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_clone_member_shm.tsv",
    benchmark:
        log_dir / "benchmarks/calculate_donor_igh_summaries/{scope}/{donor_id}.tsv"
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        feature_files = lambda wc, input: ",".join(input.features),
        min_subregion_fraction = config['params']['shm']['min_supported_subregion_fraction'],
        min_germline_occupied_v = config['params']['shm']['min_germline_occupied_v_positions'],
    log:
        log_dir / "calculate_donor_igh_summaries/{scope}/{donor_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.shm:q})"
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.shm_script:q} \
            --qc_audit_file={output.qc_audit:q} \
            --out_shm={output.shm:q} \
            --pooled_file={input.pooled:q} \
            --grouping_file={input.grouping:q} \
            --out_clone_shm={output.clone_shm:q} \
            --out_clones={output.clones:q} \
            --out_lineage_isotype={output.clone_isotype:q} \
            --scope_id={wildcards.scope:q} \
            --donor_id={wildcards.donor_id:q} \
            --out_members={output.members:q} \
            --out_clone_member_shm={output.clone_member_shm:q} \
            --out_clone_member_shm_coverage={output.clone_member_shm_coverage:q} \
            --records_file={input.records:q} \
            --manifest_file={input.manifest:q} \
            --feature_files={params.feature_files:q} \
            --membership_file={input.membership:q} \
            --min_subregion_fraction={params.min_subregion_fraction:q} \
            --min_germline_occupied_v={params.min_germline_occupied_v:q} \
            --reference_metadata_file={input.reference_metadata:q} \
            > {log:q} 2>&1
        """

# Sum member read support within each library and donor/scope clone through
# aggregate_clone_support_per_library() in diversity.R. Diversity retains every
# manifest library, including empty repertoires. The IGH diversity minimum gates
# entropy, Hill diversity, dominance and target-depth estimates; richness,
# support and curves remain available. rarefaction.targets supplies the depth.
# gene_usage.R counts each clone once, with its support, for V/J usage and for
# V and J pair usage, which collect_scope_gene_usage and
# collect_scope_vj_pair_usage join with the clonotype arm. This job needs clone
# membership, independently of SHM and selection.
rule calculate_scope_clone_diversity:
    input:
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        membership = lambda wc: _scope_membership_files(wc.scope),
        reference_metadata = rules.obtain_imgt_reference.output.metadata,
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "calculate_scope_clone_diversity.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        diversity = DIVERSITY,
        gene_calls = GENE_CALLS,
        gene_usage = GENE_USAGE,
    output:
        curves = results_dir / "scopes/{scope}/intermediate/{scope}_IGH_rarefaction_curves.tsv",
        clone_diversity = results_dir / "scopes/{scope}/intermediate/{scope}_IGH_diversity.tsv",
        group_abundance = temp(results_dir / "scopes/{scope}/intermediate/{scope}_IGH_group_abundance.tsv"),
        gene_usage = results_dir / "scopes/{scope}/intermediate/{scope}_IGH_gene_usage.tsv",
        vj_pair_usage = results_dir / "scopes/{scope}/intermediate/{scope}_IGH_vj_pair_usage.tsv",
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        membership_files = lambda wc, input: ",".join(input.membership),
        rarefaction_targets = named_values({k: v for k, v in config["params"]["rarefaction"]["targets"].items() if k == "IGH"}),
        minimum_consensus_count = named_values({k: v for k, v in config["params"]["diversity"]["min_consensus_count"].items() if k == "IGH"}),
    log:
        log_dir / "calculate_scope_clone_diversity/scopes/{scope}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --manifest_file={input.manifest:q} \
            --membership_files={params.membership_files:q} \
            --scope_id={wildcards.scope:q} \
            --reference_metadata_file={input.reference_metadata:q} \
            --rarefaction_targets={params.rarefaction_targets:q} \
            --minimum_consensus_count={params.minimum_consensus_count:q} \
            --curves_file={output.curves:q} \
            --clone_diversity_file={output.clone_diversity:q} \
            --group_abundance_file={output.group_abundance:q} \
            --gene_usage_file={output.gene_usage:q} \
            --vj_pair_usage_file={output.vj_pair_usage:q} \
            > {log:q} 2>&1
        """

# Fit one donor threshold from prepared sequence representatives with
# donor_distance_threshold.R. min_fit_distances and the fit_* settings control
# the SHazaM fit; insufficient distances or failed fits retain their status and
# diagnostics. These temporary outputs feed scope resolution even with a numeric
# configured threshold. Scope resolution retains their calibration evidence.
rule fit_donor_threshold:
    input:
        grouping = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_grouping.tsv",
        pooled = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_pooled.tsv",
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        fit_script = str(SCRIPTS_DIR / "fit_donor_threshold.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        donor_distance_threshold = DONOR_DISTANCE_THRESHOLD,
        gene_calls = GENE_CALLS,
        assembly_backbone = ASSEMBLY_BACKBONE,
        igh_clones = IGH_CLONES,
        igh_grouping = IGH_GROUPING,
        rearrangement_keys = REARRANGEMENT_KEYS,
    output:
        donor_threshold = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_donor_threshold.tsv"),
        donor_diagnostics = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_distance_diagnostics.tsv"),
    threads: 1
    params:
        fit_edge = config["params"]["lineage"]["fit_edge"],
        fit_specificity = config["params"]["lineage"]["fit_specificity"],
        fit_subsample = str(config["params"]["lineage"]["fit_subsample"]) if config["params"]["lineage"]["fit_subsample"] is not None else "null",
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        min_fit_distances = config['params']['lineage']['min_fit_distances'],
        # max_n is not passed here. The junction allowance was applied once in
        # prepare_igh_grouping and reaches this rule through the grouping file, so
        # calibration and assignment share one record gate.
    log:
        log_dir / "fit_donor_threshold/scopes/{scope}/{donor_id}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.donor_threshold:q})"
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.fit_script:q} \
            --input_file={input.pooled:q} \
            --grouping_file={input.grouping:q} \
            --output_file={output.donor_threshold:q} \
            --diagnostic_file={output.donor_diagnostics:q} \
            --scope_id={wildcards.scope:q} \
            --donor_id={wildcards.donor_id:q} \
            --minimum_n={params.min_fit_distances:q} \
            --fit_edge={params.fit_edge} \
            --fit_specificity={params.fit_specificity} \
            --fit_subsample={params.fit_subsample:q} \
            > {log:q} 2>&1
        """

# Freeze one operative threshold per scope and retain the compact donor
# calibration table used to review that decision. A number set for
# params.lineage.clonal_threshold is used as the operative cutoff. With "auto",
# the cutoff is the median of the usable donor fits. The donor fits still run
# and are still released, so the calibration evidence stays available whichever
# source applies. Where the setting is "auto" and no donor could be fitted at
# all, the cutoff comes from the fixed 0.1 in
# workflow/lib/scope_clone_threshold.R. Every donor in a scope is clustered
# with the same cutoff.
rule resolve_scope_clone_threshold:
    input:
        donor_thresholds = lambda wc: _scope_donor_threshold_files(wc.scope),
        donor_diagnostics = lambda wc: _scope_donor_diagnostic_files(wc.scope),
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        fit_script = str(SCRIPTS_DIR / "resolve_scope_clone_threshold.R"),
        runtime = RUNTIME,
        donor_distance_threshold = DONOR_DISTANCE_THRESHOLD,
        scope_clone_threshold = SCOPE_CLONE_THRESHOLD,
        table_io = TABLE_IO,
    output:
        scope_threshold = results_dir / "scopes/{scope}/clones/{scope}_clone_threshold.tsv",
        donor_calibration = results_dir / "scopes/{scope}/clones/{scope}_donor_clone_thresholds.tsv",
        distance_diagnostics = results_dir / "scopes/{scope}/clones/{scope}_clone_distance_diagnostics.tsv",
    threads: 1
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        clonal_threshold = config['params']['lineage']['clonal_threshold'],
        donor_files = lambda wc, input: ",".join(input.donor_thresholds),
        diagnostic_files = lambda wc, input: ",".join(input.donor_diagnostics),
    log:
        log_dir / "resolve_scope_clone_threshold/scopes/{scope}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.scope_threshold:q})"
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.fit_script:q} \
            --output_file={output.scope_threshold:q} \
            --donor_output_file={output.donor_calibration:q} \
            --diagnostics_output_file={output.distance_diagnostics:q} \
            --scope_id={wildcards.scope:q} \
            --configured_threshold={params.clonal_threshold:q} \
            --donor_files={params.donor_files:q} \
            --diagnostic_files={params.diagnostic_files:q} \
            > {log:q} 2>&1
        """

# Calculate library- and isotype-group-level B-cell antigen selection with
# BASELINe, reusing donor-pooled clone identity and the authoritative masked C
# calls. Each clone is collapsed separately within a library and isotype group,
# then the retained densities are convolved. Libraries from one donor remain
# repeated observations. baseline_selection.R owns eligibility, clone consensus
# and density estimation. collapse_min_freq controls the clone consensus base
# frequency. Partial sequences remain eligible; uncovered positions and incomplete
# codons are masked in the effective sequence, preserving the IMGT frame and gaps.
# This masking does not change SHM. Selection reads C calls directly and does not
# wait for SHM summaries. Empty eligible input retains every library/isotype row
# with estimate status and missing Sigma rather than inventing an estimate.
rule calculate_donor_selection:
    input:
        pooled = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_pooled.tsv",
        membership = results_dir / "scopes/{scope}/clones/{donor_id}_IGH_clone_membership.tsv",
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        isotype_calls = lambda wc: _scope_donor_isotype_call_files(wc.scope, wc.donor_id),
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        selection_script = str(SCRIPTS_DIR / "calculate_donor_selection.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        baseline_selection = BASELINE_SELECTION,
        gene_calls = GENE_CALLS,
        isotype_resolution = ISOTYPE_RESOLUTION,
        isotype_calls_module = ISOTYPE_CALLS,
        pooled_records = POOLED_RECORDS,
        qc_audit_module = QC_AUDIT,
    output:
        qc_audit = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_qc_audit_selection_inputs.tsv",
        selection = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_selection.tsv"),
    threads:
        1
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        collapse_min_freq = config['params']['selection']['collapse_min_freq'],
        isotype_call_files = lambda wc, input: ",".join(input.isotype_calls),
    log:
        log_dir / "calculate_donor_selection/scopes/{scope}/{donor_id}.log"
    shell:
        r"""
        set -euo pipefail
        # Compute one shard per donor without pooling the reported library metric.
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.selection_script:q} \
            --qc_audit_file={output.qc_audit:q} \
            --pooled_file={input.pooled:q} \
            --membership_file={input.membership:q} \
            --manifest_file={input.manifest:q} \
            --out_file={output.selection:q} \
            --scope_id={wildcards.scope:q} \
            --donor_id={wildcards.donor_id:q} \
            --threads={threads:q} \
            --collapse_min_freq={params.collapse_min_freq:q} \
            --isotype_call_files={params.isotype_call_files:q} \
            > {log:q} 2>&1
        """
