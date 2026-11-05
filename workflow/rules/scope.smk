# Collect supplied library and donor measurements for one publishable scope.
# scope_tables.py joins manifest context and checks collection contracts; each
# producer owns its eligibility, estimates and denominators. Selection, overlap
# and abundance keep separate collectors so they can run independently.
# The metric wildcard names only tables that share the library collection path.
LIBRARY_METRIC_PATTERNS = {
    "clonotype_members": ["metrics/{library_id}_clonotype_members.tsv"],
    "isotype_usage": ["metrics/{library_id}_isotype_usage.tsv"],
    "isotype_summary": ["metrics/{library_id}_isotype_summary.tsv"],
    "clonotypes": [
        "clonotypes/{library_id}_" + locus + "_clonotypes.tsv"
        for locus in LOCI if locus != "IGH"],
    "yield": ["metrics/{library_id}_yield.tsv"],
}
# Collection names each V/J entity beside its column, from the pinned reference
# metadata, so a VSEQ_ or JSEQ_ composite reads as its member genes.
LIBRARY_METRIC_LABELS = {
    "clonotypes": {"v_gene_key": "v_gene_names", "j_gene_key": "j_gene_names"},
}


def _scope_library_metric_files(wildcards):
    return [str(results_dir / "libraries" / library / pattern.format(library_id=library))
            for library in sorted(SCOPE_MEMBERS[wildcards.scope])
            for pattern in LIBRARY_METRIC_PATTERNS[wildcards.metric]]


# collect_scope_metric() preserves source rows and their units, adds scope, tissue
# and donor IDs from the manifest, leads with the identifiers and rejects
# libraries outside this scope. Valid header-only inputs remain valid empty tables.
rule collect_scope_library_metric:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        metric = _scope_library_metric_files,
        reference_metadata = lambda wc: ([rules.obtain_imgt_reference.output.metadata]
                                         if wc.metric in LIBRARY_METRIC_LABELS else []),
    output:
        metric = results_dir / "scopes/{scope}/metrics/{metric}.tsv",
    wildcard_constraints:
        metric = "|".join(LIBRARY_METRIC_PATTERNS),
    params:
        context_columns = ["scope_id", "tissue_id", "donor_id"],
        quoting = csv.QUOTE_NONE,
        label_columns = lambda wc: LIBRARY_METRIC_LABELS.get(wc.metric, {}),
    benchmark:
        log_dir / "benchmarks/collect_scope_library_metric/{scope}/{metric}.tsv"
    log:
        log_dir / "collect_scope_library_metric/{scope}/{metric}.log"
    script:
        "../scripts/collect_scope_metric.py"

# Donor jobs write shards sorted by the released table's key. Collection merges
# them holding one row per shard and stops on an unsorted or repeated key, so
# memory stays bounded for the wide sequence-member table too. A key column that
# DONOR_METRIC_ORDERS lists for its table follows the listed value order; other
# key columns sort as text.
DONOR_METRIC_KEYS = {
    "clone_shm": ["library_id", "sort_group_support", "sort_group_junction", "clone_id", "isotype_group", "v_region"],
    "shm": ["library_id", "isotype_group", "v_region"],
    "clones": ["donor_id", "sort_pooled_support", "clone_id"],
    "clone_isotype": ["library_id", "sort_group_support", "sort_group_junction", "clone_id", "isotype_level", "isotype"],
    "clone_members": ["library_id", "sort_group_support", "sort_group_junction", "clone_id", "sequence_id"],
    "clone_member_shm": ["library_id", "sort_group_support", "sort_group_junction", "clone_id", "sequence_id", "v_region"],
    "clone_member_shm_coverage": ["library_id", "sort_group_support", "sort_group_junction", "clone_id", "sequence_id", "v_region"],
}
# The SHM tables list isotype groups in the order of SHM_ISOTYPE_GROUPS and V
# regions in the order of shm_v_regions(V_REGION_LEVELS), as shm_summaries.R
# reports them.
SHM_KEY_ORDERS = {
    "isotype_group": ["all_IGH", "IgM", "IgG", "IgA", "IgG1", "IgG2", "IgG3", "IgG4",
                      "IgA1", "IgA2"],
    "v_region": ["whole_v", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3"],
}
# The clone isotype table lists class rows before subclass rows. Within a level,
# isotypes keep the order of ISOTYPE_COMPOSITION_LEVELS in isotype_metrics.R, and
# this one list holds both levels' values in that order.
CLONE_ISOTYPE_KEY_ORDERS = {
    "isotype_level": ["class", "subclass"],
    "isotype": ["IgM", "IgD", "IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgA1", "IgA2", "IgE"],
}
DONOR_METRIC_ORDERS = {"clone_shm": SHM_KEY_ORDERS, "shm": SHM_KEY_ORDERS,
                       "clone_member_shm": {"v_region": SHM_KEY_ORDERS["v_region"]},
                       "clone_isotype": CLONE_ISOTYPE_KEY_ORDERS,
                       "clone_member_shm_coverage": {"v_region": SHM_KEY_ORDERS["v_region"]}}
# The released member SHM table keeps a member's mutation count and rate only
# where shm_eligible is TRUE, the eligibility shm.tsv and clone_shm.tsv summarise.
# The donor shards keep every measured value.
DONOR_METRIC_GATES = {
    "clone_member_shm": {"n_mutations": "shm_eligible", "shm_frequency": "shm_eligible"},
}
# The released member table has these columns, in this order. Three take a
# release name for a shard column: v_gene_key and j_gene_key are the shard's
# v_gene and j_gene entities, and member_consensus_count its consensus_count.
# The shards also keep each member's sequence, alignments and qualified C call,
# its mutation measurements and the SHM eligibility and evidence ID of every V
# region, which the release checks read. The AIRR export carries the sequence,
# alignments and calls under the same sequence_id, and clone_member_shm.tsv the
# mutation measurements.
CLONE_MEMBER_RELEASE_COLUMNS = [
    "scope_id", "library_id", "tissue_id", "donor_id", "locus", "clone_id", "sequence_id",
    "clonotype_key", "v_call", "v_gene_key", "v_gene_names", "j_call", "j_gene_key",
    "j_gene_names", "junction", "junction_aa", "junction_length", "isotype_class",
    "isotype_subclass", "isotype_switch_state", "member_consensus_count"]
CLONE_MEMBER_SHARD_NAMES = {"v_gene_key": "v_gene", "j_gene_key": "j_gene",
                            "member_consensus_count": "consensus_count"}
DONOR_METRIC_COLUMNS = {
    "clone_members": {name: CLONE_MEMBER_SHARD_NAMES.get(name, name)
                      for name in CLONE_MEMBER_RELEASE_COLUMNS},
}
# The released clone SHM table leaves out the clone, isotype and region rows
# with no SHM-eligible member, which the shards keep. group_abundance.tsv keeps
# those clones' support, and shm.tsv and the headline come from the donor
# records rather than from these rows.
DONOR_METRIC_DROPPED = {"clone_shm": {"mutation_state": ["no_evaluable_members"]}}


def _scope_donor_metric_files(wildcards):
    return [str(results_dir / "scopes" / wildcards.scope / "intermediate" /
                donor / (donor + "_" + wildcards.donor_metric + ".tsv"))
            for scope, donor in SCOPE_DONOR_PAIRS if scope == wildcards.scope]


rule collect_scope_donor_metric:
    input:
        exports = lambda wc: [str(results_dir / "libraries" / library / "export" /
                                  (library + "_rearrangements.tsv.gz"))
                              for library in sorted(SCOPE_MEMBERS[wc.scope])]
                             if wc.donor_metric == "clone_member_shm_coverage" else [],
        collection_module = SCOPE_TABLES_MODULE,
        tables = _scope_donor_metric_files,
    output:
        table = results_dir / "scopes/{scope}/metrics/{donor_metric}.tsv",
    wildcard_constraints:
        donor_metric = "|".join(DONOR_METRIC_KEYS),
    params:
        key_columns = lambda wc: DONOR_METRIC_KEYS[wc.donor_metric],
        key_orders = lambda wc: DONOR_METRIC_ORDERS.get(wc.donor_metric, {}),
        gated_columns = lambda wc: DONOR_METRIC_GATES.get(wc.donor_metric, {}),
        output_columns = lambda wc: DONOR_METRIC_COLUMNS.get(wc.donor_metric, {}),
        dropped_values = lambda wc: DONOR_METRIC_DROPPED.get(wc.donor_metric, {}),
    benchmark:
        log_dir / "benchmarks/collect_scope_donor_metric/{scope}/{donor_metric}.tsv"
    log:
        log_dir / "collect_scope_donor_metric/{scope}/{donor_metric}.log"
    script:
        "../scripts/collect_scope_donor_metric.py"

# Collect one curve row per library, locus and support depth from each arm.
# Producers calculate curves from the same support vectors as their diversity
# rows; collection adds scope, donor and tissue IDs without estimating diversity.
# Rows are sorted by library and locus, keeping each curve's support_depth order.
rule collect_scope_rarefaction_curves:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        metric = lambda wc: [str(results_dir / "libraries" / library / "metrics" /
                                  (library + "_rarefaction_curves.tsv"))
                             for library in sorted(SCOPE_MEMBERS[wc.scope])] +
                            [str(results_dir / "scopes" / wc.scope / "intermediate" /
                                 (wc.scope + "_IGH_rarefaction_curves.tsv"))],
    output:
        metric = results_dir / "scopes/{scope}/metrics/rarefaction_curves.tsv",
    params:
        context_columns = ["scope_id", "donor_id", "tissue_id"],
        quoting = csv.QUOTE_MINIMAL,
        sort_columns = ["library_id", "locus"],
    benchmark:
        log_dir / "benchmarks/collect_scope_rarefaction_curves/{scope}.tsv"
    log:
        log_dir / "collect_scope_rarefaction_curves/{scope}.log"
    script:
        "../scripts/collect_scope_metric.py"

# Collect one diversity row per library and locus from both arms: the per-library
# non-IGH clonotype rows and the scope's IGH clone rows. group_unit says which.
# Rows are sorted by library and locus.
rule collect_scope_diversity:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        metric = lambda wc: [str(results_dir / "libraries" / library / "metrics" /
                                  (library + "_diversity.tsv"))
                             for library in sorted(SCOPE_MEMBERS[wc.scope])] +
                            [str(results_dir / "scopes" / wc.scope / "intermediate" /
                                 (wc.scope + "_IGH_diversity.tsv"))],
    output:
        metric = results_dir / "scopes/{scope}/metrics/diversity.tsv",
    params:
        context_columns = ["scope_id", "donor_id", "tissue_id"],
        quoting = csv.QUOTE_NONE,
        sort_columns = ["library_id", "locus"],
    benchmark:
        log_dir / "benchmarks/collect_scope_diversity/{scope}.tsv"
    log:
        log_dir / "collect_scope_diversity/{scope}.log"
    script:
        "../scripts/collect_scope_metric.py"

# Collect V and J usage per library, locus and gene type from both arms: the
# per-library non-IGH clonotype rows and the scope's IGH clone rows. group_unit
# says which unit n_groups counts. Rows are sorted by library and locus, and the
# producers' order, which keeps a gene equivalence group's rows together, is
# kept within.
rule collect_scope_gene_usage:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        metric = lambda wc: [str(results_dir / "libraries" / library / "metrics" /
                                  (library + "_clonotype_gene_usage.tsv"))
                             for library in sorted(SCOPE_MEMBERS[wc.scope])] +
                            [str(results_dir / "scopes" / wc.scope / "intermediate" /
                                 (wc.scope + "_IGH_gene_usage.tsv"))],
    output:
        metric = results_dir / "scopes/{scope}/metrics/gene_usage.tsv",
    params:
        context_columns = ["scope_id", "donor_id", "tissue_id"],
        quoting = csv.QUOTE_NONE,
        sort_columns = ["library_id", "locus"],
    benchmark:
        log_dir / "benchmarks/collect_scope_gene_usage/{scope}.tsv"
    log:
        log_dir / "collect_scope_gene_usage/{scope}.log"
    script:
        "../scripts/collect_scope_metric.py"

# Collect V and J pair usage per library and locus from the same two arms.
# group_unit says which unit n_groups counts, and a composite entity counts once
# under its own key. Rows are sorted by library and locus, and the producers'
# order, by n_groups and then the V and J keys, is kept within.
rule collect_scope_vj_pair_usage:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        metric = lambda wc: [str(results_dir / "libraries" / library / "metrics" /
                                  (library + "_clonotype_vj_pair_usage.tsv"))
                             for library in sorted(SCOPE_MEMBERS[wc.scope])] +
                            [str(results_dir / "scopes" / wc.scope / "intermediate" /
                                 (wc.scope + "_IGH_vj_pair_usage.tsv"))],
    output:
        metric = results_dir / "scopes/{scope}/metrics/vj_pair_usage.tsv",
    params:
        context_columns = ["scope_id", "donor_id", "tissue_id"],
        quoting = csv.QUOTE_NONE,
        sort_columns = ["library_id", "locus"],
    benchmark:
        log_dir / "benchmarks/collect_scope_vj_pair_usage/{scope}.tsv"
    log:
        log_dir / "collect_scope_vj_pair_usage/{scope}.log"
    script:
        "../scripts/collect_scope_metric.py"

# build_headline() writes one scientific summary row per manifest library.
# Required source rows must exist even when an estimate is unavailable; missing
# values stay missing. Project annotations belong only to report preparation.
rule build_scope_headline:
    input:
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        yields = results_dir / "scopes/{scope}/metrics/yield.tsv",
        qc = results_dir / "scopes/{scope}/metrics/qc.tsv",
        diversity = results_dir / "scopes/{scope}/metrics/diversity.tsv",
        isotype = results_dir / "scopes/{scope}/metrics/isotype_usage.tsv",
        isotype_summary = results_dir / "scopes/{scope}/metrics/isotype_summary.tsv",
        shm = results_dir / "scopes/{scope}/metrics/shm.tsv",
        headline_module = str(Path(WORKFLOW_DIR) / "workflow/lib/headline.py"),
    output:
        headline = results_dir / "scopes/{scope}/metrics/headline.tsv",
    benchmark:
        log_dir / "benchmarks/build_scope_headline/{scope}.tsv"
    log:
        log_dir / "build_scope_headline/{scope}.log"
    script:
        "../scripts/build_scope_headline.py"


# scope_qc_audit.py validates and combines supplied population counts and read
# support into the scope's working audit. Audit rows identify their stage, unit
# and locus; overlapping failures are not sequential losses. The same job writes
# qc.tsv, which divides those support counts per library and locus and sums
# grouping reasons from the library AIRR exports, and the released view
# qc_trail.tsv, whose strata follow the SHM table orders.
# This collection is downstream of its producers and never gates an analysis.
def _scope_qc_shards(wildcards):
    libraries = [
        str(results_dir / "libraries" / library / "qc" / (library + "_qc_audit_" + stage + ".tsv"))
        for library in sorted(SCOPE_MEMBERS[wildcards.scope])
        for stage in ("reannotate", "qc", "igh_grouping", "nonigh_grouping",
                      "isotype_calls", "isotype", "shm_measurement")]
    donors = [
        str(results_dir / "scopes" / wildcards.scope / "intermediate" / donor /
            (donor + "_qc_audit_" + stage + ".tsv"))
        for scope, donor in SCOPE_DONOR_PAIRS if scope == wildcards.scope
        for stage in ("clones", "shm", "selection_inputs")]
    return libraries + donors


rule collect_scope_qc_audit:
    input:
        coverage = results_dir / "scopes/{scope}/metrics/clone_member_shm_coverage.tsv",
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
        shards = _scope_qc_shards,
        diversity = [results_dir / "scopes/{scope}/metrics/diversity.tsv"],
        selection = results_dir / "scopes/{scope}/metrics/selection_pressure.tsv",
        # A library with no assembly provenance sidecar is left out here, and the
        # collector records its assembly_history as unknown_history. The check runs
        # while the DAG is built, so a sidecar written later in the same run does
        # not change that status.
        producers = lambda wc: [
            str(path) for library in sorted(SCOPE_MEMBERS[wc.scope])
            for path in [results_dir / "libraries" / library / "trust4" /
                         (library + "_assembly_provenance.json")] if path.exists()],
        exports = lambda wc: [
            str(results_dir / "libraries" / library / "export" / (library + "_rearrangements.tsv.gz"))
            for library in sorted(SCOPE_MEMBERS[wc.scope])],
        module = str(Path(WORKFLOW_DIR) / "workflow/lib/scope_qc_audit.py"),
    output:
        audit = results_dir / "scopes/{scope}/metrics/qc_audit.tsv",
        qc = results_dir / "scopes/{scope}/metrics/qc.tsv",
        trail = results_dir / "scopes/{scope}/metrics/qc_trail.tsv",
    params:
        value_orders = SHM_KEY_ORDERS,
    benchmark:
        log_dir / "benchmarks/collect_scope_qc_audit/{scope}.tsv"
    log:
        log_dir / "collect_scope_qc_audit/{scope}.log"
    script:
        "../scripts/collect_scope_qc_audit.py"


def _scope_donor_nonigh_group_files(scope, donor, locus):
    # Non-IGH overlap consumes the authoritative rearrangement-group tables.
    # IGH is deliberately absent: its scope-specific clone membership already
    # contains sample, clone, and read-support fields and is its sole substrate.
    tmpl = "libraries/{}/clonotypes/{}_{}_clonotypes.tsv"
    return [str(results_dir / tmpl.format(g, g, locus))
            for (t, g) in scope_donor_to_libraries[(scope, donor)]]

def _scope_donor_airr_files(scope, donor):
    # In the same library order as the clonotype tables, which name each
    # representative record whose sequence these tables hold.
    return [str(results_dir / "libraries/{}/airr/{}_airr_qc_pass.tsv".format(g, g))
            for (t, g) in scope_donor_to_libraries[(scope, donor)]]

def _scope_selection_files(scope):
    # Each donor shard contains library-by-isotype BASELINe summaries.
    return [str(results_dir / "scopes/{}/intermediate/{}/{}_IGH_selection.tsv".format(scope, d, d))
            for d in _scope_donors(scope)]

def _scope_nonigh_group_abundance_files(scope):
    # The non-IGH half of the clone-consensus_count table is scope-invariant, so it is
    # gathered per library rather than recomputed; the IGH half arrives as one
    # scope-level shard from calculate_scope_clone_diversity. Collection merges
    # the two halves, so this list stays in library order.
    return [str(results_dir / "libraries/{}/metrics/{}_group_abundance.tsv".format(g, g))
            for g in sorted(SCOPE_MEMBERS[scope])]

def _scope_donor_overlap_pair_files(scope):
    return [str(results_dir / "scopes/{}/intermediate/{}/{}_donor_overlap_pairs.tsv".format(scope, d, d))
            for (s, d) in SCOPE_MULTI_PAIRS if s == scope]

def _scope_donor_overlap_files(scope):
    return [str(results_dir / "scopes/{}/intermediate/{}/{}_donor_overlap_summary.tsv".format(scope, d, d))
            for (s, d) in SCOPE_MULTI_PAIRS if s == scope]

def _scope_shared_clone_files(scope):
    return [str(results_dir / "scopes/{}/intermediate/{}/{}_shared_groups.tsv".format(scope, d, d))
            for (s, d) in SCOPE_MULTI_PAIRS if s == scope]


# Calculate pairwise overlap among a donor's libraries within one scope, and list
# the clones it is counting. IGH consumes this scope's IGH clone membership;
# non-IGH loci use the authoritative nucleotide rearrangement-group key. Only
# donors with at least MIN_OVERLAP_SAMPLES libraries in the scope enter this rule,
# so a donor with one library has no row in any of the three output tables.
# calculate_donor_overlap.R builds positive-support sharing units;
# overlap_metrics.R calculates pair and donor summaries. Support fractions use
# each endpoint's complete locus support. params.overlap.min_consensus_count flags
# shallow pairs without masking values. Empty endpoints make comparisons
# unavailable.

rule calculate_donor_overlap:
    input:
        unpack(lambda wc: {
            locus.lower(): _scope_donor_nonigh_group_files(
                wc.scope, wc.donor_id, locus
            )
            for locus in config['params']['overlap']['loci']
            if locus != "IGH"
        }),
        membership = results_dir / "scopes/{scope}/clones/{donor_id}_IGH_clone_membership.tsv",
        # Lineage records supply junctions and calls for shared-clone descriptions,
        # and sequences that decide between tied observations. The overlap
        # arithmetic uses membership support.
        records = results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_IGH_lineage_records.tsv",
        # Common-QC AIRR supplies the sequence of each clonotype's representative,
        # which decides between equally supported observations of a clonotype.
        airr_qc_pass = lambda wc: _scope_donor_airr_files(wc.scope, wc.donor_id),
        analysis_image = ANALYSIS_IMAGE_RECEIPT,
        container_exec = CONTAINER_EXEC,
        script = str(SCRIPTS_DIR / "calculate_donor_overlap.R"),
        runtime = RUNTIME,
        table_io = TABLE_IO,
        gene_calls = GENE_CALLS,
        overlap_metrics = OVERLAP_METRICS,
        observed_representatives = str(LIB_DIR / "observed_representatives.R"),
        diversity = DIVERSITY,
        gene_usage = GENE_USAGE,
        pooled_records = POOLED_RECORDS,
        reference_metadata = rules.obtain_imgt_reference.output.metadata,
    output:
        donor_overlap_pairs = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_donor_overlap_pairs.tsv"),
        donor_overlap_summary       = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_donor_overlap_summary.tsv"),
        shared_groups       = temp(results_dir / "scopes/{scope}/intermediate/{donor_id}/{donor_id}_shared_groups.tsv"),
    params:
        docker_image = config['deps']['immcantation_docker'],
        analysis_sif = str(ANALYSIS_SIF),
        tissue_ids = lambda wc: ",".join(t for (t, g) in scope_donor_to_libraries[(wc.scope, wc.donor_id)]),
        library_ids = lambda wc: ",".join(g for (t, g) in scope_donor_to_libraries[(wc.scope, wc.donor_id)]),
        airr_files = lambda wc, input: ",".join(input.airr_qc_pass),
        loci = ",".join(config['params']['overlap']['loci']),
        minimum_consensus_count = named_values({
            locus: config['params']['overlap']['min_consensus_count'].get(locus, 0)
            for locus in config['params']['overlap']['loci']
        }),
        # One --<locus>_files option per configured non-IGH locus, built from the
        # same helper as the unpacked inputs so an unconfigured locus passes no
        # option. Each value is quoted here because the option list is
        # substituted unquoted below. The option name avoids an f-string because
        # Snakemake 7.32.4 inserts spaces into f-strings in .smk files under
        # Python 3.12.
        group_file_options = lambda wc: " ".join(
            "--" + locus.lower() + "_files=" + shlex.quote(",".join(
                _scope_donor_nonigh_group_files(wc.scope, wc.donor_id, locus)))
            for locus in config['params']['overlap']['loci'] if locus != "IGH"),
    log:
        log_dir / "calculate_donor_overlap/scopes/{scope}/{donor_id}.log"
    shell:
        r"""
        set -euo pipefail
        {input.container_exec:q} {params.analysis_sif:q} {params.docker_image:q} \
            Rscript {input.script:q} \
            --membership_file={input.membership:q} \
            --records_file={input.records:q} \
            --scope_id={wildcards.scope:q} \
            --donor_id={wildcards.donor_id:q} \
            --library_ids={params.library_ids:q} \
            --tissue_ids={params.tissue_ids:q} \
            --loci={params.loci:q} \
            --minimum_consensus_count={params.minimum_consensus_count:q} \
            {params.group_file_options} \
            --airr_files={params.airr_files:q} \
            --donor_overlap_pairs_file={output.donor_overlap_pairs:q} \
            --donor_overlap_file={output.donor_overlap_summary:q} \
            --shared_groups_file={output.shared_groups:q} \
            --reference_metadata_file={input.reference_metadata:q} \
            > {log:q} 2>&1
        """

# The pair table, the donor summary and the shared-clone table are collected
# together because one donor calculation writes all three from the same
# repertoire pairs. collect_scope_overlap() writes headers when the scope has no
# eligible donors.
rule collect_scope_overlap:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        donor_overlap_pairs = lambda wc: _scope_donor_overlap_pair_files(wc.scope),
        donor_overlap_summary       = lambda wc: _scope_donor_overlap_files(wc.scope),
        shared_groups       = lambda wc: _scope_shared_clone_files(wc.scope),
    output:
        donor_overlap_pairs = results_dir / "scopes/{scope}/metrics/donor_overlap_pairs.tsv",
        donor_overlap_summary       = results_dir / "scopes/{scope}/metrics/donor_overlap_summary.tsv",
        shared_groups       = results_dir / "scopes/{scope}/metrics/shared_groups.tsv",
    params:
        scope_id = lambda wc: wc.scope,
    log:
        log_dir / "collect_scope_overlap/scopes/{scope}.log"
    script:
        "../scripts/collect_scope_overlap.py"

# collect_scope_selection() combines supplied library-by-isotype BASELINe
# summaries in library order. Donor selection jobs own clone consensus and
# usable-site decisions. This collection needs neither overlap nor non-IGH tables.
rule collect_scope_selection:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        selection = lambda wc: _scope_selection_files(wc.scope),
    output:
        selection = results_dir / "scopes/{scope}/metrics/selection_pressure.tsv",
    params:
        scope_id = lambda wc: wc.scope,
    log:
        log_dir / "collect_scope_selection/scopes/{scope}.log"
    script:
        "../scripts/collect_scope_selection.py"

# Put both arms' clone-size vectors in one table: for every library and locus
# that receives a diversity row, one row per clone with its read support.
# n_members counts the distinct retained member sequences of each IGH clone or
# clonotype. The released diversity indices use these support vectors, so a
# reader can reproduce them or compute an estimator this pipeline does not
# publish. Nothing is recomputed here. Each
# arm wrote its vector from the stage that already held it. Rows are sorted by
# library, then locus, keeping each producer's row order within a locus, which
# is group_rank order.
rule collect_scope_group_abundance:
    input:
        collection_module = SCOPE_TABLES_MODULE,
        nonigh   = lambda wc: _scope_nonigh_group_abundance_files(wc.scope),
        igh      = results_dir / "scopes/{scope}/intermediate/{scope}_IGH_group_abundance.tsv",
        manifest = results_dir / "scopes/{scope}/manifest/scope_manifest.tsv",
    output:
        group_abundance = results_dir / "scopes/{scope}/metrics/group_abundance.tsv",
    params:
        scope_id = lambda wc: wc.scope,
    log:
        log_dir / "collect_scope_group_abundance/scopes/{scope}.log"
    script:
        "../scripts/collect_scope_group_abundance.py"


rule write_scope_intermediate_readme:
    input:
        template = str(Path(WORKFLOW_DIR) / "workflow/release/intermediate_README.md"),
    output:
        readme = results_dir / "scopes/{scope}/intermediate/README.md",
    shell:
        "cp {input.template:q} {output.readme:q}.tmp && mv {output.readme:q}.tmp {output.readme:q}"
