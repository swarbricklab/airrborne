# QC matrix

This guide lists the quality-control (QC) and evidence requirements applied by the pipeline, their effects and the code that implements them. It covers record QC, analysis eligibility and estimate availability. The counts for one run are in `metrics/detail/qc_trail.tsv`, described under [QC trail table](#qc-trail-table). The [glossary](glossary.md) defines exclusion, eligibility, call masking and grouping, and the [interpretation guide](interpretation.md) explains how to read the resulting measurements.

## Path of a record

1. TRUST4 assembles receptor sequences from reads and reports their CDR3 calls. The pipeline removes calls that TRUST4 imputed.
2. IgBLAST annotates the remaining records on the seven canonical loci. Identical sequences keep their own record and assembly identities.
3. Every locus uses the same six common-QC criteria.

Records passing common QC then supply several analyses:

- Grouping selects one representative per IGH, IGK or IGL contig, then merges identical full sequences within each library and locus. TCR loci retain their eligible variants before exact-sequence merging. Each eligible source record contributes its support once, and the source records remain available.
- IGH clones are inferred within each donor and scope. The other loci form single-chain clonotypes within each library.
- Isotype composition uses common-QC IGH records, including those outside grouping.
- SHM and selection use assigned IGH members, with separate evidence requirements for each analysis.
- Diversity and pair comparisons use grouped receptor units.

## Requirements

Owner paths are relative to `workflow/`, except `Snakefile` at the repository root. A requirement written as `section.name=value` is a setting under `params` in `config/airrborne/config.yaml`, shown at its pipeline default. "Recorded in" names the `qc_trail.tsv` stage or the status column that shows the decision in a release.

### Membership and assembly

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| All | `usage=fail_qc` in the membership CSV | Leaves the library out of that project's scope. The `cohort` scope is built from the samplesheet and still includes it. | `Snakefile` | Scope manifest |
| All | `trust4.contig_min_cov=0` | TRUST4's own contig coverage filter, disabled at zero. A reused assembly with no producer history has no known setting. | `rules/assembly.smk` | `<library_id>_assembly_provenance.json`, `assembly_producers` in `provenance.json` |

### Parsing and annotation

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| All | CDR3 score within 1e-6 of TRUST4's imputation score 0.01 | Removes the record before annotation, because TRUST4 imputed that CDR3. | `lib/trust4_parsing.R::flag_imputed_cdr3_records` | `observed_cdr3` |
| All | TRUST4 locus is IGH, IGK, IGL, TRA, TRB, TRD or TRG | Only these records go to IgBLAST. A record on any other TRUST4 locus is not annotated. | `scripts/run_igblast_all_loci.sh`, counted by `scripts/parse_trust4_airr.R` | `igblast_input` |
| All | IgBLAST and MakeDb return a record | A record they cannot annotate does not continue. This includes a record whose first-listed J allele has no IgBLAST auxiliary data entry, as the [reference explanation](interpretation.md#germline-reference) describes. | `scripts/build_reannotated_airr.R` | `igblast_annotated` |
| TRA, TRD | A V call from the shared TRA/TRD V pool | Resolves the record to TRA or TRD by its J call, when every V candidate belongs to that pool and every present D or C call agrees with the J locus. Raw calls are kept. | `lib/reannotated_airr.R::resolve_shared_trav_dv_locus` | Annotated `locus` |
| All | IgBLAST's locus matches TRUST4's locus | A diagnostic comparison that excludes nothing. | `scripts/build_reannotated_airr.R` | `locus_agreement` |

### Annotation and common QC

A record is retained only when it passes all six criteria. Each criterion sees the same input, so a record can fail several at once.

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| All | Every called V and J gene agrees with the assigned locus | An absent V or J call fails. A call naming several genes of the locus passes, so this criterion does not require one resolved gene. | `lib/record_qc.R::airr_qc_flags` | `qc_locus_match` |
| All | Any present D call agrees with the locus | A missing D call passes. A D call on IGK, IGL, TRA or TRG fails, because those loci have no D genes. | `lib/record_qc.R::airr_qc_flags` | `qc_d_locus_match` |
| All | Aligned N fraction <= `airr_qc.max_alignment_n_fraction=0.1` | The denominator is aligned positions excluding gaps, and it must be positive. At 1 the ceiling is removed and the positive denominator still applies. | `lib/record_qc.R::airr_qc_flags` | `qc_low_n` |
| All | IgBLAST and MakeDb `productive` is true | Applies at every locus. Unaligned terminal bases are outside the stop-codon assessment. | `lib/record_qc.R::airr_qc_flags` | `qc_productive` |
| All | Junction length is finite, nonnegative and divisible by three | A zero-length junction passes here. Grouping applies its own nonempty-junction requirement. | `lib/record_qc.R::airr_qc_flags` | `qc_junction_mod3` |
| All | Read support is positive | Malformed support stops the run as an input-contract error. | `lib/record_qc.R::airr_qc_flags` | `qc_positive_support` |

Common QC does not test the AIRR `vj_in_frame` flag, impose a whole-alignment length floor or apply any diversity or SHM requirement. The retained population is recorded as `common_qc`.

### Grouping

Grouping combines eligible records into members, clones and clonotypes. The AIRR export retains every common-QC record with its original support, including records outside grouping.

A V or J assignment must have a complete candidate list. If equally best-scoring candidates reach the annotation cap, their resolution remains uncertain. The default caps are 50 V alignments and 20 J alignments per query. Multi-gene equivalence also requires complete, unambiguous reference sequences for every named allele. Composite entities join only themselves.

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| IGH | Nonempty junction of A, C, G, T, N and ? only, with at most `lineage.max_n=0` N or ? symbols, and one resolved V and one resolved J germline entity | Gates threshold fitting and clone assignment. A present C call must have a complete candidate list and agree with the resolved locus; an absent C call passes. SHM coverage is not required. | `lib/igh_grouping.R::prepare_igh_grouping` | `grouping_eligible` |
| IGK, IGL, TCR | Nonempty A, C, G, T junction with one resolved V and one resolved J germline entity | Defines the clonotype key. A complete candidate list naming alleles of one exact reference sequence class qualifies. Present C calls must have complete candidate lists and agree with the resolved locus; absent C calls pass. | `lib/rearrangement_groups.R::build_rearrangement_group_table` | `grouping_eligible` |
| IGH, IGK, IGL | The dominant common-QC variant of one contig is chosen before grouping eligibility | Ranks by original record support, then an A/C/G/T-only junction with single-gene V and J calls, V identity, usable alignment length, junction, full sequence, V, J and C calls, and sequence ID. If the dominant variant fails grouping, no minor sibling replaces it. Eligible source support transfers once to a retained representative, and source rows and links are kept. TCR loci keep every eligible variant. | `lib/assembly_backbone.R::prepare_assembly_mapping` | `contig_collapse` |
| All | Identical full DNA among contig representatives becomes one member within each library and locus | Analytical annotations must agree, or the run stops. The member representative is chosen by the same ranking on each contig representative's own original support. Each eligible source record contributes its support once. Source contig identities and both representative links remain. | `lib/assembly_backbone.R::prepare_assembly_mapping` | `exact_sequence_collapse` |
| IGH | A donor fit needs `lineage.min_fit_distances=500` finite distances, at least two distinct values and a fitted cutoff from 0 to 1 | Decides whether the donor's fit informs the scope threshold. A donor that cannot be fitted is still clustered and keeps its whole repertoire. | `lib/donor_distance_threshold.R` | `fit_status` in `igh_clone_threshold/donor_fits.tsv` |
| IGH | `lineage.clonal_threshold=auto` takes the median usable donor fit and falls back to 0.1 when no fit is usable. A configured number is used instead. | One threshold for the whole scope. SCOPer `hierarchicalClones` clusters the pooled donor representatives with single linkage on junction Hamming distance normalised by length, within partitions sharing exact V entity, J entity and junction length. A composite entity joins only itself in fitting and assignment. The clone label propagates to each representative's source records. | `lib/scope_clone_threshold.R`, `scripts/assign_scope_donor_clones.R`, `lib/igh_clones.R` | `threshold_source`, `receptor_units` |
| IGK, IGL, TCR | Exact nucleotide junction with the resolved V and J entities, within one library and locus | Combines observations and their support. It is not a sequence-quality test. | `lib/rearrangement_groups.R` | `receptor_units` |

### Isotype

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| IGH | A present constant call compatible with IGH, with a complete candidate list and finite informative C alignment positions >= `airr_qc.min_informative_c_alignment_positions=40` | Masks `qualified_c_call` only, and the record stays where it was. At zero, presence, locus compatibility, candidate completeness and finite evidence still apply. Uses every common-QC IGH record, including grouping-ineligible ones. | `lib/isotype_calls.R::c_call_qc_flags` | `c_call_present`, `c_call_locus_match`, `c_call_candidates_complete`, `c_call_enough_evidence`, `isotype_callable` |
| IGH | The qualified call resolves to one class or one subclass | Class and subclass have separate callable denominators, so their rows overlap and must not be summed. A subclass tie can still resolve at class level. A C candidate list at its annotation cap cannot establish resolution. SHM and selection isotype strata need the relevant resolution, and all-IGH summaries do not. | `lib/isotype_resolution.R`, `lib/isotype_metrics.R` | `class_callable`, `subclass_callable` |

### Somatic hypermutation

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| IGH | `shm.min_position_votes=2` agreeing votes for a supported IMGT V position | Defines the supported positions SHM coverage is measured against. At zero the vote floor is removed, and unmapped positions are still missing. Common QC and clone membership are unchanged. | `lib/shm_regions.R`, `lib/igh_features.R` | `shm_measured`, `shm_supported_positions` |
| IGH | For one region: a finite rate, a positive germline-occupied denominator and a supported fraction >= `shm.min_supported_subregion_fraction=1.0` | Tests each retained contig representative carrying the member's exact full DNA on its own. At zero the finite rate and positive denominator are still required. | `lib/shm_summaries.R::shm_region_eligible` | `shm_eligible_members` |
| IGH | For whole V: one source passes all five subregions and has >= `shm.min_germline_occupied_v_positions=270` germline-occupied V positions | Applies to whole V only. At zero the position floor is removed. Coverage from different sources cannot be combined. | `lib/shm_summaries.R` | `shm_eligible_members` |
| IGH | At least one same-sequence contig representative passes the region's rule | Selects the passing source with the most supported regional positions, then the greatest original record support, then the stable sequence ID, and records its ID in the released coverage audit and donor member table. With no passing source, SHM is unavailable for that member and region. Votes are never pooled. | `lib/shm_summaries.R` | `evidence_sequence_id` in `metrics/detail/igh/clone_member_shm_coverage.tsv`, `n_clones_eligible` in `shm.tsv` |

Each clone contributes the median of its eligible members to the library SHM summary.

### Selection

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| IGH | An assigned clone, a sequence alignment and a D-masked germline alignment | No SHM coverage or full-length requirement. Incomplete effective-sequence codons are masked before observed and expected selection calculations. The all-IGH group and each qualified isotype group use their own inputs. | `lib/baseline_selection.R::selection_record_eligible` | `selection_input` |
| IGH | `selection.collapse_min_freq=0.6` base frequency, and a usable regional mutation density | The frequency must be in `(0, 1]`. It defines one effective sequence per clone within a library and isotype group. Non-covering members remain in the consensus denominator. A missing density gives a missing estimate with its status and removes nothing from other analyses. | `lib/baseline_selection.R` | Selection `n_clone_consensus_sequences`, `n_clone_consensus_sequences_used_<region>` and status columns |

### Diversity and pairs

| Loci | Requirement | Effect and conditions | Owner | Recorded in |
|---|---|---|---|---|
| All seven loci | Per-locus `diversity.min_consensus_count=0` | Below a raised minimum, Shannon entropy, Hill q = 1 and q = 2, dominance and target-depth values are unavailable. Richness, support and diagnostic curves stay available. | `lib/diversity.R`, `scripts/calculate_library_clonotype_diversity.R`, `scripts/calculate_scope_clone_diversity.R` | `diversity_status` |
| Targeted loci | `rarefaction.targets={IGH: 1000, TRB: 400}`, with observed support reaching the target and meeting the diversity minimum | An omitted locus has no target value and no diagnostic curve. Falling short of a target removes no group and leaves observed diversity available. | `lib/diversity.R` | `rarefaction_status` |
| `overlap.loci`, by default IGH and TRB | Two of the donor's libraries in the scope, and defined estimator denominators | Empty endpoints have unavailable comparisons. IGH compares clones and the other loci compare clonotypes. | `rules/scope.smk`, `scripts/calculate_donor_overlap.R`, `lib/overlap_metrics.R`, `lib/diversity.R` | `pair_status` |
| `overlap.loci` | Per-locus `overlap.min_consensus_count=0` at both endpoints | Controls default pair report inclusion only. Zero applies no minimum. Values below the minimum stay released, and the setting changes no clone membership or diversity eligibility. | `lib/overlap_metrics.R`, report pair views | `pair_status` |

## Configured values

Each `section.name=value` above is a key under `params` in `config/airrborne/config.yaml`. The comment beside each key gives its pipeline default and, where they exist, its disabling value and the tool's own default. `workflow/schemas/config.schema.yaml` holds the valid bounds. A release records the values it used in `effective_config.yaml`.

Some settings have no disabling value. `overlap.loci` cannot be empty, so removing a locus removes its rows. `selection.collapse_min_freq` rejects zero.

## Input-contract errors

The following input problems stop the run and need correction before analysis can continue. They are reported as errors outside the biological exclusion counts:

- a duplicate identifier
- unreadable or malformed read support
- contradictory analytical annotations for identical full sequences
- a V, D or J CIGAR whose query coordinates contradict the sequence

Reannotation converts V, D and J CIGARs from the pinned Change-O convention to AIRR query and reference semantics and checks their query coordinates. Publication checks the staged files again. Records with productive indels remain subject to the same criteria as other records.

## QC trail table

`metrics/detail/qc_trail.tsv` records the populations used by each analysis, with input and output counts for each library. Every count has a stated unit, and read support is included where a support denominator is defined.

Estimate availability is recorded beside the estimates in the metric tables; pair availability is in the overlap tables. `metrics/qc.tsv` summarises the trail as fractions of annotated, common-QC, grouping-eligible and callable support for each library and locus. It also includes grouping-reason fractions calculated from the AIRR exports.

### Columns

- The key is scope, library, locus, stage, isotype group and V region.
- `stage_order` numbers the stages from 1 to 29 in the order of the table below. It orders display only and does not control execution.
- `locus_source` is `trust4` where `locus` is the locus TRUST4 assigned, `igblast` where it is IgBLAST's and empty for counts over the whole assembly.
- `kind` identifies how to read the population change, as listed below.
- `n_units_in` and `n_units_out` count the stated `input_unit` and `output_unit`. The [audit unit vocabulary](glossary.md#audit-units) defines each token. A grouping row can change the unit.
- `consensus_count_in` and `consensus_count_out` are empty when no support denominator is defined, such as for positions.
- A blank count is unavailable or does not apply. An observed empty population has zero counts.

| `kind` | Meaning |
|---|---|
| `checkpoint` | Counts one population. |
| `criterion` | Counts the input and records meeting one of several criteria tested on that input. |
| `exclusion` | Counts the records continuing from an input population; the others take no further part. |
| `eligibility` | Counts the subset used by one analysis; other records remain in their population. |
| `masking` | Counts records whose constant-region call is retained; other records remain with the call blanked. |
| `grouping` | Counts source units and the units receiving their support. |
| `diagnostic` | Records a comparison with no exclusions. |

`decision_status` is `applied` or `disabled`:

- `disabled` marks three rows whose numeric requirement sits at its disabling value: `qc_low_n` at an N fraction of 1, `c_call_enough_evidence` at zero C positions and `shm_supported_positions` at zero votes. The remaining conditions of each still apply.
- The two SHM coverage floors at zero leave `shm_eligible_members` marked `applied`.
- The TRUST4 coverage filter, the diversity minima and the pair-support minima have no trail row. Their values are in `effective_config.yaml`, and the minima's effects are in the metric tables' status columns.
- `assembly_producers` in `provenance.json` lists the libraries whose reused assemblies have no producer record.

### Stages

| Order | Stage | Input and output population | Context |
|---|---|---|---|
| 1 | `trust4_contigs` | Distinct contigs with a CDR3 call TRUST4 reported. | `locus=all`; no support denominator. |
| 2 | `trust4_cdr3_calls` | CDR3 calls TRUST4 reported. | `locus=all`; no support denominator. |
| 3 | `trust4_records` | TRUST4's AIRR records on the seven canonical loci and their source support. | TRUST4's locus. Matches `yield.tsv`. |
| 4 | `observed_cdr3` | TRUST4's records to the records whose CDR3 TRUST4 did not impute. | TRUST4's locus. |
| 5 | `igblast_input` | Parsed records to records routed to annotation. | `locus=all` counts every parsed record, including other TRUST4 loci. |
| 6 | `igblast_annotated` | Routed records to annotated records. | TRUST4's locus. |
| 7 | `locus_agreement` | Annotated records to records whose IgBLAST locus equals their TRUST4 one. | IgBLAST's locus; diagnostic. |
| 8 to 13 | `qc_*` criteria | Annotated records to records passing that criterion. | IgBLAST's locus; failures overlap. |
| 14 | `common_qc` | Annotated records to the common-QC population. | IgBLAST's locus. |
| 15 | `grouping_eligible` | Common-QC records to grouping-eligible records. | IgBLAST's locus. |
| 16 | `contig_collapse` | Eligible records to one representative per contig. | IGH, IGK or IGL; support conserved. |
| 17 | `exact_sequence_collapse` | Contig representatives to one member per exact full sequence within a library and locus. | Support conserved. |
| 18 | `receptor_units` | Members to the IGH clones observed in the library, or to clonotypes at the other loci. | The seven loci; support conserved. |
| 19 to 22 | `c_call_*` | IGH records to records meeting each qualification condition. | Call masking keeps source rows. |
| 23 | `isotype_callable` | IGH records to records whose constant-region call meets every qualification condition. | Call masking keeps source rows. |
| 24, 25 | `class_callable`, `subclass_callable` | IGH records to records resolved at each level. | Each level has its own denominator. |
| 26 | `shm_measured` | IGH records to records with a positive informative-position count and a finite regional mutation rate. | Each of the five V subregions, with no whole-V row; source support. Zero informative positions means no measured rate. |
| 27 | `shm_supported_positions` | Informative positions to positions meeting the vote minimum. | Each of the five V subregions; no support denominator. |
| 28 | `shm_eligible_members` | Clone members in the isotype group to members with a passing source for the region. | Isotype group and V region. |
| 29 | `selection_input` | Assigned members to members eligible for selection before clone consensus. | Isotype group; transferred member support. Grouping and contig losses appear upstream. |

An `all` row built from a locus partition totals the seven canonical loci. Only the scope collector adds totals spanning IGH and non-IGH grouping. Contig collapse has no TCR or `all` row, and `receptor_units` has no `all` row. A clone observed in two libraries counts once in each library's `receptor_units` row.

### Reading the counts

- **Failures overlap.** The common-QC criteria see the same input records, so their counts do not subtract from one another and do not add up to the records lost. Read `common_qc` for the population that passed.
- **Read status and missing counts separately.** A disabled numeric requirement can leave other conditions active, and those conditions can still exclude records. Blank counts are unavailable or inapplicable. Zero records an observed empty population.
- **Collapse conserves support.** Contig and exact-sequence collapse reduce rows and carry the read support over. Nothing is excluded. The AIRR export keeps every common-QC row with its own `consensus_count`, and `member_consensus_count` carries the transferred total on the member representative alone, so the two columns must not be summed together.
