# Interpreting the results

These results describe receptor sequences recovered from bulk RNA-seq. Read support depends on sequencing depth, assembly and transcript abundance. Cell counts, molecule counts and receptor-chain pairings are unavailable from this analysis. Limited recovery can leave a receptor undetected in a library.

Start with `metrics/headline.tsv` in the scope's `release/` folder. It has one row per library in the scope manifest. Read each measurement alongside its support count and availability status. The report displays the measurements from the released tables.

The headline covers receptor yields, IGH and TRB diversity and dominance, antibody-class fractions and somatic hypermutation (SHM). Sequencing depth and recovery affect these values. Rarefied diversity compares libraries at a common depth of recovered read support. The table for each topic contains further measurements, such as Shannon entropy and clonality.

For an unavailable estimate, check the corresponding status column. `diversity_status` in `metrics/diversity.tsv` explains missing diversity and dominance values, including top-5 fractions below a configured support minimum. `yield_status` in `metrics/yield.tsv` explains yield availability. QC fractions are in `metrics/qc.tsv`.

All table paths below are relative to the scope's `release/` folder unless stated otherwise.

Use the [glossary](glossary.md) for terms, the [QC matrix](qc_matrix.md) for requirements and the [data dictionary](data_dictionary.csv) for column definitions. Each release's `README.md` lists table paths and join keys.

## Find the table for a question

| Question | Tables | Section |
|---|---|---|
| How much receptor evidence did each library yield? | `metrics/headline.tsv`, `metrics/yield.tsv`, `metrics/qc.tsv`, `metrics/detail/qc_trail.tsv` | [Receptor recovery and QC](#receptor-recovery-and-qc) |
| How diverse is each repertoire? | `metrics/diversity.tsv`, `metrics/rarefaction_curves.tsv`, `metrics/detail/group_abundance.tsv` | [Diversity and dominance](#diversity-and-dominance) |
| What is the IgA, IgG and IgM composition? | `metrics/isotype_usage.tsv`, `metrics/detail/igh/clone_isotype.tsv`, and `class_iga_fraction`, `class_igg_fraction` and `log2_igg_iga` in the headline | [Antibody composition](#antibody-composition) |
| Which clones appear in two tissues or visits of one donor, and how much do they change? | `metrics/detail/shared_groups.tsv`, `metrics/donor_overlap_pairs.tsv`, `metrics/donor_overlap_summary.tsv`, `metrics/detail/group_abundance.tsv` | [Sharing between libraries](#sharing-between-libraries), including the endpoint order |
| Which V and J genes are used? | `metrics/gene_usage.tsv`, one row per library, locus, gene type and V or J entity, and `metrics/vj_pair_usage.tsv` for V and J pairs | [V and J gene usage](#v-and-j-gene-usage) |
| How mutated and selected are the antibodies? | `metrics/shm.tsv`, `metrics/selection_pressure.tsv`, `metrics/detail/igh/clone_shm.tsv`, `metrics/detail/igh/clone_member_shm.tsv`, `metrics/detail/igh/clone_member_shm_coverage.tsv` | [Somatic hypermutation](#somatic-hypermutation) and [Selection](#selection) |
| How do donor groups, such as responders and non-responders, compare? | Join your own metadata to these tables by `tissue_id` or `library_id` downstream | [Project metadata](#project-metadata) |

The unit being counted depends on the measurement: richness counts receptor groups, SHM eligibility counts members or clones, and SHM means average clone medians. Support-based measurements use recovered reads.

## Known limitations

The following limitations affect comparisons between libraries. Check the counts and status fields in each release to see where they apply.

- **Observed diversity depends on support and composition.** A library with more receptor reads need not have greater inverse Simpson diversity. Adding support to one dominant group can lower it. Shallow recovery can miss groups and inflate observed dominance, so read `top_5_support_fraction` and `top_10_support_fraction` alongside total support. Compare libraries on the rarefied values in `diversity.tsv` and read `rarefaction_status` alongside observed support. Low recovered support limits comparisons. No diversity minimum applies by default; `params.diversity.min_consensus_count` sets one.
- **Ambiguous V and J calls stay out of grouping.** A record whose V or J call names genes with different reference sequences is not grouping-eligible. Genes with near-identical paralogues, such as TRBV12-3 and TRBV12-4, lose more support this way, so their usage and the diversity and SHM of their clones rest on less support. Some J alleles differ only near the 5' end of the J-REGION and tie when the J is trimmed there, whatever the V coverage. The AIRR export keeps each excluded record with its `v_call`, `j_call` and `grouping_status`, so the loss can be measured per gene.
- **SHM eligibility rises with support.** Whole-V SHM uses a member only when one assembled representative covers all five V subregions with supported positions. Members with more reads are more often fully covered, so library SHM summaries describe the better-supported members. Compare libraries alongside `igh_shm_eligible_fraction_of_measured_members_whole_v` and the grouped-member coverage fractions in `qc.tsv`. The [SHM explanation](#somatic-hypermutation) also defines a downstream eligible-clone fraction.
- **The clone threshold depends on donor depth.** A donor's fitted threshold depends on how many nearest-neighbour distances it has, and the scope threshold is the median of the donors whose `fit_status` is `fitted`. `igh_clone_threshold/donor_fits.tsv` releases every donor fit, so each donor's own fit can be compared with the scope value.
- **Selection is exploratory.** BASELINe runs on partial sequences from bulk RNA-seq without a coverage gate, and `estimated` means only that a value exists. [Selection](#selection) gives the details.

## Units of observation

`library_id` identifies a sequencing library, `tissue_id` a sampled tissue and `donor_id` its donor. Several libraries can describe one tissue. Keeping those rows separate does not make them independent biological replicates. Supplied metadata labels the observations in the report, but it does not choose replicates or adjust for study covariates.

`consensus_count` is the recovered read support of a reported sequence or grouped receptor unit. Each table states its population. TRUST4-reported support, common-QC support and grouping-eligible support can differ, so they do not share a denominator.

TRUST4 itself can already combine contigs sharing its reported V, D, J, normalised C and CDR3 calls, and it assigns that support to a representative assembly. A record's support therefore need not equal the support of its single contig in `cdr3.out`. This happens before the pipeline's own collapses. The pinned implementation is [TRUST4 1.1.10 trust-simplerep.pl](https://github.com/liulab-dfci/TRUST4/blob/v1.1.10/trust-simplerep.pl).

## Clones and clonotypes

Grouping turns common-QC records into the units that diversity, gene usage and sharing count:

- An IGH clone is a lineage inferred within one donor and scope. Its members share V, J and junction length, and their junctions can differ through somatic hypermutation within one distance threshold.
- A single-chain clonotype, for the other six loci, groups retained members within one library and locus that share an exact nucleotide junction and the same resolved V and J entities. Its support follows the contig-representative transfer rule.

Before either grouping, duplicate records are reduced so that each eligible record's read support counts once. A member is one unique retained sequence within a library and locus. [Representatives before grouping](#representatives-before-grouping) describes the two reduction steps.

### IGH clones

IGH clone inference combines the members from one donor's libraries within a scope. The same members supply threshold fitting and clustering. Every donor in the scope uses one shared threshold, and the resulting clones are summarised within each contributing library.

Each scope uses only its own libraries for fitting and assignment. A library in two scopes can therefore have different IGH clone assignments in each. For a combined analysis, use the cohort release, where all libraries were analysed together.

Single-linkage clustering joins a member to a clone when its junction is within the threshold of at least one existing member. A chain of similar junctions can therefore have endpoints farther apart than the threshold. Every link must still share the same V and J germline entities and junction length.

`igh_clone_threshold/scope_threshold.tsv` records the threshold used. Its `threshold_source` identifies how it was chosen: a configured value, the median of usable donor fits or the fallback of 0.1 when no fit was usable. `n_donors` counts all donors in the scope, and `n_donors_fitted` counts those with a usable fit.

The donor evidence is in `igh_clone_threshold/donor_fits.tsv` and `igh_clone_threshold/donor_distance_histograms.tsv`, whichever threshold source applies. A fitted donor's row includes the mixture weights, shape, scale, mean and standard deviation. These values can be used to recompute component densities at each bin midpoint and compare them with `observed_density`. `n_distances_fitted` counts the distances used in fitting; `params.lineage.fit_subsample` can limit this below `n_distances`.

The fit uses a gamma component for near distances and a normal component for far distances. `params.lineage.fit_specificity` defaults to 0.995. SHazaM selects a cutoff targeting that fraction of the far component above the cutoff, after restricting the component to the fitted distances' minimum and maximum. This quantity calibrates the fitted model. The fraction of correctly separated biological clones remains unknown. The [pinned SHazaM fitting code](https://github.com/immcantation/shazam/blob/1.3.2/R/DistToNearest.R) defines that normalization and the cutoff search.

A donor fit needs at least `params.lineage.min_fit_distances` nearest-neighbour distances, 500 by default. A sequence has a distance only when the donor has another sequence with the same V and J germline entities and junction length but a different junction, so a donor can have many sequences and few distances.

`fit_status` records why a donor was not fitted: too few distances, distances that are all equal, a failed fit or a cutoff outside 0 to 1. `n_members` and `n_libraries` count the sequences and libraries that entered the fit. A donor whose fit is not used is still clustered at the scope threshold.

Clone IDs encode `scope_id~donor_id~cluster`, and scope and donor names cannot contain the separator. An ID has meaning within its inference run, scope and donor. It does not identify a clone across donors, or across scopes that resolved their own thresholds.

Members of one clone can carry different reconstructed germline alignments in the exported AIRR file. A method that needs one germline per clone has to select clones whose members agree, or reconstruct its own germline.

### Single-chain clonotypes

For IGK, IGL, TRA, TRB, TRD and TRG, a single-chain clonotype has an exact nucleotide junction and one assigned V and one assigned J germline entity, within one library. A clonotype can contain different full sequences sharing that key. An IGK or IGL clonotype cannot establish B-cell clonality without its heavy chain. The junction must match exactly, so somatic hypermutation in it can split one B-cell lineage into several light-chain clonotypes. Unrelated B cells with an identical light-chain junction and the same V and J entities count as one clonotype.

### Representatives before grouping

Non-imputed TRUST4 records keep their source IDs, support and assembly links through annotation and common QC. Grouping then reduces them in two steps:

1. For IGH, IGK and IGL, the dominant common-QC variant of one TRUST4 contig at one locus is chosen before grouping eligibility is tested. Ranking uses original record support, then an unambiguous call, V identity, usable alignment length, the junction, the full sequence and the V, J and C calls, and last the sequence ID. An ineligible dominant variant cannot be replaced by a minor eligible sibling. TCR records keep their individual eligible junction observations at this step.
2. Identical full DNA among the retained representatives merges into one member within each library and locus. Their analytical annotations must agree.

A member can therefore represent several original contigs. Each eligible source record contributes its support once through these two links. For light chains, an eligible sibling's junction can differ from the selected representative's junction. In an illustrative contig with eligible variant support 10 and 2, the dominant member can receive 12 support units. Those 12 units include support transferred from the sibling variant, whose junction can differ from the representative's.

Identical sequences in different libraries remain separate members with separate support. A contig representative is an analytical unit. Alternative junctions on its contig need not describe one biological rearrangement.

The pinned TRUST4 report script already filters TCR siblings at a default support fraction of 0.05 relative to the representative CDR3. The pipeline retains the surviving TCR variants.

TRUST4 can number the CDR3 variants of one contig differently in two runs of the same reads, so ties are broken on content and the sequence ID decides only between identical records. Released `sequence_id` and `clone_id` values still follow that numbering and SCOPer's labels. Compare runs using sequences and clone membership.

### Grouping in the AIRR export

The exported AIRR file keeps every common-QC row with its original support. `grouping_eligible`, `contig_id`, `contig_representative_sequence_id`, `member_representative_sequence_id` and the two representative flags, `is_contig_representative` and `is_member_representative`, record both grouping steps. Where no identical sequences merge within a library and locus, the contig and member links coincide. `member_consensus_count` carries the transferred support on member representatives and is zero on every other row, while `consensus_count` stays unchanged. Do not sum the two columns.

To rebuild the member tables from AIRR files, select rows where `is_member_representative` is T. Weight them by `member_consensus_count` and group them by `clone_id`. In this column, IGH records carry the scope's clone ID and other loci carry their clonotype ID.

For isotype composition, use every common-QC IGH row, including records outside grouping. The qualified constant-region call is `c_call` where `c_call_callable` is T. `isotype_class` and `isotype_subclass` give the resolved labels. These columns supply the isotype information in the export; there is no separate qualified-call column.

`n_contigs` counts original contributing TRUST4 contigs. `n_members` counts retained unique sequences within a library, for both IGH clones and other-chain clonotypes.

### Identifiers

| Identifier | Names | Unique within | Tables |
|---|---|---|---|
| `library_id` | One sequencing library. | The release. | Every per-library table. |
| `sequence_id` | One AIRR record, written `<library_id>~<TRUST4 record ID>`. | The release. | `clone_members.tsv`, `clone_member_shm.tsv`, `clonotype_members.tsv` and the AIRR export. |
| `contig_id` | One TRUST4 contig, written `<library_id>~<contig ID>`. | The release. | The AIRR export. |
| `clone_id` | One IGH clone, written `<scope_id>~<donor_id>~<number>`. | The scope. | `clones.tsv`, `clone_members.tsv`, `clone_member_shm.tsv`, `clone_isotype.tsv`, `clone_shm.tsv`, and AIRR `clone_id` on IGH member representatives. |
| `clonotype_id` | One clonotype in one library and locus. | The release. | `clonotypes.tsv`, `clonotype_members.tsv`, and AIRR `clone_id` on non-IGH member representatives. |
| `clonotype_key` | An exact junction with one V and one J entity, written `<junction>~<V entity>~<J entity>`. | Comparable across libraries annotated with the same reference. | `clonotypes.tsv`, and `clone_members.tsv` for each IGH member's own key. |
| `group_key` | The group a row counts: `clone_id` at IGH and `clonotype_key` at the other loci. | A library and locus in `group_abundance.tsv`, and a donor and locus in `shared_groups.tsv`, where it repeats on each of the donor's library rows. | `group_abundance.tsv` and `shared_groups.tsv`. |
| `gene_key` | One V or J germline entity: a gene, or a `VSEQ_` or `JSEQ_` composite. | The release. | `gene_usage.tsv`, and `v_gene_key` and `j_gene_key` in the other tables. |

`sequence_id`, `contig_id` and `clone_id` follow TRUST4's record numbering and SCOPer's labels, so they are specific to one run, as [representatives before grouping](#representatives-before-grouping) describes.

### Assess a threshold fit

`fit_status=fitted` means the donor met the computational requirements: at least `lineage.min_fit_distances` values, at least two distinct distances and a finite threshold in [0, 1]. Biological clone assignments still need assessment against the sequence evidence.

The default minimum of 500 distances is a pipeline policy. The tests check that it is applied; its suitability depends on the repertoire. Members sharing a junction can repeat distances, and identical full sequences in different libraries remain separate members. The distance count therefore includes dependent observations.

Inspect `igh_clone_threshold/donor_distance_histograms.tsv` with `donor_fits.tsv`, or the report's calibration plots. A fitted curve with separated components and a cutoff between their modes is easier to interpret than a boundary cutoff or strongly overlapping components. Neither appearance is proof of a correct clone partition. SHazaM recommends visual inspection and discusses roughly 1000 distances as empirical context for assessing fits. The pipeline retains its configured minimum. [SHazaM threshold guidance](https://shazam.readthedocs.io/en/stable/topics/findThreshold/) gives the method's qualifications.

The public fixture demonstrates fallback: no donor has enough distances, so `scope_threshold.tsv` records `automatic_fallback` and 0.1. The constructed 600-record threshold integration test demonstrates a successful fit and its handoff to assignment.

Automatic mode uses the median of usable donor fits; if none is usable, it uses 0.1. These choices prevent one donor's library count from directly weighting the scope threshold. They do not remove depth or sampling effects on the donor fits. A configured numeric threshold overrides the operative value while retaining calibration evidence.

To assess sensitivity, run plausible thresholds in separate results directories with the same inputs and eligibility settings. Compare member co-assignment, clone counts, support distributions and downstream summaries. Clone numbers are run-specific labels, so compare memberships by qualified `sequence_id`. Do not select a threshold solely because it strengthens a downstream association.

## Receptor recovery and QC

`yield.tsv` reports each library's input depth and TRUST4-reported support for the seven loci before pipeline QC. The trail stage `trust4_records` counts the same rows, so the two tables agree at each locus. For each locus:

```text
<locus>_yield_per_million_read_pairs = <locus>_trust4_consensus_count / input_read_pairs * 1,000,000
```

`input_read_pairs` counts each R1/R2 pair once. Positive depth with zero TRUST4-reported support gives zero yield. Zero depth gives an unavailable yield with `yield_status=no_input_read_pairs`.

Yield includes TRUST4-reported records with incomplete V(D)J sequences, imputed CDR3s and productive or non-productive annotations, at TRUST4's assigned locus. Partial CDR3 calls are excluded by TRUST4's AIRR conversion. A full CDR3 can still occur in an incomplete receptor sequence.

TRUST4 selects which assemblies and CDR3 calls to report. Yield therefore describes recovered receptor signal, with the sequencing depth and expression context needed for interpretation. It cannot count every receptor-containing fragment or distinct read.

The [QC matrix](qc_matrix.md) lists each requirement after that point, and `detail/qc_trail.tsv` records the populations each analysis used. The [reading rules](qc_matrix.md#reading-the-counts) explain overlapping criteria, missing counts and support transfer during grouping.

Productivity and stop-codon annotations describe the aligned V(D)J sequence. Unaligned terminal bases are outside the stop assessment, so passing common QC does not establish that the complete assembled transcript encodes a functional receptor. [AIRR defines the annotation fields](https://docs.airr-community.org/en/v2.0.0/datarep/rearrangements.html), and [IgBLAST uses local alignment that can omit changed terminal bases](https://www.ncbi.nlm.nih.gov/igblast/).

`qc.tsv` has one row per library and IgBLAST locus, including empty populations. It reports annotated, common-QC and grouping-eligible read support, with fractions showing how much support passes each stage.

The denominators depend on the fraction:

- Non-productive and common-QC fractions divide by annotated support.
- Grouping-eligible and grouping-reason fractions divide by common-QC support.
- IGH class- and subclass-callable fractions also divide by common-QC support.

The five grouping-reason fractions describe unresolved V calls, unresolved J calls, ineligible junctions, ineligible constant-region calls and contigs whose dominant variant fails grouping. A record can fail for several reasons, so the fractions can overlap. A fraction is unavailable at a zero denominator. Both contig and exact-sequence collapse conserve grouping-eligible support.

### Germline reference

The pipeline supports human receptor loci only. It curates `human_IMGT+C.fa` from the pinned Immcantation image's IMGT FASTAs, keeping the first record for each allele name. IgBLAST's databases and the MakeDb/CreateGermlines reference come from that same image. The reference audit checks their agreement before publication. `provenance.json` records the reference digest, IMGT download date and audit outcome.

IgBLAST's auxiliary data file supplies each J allele's coding frame and CDR3 end. Without an entry, a first-listed J hit has no proper frame or junction annotation, as [NCBI's setup guide](https://ncbi.github.io/igblast/cook/How-to-set-up.html) describes. MakeDb then rejects the record. Strict MakeDb also excludes records with inconsistent junction or alignment information. Raw annotation and failure evidence stay in each library's `airr/` directory and annotation log. These losses do not change the yield, which is measured before annotation.

A record can pass common QC yet lack a reconstructed germline. SHM and selection evidence is then unavailable. The pipeline retains the record and its original D call; it does not erase reverse-strand D calls to force reconstruction.

The pipeline supplements the image's NCBI file using two rules, in order:

1. An allele with an identical sequence to an allele that has an entry copies that entry.
2. Otherwise, a same-gene allele of equal length supplies the entry when its conserved codon and frame are preserved.

The build checks these rules against the entries already present. The pinned image gains entries for IGHJ5\*04, IGKJ4\*03, IGLJ2A\*01, TRAJ41\*02, TRAJ57\*02, TRDJ4\*02 and TRGJP2\*02. The three IMGT pseudogenes TRAJ51\*01, TRAJ55\*01 and TRAJ60\*01 have no qualifying source, so they remain without entries and records assigned to them are not annotated. The [annotation tests](../test/data/annotation/README.md) check constructed queries with the actual tools. The [reference receipt](maintenance.md#audit-the-germline-reference) records every copied entry and remaining gap.

## Diversity and dominance

`diversity.tsv` contains one row per library at each of the seven loci. Its `group_unit` column identifies IGH clones or [single-chain clonotypes](#single-chain-clonotypes).

The support used for these measurements is in `group_abundance.tsv`. Each row describes an observed group in one library, identified by `locus`, `group_unit` and `group_key`. Rows are ordered by `group_rank`. Tied groups share the smaller rank, so selecting the top k ranks can return more than k groups.

For a locus in `overlap.loci`, `shared_groups.tsv` lists the groups detected in at least two libraries from the same donor.

Let `p` be each group's fraction of the library's support at that locus. The observed measures are:

| Measure | Meaning |
|---|---|
| `richness` | Number of observed groups; Hill q = 0. |
| `shannon_entropy` | `-sum(p * log(p))`, using natural logarithms, in nats. |
| `hill_shannon` | `exp(shannon_entropy)`; Hill q = 1. |
| `inverse_simpson` | `1 / sum(p^2)`; Hill q = 2. |
| `clonality` | `1 - shannon_entropy / log(richness)`, defined only with at least two groups. |
| `top_1_support_fraction`, `top_5_support_fraction`, `top_10_support_fraction` | Fraction of support in up to the stated number of largest groups. |

Hill values are the effective number of equally supported groups. Higher orders give more weight to common groups, and none of them counts lymphocytes.

- A top-k fraction equals one when a nonempty repertoire has k or fewer groups.
- A repertoire with one observed group has Hill q = 1 and q = 2 equal to one, and unavailable clonality.
- An empty repertoire keeps zero support and richness. Its Shannon entropy, q = 1, q = 2, clonality and support fractions are unavailable.

The per-locus `params.diversity.min_consensus_count` defaults to zero. Below a positive minimum, Shannon entropy, Hill q = 1 and q = 2, dominance and target-depth values are unavailable. Observed `richness`, support and diagnostic curves stay available. Use `diversity_status` and `rarefaction_status` to find the reason for a missing estimate.

### Support rarefaction

The pipeline interpolates within the observed support vector without replacement. At depth `m`, it reports expected richness for q = 0, `exp(E[H])` for q = 1 and the reciprocal of expected Simpson concentration for q = 2. The q = 1 and q = 2 values are therefore not averages of diversity over repeated draws. The calculation is analytical and simulates no draws.

The production targets are 1000 [support units](glossary.md#assembly-and-annotation) for IGH and 400 for TRB. Nothing is extrapolated above observed support:

- A repertoire meeting its diversity minimum with total support below the target has `rarefaction_status=below_target_support` and no target value.
- An empty or below-minimum repertoire reports that reason instead.
- A locus with no target has `not_configured` and no curve.

`rarefaction_curves.tsv` holds q = 0, 1 and 2 on a bounded depth grid, including the full-depth endpoint and any reachable target. An empty repertoire has a depth-zero endpoint with zero richness and unavailable q = 1 and q = 2. Target values and curve points use the same calculation, and the endpoint equals the observed Hill values.

The curves describe subsampling of recovered support. They do not predict what TRUST4 would assemble from fewer FASTQ reads, and the pipeline does not downsample or reassemble FASTQs.

`inverse_simpson_asymptote` describes a model-based plateau for q = 2, calculated from corrected Simpson concentration. Read it alongside `asymptote_status`: all-singleton support has no finite asymptote. This is a diagnostic for the recovered support distribution. Sequencing-depth requirements need separate assessment.

For asymptote `A`, the reciprocal of expected Simpson concentration at support depth `m` is `m / (1 + (m - 1) / A)`. For `A > 1`, the continuous depth reaching 95% of `A` is `19 (A - 1)`, and the first positive integer depth at or above it is `max(1, ceiling(19 * (A - 1)))`. When `A = 1`, diversity is already one at every positive depth. Depth zero has no q = 2 value. These calculated depths are separate from the released curves, which stop at observed support.

Every value in this section can be recomputed from the support vector in `group_abundance.tsv`. The corresponding diversity row records the minimum and target. A depth-standardised comparison still inherits assembly, expression and eligibility differences between libraries.

## Antibody composition

`isotype_usage.tsv` reports subclass and class composition from qualified constant-region calls. Each level has its own callable-support denominator. A tie between IgG subclasses can resolve to IgG while remaining unavailable at subclass level.

IgM, IgD and IgE appear at both levels with the same support as a reporting convention; this does not assign them recognised human subclasses. Their two fractions differ only by denominator. Do not add subclass percentages and compare the sum with a class percentage.

Each row's `callable_consensus_count` is the support with a call at its `isotype_level`, the denominator of `support_fraction`. The headline reports the class fractions with `class_callable_consensus_count`, which separates missing isotype evidence from a low observed class fraction, and `log2_igg_iga`. `class_callable_fraction` and `subclass_callable_fraction` in `qc.tsv` give the callable share of common-QC IGH support. Composition does not need evaluable SHM or an assigned IGH clone.

A qualified C call can still span unresolved classes or subclasses. For example, discrimination between the named reference alleles IGHA1*01 and IGHA2*01, or IGHG1*01 and IGHG3*04, depends on which distinguishing bases were observed. Those allele pairs do not define universal subclass positions.

`clone_isotype.tsv` holds each clone's members and support per isotype within one library, for a downstream clone-weighted view. Its support sums grouping-eligible members, and `isotype_usage.tsv` sums every common-QC record's own support, so the two tables differ. Assessing treatment-associated class switching requires evidence beyond these support changes.

Class-switch recombination changes the heavy-chain constant region while retaining the rearranged V(D)J segment. Switching can occur directly from IgM or sequentially through another class or subclass. The human lineage study by [Horns et al. (2016), eLife](https://doi.org/10.7554/eLife.16578) describes direct and sequential pathways. This pipeline does not reconstruct switching events or their direction from its isotype fractions.

The log ratio is `log2(IgG support / IgA support)`. Both supports must be positive. The value is unavailable when class-callable support is zero or either the IgG or IgA fraction is zero. There is no pseudocount. Zero observed IgA support does not establish biological absence of IgA, and transcript-support composition is not a B-cell proportion.

## Somatic hypermutation

`shm.tsv` leads with all-IGH whole-V measurements. It also has rows for IgM, IgG, IgA, IgG1 to IgG4, IgA1 and IgA2, and for each of the five V subregions. Subclass rows need not sum to their class. IgD and IgE have no stratum of their own and stay in all-IGH. A missing measurement indicates insufficient evaluable evidence. Preserve it as missing when importing or summarising the table.

The [whole-V definition](glossary.md#isotype-and-somatic-hypermutation) includes FWR1, CDR1, FWR2, CDR2 and FWR3, excluding CDR3 and FWR4. Coverage uses the record's reconstructed germline occupancy. IMGT numbering gaps do not count, and legitimate V-end trimming and reference completeness affect that occupancy. It is not the length of a complete reference allele.

A member's coverage is assessed separately on every retained contig representative with its exact full DNA. For each region, the pipeline first applies the coverage rule to each source. It then chooses among passing sources by supported positions, original record support and stable sequence ID. Whole V needs one source to pass all five subregions and the occupancy floor. The pipeline never pools vote matrices or combines partial coverage from different sources. Identical sequences therefore contribute once while keeping their distinct assembly evidence.

`clone_member_shm.tsv` has one row per member and V region with a measured rate. It covers whole V and the five subregions. An absent row means no measured rate; within a library, isotype group and region, the row count equals `n_members_measured` in `shm.tsv`.

For rows where `shm_eligible` is FALSE, `n_mutations` and `shm_frequency` are empty. The remaining measurements supply the summaries in `shm.tsv` and `clone_shm.tsv`. Join to `clone_members.tsv` for each member's isotype to rebuild those tables.

`n_informative_positions` is the rate denominator. Coverage compares `n_vote_supported_positions` with `n_germline_occupied_positions`. Whole V also requires one source to pass all five subregions, so its summed counts alone are insufficient to determine eligibility.

`metrics/detail/igh/clone_member_shm_coverage.tsv` holds every grouped IGH member in each of six regions: whole V, FWR1, CDR1, FWR2, CDR2 and FWR3. `measurement_status` distinguishes `unmeasured`, `measured_ineligible` and `eligible` rows. A member-region with no measured rate has an `unmeasured` row; a missing member-region row is not its normal encoding.

Whole V overlaps the five subregions, so count one selected region or distinct members instead of summing region rows. The sparse `clone_member_shm.tsv` omits unmeasured regions.

The complete audit's `evidence_sequence_id` identifies the selected qualifying source for each member and region. It joins to the source AIRR record and is empty when no source qualifies. `germline_reconstruction_status` records reconstruction availability separately from coverage. The pinned Change-O parser cannot reconstruct affected reverse-D cases, whose SHM remains unavailable. `n_members_with_several_contig_sources` in `airr/validation.json` counts members with a choice of source.

Vote attribution can also be unavailable when an annotated record is opposite to its original contig orientation. The mapper leaves all V-position votes unavailable in that case, even if reads covered the region. This conservative mapping limitation is separate from missing reconstruction, an unmeasured rate and insufficient supported coverage. Changing the exported AIRR orientation flag does not change this internal vote mapping.

`junction_vote_status=opposite_contig_orientation` and `imgt_v_positions_mapped=0` identify the affected common-QC IGH source records in the originating run's `*_igh_consensus_votes.tsv` files. From the checkout, count them with `python examples/audit_vote_mapping.py results/airrborne`. The diagnostic counts source records. Another source with the same sequence can still supply eligible SHM evidence for the member. The portable release omits these vote tables, so its coverage status alone cannot distinguish every reason for failed attribution.

The AIRR alignments support recalculation of mutation counts with the pinned SHazaM method, including over partial coverage. Join the selected evidence ID when checking an eligible region's source. The position counts in `clone_member_shm.tsv` and `params.shm` reproduce each subregion's coverage comparison. The complete audit supplies the chosen source and status, but no position votes.

Rechecking those votes requires the unreleased source features in `libraries/<library_id>/metrics/<library_id>_igh_features.tsv` and `libraries/<library_id>/airr/<library_id>_igh_consensus_votes.tsv` in the results tree. Donor working tables under `scopes/<scope>/intermediate/<donor_id>/` also retain `shm_evidence_sequence_id_<region>` and unmasked measured values.

`mean_clone_median_shm` is the arithmetic mean of library-local clone medians, and `median_clone_median_shm` is their median. Each clone contributes the median SHM frequency of its eligible members in the stated region and isotype group. These medians receive equal weight, regardless of the clone's read support or member count. Read `n_clones_eligible`, `n_members_eligible` and `n_members_measured` beside it.

SHM eligibility fractions use different populations:

- The headline's sole SHM eligibility fraction, `igh_shm_eligible_fraction_of_measured_members_whole_v`, divides `n_members_eligible` by `n_members_measured` in the all-IGH whole-V row of `shm.tsv`. Unmeasured members are outside this denominator.
- In `qc.tsv`, `shm_eligible_fraction_of_grouped_members_whole_v` and `shm_eligible_fraction_of_grouped_members_any_subregion` divide by `n_grouped_igh_members`, including unmeasured members. Any-subregion counts a member once when at least one of the five subregions qualifies. Non-IGH rows have no applicable SHM fraction.
- A downstream eligible-clone fraction can be derived as `igh_n_shm_eligible_clones_whole_v / igh_richness`. Its numerator is `n_clones_eligible` in the all-IGH whole-V row of `shm.tsv`. The headline exports both counts, but does not export this ratio.
- A zero or unavailable denominator gives an unavailable fraction. Zero eligibility in a nonempty population gives zero.

`fraction_clones_zero_median_shm` divides by `n_clones_eligible`, the clones with at least one eligible member, as `mean_clone_median_shm` does. Its numerator counts clones whose members' median scored mutation count is zero. A clone with both zero and positive members can have a zero median, so this is not the fraction of fully unmutated clones.

`clone_shm.tsv` keeps library-local support and an observed state for each clone, isotype group and region with at least one evaluable member:

| `mutation_state` | Eligible observed members |
|---|---|
| `all_zero` | Every evaluable member has zero scored mutations. |
| `mixed` | Both zero and positive scored mutation counts occur. |
| `all_positive` | Every evaluable member has a positive scored mutation count. |

A clone with no evaluable member in a region has no row for it, and `group_abundance.tsv` keeps its support, so an unmeasured dominant clone remains visible. These states describe observed, evaluable members in one region. They do not establish exact germline identity or the mutation history of every cell in a biological clone. The library mean and zero-median fraction do not use the diversity minimum.

The SHM numerator is the sum of replacement and silent substitutions scored by SHazaM in germline codon context. It differs from a raw mismatch count or `1 - v_identity`. The pinned method excludes substitutions it cannot score, including a substitution that would create a stop codon when considered alone in the germline codon. It also excludes mutations in codons containing N or alignment gaps.

The denominator counts positions with a non-N, non-gap base in both sequences, including other IUPAC ambiguity codes. An incomplete codon can contribute positions to the denominator without a scored mutation. Zero scored mutations can therefore occur between sequences that differ. [SHazaM's mutation-count definition](https://shazam.readthedocs.io/en/stable/topics/calcObservedMutations/) describes these rules.

IMGT alignment removes insertions from this calculation. Deletions can reduce supported coverage and affect eligibility. Antibody SHM can include genuine insertions and deletions, so investigate an indel using its read support, frame and alignment context.

Germline reconstruction uses the first-listed reference allele before grouping. Within-gene ties remain possible. This choice is not donor genotyping, and an unrepresented inherited variant can be scored as a difference from germline. Missing reconstruction follows the [reference evidence policy](#germline-reference).

SHM depends on reference assignment, assembly fidelity and coverage. Compare libraries alongside their region-specific evidence and technical context. Changes in these measurements can suggest hypotheses about antibody maturation, but they cannot establish its mechanism or direction on their own.

## Selection

`selection_pressure.tsv` keeps separate BASELINe CDR and framework estimates, intervals and availability states. One row is one library within one isotype group: all IGH, IgM, IgG or IgA. The CDR estimate covers CDR1 and CDR2, and the framework estimate covers FWR1 to FWR3. CDR3 and FWR4 lie outside the estimated region.

Each clone contributes one effective sequence within each library and isotype group it appears in. Each member contributes once to that sequence, whatever its read support or number of contigs. `selection.collapse_min_freq` sets the base frequency the effective sequence needs. It defaults to 0.6, must lie in `(0, 1]` and has no disabled value. Zero is rejected because the pinned SHazaM method can then select unobserved bases at every position where member sequences differ. Positive frequencies follow the [SHazaM 1.3.2 collapseClones frequency comparison and character order](https://github.com/immcantation/shazam/blob/1.3.2/R/MutationProfiling.R).

The pinned `collapseClones` call builds a separate effective germline with `mostCommon` at each position. It does not choose one member's whole germline.

Both consensuses use `includeAmbiguous=FALSE` and `breakTiesStochastic=FALSE`. Ties follow the pinned character order A, T, G, C, N, dot and dash. No random base is drawn and no new IUPAC ambiguity code is introduced; unsupported input ambiguity codes stop this consensus calculation.

Each member's observed and germline alignments must have equal length. Across members, the consensus ends where no strict majority has sequence remaining, bounded by the region definition. Shorter members remain in the frequency denominator. These are [SHazaM 1.3.2 consensus rules](https://github.com/immcantation/shazam/blob/1.3.2/R/MutationProfiling.R), followed by the pipeline's observed-sequence masking described below.

BASELINe Sigma describes replacement-mutation enrichment or depletion relative to a context-dependent mutation model. It does not measure mutation burden, affinity or a treatment effect.

Before estimation, the pipeline masks gaps opposite germline nucleotides and incomplete codons in the effective sequence. The same masked sequence supplies observed mutations and expected opportunities. IMGT frame and shared numbering gaps stay aligned. This selection correction does not change the SHM endpoint.

Partial sequences remain eligible. Non-covering members enter the consensus denominator, so a base must meet the configured frequency across all members to survive collapse. Only usable positions in the effective sequence contribute to selection. No full-coverage gate or minimum usable-density cutoff is applied. Missing germline follows the [reference evidence policy](#germline-reference).

Treat the estimates as exploratory model-based evidence, read with their contributing clone counts and usable-density status. `sigma_cdr_status` and `sigma_fwr_status` are `estimated` when a finite Sigma with usable densities exists, which confirms that an estimate is available. Assess its reliability using the evidence below. A missing estimate is not evidence of neutral selection, and a CDR estimate does not substitute for a framework estimate.

Sparse evidence can produce a high posterior mean with a wide interval. Read each region's clone counts and 95% interval together. Comparisons require comparable mutation load and coverage, and a causal interpretation requires additional evidence. No universal neutral band or constant offset applies.

Of a row's `n_clone_consensus_sequences`, `n_clone_consensus_sequences_with_positions_<region>` have at least one covered position in the region, where the masked sequence and its germline both hold A, C, G or T. `n_clone_consensus_sequences_used_<region>` gave a selection density. A clone gives no density for a region when the focused test counts no mutation for it, that is no replacement or silent mutation in the region and no silent mutation in the other region, or when the region has no expected replacement mutation, as in a region the clone does not cover.

## V and J gene usage

`gene_usage.tsv` counts V and J entities at each locus's receptor unit: `group_unit` is clone at IGH and clonotype elsewhere. `n_groups` counts each unit once whatever its read support, and `consensus_count` sums the units' read support. An IGH clone's gene is the one entity its members vote for, and a clonotype's gene is part of its key, so each unit has one gene of each type. `group_fraction` and `support_fraction` each sum to one over each library, locus and gene type. A tied call or excluded assignment is not an extra gene observation.

The two fractions measure different things. `group_fraction` is the share of distinct clones that use an entity, and `support_fraction` is the share of recovered receptor transcript they carry. Support scales with each cell's receptor transcript level as well as with the number of cells. Antibody-secreting cells carry far more immunoglobulin transcript than other B cells, so at IGH `support_fraction` weights the genes of their clones more heavily.

Read usage by `gene_equivalence_group`. `gene_key` is a gene name or a composite entity, and the two levels can describe one gene:

- A call whose tied alleles all belong to one gene counts under that gene, however many alleles it names.
- A call whose tied alleles belong to two or more genes counts only when every named allele has the same complete reference sequence. Its entity is that sequence class, written as `VSEQ_` or `JSEQ_` and a SHA-256 digest of the sequence. `gene_names` lists its genes, and `identical_alleles` lists the alleles that share the sequence, so they can be checked against the reference.

A gene can appear under its own name and in a composite, depending on the observed calls. A complete call naming only that gene resolves at gene level, even when its reference sequence is also present under another gene. A complete multi-gene call resolves to a composite only when all named alleles share one complete sequence. Reference sharing alone does not determine which rows a library contains.

In the pinned human reference, IGKV1-39\*01 and IGKV1D-39\*01 share one sequence and IGKV1-39\*02 differs. A complete call naming the two shared alleles resolves to their composite; a call naming only IGKV1-39 alleles resolves to IGKV1-39. Likewise, IGLJ2\*01, IGLJ2A\*01 and IGLJ3\*01 share a sequence while IGLJ3\*02 differs. These identities support multi-gene resolution but do not prohibit single-gene entities.

`gene_equivalence_group` connects genes that share a complete reference sequence, directly or through another gene. Each analytical entity belongs to one group, and the group's rows are adjacent in `gene_usage.tsv`. This use of "group" is separate from IMGT gene groups such as IGHV.

Summing by `gene_equivalence_group` preserves total `n_groups` and `consensus_count` within a library, locus and gene type. The aggregated `group_fraction` and `support_fraction` each sum to one for a nonempty repertoire. Genes connected through different shared alleles may share no single sequence across the whole group.

These groups organise the usage table. Clone inference and sharing still use the individual entity in `gene_key`, and `gene_equivalence_group` appears only in `gene_usage.tsv`. Use `gene_key` for joins, because two composite entities can list the same genes.

Distinguishing genes within a group can be difficult. Some alleles differ from a shared allele by one or two nucleotides, and the pipeline does not infer donor genotypes. A somatic mutation or assembly error can therefore move a record between a gene's own row and its composite.

`vj_pair_usage.tsv` counts the same units by their V and J entities together, which the separate V and J rows of `gene_usage.tsv` cannot recover. A composite entity counts once under its own key. Summing its rows over one segment gives the other segment's rows in `gene_usage.tsv`.

### Why grouping needs one V and one J entity

Grouping requires one resolved V entity and one resolved J entity. This gives up sensitivity for completeness and keeps four steps from depending on an unresolved call:

- **Usage.** Each unit counts once, and a composite lists only genes that share its sequence.
- **Clone inference and sharing.** IGH clone inference compares junctions within one V entity, one J entity and one junction length, so a tied call cannot link two partitions. Clonotype keys, which sharing compares at the other loci, carry one entity of each type.
- **Record completeness.** A short V alignment distinguishes fewer alleles, so the requirement tends to remove records with less V coverage. It is not a coverage threshold: a short alignment resolves when one allele scores best.
- **SHM.** Grouping resolution does not select a donor-specific germline allele. A call naming several alleles of one gene resolves at gene level. Germline reconstruction already used the first-listed allele before grouping, as [Somatic hypermutation](#somatic-hypermutation) describes.

Some J alleles differ only near the 5' end of the J-REGION, inside the junction. A J trimmed there can tie two sequence classes and leave grouping whatever its V coverage. [Known limitations](#known-limitations) describes the support this costs.

The original `v_call` and `j_call` stay beside the resolved analytical entities. Single-gene calls resolve at gene level without erasing duplicated-location D markers. A multi-gene call resolves only when every named allele belongs to one exact sequence class derived from complete, unambiguous V-REGION or J-REGION sequences in the pinned snapshot. Partial entries, unknown bases and potentially truncated candidate lists cannot establish identity. The reference receipt records class membership.

A composite entity joins only itself. It cannot join a single member gene or a partly overlapping composite. Calls that remain ambiguous stay outside grouping. An identical observed junction under conflicting V/J entities is a diagnostic observation; it does not establish a common clone or justify attachment.

Reference functionality is allele-specific context, separate from observed rearrangement productivity. The pinned reference metadata records each allele's functionality and orphon status; the AIRR export leaves them out, since they follow from the calls. A functional allele elsewhere in a gene does not establish that the observed receptor functions. Leader, splicing and regulatory defects can lie outside the analysed sequence. Reference P/ORF entries are retained.

## Sharing between libraries

`donor_overlap_pairs.tsv` compares pairs of libraries from one donor within one scope, at each locus in `params.overlap.loci`. A donor with more than two libraries contributes every unordered pair. Morisita-Horn is the primary support-weighted similarity and emphasises dominant groups. Shared inferred IGH clones can contain different observed junctions through SHM; identical junctions or two members per endpoint are not required. The table also reports each endpoint's support in shared groups, and that support as a fraction of the endpoint's total. The two fractions can differ substantially when the libraries recovered different amounts of support.

`donor_overlap_summary.tsv` has one row per donor with at least two libraries in the scope, for each locus in `params.overlap.loci`. It carries the union and occupancy counts, the number of library pairs and of pairs with an available comparison (`n_pairs_total`, `n_pairs_available`), and two donor-level mean similarities. Each mean averages all of the donor's pairs at that locus. It is empty unless every in-scope library of the donor has nonempty support there, and it ignores the pair-support minimum.

`shared_groups.tsv` describes groups observed across libraries, and `group_abundance.tsv` keeps the complete vectors for pair plots. A group absent from one nonempty library's vector is displayed at zero. This records an undetected group; biological loss requires additional evidence.

Relative support is unavailable when a library has no retained support. Inspect support and diversity alongside similarity, because depth sensitivity remains possible even for an abundance-weighted index. Bray-Curtis on proportions and Morisita-Horn do not equalise detection depth or remove sampling uncertainty. Pairs that share a donor or an endpoint are not statistically independent.

### Redistribution and differences between endpoints

The pair table defines `p` at `library_id_1` and `q` at `library_id_2`. These endpoint identifiers are ordered lexically by `(tissue_id, library_id)` and do not establish visit order. With tissues named `P01_post` and `P01_pre`, endpoint 1 is the later visit. Every difference compares endpoint 2 with endpoint 1.

The report orients pairs by `timepoint_order` when it is supplied and differs between the two tissues; code that reads the table must orient them itself. Every proportion uses the full retained support at that endpoint and locus, including groups detected only there. IGH groups are donor and scope clones, and the other loci use single-chain clonotypes.

| Measurement | Interpretation |
|---|---|
| `bray_curtis_proportion` | `0.5 * sum(abs(q - p))` over the complete union, computed on within-library proportions, so it equals the total variation distance. It is a dissimilarity, where `jaccard` and `morisita_horn` are similarities. Zero means identical proportions and one means no groups in common. |
| `n_shared_higher_in_2`, `n_shared_higher_in_1`, `n_shared_equal` | Shared groups whose unrounded proportion is higher at endpoint 2, higher at endpoint 1 or equal. Every shared group has one vote. |

The report derives four further measurements from the pair table and `diversity.tsv`:

| Measurement | Interpretation |
|---|---|
| Simpson concentration difference | `1 / inverse_simpson` at endpoint 2 minus the same at endpoint 1, that is `sum(q^2) - sum(p^2)`. Positive means endpoint 2 is more concentrated than endpoint 1. Both endpoints must meet the diversity minimum. |
| Shared-group fractions | Each of the three counts above divided by `n_groups_shared`. |
| Shared-group balance | `(n_shared_higher_in_2 - n_shared_higher_in_1) / n_groups_shared`. Positive means more shared groups have a higher proportion at endpoint 2. |
| Support only at one endpoint | `1 - shared_support_fraction_k`: endpoint k support share in groups not detected at the other endpoint. |

Bray-Curtis on proportions is total variation distance. The sum of all positive differences equals the sum of all negative ones, because each endpoint sums to one. A value of 0.20 describes redistribution of 20% of recovered support across clone identities. It does not establish proliferation of 20% of cells. See the [SciPy Bray-Curtis definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.braycurtis.html).

The count-based fractions describe how widespread differences are among shared groups. Tiny differences count equally with large ones, and no effect-size tolerance or significance test applies. One dominant clone with a much higher share at endpoint 2 can give most other shared clones a lower share there. The balance describes relative differences among shared groups. Expansion or contraction of the whole cell population requires separate evidence.

Read `n_shared_higher_in_2`, `n_shared_higher_in_1`, `n_shared_equal` and `n_groups_shared` beside their fractions:

- With two nonempty repertoires and no shared groups, the counts are zero and the fractions and balance are unavailable.
- With an empty endpoint, every comparative measurement and difference count is unavailable.
- A group absent from a nonempty repertoire has observed zero support, while an empty repertoire has no proportion denominator.

Detection turnover does not establish biological emergence or extinction.

### Support and report inclusion

`params.overlap.min_consensus_count` sets a minimum both endpoints must meet at each locus. It defaults to zero, including for an omitted locus. The pair table records `min_consensus_count_per_library` and `pair_status`. Values below the minimum stay released and inspectable, and default report summaries exclude them. An `available` status means the configured support rule passed. Sampling uncertainty still needs assessment. The diversity minimum separately controls the report's Simpson concentration difference.

The minimum applies to the total support of each library in the pair. Individual shared clones can have lower support. Low-support sharing remains included under the configured rule. Low support makes proportions coarse and sensitive to sampling. A small sharing observation alone does not establish contamination or sampling noise. Total support also does not establish adequate support for every shared clone. Inspect both endpoint totals, their clone counts and the shared-clone count.

In the report's Sharing between libraries chapter:

- The smaller endpoint support is plotted against the selected pair score.
- The support filter can be changed for exploration. The active value is shown, and reset restores the configured minimum. Changing it does not change released values or flags, so record the configured rule used for publication separately.
- Open circles mark pairs below the configured minimum.
- Inclusion counts give candidate and plotted pairs and distinct donors for each displayed group.

For a study, inspect support distributions and score sensitivity before fixing publication minima, then report donor retention by group. Rarefaction targets are separate settings and do not establish reliable longitudinal cutoffs. These plots and measurements are descriptive. Donor-group tests remain a downstream analysis with an explicit repeated-measures design.

### Visit order

Optional `timepoint` and `timepoint_order` metadata supply visit labels and order. The report orients directional measurements by distinct visit orders. It swaps directional fractions and detection shares, and negates signed changes, when the stored endpoint order is reversed.

Missing visits do not exclude a donor's available pairs. Pairs without an order stay in the release and in the named A-to-B pair inspector, but have no grouped directional value. Response categories remain separate from visit order, and choosing an order does not establish that different specimens are anatomically comparable.

The pair inspector shows the complete clone union, the original denominators and the shared-group change counts. When "Abundance axes" is set to share of library read support, it also draws an equality line.

## The report

The report's chapters are receptor recovery, diversity and dominance, isotypes and mutation, and sharing between libraries, followed by methods and sources. Diagnostic panels expand in place, and navigation keeps grouping, filters, selections and scale choices. The running guide explains how to [open the report](running.md#open-the-report).

Groups and colours use supplied metadata, separately from timepoint panels. The labels come from the supplied metadata. Types are inferred from the source values, and numeric fields use continuous colour with fixed scope-wide limits and grey for missing values.

Distributions show one value per library, including when libraries are grouped by clinical subtype, and repeated libraries from one tissue or donor remain separate observations:

- A violin describes the plotted values when a group has at least five observations and three distinct values. Sparse and constant groups show points alone.
- A dark tick marks the median, calculated in original units.
- On log axes the violin and median use positive values only. Exact zeros are shown separately and remain selectable.
- Group labels count plotted values and represented donors.
- The report shows eligible and total library, tissue and donor counts, keeps missing measurements and applies no report-level minimum.

Selections link across chapters:

- The donor search and selectors work in every chapter, including reports without annotations.
- Selecting a pair highlights its two endpoint libraries. "Select all libraries from these donors" extends the selection across visits and tissues, and pairs then follow their selected endpoints.
- Visibility filters stay active, and the report counts selected libraries outside them.
- Clicking the sole selected point again, double-clicking the plot, pressing Escape or choosing Clear selection deselects without resetting the display. Clearing also turns off "Selection only" and closes the inspected pair.

## The release

The README [results section](../README.md#read-the-results) identifies the entry files. The release README lists what each row represents, the join keys, loci and inclusion rules, including threshold file names and AIRR and member joins. Read large detail tables in chunks and keep identifiers as strings. The release's `data_dictionary.csv` defines each released column's unit, denominator and missing-value rule.

### Join tables

Join library results on `scope_id, library_id`. Join IGH abundance to donor clones on `scope_id, group_key=clone_id`, and non-IGH abundance to clonotypes on `scope_id, library_id, locus, group_key=clonotype_key`. Join member-region rows to members on `scope_id, sequence_id`. Include a release identifier when joining across frozen outputs; a scope name alone does not identify a revision.

`representative_sequence_id` in the donor clone table identifies the highest-support observed member, also used for sharing descriptors. Join to AIRR using its `representative_library_id`, which can differ from a sharing row's library. Inferred ancestors and clone consensuses are absent from this table.

AIRR source `consensus_count` and transferred `member_consensus_count` describe different populations. Sum member support once per member. A source-to-member join can repeat that total on several rows, so deduplicate members before summing it.

### AIRR files

The release holds one compressed AIRR rearrangement file per library under `airr/`. `sequence_id` is qualified as `<library_id>~<source_sequence_id>`, and `source_sequence_id` keeps the original library-local ID. `contig_id` and the contig and member representative links use the same qualified IDs. Each library's export is shared between scopes, and the release copy adds this scope's `clone_id`.

At every locus `clone_id` sits only on member-representative records. It holds the scope's IGH clone at IGH and the `clonotype_id` at the other loci. Contig siblings and records outside grouping keep an empty `clone_id`, because a sibling's junction can differ from its member's. The release README gives the join keys to `clone_members.tsv` and `clonotype_members.tsv`.

Released V and J CIGARs and supported D CIGARs use AIRR query and reference semantics: I consumes query positions and D consumes reference positions. Publication checks their query and reference intervals, selected alleles and bases against the stored sequence and pinned reference. `cigar_query_spans_valid` and `cigar_reference_spans_and_bases_valid` in `airr/validation.json` record those semantic checks. They are separate from AIRR structural validation and metadata completeness.

`v_germline_start` and `v_germline_end` are one-based inclusive positions on the selected ungapped reference allele, derived from `v_cigar`. The native padded IMGT positions remain in `changeo_v_germline_start` and `changeo_v_germline_end`.

Standard D CIGARs and germline endpoints are empty where the selected hit cannot be represented, with `d_alignment_status` explaining availability. Reverse-D hits retain their selected native evidence and checked query coordinates. Validation of a short D match does not establish a unique biological D assignment or restore an unavailable germline. Raw MakeDb rows stay in each library's `airr/*_makedb_pass.tsv.gz` outside the release.

**Restore the native Change-O coordinate view.** Run the converter from the pipeline checkout before reusing exported coordinates with the pinned Change-O tools:

```bash
python workflow/scripts/restore_changeo_coordinates.py \
  --input rearrangements.tsv.gz --output changeo_input.tsv \
  --unavailable reconstruction_unavailable.tsv
```

The converter restores native V/D coordinates and the native D CIGAR only for rows with `germline_reconstruction_status=available`. It writes every other row separately and leaves sequences and alignments unchanged. `CreateGermlines.py` requires the matching pinned, IMGT-gapped reference. A lineage analysis needing one germline per clone still requires an explicit clone-level reconstruction step.

`fwr1`, `cdr1`, `fwr2`, `cdr2`, `fwr3`, `cdr3` and `fwr4` divide `sequence_alignment` into its IMGT regions. MakeDb writes them from the gapped alignment, so they keep its gap (`.`) and deletion (`-`) characters and can contain N. A `.` marks an IMGT numbering gap or a position before the first sequenced base.

A V-domain region with no sequenced base is empty, as the FWR1 of a record whose sequence starts in FWR3 is. Joining the regions in order therefore spells `sequence_alignment` without the gap characters of the empty regions. Remove the gap characters before translating a region or measuring its length.

`sequence` is stored in V orientation, and standard `rev_comp` is F. Query coordinates, CIGARs and alignments index this stored sequence directly. `changeo_rev_comp` preserves MakeDb's historical flag, which is T when it reverse-complemented the original input query. Do not reverse the exported sequence because of this custom field. This convention applies from output schema 2.0.0; schema 1.0.0 retained the historical flag in `rev_comp` and requires this correction before standards-aware reuse.

`c_call` names IGHG4 reference sequences as `IGHG4A*01` and `IGHG4A*02`, as the IMGT reference does. Isotype analyses count them as IgG4. See the [AIRR 2.0.0 CIGAR specification](https://docs.airr-community.org/en/v2.0.0/datarep/rearrangements.html#cigar-specification).

`isotype_class` and `isotype_subclass` name the isotype that a qualified constant-region call resolves to, as the isotype tables and `clone_members.tsv` count it. Each is empty where the call fails qualification or does not resolve at its level, and on every non-IGH record.

`airr/Repertoire.json` links the scope's libraries, tissues, donors, processing identities and files. Its primary `data_processing_id` identifies shared library annotation. Clone assignments are added for the enclosing `scope_id`; use that field and the scope provenance to identify clone inference.

The official AIRR 2.0.0 schema is pinned separately from the AIRR packages in the analysis image. `airr/validation.json` separates structural validation from `strict_metadata_gaps`. The strict validator also checks omitted optional objects declared nonnullable by the schema, including sequencing files and PCR targets. These gaps do not mean the pipeline accepts those metadata as inputs.

`semantic_checks` in the same file states the V germline coordinate convention and counts, per locus, member representatives, members with several contig sources, members whose support exceeds their source record, blank germlines and noncanonical junction peptides.

Known bulk RNA and unpaired context is filled in. Required unavailable preparation, clinical and study values stay null, with the reason stated once in `Info.description`. Optional unavailable fields, including sequencing files and repertoire type, are omitted. FASTQ read counts do not establish reads passing an upstream QC filter, and sequence completeness does not establish the library preparation protocol. Neither is inferred. The release does not claim complete MiAIRR reporting.

The AIRR `quality_thresholds` field describes base-quality filtering. `data_processing_protocols` describes annotation QC, while `collapsing_method` states that exported common-QC rows remain separate. Primer processing and unavailable assembly history have explicit descriptions. Project annotations and timepoints remain exclusive to the report.

### Provenance

`provenance.json` records the pipeline source, image identities, measured tool versions, reference digests, IMGT/GENE-DB download date and reference audit outcome. Detailed reference sequence classes are in `shared/provenance/imgt_reference.json` in the original results tree.

Assembly producer records preserve the settings used to create each assembly. These can differ from the analysis configuration when assemblies are reused. Without a producer record, the assembly's history remains unknown.

`effective_config.yaml` records the scientific settings and seed. Its `fixed_methods` section records choices fixed in code, such as clone linkage and the selection model. Paths under `/usr/local/share` refer to files inside the named analysis image. The generated methods text combines these settings and provenance with the scope manifest and chosen threshold.

### Project metadata

Project metadata appears only in the report, under its original column names. To analyse annotations with the scientific tables, join the original metadata to `metrics/headline.tsv` on `tissue_id` or `library_id`:

- Check that the metadata key is unique. A tissue key can annotate several library rows.
- Keep every headline row, and leave unmatched annotations missing.
- The report's internal `meta_` prefix does not rename columns in the source file.
