# Glossary

This glossary defines the concepts used across the pipeline's tables, report and guides. The [data dictionary](data_dictionary.csv) defines each released column with its unit and missing-value rule. The [QC matrix](qc_matrix.md) lists each requirement a record meets, and the [interpretation guide](interpretation.md) explains how to compare the measurements.

## Samples and scopes

| Term | Meaning |
|---|---|
| Library | One sequenced preparation, identified by `library_id`. Assembly, annotation and most measurements are calculated separately for each library. |
| Tissue | One specimen, such as a biopsy or blood sample, identified by `tissue_id`. A donor's tumour and blood samples have separate IDs. Several libraries can come from one specimen and remain separate observations. Tissue type is a separate metadata field. |
| Donor | The person a tissue came from, identified by `donor_id`. IGH clones are inferred within one donor. |
| Scope | One analysed and separately released set of libraries, identified by `scope_id`. Each project in the membership CSV is a scope. The reserved `cohort` scope holds every samplesheet library. |
| Scope manifest | `manifest/scope_manifest.tsv`, the authoritative list of a scope's libraries with their tissues and donors. |
| Identifier | Library, tissue, donor and scope IDs are strings, so `0460` and a literal `NA` remain labels. Released numeric columns use an empty cell for an unavailable value. |
| Timepoint | An optional supplied visit label (`timepoint`) and numeric visit order (`timepoint_order`) for a tissue. The order does not imply an elapsed time. |
| Report annotation | A tissue or library metadata column selected for one scope's report. The report stores it as an internal `meta_` field labelled with the original header. Scientific tables and AIRR files do not carry it. |

## Assembly and annotation

| Term | Meaning |
|---|---|
| TRUST4 | The assembler that finds receptor reads in bulk RNA-seq and builds receptor sequences from them. It is the first step of the pipeline. |
| IMGT | The international ImMunoGeneTics information system. Its germline gene database supplies the reference sequences, gene names and region numbering the pipeline uses. |
| IgBLAST | The aligner that names the V, D, J and constant genes of each assembled sequence against the human IMGT germline reference. Change-O, SHazaM and SCOPer from the Immcantation suite then annotate, measure mutation and group clones. |
| Locus | One of IGH, IGK, IGL, TRA, TRB, TRD and TRG. |
| Chain | The receptor polypeptide one locus encodes. IGH is the B-cell heavy chain, IGK and IGL are its light chains, and the TR loci encode T-cell receptor chains. Bulk RNA-seq does not pair chains, so one record describes one chain. |
| BCR | B-cell receptor. The pipeline analyses its IGH, IGK and IGL chains separately. |
| TCR | T-cell receptor. The pipeline analyses its TRA, TRB, TRD and TRG chains separately. |
| V, D and J | Variable, diversity and joining gene segments. Rearrangement combines these segments to form a receptor sequence. Some loci use V and J without D. |
| CDR | Complementarity-determining region. The region labels are CDR1, CDR2 and CDR3. |
| CDR3 | The third complementarity-determining region, spanning the V-J or V-D-J join according to the locus. The AIRR junction is CDR3 with its two conserved flanking codons. |
| FWR | Framework region. The annotated regions are FWR1 to FWR4. |
| Contig | A nucleotide sequence assembled from reads by TRUST4, identified by `contig_id` in the AIRR export. One contig can have several alternative CDR3 calls. Contig counts describe assembled sequences. |
| Assembly variant | One of the alternative CDR3 calls TRUST4 reported for one contig. |
| Imputed CDR3 | A CDR3 inferred by TRUST4, marked by the score 0.01 in `cdr3.out`. The pipeline removes these records before annotation. |
| Record | One row of AIRR rearrangement data: one assembled sequence with its annotation. `sequence_id` identifies it within its library. Pooled and exported tables qualify it as `<library_id>~<source_sequence_id>`. |
| Rearrangement | The annotated V(D)J sequence a record describes, in the sense of the AIRR Rearrangement schema. |
| Junction | The annotated nucleotide junction (`junction`). Non-IGH clonotypes require an exact match. IGH clone inference compares junctions by distance. |
| Read support | The recovered read support recorded in `consensus_count`. It depends on sequencing depth, assembly and transcript abundance. Each table defines the records or groups whose support it counts. Cell and molecule counts are unavailable from this analysis. |
| Support unit | One unit of read support. Rarefaction targets, diversity minima and pair-support minima are counts of support units, so a target of 1000 means a summed `consensus_count` of 1000. |
| Input read pairs | Paired sequencing records in the input FASTQs (`input_read_pairs`). The two mates of a pair form one unit. |
| TRUST4-reported support | The read support of TRUST4's full-CDR3 AIRR records at the locus TRUST4 assigned, before pipeline QC, released as `trust4_consensus_count` and its per-locus columns. |
| Receptor yield | TRUST4-reported support at one locus per million input read pairs. |
| Productivity | The IgBLAST and MakeDb prediction (`productive`) for the aligned V(D)J sequence. Unaligned terminal bases are outside its stop-codon assessment, so a true value does not establish that the whole assembled transcript is functional. |
| Original gene calls | `v_call` and `j_call` as annotated, which can list tied genes or alleles. The pipeline keeps them beside its normalised `v_gene_key` and `j_gene_key`. |
| Germline entity | The analytical V or J assignment. A call naming alleles of one gene resolves to that gene, however many alleles it names. A call naming alleles of several genes resolves to their complete-reference sequence class only when every named allele belongs to it. `v_gene_key`, `j_gene_key` and `gene_key` hold it. |
| Composite entity | A germline entity for alleles of two or more genes that share one complete reference sequence, so the data cannot say which gene was used. It appears as `VSEQ_` or `JSEQ_` and a SHA-256 digest of the sequence. A `_names` column beside it lists its genes, and `identical_alleles` in `gene_usage.tsv` lists its alleles. It joins only itself and is counted apart from its member genes. |
| Gene equivalence group | Genes connected by shared complete reference sequences, directly or through another gene. `gene_equivalence_group` in `gene_usage.tsv` lists their sorted names. Each gene and composite entity belongs to one group; a gene sharing no sequence forms its own group. This classification concerns reference sequence identity. Genomic linkage is outside its definition. |
| Assigned group | The `assigned_group` value of `v_gene_status` or `j_gene_status`. Every named allele in a multi-gene call shares one complete, unambiguous V-REGION or J-REGION reference sequence, so the assignment uses a composite entity. Partial entries and truncated candidate lists cannot establish this identity. Genomic and functional equivalence require separate evidence. |
| Qualified constant-region call | The IGH constant-region call that passes isotype qualification: `c_call` on a record whose `c_call_callable` is T. `isotype_class` and `isotype_subclass` name its isotype. Pipeline tables carry it as `qualified_c_call`, which is empty when the call is masked. |

## QC and eligibility

| Term | Meaning |
|---|---|
| Common QC | The checks applied to annotated records at every locus. A record must pass all six criteria in [Annotation and common QC](qc_matrix.md#annotation-and-common-qc). |
| Exclusion | A record is dropped and takes no further part. |
| Eligibility | Whether a record has the evidence required for a particular analysis. An ineligible record remains available to other analyses whose requirements it meets. Grouping, isotype, SHM and selection each have their own requirements. |
| Call masking | Blanking one field, such as the qualified constant-region call, while the record stays in every population it belonged to. |
| Grouping | Combining records into a shared unit such as a clone or a clonotype. The records and their support remain, and the counted unit changes. |
| Grouping eligibility | Whether a common-QC record meets its locus's grouping requirements (`grouping_eligible`). |
| Grouping-eligible fraction | `grouping_eligible_fraction_of_common_qc` in `qc.tsv`: grouping-eligible read support divided by common-QC read support, for one library and IgBLAST locus. It is unavailable when common-QC support is zero. |
| Disabled requirement | A numeric requirement set to its disabling value. Any other condition of the same requirement still applies. |
| Unavailable value | An empty measurement whose status or denominator explains why it is unavailable. Preserve it as missing when importing or summarising tables. |
| Input-contract error | A malformed input or contradictory annotation. It stops the run and is not counted as an exclusion. |

## Audit units

The `input_unit` and `output_unit` fields of `metrics/detail/qc_trail.tsv` use these tokens. Counts describe the stated population; overlapping rows are not sequential losses.

| Token | Counted unit |
|---|---|
| `contig` | One TRUST4 contig. The grouping stages count it by the record that represents it: its dominant common-QC variant at IGH, IGK and IGL, and each eligible variant at the TR loci, which keep them all. |
| `cdr3_call` | One CDR3 call TRUST4 reported. |
| `airr_record` | One AIRR rearrangement record. |
| `member` | One distinct full nucleotide sequence within a library and locus, as defined under [Grouping and clones](#grouping-and-clones). |
| `clone` | One IGH clone observed in one library. |
| `clonotype` | One non-IGH single-chain clonotype. |
| `position` | One aligned position in one record. |

## Grouping and clones

| Term | Meaning |
|---|---|
| Contig representative | The one record that stands for a contig in grouping. The AIRR export links each eligible record to it in `contig_representative_sequence_id` and flags it with `is_contig_representative`. At IGH, IGK and IGL it is the dominant common-QC variant of the contig at that locus, chosen by original record support and fixed tie-breaks before grouping eligibility. An ineligible dominant variant is not replaced by an eligible sibling. Every eligible TCR record represents itself. Siblings need not be one biological rearrangement. |
| Exact sequence collapse | Merging identical full nucleotide sequences among contig representatives, within one library and locus. It keeps each source record and its links. |
| Member | One distinct full nucleotide sequence within a library and locus, and the unit IGH clone inference clusters. It is what remains after contig representatives are chosen and exact sequences are collapsed, so one member can represent several contigs with the same sequence. `member_representative_sequence_id` in the AIRR export names its representative record. Internal grouping tables call it the inference representative. |
| Member support | `member_consensus_count` in the AIRR export, the eligible source support transferred once to a member representative. Every other row carries zero and keeps its own support in `consensus_count`. The members tables carry it as `member_consensus_count`. |
| IGH clone | An inferred heavy-chain clonal lineage within one donor and scope. `clone_id` encodes it as `scope_id~donor_id~cluster`, and the label has meaning only within that inference run. |
| Single-chain clonotype | A non-IGH group of retained members sharing an exact nucleotide junction and the same resolved V and J entities within one library and locus. Support follows contig-representative transfer, so an eligible light-chain sibling can contribute even when its own junction differs. It does not pair chains. |
| Receptor unit | The grouped unit a support vector or sharing table counts: an IGH clone at IGH and a single-chain clonotype at the other loci. `group_unit` names the kind and `group_key` identifies the unit. Read both with scope, donor and locus. |
| Clone threshold | The one junction-distance cutoff at which every donor in a scope is clustered, as a Hamming distance normalised by junction length. |
| Threshold source | How the operative clone threshold was obtained (`threshold_source`): `configured`, `automatic_fit` (the median of the usable donor fits) or `automatic_fallback` (0.1, when no donor fit is usable). |
| Donor fit | One donor's own fitted cutoff and distance histogram. Fits run and are released whichever threshold source applies. A donor that cannot be fitted is still clustered. |
| Shared group | A receptor unit detected in both libraries of a pair. `shared_groups.tsv` lists a unit detected in at least two of a donor's libraries. Limited recovery can leave the group undetected in other libraries. |

## Diversity

| Term | Meaning |
|---|---|
| Support vector | The read support of each receptor unit in one library and locus, released in `group_abundance.tsv`. Diversity and pair plots read it. |
| Repertoire | A collection of receptor sequences. The diversity tables describe one library at one locus. |
| Observed diversity | Diversity calculated from the full support vector. |
| Hill diversity | The effective number of equally supported groups. Orders q = 0, 1 and 2 give progressively more weight to common groups, and are released as `richness`, `hill_shannon` and `inverse_simpson`. |
| Support rarefaction | Analytical interpolation without replacement within recovered support. It does not estimate what TRUST4 would recover from fewer reads. |
| Rarefaction target | The support depth at which diversity is compared across libraries. The value is unavailable above observed support. |
| Diversity minimum | The per-locus `diversity.min_consensus_count`. Below it, Shannon entropy, Hill q = 1 and q = 2, dominance and target-depth values are unavailable. Observed richness, read support and diagnostic curves remain available. |
| Top-k support fraction | Support in the k largest groups, or in all groups when there are fewer than k, divided by all support at that locus. |
| Clonality | One minus Shannon evenness. It is unavailable with fewer than two observed groups. |
| q = 2 asymptote | A model-based plateau `inverse_simpson_asymptote` from corrected Simpson concentration. It is separate from the observed curve, is not gated by the diversity minimum and has its own `asymptote_status`. |

## Isotype and somatic hypermutation

| Term | Meaning |
|---|---|
| Isotype class | The major IGH class: IgM, IgD, IgG, IgA or IgE. |
| Isotype subclass | A finer constant-region category, such as IgG1 or IgA2. Its callable population can differ from the class's. |
| Callable support | Support whose qualified call resolves to a category at the stated level, class or subclass. |
| Class switching | Recombination that changes the antibody heavy-chain constant region while retaining the rearranged V(D)J segment. It can occur directly from IgM or sequentially between switched classes or subclasses. The bulk RNA-seq results describe constant-region calls. Reconstructing switching events requires additional analysis. See [antibody composition](interpretation.md#antibody-composition). |
| IGHG4A | A duplicated copy of IGHG4 that some haplotypes carry. The pinned IMGT reference names it as a separate gene, with alleles `IGHG4A*01` and `IGHG4A*02`. It encodes IgG4, and isotype analyses count it as IgG4. |
| SHM | Somatic hypermutation, scored against the assigned D-masked germline over an evaluable region. |
| Whole V | The five IMGT V subregions FWR1, CDR1, FWR2, CDR2 and FWR3, excluding CDR3 and FWR4. |
| Germline-occupied position | A position occupied by a base that the pinned SHazaM method counts in the record's reconstructed germline, including IUPAC ambiguity codes other than N. N and numbering gaps are excluded. Reconstruction reflects rearrangement boundaries and reference completeness; occupancy is the SHM coverage denominator. |
| Informative position | An IMGT V position with a counted base in both the record and its germline. SHazaM includes IUPAC ambiguity codes other than N and excludes N, dot and dash. These positions form the mutation-frequency denominator; an incomplete codon can still prevent scoring a mutation at a counted position. |
| Supported position | An informative position whose consensus base has at least `shm.min_position_votes` agreeing votes. A vote is not a read. |
| SHM evidence source | The retained contig representative that supplies a member's SHM evidence for one region (`evidence_sequence_id` in `metrics/detail/igh/clone_member_shm_coverage.tsv`, also retained as `shm_evidence_sequence_id_<region>` in donor member tables under `scopes/<scope>/intermediate/`). Each representative with the member's exact sequence is tested alone, and votes from different sources are never pooled. With no passing source, the ID is empty and the coverage audit's `shm_eligible` is false. |
| SHM eligibility fractions | The headline divides whole-V eligible members by measured members. QC divides whole-V or any-subregion eligible members by all grouped IGH members. An eligible-clone fraction can be derived from whole-V eligible clones and observed library-local IGH richness. These fractions do not count cells. |
| Clone SHM mean | `mean_clone_median_shm`, the arithmetic mean of library-local clone medians, with each clone weighted equally. |
| Zero-median clone fraction | The fraction of clones with an eligible member whose members' median scored mutation count is zero. A clone with both zero and positive members can have a zero median. |
| Mutation state | `mutation_state`, the observed state of a clone's eligible members in one region and isotype group: `all_zero`, `mixed` or `all_positive`. |
| BASELINe | The SHazaM model-based estimate of selection on mutations, reported separately for CDR1 with CDR2 and for FWR1 to FWR3. CDR3 and FWR4 lie outside the estimated region. |
| Effective clone sequence | One collapsed sequence per clone within one library and isotype group. Each member enters the consensus denominator, including members that do not cover a position. BASELINe convolves the usable densities from these sequences. |

## Sharing between libraries

| Term | Meaning |
|---|---|
| Jaccard | Groups detected in both libraries divided by groups detected in either. Read support does not enter. |
| Morisita-Horn | A support-weighted similarity between two libraries that emphasises common groups. |
| Bray-Curtis redistribution | Half the summed absolute difference in support proportion over the complete group union. It equals total variation distance. |
| Simpson concentration difference | The sum of squared support proportions at endpoint 2 minus the same sum at endpoint 1, which the report derives from each endpoint's `inverse_simpson`. A positive value means endpoint 2 is more dominated by its largest groups. |
| Shared-group difference fraction | The number of shared groups whose proportion is higher at endpoint 2, higher at endpoint 1 or equal, divided by all shared groups. Every difference counts once, however small. |
| Shared-group difference balance | Shared groups higher at endpoint 2 minus shared groups higher at endpoint 1, divided by all shared groups. |
| Support share only in one endpoint | The support share at one endpoint in groups not detected at the other. |
| Pair-support minimum | The per-locus `overlap.min_consensus_count` both endpoints must meet for default report inclusion. Pairs below it stay released and are flagged in `pair_status`. |

## Software and execution

| Term | Meaning |
|---|---|
| Snakemake | The software used to organise the workflow. It determines which jobs are needed from their inputs and outputs, then runs them locally or submits them to a scheduler. |
| `deps`, `outs`, `params` | The three parts of the configuration: inputs and images, output directories, and analysis settings. |
| Rule and job | A rule is one step of the pipeline, such as assembly or annotation. A job is that step applied to one library, donor or scope. |
| Dry run | `run_module.sh -n`, which lists the planned jobs without running them. |
| Container image | Packaged software run by a compatible container runtime. The host still supplies the kernel and a compatible CPU architecture or emulation. The pipeline uses one image for TRUST4 and one for IgBLAST, R and the Immcantation packages. See [platform evidence](running.md#platform-and-runtime-evidence). |
| SIF | The single-file form of a container image that Apptainer and Singularity run. Setup pulls one SIF per image into `deps.resources`. |
| Image digest | The `sha256:` checksum, or OCI digest, that names one exact build of a container image. The configuration pins each image by digest, so every pull returns the same software. |
| Identity receipt | A JSON file beside each SIF in `deps.resources`. It records the configured image, its digest, the SIF's checksum and size, and the runtime that pulled it. A launch that finds the receipt reuses the SIF, and release provenance reads the image identity from it. |
| Container runtime | The program that runs images: Apptainer, Singularity or Docker. `AIRRBORNE_CONTAINER_RUNTIME` names one explicitly. |
| Bind path | A host directory made visible inside a container. The launch directory is always bound; `AIRRBORNE_CONTAINER_BINDS` adds others. |
| Conda environment | A directory of installed software. Setup creates the `snakemake_9.27.1` environment for the launcher, the release environment and, with `--tests`, the test environment. |
| Release environment | The Conda environment, defined by `envs/release.yaml`, in which the AIRR export, release and report steps run. Assembly, annotation, clone inference and the other analyses run in the container images instead. |
| Site profile and workflow profile | Snakemake settings. A site profile submits jobs to a scheduler, such as `profiles/global/nci` for PBS on Gadi. The workflow profile `profiles/workflow` requests cores, memory and time per rule, and `profiles/local` does the same for runs on one host. |
| Local execution | Running every job on the launching host or inside one allocation, selected with `--profile none` where a site profile would otherwise apply. |
| Post-TRUST4 run | `--post-trust4`, a run that reuses the retained assemblies and depth counts of an earlier run and skips FASTQ discovery and assembly. |
| Retained assemblies | The TRUST4 outputs kept under `libraries/<library_id>/trust4/` after a run, which a post-TRUST4 run reads. |
| Lock | Snakemake's record that a run is active in a directory. `run_module.sh -u` removes a stale lock after a run was interrupted. |
| Hard link | A second name for the same file. Released tables can be hard links to working tables, so editing a release in place also changes the working copy. |
| Submodule | A Git repository checked out inside another. `profiles/global` is one, initialised with `git submodule update --init --recursive`. |
| Scheduler | The cluster software that queues and starts jobs: PBS Pro on Gadi, Slurm on many other clusters. |
| Allocation | The cores, memory and time a scheduler grants one job. Inside one, the launcher runs every pipeline job within those limits. |
| Login node and compute node | A login node is the shared host you connect to; it limits long processes and usually has network access. Compute nodes run scheduled jobs and, on many clusters, have no network access. |
| `copyq` | The Gadi queue whose jobs have network access. Setup runs there. |
| Service unit (SU) | Gadi's charging unit. Each job is charged for its walltime at a rate set by the larger of its core count and its share of the node's memory. |
| `PROJECT` | On Gadi, the environment variable naming the project charged for jobs. The bundled profile also uses it for the default `/g/data` and `/scratch` paths. |

## Release and provenance

| Term | Meaning |
|---|---|
| AIRR | The Adaptive Immune Receptor Repertoire community standards. The release uses AIRR Rearrangement fields where their meanings match and declares its extensions. |
| MiAIRR | The minimum information about an AIRR study. Structural validation of the release does not establish that this information is available. |
| Repertoire ID | `repertoire_id`, the AIRR identifier linking a library's rearrangements to its Repertoire metadata. |
| Processing identity | `data_processing_id`, the AIRR identifier for how a library was analysed. It is derived from the library analysis code, its settings, the analysis image, the germline reference and auxiliary data files and the AIRR schema, so a change to any of them gives a new value. |
| Effective settings | `effective_config.yaml`, the scientific settings and seed a release was produced with, with the hard-coded choices under `fixed_methods`. |
| Producer provenance | A record written by the job that produced an output. A reused assembly without one has unknown history. |
| Reference receipt | `shared/provenance/imgt_reference.json`, written after the curated germline reference is checked against the analysis image. The [maintenance guide](maintenance.md#audit-the-germline-reference) lists its fields. |
