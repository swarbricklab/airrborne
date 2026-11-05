# Immune repertoire: {{scope}}

This release contains bulk RNA-seq immune receptor results for {{n_libraries}} libraries from {{n_tissues}} tissues and {{n_donors}} donors. {{donor_repeats}}

## Start here

- Open `report/index.html` in a browser to explore the results offline.
- Start table-based analysis with `metrics/headline.tsv`. It has one row per library covering input depth, receptor yields, IGH/TRB diversity, antibody classes and somatic hypermutation (SHM).
- Check `metrics/qc.tsv` for the fraction of annotated read support that passes quality control (QC) and the fraction available for grouping at each locus.
- Use `data_dictionary.csv` for column definitions across {{n_tables}} table schemas, including units, denominators and missing values.
- For downstream analysis, join metadata from its original source by unique `library_id` or `tissue_id`, keeping libraries without a match. Project metadata is included only in the report.

The results describe receptor sequences recovered from bulk RNA-seq. `consensus_count` measures read support for a sequence. Cell counts, molecule counts and receptor-chain pairings are unavailable from this analysis, and low recovery can leave a receptor undetected.

Read estimates alongside their support counts and status columns. The [interpretation guide]({{guide_url}}#known-limitations) explains how these limits affect comparisons:

- Observed diversity and dominance depend on both the amount of recovered support and how it is distributed across groups. Compare libraries on rarefied values and read their status. {{target_summary}}
- Low recovered support limits comparison. Read observed support, the configured diversity minimum and `rarefaction_status`. No diversity minimum is applied by default.
- Whole-V SHM uses only members that pass the evidence requirements in all five analysed V subregions. Report eligibility alongside SHM and check how the eligible population affects the result. The methods below give this run's support and occupancy settings.
- Selection estimates are exploratory.

Keep the whole `release/` folder together when copying it. Save edited tables separately, because released tables can be hard links to the original results. To rerun from assemblies, also retain the native assembly and depth files listed in the [running guide]({{running_url}}#reuse-retained-assemblies).

## Receptors, clones and clonotypes

A library is one sequenced preparation. Several libraries can come from the same tissue or donor. Each sequence record describes one receptor chain, without identifying its partner in a cell.

The seven locus names follow [IMGT nomenclature](https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGTnomenclature.php):

| Locus | Receptor chain |
|---|---|
| IGH | Immunoglobulin heavy chain |
| IGK, IGL | Immunoglobulin kappa and lambda light chains |
| TRA, TRB, TRD, TRG | T-cell receptor alpha, beta, delta and gamma chains |

The variable (V), diversity (D) and joining (J) gene segments form the receptor's variable region. Only IGH, TRB and TRD use a D segment. A gene call assigns a sequence to a reference gene or allele. When annotation cannot choose between candidates, the call retains the ties.

The junction includes the third complementarity-determining region (CDR3) and its two flanking anchor residues. The [AIRR definitions](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#junction-versus-cdr3) explain the distinction.

The pipeline groups related sequences in two ways:

- An IGH clone is a heavy-chain lineage inferred within one donor, pooling the donor's libraries in this scope. `clone_id` names it.
- At IGK, IGL, TRA, TRB, TRD and TRG, a clonotype groups retained members within one library by an exact nucleotide junction and the same resolved V and J entities. `clonotype_id` identifies it.

Tables containing both clones and clonotypes use `group_unit` to distinguish them. A member represents one distinct full sequence within a library and locus. Its support can include eligible source records assigned to that representative. For light chains, these can include records from the same contig whose junction differs from the representative's. See [records and read support](#records-and-read-support) before summing counts.

### Ambiguous gene calls

Unresolved multi-gene V or J calls stay out of grouping. Allele ties within one gene can qualify; multi-gene calls require complete candidate lists with identical complete reference sequences. Near-identical paralogues can therefore lose grouping support.

A call naming several alleles of one gene counts under that gene. When alleles from different genes share one complete reference sequence, they count under a composite entity. Its identifier starts with `VSEQ_` or `JSEQ_`, followed by a digest of the sequence; the `_names` columns list the genes.

A gene can appear under its own name, inside a composite entity or both. To combine these entries, use `gene_equivalence_group` in `metrics/gene_usage.tsv`.

## Files

Paths below are relative to this release folder. Start with `metrics/headline.tsv` for library summaries. The other `metrics/` tables contain the full measurements for each topic, and `metrics/detail/` contains clone, clonotype and member records.

The Key column lists the columns that together identify a row. It omits `scope_id`, which is constant within this release. Individual key columns can repeat.

| File | One row per | Key | Loci | Contents |
|---|---|---|---|---|
| `manifest/scope_manifest.tsv` | library | `library_id` | all | Donor and tissue IDs of each included library. |
| `metrics/headline.tsv` | library | `library_id` | all | Yields at all seven loci and selected IGH and TRB summaries. |
| `metrics/yield.tsv` | library | `library_id` | all | Input read pairs, bases and FASTQ file pairs, TRUST4-reported support and yields. |
| `metrics/qc.tsv` | library and locus | `library_id`, `locus` | all | IgBLAST-annotated, common-QC and grouping-eligible support, with QC and grouping-reason fractions. |
| `metrics/diversity.tsv` | library and locus | `library_id`, `locus` | all | Richness, Hill diversity, dominance and rarefied values of IGH clones and clonotypes. |
| `metrics/rarefaction_curves.tsv` | library, locus and support depth | `library_id`, `locus`, `support_depth` | loci with a rarefaction target | Interpolated diversity at each depth up to observed support, for loci with a configured target. Other loci have no curve rows, whatever their repertoire. An empty repertoire has a depth-zero endpoint. |
| `metrics/isotype_usage.tsv` | library, isotype level and isotype | `library_id`, `isotype_level`, `isotype` | IGH | Class and subclass support and fractions of callable support, including zero counts. |
| `metrics/shm.tsv` | library, isotype group and V region | `library_id`, `isotype_group`, `v_region` | IGH | Member and clone SHM summaries with their eligible counts. |
| `metrics/selection_pressure.tsv` | library and isotype group | `library_id`, `isotype_group` | IGH | BASELINe selection estimates for the CDRs and FWRs, with the clone counts behind them. |
| `metrics/gene_usage.tsv` | library, locus, segment and gene entity | `library_id`, `locus`, `gene_type`, `gene_key` | all | Separate V and J usage as shares of groups and of read support. |
| `metrics/vj_pair_usage.tsv` | library, locus and V and J entity pair | `library_id`, `locus`, `v_gene_key`, `j_gene_key` | all | Observed V and J pairs as shares of groups and of read support. |
| `metrics/donor_overlap_pairs.tsv` | donor, locus and library pair | `donor_id`, `locus`, `library_id_1`, `library_id_2` | overlap loci | Shared groups and overlap indices for every pair of one donor's libraries. |
| `metrics/donor_overlap_summary.tsv` | donor and locus | `donor_id`, `locus` | overlap loci | Union and sharing counts and mean indices, for donors with at least two libraries. |
| `metrics/detail/igh/clones.tsv` | IGH clone | `clone_id` | IGH | V and J entities, junction length, size and isotype classes of each clone across the donor's libraries. |
| `metrics/detail/igh/clone_members.tsv` | IGH member | `sequence_id` | IGH | Each member's clone, calls, junction, isotype and support. The AIRR export holds its sequence and alignments under the same `sequence_id`. |
| `metrics/detail/igh/clone_member_shm.tsv` | IGH member and V region | `sequence_id`, `v_region` | IGH | Each member's SHM eligibility, position counts, mutation count and rate, where it has a measured rate. |
| `metrics/detail/igh/clone_member_shm_coverage.tsv` | IGH member and V region | `sequence_id`, `v_region` | IGH | Complete coverage audit, including unmeasured members and reconstruction availability. |
| `metrics/detail/igh/clone_isotype.tsv` | library, IGH clone, isotype level and isotype | `library_id`, `clone_id`, `isotype_level`, `isotype` | IGH | Members and read support of each isotype a clone shows within one library. These sum grouping-eligible member support, so they do not sum to `metrics/isotype_usage.tsv`. |
| `metrics/detail/igh/clone_shm.tsv` | library, IGH clone, isotype group and V region | `library_id`, `clone_id`, `isotype_group`, `v_region` | IGH | Clone SHM for each subset that has an SHM-eligible member. Counts include the subset's other members. |
| `metrics/detail/non_igh/clonotypes.tsv` | clonotype | `library_id`, `locus`, `clonotype_id` | non-IGH | Junction, V and J entities, support and representative record of each clonotype. |
| `metrics/detail/non_igh/clonotype_members.tsv` | clonotype member | `sequence_id` | non-IGH | Each member's clonotype and support. |
| `metrics/detail/group_abundance.tsv` | library, locus and clone or clonotype | `library_id`, `locus`, `group_key` | all | Support and rank of each group with positive support. |
| `metrics/detail/shared_groups.tsv` | donor, locus, group and library | `donor_id`, `locus`, `group_key`, `library_id` | overlap loci | Groups detected in at least two of a donor's libraries, with a row for each of the donor's libraries. |
| `metrics/detail/qc_trail.tsv` | library, locus, stage, isotype group and V region | `library_id`, `locus`, `stage`, `isotype_group`, `v_region` | all | Units and read support entering and leaving each assembly, annotation, QC and analysis stage. Criterion failures can overlap. |
| `igh_clone_threshold/scope_threshold.tsv` | scope | one row | IGH | The operative IGH clone threshold, its source and the donor counts. |
| `igh_clone_threshold/donor_fits.tsv` | donor | `donor_id` | IGH | Each donor's fit, its inputs, `fit_status` and the fitted mixture's parameters, including donors whose fit is not used. |
| `igh_clone_threshold/donor_distance_histograms.tsv` | donor and distance bin | `donor_id`, `distance_bin_start` | IGH | Nearest-neighbour distance histograms and fitted densities. Donors without distances have zero observed counts and density; fitted component densities are unavailable without a usable fit. |
| `airr/*_rearrangements.tsv.gz` | AIRR record | `sequence_id` | all | Every common-QC record of one library, in the AIRR Rearrangement format. |
| `manifest/release_files.tsv` | released file | `path` | all | Size and SHA-256 hash of every file except itself and `.complete.json`. |

{{scope_assignment}}The scope manifest defines this release's population.

### Read the data dictionary

`data_dictionary.csv` defines the columns, units and missing values for each table. A table's dictionary ID is its filename without `.tsv`; AIRR files use `exported_AIRR`.

| Dictionary field | Meaning |
|---|---|
| `is_primary_key` | Columns marked TRUE together identify a row. |
| `nullable` | FALSE means the column has no empty cells. |
| `references_table`, `references_column` | The table and column containing the referenced values. |
| `references_when` | The rows to which a reference applies. `\|` separates alternatives. |
| `airr_field`, `relationship` | How the column maps to AIRR fields. |
| `source` | The table, column, setting or tool from which a value comes. |
| `source_level` | Whether an identifier describes a library, tissue or donor. |

Some references name intermediate pipeline tables that are outside this release.

## Read the measurements

### Records and read support

Each analysis uses a particular set of records. Isotype composition and grouping apply different requirements, so callable IGH support can exceed grouping-eligible support in the same library.

| Record set | Definition and use |
|---|---|
| TRUST4-reported records | Full-CDR3 records before pipeline QC, including imputed and non-productive records. These supply receptor yield in `yield.tsv`. |
| IgBLAST-annotated records | Non-imputed records successfully annotated at the assigned locus. Their support is the denominator of `nonproductive_fraction` and `common_qc_fraction_of_annotated` in `qc.tsv`. |
| Common-QC records | Annotated records passing the shared checks in the [QC matrix]({{qc_matrix_url}}). The AIRR export includes all of them. Their support is the denominator of grouping, grouping-reason and isotype-callability fractions in `qc.tsv`. |
| Callable IGH support | Common-QC support with a constant-region call resolving to one class or subclass. `isotype_usage.tsv` uses the callable support at the stated resolution. |
| Grouping-eligible records | Common-QC records meeting the locus's junction and V/J requirements. Their support contributes once to clones or clonotypes through members. Grouping failures can overlap. |

Diversity, abundance, overlap and V/J usage describe the groups built from eligible members. SHM uses IGH members that pass the coverage rule for the region being measured. Library SHM summaries give each clone median equal weight.

Selection estimates use IGH clone members with observed and reconstructed germline alignments. They do not require the SHM coverage rule.

Annotated counts reflect several steps: exclusion of imputed inputs, removal of MakeDb failures and locus assignment by IgBLAST.

A record with a full CDR3 can still have an incomplete V(D)J sequence. TRUST4's AIRR conversion excludes partial CDR3 calls before yield is measured.

### SHM eligibility fractions

SHM eligibility fractions use different denominators:

- In `qc.tsv`, coverage fractions divide by `n_grouped_igh_members`, including unmeasured members.
- In the headline, `igh_shm_eligible_fraction_of_measured_members_whole_v` divides by measured members only.
- To calculate the eligible-clone fraction, divide `igh_n_shm_eligible_clones_whole_v` by `igh_richness`. This fraction is unavailable at zero richness and is not exported as a headline column.

### Missing values

An empty cell means missing, unavailable or inapplicable. Read the status column to find the reason, and preserve empty cells when importing the tables. Zero has the meaning given in the data dictionary. The dictionary uses `empty` in an `options` list to mean a blank cell.

Keep identifiers as strings so that leading zeroes and literal `NA` IDs survive import. The [import examples](#import-identifiers-safely) show how.

Diversity and dominance estimates require the configured minimum support. Richness, support and diagnostics remain available below that minimum. A headline top-5 fraction can therefore be empty even when richness is nonzero. Check `diversity_status` in `metrics/diversity.tsv` for the reason.

`log2_igg_iga` is unavailable when class-callable support is zero or when either IgG or IgA has zero support. Class and subclass compositions use separate callable-support denominators.

### Diversity and rarefaction

Shannon entropy uses natural logarithms and is measured in nats. Hill diversity expresses the number of equally supported groups that would give the same diversity. Support fractions have no units.

Rarefaction estimates diversity in a smaller draw of recovered support, without replacement. It compares libraries at a common depth of recovered read support. Differences in FASTQ depth, cell diversity and receptor expression per cell remain unadjusted. Curves require a configured target and stop at observed support. The inverse Simpson asymptote does not prescribe a sequencing depth.

The rarefied Shannon-derived value is `exp(E[H])`, and the inverse Simpson value is `1 / E[sum(p^2)]`. Here, the expectation averages over possible draws and `p` is each group's fraction within a draw. Neither value is the average of the corresponding Hill diversity across draws.

### IGH clone thresholds

All donors in a scope use the same IGH clustering threshold, including donors whose own fit was unusable. The threshold measures differences between junction sequences. Clustering uses single linkage, so members of a clone can be further apart than the threshold if other members connect them.

For each sequence, SHazaM `distToNearest` finds the nearest sequence from the same donor with matching V/J germline entities and junction length, but a different junction. It records the fraction of differing junction positions. A sequence with no such neighbour has no distance.

With `params.lineage.clonal_threshold: auto`, the pipeline takes the median of usable donor thresholds. If no donor fit is usable, it uses 0.1. A donor needs at least `params.lineage.min_fit_distances` distances for fitting, with a default of 500. Read `fit_status` in `igh_clone_threshold/donor_fits.tsv` for each donor's outcome and `igh_clone_threshold/scope_threshold.tsv` for the chosen threshold and its source.

SHazaM `findThreshold` models smaller distances with a gamma component and larger distances with a normal component representing unrelated sequences. It restricts the normal component to the observed fitting range, then places the cutoff so that the fraction `params.lineage.fit_specificity` lies above it. This specificity describes the fitted model within the observed distance range. The fit does not measure the probability that two sequences belong to different biological clones. `effective_config.yaml` records the settings used.

## Import identifiers safely

Python's CSV reader keeps cells as strings, including literal `NA` and leading zeroes. Convert only named numeric columns after import:

```python
import csv
with open("metrics/headline.tsv", newline="") as handle:
    headline = list(csv.DictReader(handle, delimiter="\t"))
input_pairs = [None if row["input_read_pairs"] == "" else int(row["input_read_pairs"])
               for row in headline]
```

In R, disable the default `NA` missing token while retaining numeric blanks:

```r
headline <- read.delim("metrics/headline.tsv", colClasses = "character",
                       na.strings = "", quote = "", comment.char = "", check.names = FALSE)
headline$input_read_pairs <- as.numeric(headline$input_read_pairs)
```

## Join and compare tables

Use tables from one jointly analysed release for combined analysis. Scopes can infer different IGH clones for the same library, and concatenating overlapping releases duplicates observations. The [identifier guide]({{guide_url}}#identifiers) explains what each identifier names and where it is unique.

### Join keys

| Join | Keys |
|---|---|
| Library tables to the manifest or headline | `scope_id, library_id` |
| IGH abundance or shared groups to clones | `scope_id, group_key=clone_id` |
| Non-IGH abundance to clonotypes | `scope_id, library_id, locus, group_key=clonotype_key` |
| Non-IGH members to clonotypes | `scope_id, clonotype_id` |
| Member-region tables to members | `scope_id, sequence_id` |
| Members to their library's AIRR export | `sequence_id` |

Retain donor and locus context, and keep `group_unit` when combining clones and clonotypes. Non-IGH `clonotype_id` identifies one clonotype in one library and locus and is unique within the release. When comparing frozen outputs, include a release identifier because a scope name alone does not identify a revision.

AIRR source records link to their member through `member_representative_sequence_id`. Only grouping-eligible records have this link; the AIRR export also includes common-QC records excluded from grouping. Each member representative's AIRR `clone_id` identifies its clone or clonotype in this scope.

### Avoid counting support twice

`member_consensus_count` in the member tables equals `member_consensus_count` on the AIRR representative row. AIRR `consensus_count` retains each source record's original support. Choose either eligible source support or member support for a total. Joining sources to members repeats member totals, so each member must contribute only once.

Abundance `n_members` counts the retained unique sequences of each IGH clone or clonotype within a library. Isotype and region strata overlap. Filter to the intended stratum before summing counts.

### Shared groups and representatives

Shared-group descriptors repeat across library rows, including rows for non-detections. `v_gene_key` and `j_gene_key` identify the analytical gene entities; adjacent columns give readable labels. Use the identifier columns when joining tables.

In the donor clone table, `representative_sequence_id` identifies the highest-support observed member. This member also supplies the sharing descriptors. For an AIRR join, use its `representative_library_id`, which can differ from the sharing row's library. Inferred ancestors and clone consensuses are absent from this table.

Library clone details are sorted by abundance, with tied ranks such as 1, 2, 2, 4. Selecting rank <= 5 can therefore return more than five groups. Donor clones are sorted by pooled support, and shared groups keep their library rows together.

### Pair comparisons

Pair-table differences compare `library_id_2` with `library_id_1`. This order does not imply chronology. Pairs below the support minimum retain their values but carry a status flag. Donor mean statuses identify empty libraries or low support; a below-minimum pair flag does not suppress an otherwise defined donor mean.

### SHM coverage and measurements

`metrics/detail/igh/clone_member_shm_coverage.tsv` has a row for every grouped IGH member in each of six regions: whole V, FWR1, CDR1, FWR2, CDR2 and FWR3. `measurement_status` is `unmeasured`, `measured_ineligible` or `eligible`. A member with no measured rate still has a row marked `unmeasured`.

Whole V overlaps the subregions. Select one region or count distinct members instead of summing all region rows. `evidence_sequence_id` identifies the qualifying source chosen for a region and is empty when none qualifies. `germline_reconstruction_status` records reconstruction availability separately from coverage.

`metrics/detail/igh/clone_member_shm.tsv` contains only member-region rows with a measured rate. Where `shm_eligible` is FALSE, `n_mutations`, `shm_frequency` and `median_trust4_consensus_base_votes` are empty, but position counts remain. An absent row means no measured rate. Missing rates never mean zero mutations.

Join this table to `clone_members.tsv` for isotype to rebuild `shm.tsv` and `clone_shm.tsv`. Clone SHM omits subsets with no eligible member. To inspect the source selected for an eligible region, use `evidence_sequence_id` from the complete coverage audit.

The AIRR alignments support mutation reconstruction with the pinned SHazaM method, including partial coverage. Position counts and the `params.shm` settings in `effective_config.yaml` reproduce subregion eligibility decisions. The release does not contain position votes. To recheck the votes behind `n_vote_supported_positions`, retain these files from the original results tree:

- `libraries/<library_id>/metrics/<library_id>_igh_features.tsv`
- `libraries/<library_id>/airr/<library_id>_igh_consensus_votes.tsv`

## Reuse the tables

Before using a released table in another tool, select the records and support measure needed for the analysis. Check which fields the tool reads by default.

**Select records and support.**

- Each AIRR row is one common-QC source record, and `consensus_count` is its own support.
- Members are the rows where `is_member_representative` is T. Only they carry a `clone_id`, so selecting on `clone_id` changes the population from source records to members.
- `member_consensus_count` holds the support transferred to a member. It can exceed the same row's `consensus_count`, so filtering to members and reading `consensus_count` undercounts the grouped repertoire.
- Subset one locus. IGH clones and clonotypes are different units.
- Keep a donor's libraries in one data split when training a model. A library ID is not a pairing key.

**Check the fields a tool reads by default.**

The stored `sequence` is already in V orientation, so standard `rev_comp` is F. Query coordinates and CIGARs refer directly to that sequence. `changeo_rev_comp` records MakeDb's transformation of the original query.

This convention applies to output schema 2.0.0. In schema 1.0.0, `rev_comp` incorrectly contained the historical transformation flag. Before reusing those earlier files with standards-aware tools, preserve that flag separately and set `rev_comp` to F.

- Standard V germline coordinates are one-based inclusive positions on the ungapped selected allele and agree with `v_cigar`. `changeo_v_germline_start/end` preserve padded IMGT reconstruction coordinates. Restore those native fields explicitly before using Change-O `CreateGermlines`; the standard export is not its native input.
- Reverse-D hits keep their call, validated query coordinates and selected native hit evidence. Their standard D CIGAR and germline endpoints are empty, with `d_alignment_status=unavailable_reverse_orientation`. Reconstruction remains unavailable for affected IG records. Tool-reported short reverse matches do not establish a biological D inversion.
- `c_call` is the raw call. Use it where `c_call_callable` is T, or read `isotype_class` and `isotype_subclass`, which are empty where a call was masked or does not resolve. A callable call need not name one class.
- `germline_alignment` and `germline_alignment_d_mask` describe each record's own alignment and can differ within a clone. This release holds no clone-level germline, so a lineage tool that needs one germline per clone needs a separate reconstruction step. Germlines are reconstructed for IG loci only, so `germline_alignment_d_mask` is empty at TCR loci.
- Common QC does not require a canonical junction. Some `junction_aa` values do not start with C or end with F or W, and `semantic_checks` counts them. No record is filtered on this.

The converter in the matching pipeline checkout restores native coordinates only for `germline_reconstruction_status=available` rows and writes all other rows separately:

```bash
python workflow/scripts/restore_changeo_coordinates.py \
  --input rearrangements.tsv.gz --output changeo_input.tsv \
  --unavailable reconstruction_unavailable.tsv
```

Sequences and alignments stay unchanged. `CreateGermlines.py` requires the matching pinned, IMGT-gapped reference and an explicit clone-level reconstruction step when a tool needs one germline per clone.

**Junction properties.** Pass amino-acid columns to Alakazam [`aminoAcidProperties`](https://alakazam.readthedocs.io/en/stable/topics/aminoAcidProperties/) with `nt = FALSE`. `cdr3_aa` is `junction_aa` without its first and last residues, so leave `trim` at its default of FALSE for it.

```r
aminoAcidProperties(members, seq = "junction_aa", nt = FALSE)
aminoAcidProperties(members, seq = "cdr3_aa", nt = FALSE)
```

**Repertoire Dissimilarity Index.** RDI `calcRDI` takes a count matrix with one column per repertoire. For individual genes, build it from `metrics/gene_usage.tsv` for one locus and `gene_type`. For V/J pairs, use `metrics/vj_pair_usage.tsv`. Choose `n_groups` to weight groups equally or `consensus_count` to weight by read support.

To combine equivalent gene entities, sum gene usage by `gene_equivalence_group`. The V/J pair table does not have that column, so look it up in gene usage:

- For V, filter to `gene_type=V` and join `v_gene_key=gene_key` with `scope_id, library_id, locus`.
- For J, filter to `gene_type=J` and join `j_gene_key=gene_key` with the same context.

Require unique lookup keys. Each join must preserve the pair-table row count and total support. Never assign a composite entity's full count to each of its genes.

Fill absent features with zero only for libraries with positive `grouped_consensus_count` at that locus in `metrics/diversity.tsv`. A nonempty repertoire below the support minimum remains distinct from an empty repertoire.

**sumrep.** These columns are the inputs its data specification names:

| sumrep input | In this release |
|---|---|
| `sequence_alignment`, `germline_alignment`, `v_call`, `d_call`, `j_call`, `junction`, `junction_aa`, `vj_in_frame`, `np1_length`, `np2_length`, `clone_id` | Present under these names. |
| `v_3p_del`, `d_5p_del`, `d_3p_del`, `j_5p_del`, `vd_insertion`, `dj_insertion`, `vj_insertion` | Absent. Deriving them needs the reference, and partial sequences limit deletion lengths. |
| Lineage trees | Absent. |

**Single-chain encoders.** `junction_aa` keeps the two anchor residues and `cdr3_aa` removes them, so match the model's convention. Keep `sequence_id` when deduplicating peptides, so that each embedding maps back to its observations.

## Validation and provenance

`airr/Repertoire.json` links libraries, donors, tissues and processing records. Unavailable required metadata stays null; optional unavailable metadata is omitted. Study, clinical, preparation and sequencing-run metadata are not supplied.

`airr/validation.json` separates structural and linkage checks from `strict_metadata_gaps`, including omitted nonnullable optional objects. Structural validation does not establish full MiAIRR compliance. Its `semantic_checks` state the V germline coordinate convention and count, per locus, the members, blank germlines and noncanonical junction peptides.

`provenance.json` records the pipeline source, tool versions, reference identities and available assembly history. Missing assembly history is recorded as unknown. `effective_config.yaml` records the scientific settings and, under `fixed_methods`, the choices fixed in code. Reused assemblies can have different producer settings.

Reference paths under `/usr/local/share` refer to files inside the analysis image named by `images.immcantation_docker`. To inspect the germlines, run that image:

```bash
image=$(awk '/immcantation_docker:/ {print $2}' effective_config.yaml)
# Linux with Singularity; Apptainer accepts the same arguments.
singularity exec "docker://${image}" ls /usr/local/share/germlines/imgt/human/vdj
# Linux or macOS with Docker running.
docker run --rm --platform linux/amd64 "$image" ls /usr/local/share/germlines/imgt/human/vdj
```

`output_schema_version` in `provenance.json` and `.complete.json` identifies the output format. The source revision, settings and file manifest identify this particular result.

`.complete.json` records the library, table-schema and file counts, excluding itself from the file count. It also records the SHA-256 hash of `manifest/release_files.tsv`. The manifest lists the sizes and hashes of all other files except itself and the receipt. Run this check from the release directory to verify those files:

```bash
python3 - <<'PYCODE'
import csv
import hashlib
from pathlib import Path
with Path("manifest/release_files.tsv").open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        path = Path(row["path"])
        digest = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
        assert path.stat().st_size == int(row["n_bytes"]), path
        assert digest.hexdigest() == row["sha256"], path
print("Verified released file sizes and SHA-256 hashes")
PYCODE
```

If the pipeline code was modified locally, `source.state_sha256` identifies the source evidence under `shared/provenance/source_evidence/` in the original results tree. This evidence is outside the release folder.

The [interpretation guide]({{guide_url}}) links to the documentation at this run's `source.revision`; local changes need the source evidence.

## Methods as run

{{methods}}

## Citations

Cite the pipeline, [airrborne](https://github.com/swarbricklab/airrborne), at the source revision above, and cite the tools it ran. The methods above and `provenance.json` give the tool versions. `airr/validation.json` gives the version of the release validator, and the `airr-python` entry in `provenance.json` describes the analysis image's package.

| Tool | Use | Citation |
|---|---|---|
| TRUST4 | Receptor assembly | Song et al. 2021, Nat Methods, [doi:10.1038/s41592-021-01142-2](https://doi.org/10.1038/s41592-021-01142-2) |
| IMGT/GENE-DB | Germline reference | Lefranc et al. 2015, Nucleic Acids Res, [doi:10.1093/nar/gku1056](https://doi.org/10.1093/nar/gku1056) |
| IgBLAST | V(D)J and C annotation | Ye et al. 2013, Nucleic Acids Res, [doi:10.1093/nar/gkt382](https://doi.org/10.1093/nar/gkt382) |
| Change-O and SHazaM | AIRR conversion, germline reconstruction, nearest-neighbour distances and mutation counts | Gupta et al. 2015, Bioinformatics, [doi:10.1093/bioinformatics/btv359](https://doi.org/10.1093/bioinformatics/btv359) |
| SHazaM `findThreshold` | Clone threshold fitting | Nouri and Kleinstein 2018, Front Immunol, [doi:10.3389/fimmu.2018.01687](https://doi.org/10.3389/fimmu.2018.01687) |
| SCOPer `hierarchicalClones` | IGH clone assignment | Gupta et al. 2017, J Immunol, [doi:10.4049/jimmunol.1601850](https://doi.org/10.4049/jimmunol.1601850) |
| BASELINe (SHazaM) | Selection pressure | Yaari et al. 2012, Nucleic Acids Res, [doi:10.1093/nar/gks457](https://doi.org/10.1093/nar/gks457); S5F model, Yaari et al. 2013, Front Immunol, [doi:10.3389/fimmu.2013.00358](https://doi.org/10.3389/fimmu.2013.00358) |
| Hill-number rarefaction | Diversity at a common support depth | Chao et al. 2014, Ecol Monogr, [doi:10.1890/13-0133.1](https://doi.org/10.1890/13-0133.1) |
| AIRR Python library | Release validation | Vander Heiden et al. 2018, Front Immunol, [doi:10.3389/fimmu.2018.02206](https://doi.org/10.3389/fimmu.2018.02206) |
