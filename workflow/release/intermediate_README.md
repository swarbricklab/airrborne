# Intermediate evidence

These intermediate files support debugging and checks of how the results were produced. For routine analysis or sharing, use the tables in the scope's `release/` directory.

Intermediate donor files contain more measurements than the released tables:

- Donor member files retain alignments, regional SHM measurements and selected evidence IDs.
- Member SHM files retain measured mutation counts and rates when coverage fails. Those values are masked in the release.
- Clone SHM files retain subsets with no eligible members. Those subsets are omitted from the release.

The complete member coverage audit includes all six regions for every grouped member.

The `sort_*` columns control the order of intermediate donor rows and are removed during collection. Because the intermediate and released tables include different records, check their populations before comparing or joining them. Native MakeDb and IgBLAST evidence is kept in the library analysis directories.

## Locate intermediate evidence

Paths are relative to the originating results tree, outside this folder unless shown as a donor-local filename. Each row population must be read with its scope and library context.

| File pattern | Row unit and producer | Difference from the release and regeneration |
|---|---|---|
| `libraries/<library>/airr/*_makedb_pass.tsv.gz` | Native MakeDb record, from reannotation. | Native coordinates and annotations before the export adapter; regenerate with pinned annotation inputs/reference. |
| `libraries/<library>/airr/*_igh_consensus_votes.tsv` | Common-QC IGH source record, from `build_igh_consensus_votes.R`. | Position vote vectors and mapping reasons are not shipped; requires native final.out and annotation. |
| `libraries/<library>/metrics/*_igh_features.tsv` | IGH source record, from `calculate_igh_features.R`. | Source mutation and coverage measurements precede member evidence selection; requires votes and reconstructed germlines. |
| `scopes/<scope>/intermediate/<donor>/<donor>_clone_members.tsv` | Retained IGH inference member, from donor assignment. | Keeps alignments, unmasked measurements and regional evidence IDs; requires donor pool and scope threshold. |
| Member and clone SHM files in the donor directory | Member-region or clone/library/isotype/region, from SHM summaries. | Measured ineligible rates and no-eligible-member clone subsets can remain here; collection masks or omits them in release tables. |
| `scopes/<scope>/intermediate/report_data/` | Prepared report payloads, from report preparation. | Display metadata and loading layout differ from scientific tables; regenerate from authoritative metrics and metadata. |

Keep the native input bundle and source evidence if a later audit or rerun is required. The retention matrix in `docs/maintenance.md` in the matching checkout gives the complete task-specific requirements.
