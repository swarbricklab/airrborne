# Constructed annotation panel

`segments.fasta` contains 17 public IMGT germline alleles from the pipeline's pinned human reference, whose SHA-256 is `085536aa1a10e8492de129f3420f36713d71992b3e48d426e82fa166fd8ea50e`. The source is Immcantation Suite 4.8.0, OCI digest `44e228a45b047bdb352363ab20a670336b003e98f3d1bac505f78af41a91ab67`.

These bases derive from IMGT, the international ImMunoGeneTics information system, Montpellier, France, under CC BY 4.0. See the public-source attribution and reference preparation command in [the public data guide](../README.md).

`joins.tsv` describes seven deliberately simple V-J sequences. `v_end` is the exclusive end in the gapped reference; `j_start` is a zero-based inclusive coordinate.

Construction keeps V through the conserved IMGT Cys104, removes reference dots, inserts five alanine codons and appends J from its conserved Phe/Trp-Gly-X-Gly motif. There is no constant segment or mutation. The test assigns fixed constructed support counts. IGH, TRB and TRD deliberately omit D; a missing D call is allowed by common QC. These sequences test annotation and routing. Recombination frequencies and assembler performance are outside their scope.

The integration test runs the real pinned IgBLAST, MakeDb and CreateGermlines through the production wrapper. It adds a sequence with no receptor alignment, a noncanonical native locus, and a TRD sequence labelled TRA to distinguish annotation loss, routing loss and diagnostic locus reassignment. Separate empty inputs exercise skipped tool calls.

The panel is separate from the four-library FASTQ fixture. The returned panel also passes through common QC, including the reassigned TRD record, so the diagnostic cannot silently become an exclusion.

The indel cases change the IGH join: a three-base insertion and deletion at V position 60 and in the J segment, three N bases at position 61 and the reversed sequence. They run through common QC, consensus-vote mapping and SHM measurement.

The constructed native matrices give ten agreeing votes to each called base and zero votes at N bases. The reverse-oriented record receives no mapped votes. The test counts informative positions independently from the aligned strings.

The auxiliary-data test builds the pipeline reference and the supplemented `human_gl.aux` from the image with `obtain_imgt_reference.py`, constructs one V-J query for each J allele NCBI's file lacks and one for the allele whose entry the pipeline copies, and checks each junction against the shipped entry. It uses the V alleles above through Cys104, five alanine codons padded into the J frame and the whole J.

The shared-V cases use `TRAV14/DV4*01` through Cys104, the same five alanine codons and either `TRAJ1*01` or `TRDJ1*01` from the coordinates in `joins.tsv`. One delta sequence also appends `TRDC*01`. Both delta sequences must retain their support through the annotation audit, common QC and exact grouping; the alpha sequence remains TRA.

Separate joins pair `TRDV1*01` with `TRAJ1*01` and `TRAV1-1*01` with `TRDJ1*01`, testing the shared V pool beyond dual-named genes. Unit cases retain conflicting segment calls for rejection by common QC.

The reference case builds the original allele metadata sidecar and checks partial-reference exclusions, named-allele equivalence and the digest of the pinned sequence-class memberships. These cases use the same image as annotation.

The orientation case appends the first 270 bases of `IGHG1*01` to the constructed IGH sequence and annotates both it and its reverse complement. Both must retain 270 informative C positions and a callable isotype. These cases exercise annotation directly and do not establish the frequency of reverse-oriented contigs in native TRUST4 output.

The allele FASTA and join table are committed inputs. The annotation integration suite constructs the query sequences from the declared joins and writes transient native vote cases in its temporary run directory. Tests do not rewrite either file.

Deliberate changes require review against the pinned reference, updated data checksums and execution of the actual annotation tools. Expected loci, junctions and V/J genes are derived independently from the declared joins.
