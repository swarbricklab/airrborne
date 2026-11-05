#!/usr/bin/env python3
"""Curate TRUST4's V/D/J/C reference and IgBLAST's auxiliary data from the analysis image.

The obtain_imgt_reference rule runs this script in the analysis image. It reads
the vdj and constant FASTAs of Immcantation's IMGT/GENE-DB download, keeps the
first record for each allele name and writes one FASTA with bare allele
headers. It also writes IgBLAST's auxiliary data file: NCBI's file from the
image, plus an entry for each J allele that file lacks, copied from an allele
with an identical sequence or from a same-gene allele of equal length. An
allele neither rule covers stays without an entry.
workflow/lib/germline_reference.py documents both rules, and
audit_imgt_reference.py checks the written files afterwards.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.germline_reference import build_auxiliary_data, build_reference, render_reference_metadata, write_reference


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--germline-dir", type=Path, required=True,
                        help="the image's human/vdj directory, one IMGT FASTA per V, D and J locus")
    parser.add_argument("--constant-dir", type=Path, required=True,
                        help="the image's human/constant directory, one IMGT FASTA per constant locus")
    parser.add_argument("--auxiliary-data", type=Path, required=True,
                        help="NCBI's auxiliary data file in the image, optional_file/human_gl.aux")
    parser.add_argument("--output", type=Path, required=True, help="the reference FASTA to write")
    parser.add_argument("--auxiliary-output", type=Path, required=True,
                        help="the supplemented auxiliary data file to write")
    parser.add_argument("--metadata-output", type=Path, required=True,
                        help="The original allele metadata and sequence identity classes to write.")
    args = parser.parse_args()
    try:
        text, curation = build_reference(args.germline_dir, args.constant_dir)
        auxiliary, report = build_auxiliary_data(args.germline_dir, args.auxiliary_data)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
    write_reference(text, args.output)
    write_reference(render_reference_metadata(curation["allele_metadata"]), args.metadata_output)
    write_reference(auxiliary, args.auxiliary_output)
    print(f"Wrote {curation['alleles_written']} alleles from {curation['records_read']} records in "
          f"{len(curation['germline_files'])} files to {args.output}; skipped "
          f"{len(curation['duplicate_records_skipped'])} records repeating an allele name")
    added = ", ".join(entry["allele"] for entry in report["supplemented"]) or "none"
    without = ", ".join(entry["allele"] for entry in report["j_alleles_without_entry"]) or "none"
    print(f"Wrote {args.auxiliary_output} from {report['shipped_entries']} shipped entries after checking "
          f"{report['sibling_pairs_checked']} same-gene pairs; added entries for {added}; "
          f"J alleles still without an entry: {without}")
