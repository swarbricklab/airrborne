#!/usr/bin/env python3
"""Check the curated reference and auxiliary data against the analysis image and write the receipt.

The audit_imgt_reference rule runs this script in the analysis image. It
builds the reference and the auxiliary data again from the image's files and
requires both written files to match byte for byte. It then reads every
IgBLAST database with blastdbcmd and requires each to hold the same alleles
and sequences as the reference and as the germline files it was built from.
The receipt records the image identity, the IMGT/GENE-DB download date, the
files and records the curation read, every database comparison, the
auxiliary entries added and the J alleles still without one. Scope provenance
reads the receipt.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from workflow.lib.germline_reference import audit_reference


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference", type=Path, required=True, help="the curated reference FASTA")
    parser.add_argument("--auxiliary", type=Path, required=True,
                        help="the supplemented auxiliary data file written beside the reference")
    parser.add_argument("--germline-dir", type=Path, required=True,
                        help="the image's human/vdj directory, one IMGT FASTA per V, D and J locus")
    parser.add_argument("--constant-dir", type=Path, required=True,
                        help="the image's human/constant directory, one IMGT FASTA per constant locus")
    parser.add_argument("--snapshot", type=Path, required=True,
                        help="the image's IMGT.yaml, recording the download source and date")
    parser.add_argument("--igblast-db", type=Path, required=True,
                        help="params.igblast.igblast_db, holding database/ and optional_file/")
    parser.add_argument("--image-receipt", type=Path, required=True,
                        help="the analysis image receipt obtain_image.py wrote")
    parser.add_argument("--receipt", type=Path, required=True, help="the reference receipt to write")
    parser.add_argument("--metadata", type=Path, required=True, help="The curated allele metadata table.")
    args = parser.parse_args()
    try:
        image = json.loads(args.image_receipt.read_text())
        receipt = audit_reference(reference=args.reference, auxiliary=args.auxiliary, metadata=args.metadata,
                                  germline_dir=args.germline_dir, constant_dir=args.constant_dir,
                                  snapshot=args.snapshot, igblast_db=args.igblast_db, image=image)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.receipt.with_name(args.receipt.name + ".tmp")
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(args.receipt)
    databases = receipt["igblast_databases"]
    auxiliary = receipt["auxiliary_data"]
    print(f"{receipt['reference']['name']} sha256 {receipt['reference']['sha256']} matches a fresh curation "
          f"of the IMGT/GENE-DB download dated {receipt['imgt_snapshot']['date']}; "
          f"{len(databases)} IgBLAST databases hold the same {sum(entry['alleles'] for entry in databases.values())} "
          f"alleles; {auxiliary['written']['name']} matches a fresh supplement with "
          f"{len(auxiliary['supplemented'])} added entries; {len(auxiliary['j_alleles_without_entry'])} J alleles "
          f"have no auxiliary data entry; receipt {args.receipt}")
