"""Write one library's AIRR export with its members' clone assignments in one scope.

The IGH clone_id comes from the donor's clone-member shard and the clonotype_id of
the other loci from the library's clonotype members, each keyed by the qualified
sequence_id of a member representative. annotate_member_assignments() in
airr_exports.py adds clone_id and copies every other cell unchanged.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.input.export_module).resolve().parents[2]))
from workflow.lib.airr_exports import annotate_member_assignments

csv.field_size_limit(sys.maxsize)
library = str(snakemake.wildcards.library_id)
assignments = {}
for path, column in ((snakemake.input.clone_members, "clone_id"),
                     (snakemake.input.clonotype_members, "clonotype_id")):
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["library_id"] != library:
                continue
            key = row["sequence_id"]
            if key in assignments or not row[column]:
                raise ValueError(f"{key} has a repeated or empty member assignment in {path}")
            assignments[key] = row[column]
count = annotate_member_assignments(Path(snakemake.input.airr), Path(snakemake.output.airr),
                                    assignments)
Path(snakemake.log[0]).write_text(f"Assigned {len(assignments)} members among {count} records\n")
