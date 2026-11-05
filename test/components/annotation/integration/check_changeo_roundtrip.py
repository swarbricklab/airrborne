"""Reconstruct full and partial indel inputs after restoring exported coordinates.

The annotation integration case supplies real pinned-tool outputs. Both native
MakeDb rows and restored export rows go through the same CreateGermlines binary.
"""
import csv
import gzip
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from workflow.lib.airr_exports import EXPORT_NAMES, GERMLINE_FIELDS
from workflow.lib.alignment_export import add_alignment_evidence, normalise_alignment, read_reference, validate_alignment
from workflow.scripts.restore_changeo_coordinates import restore


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def check(directory):
    native = {row['sequence_id']: row for row in read(directory / 'makedb_pass.tsv.gz')}
    annotations = {identity: {} for identity in native}
    add_alignment_evidence(annotations, [], directory / 'igblast_hits.fmt7.gz')
    reference = read_reference(directory / 'reference/human_IMGT+C.fa')
    exported = []
    for source in read(directory / 'annotated.tsv'):
        row = {EXPORT_NAMES.get(key, key): value for key, value in source.items()}
        row.update({key: native[row['sequence_id']][key] for key in GERMLINE_FIELDS})
        row.update(annotations[row['sequence_id']])
        row['changeo_d_cigar'] = native[row['sequence_id']]['d_cigar']
        normalise_alignment(row)
        validate_alignment(row, reference)
        exported.append(row)
    write(directory / 'export.tsv', exported)
    supported, _ = restore(directory / 'export.tsv', directory / 'restored.tsv', directory / 'unavailable.tsv')
    assert supported == len(exported), 'Indel panel unexpectedly lost reconstructable inputs'
    write(directory / 'native.tsv', list(native.values()))
    for name in ('native', 'restored'):
        subprocess.run(['CreateGermlines.py', '-d', str(directory / (name + '.tsv')), '-g', 'dmask',
                        '-r', '/usr/local/share/germlines/imgt/human/vdj', '--format', 'airr', '--failed',
                        '--outdir', str(directory), '--outname', name], check=True)
    original = {row['sequence_id']: row['germline_alignment_d_mask'] for row in read(directory / 'native_germ-pass.tsv')}
    restored = {row['sequence_id']: row['germline_alignment_d_mask'] for row in read(directory / 'restored_germ-pass.tsv')}
    assert len(original) == len(exported) and original == restored, 'Native-coordinate round trip changed reconstructed germlines'


if __name__ == '__main__':
    check(Path(sys.argv[1]))
