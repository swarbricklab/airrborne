"""Native coordinate restoration preserves sequence evidence and unavailable rows."""
import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from workflow.scripts.restore_changeo_coordinates import COORDINATES, restore


class ChangeoViewTests(unittest.TestCase):
    def write_supported_input(self, source):
        row = {"germline_reconstruction_status": "available",
               **{field: "1" for field in COORDINATES},
               **{"changeo_" + field: "2" for field in COORDINATES}}
        with source.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=row, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)

    def test_failed_second_replacement_restores_previous_output_pair(self):
        for existing in ((False, False), (True, False), (False, True), (True, True)):
            with self.subTest(existing=existing), tempfile.TemporaryDirectory() as directory:
                source, output, other = [Path(directory) / name for name in
                                         ("export.tsv", "native.tsv", "unavailable.tsv")]
                self.write_supported_input(source)
                original = source.read_bytes()
                for destination, present in zip((output, other), existing):
                    if present:
                        destination.write_text("previous " + destination.name)
                replace = Path.replace

                def fail_second_replacement(path, destination):
                    if destination == other and ".backup." not in path.name:
                        raise OSError("injected second output failure")
                    return replace(path, destination)

                with patch.object(Path, "replace", fail_second_replacement):
                    with self.assertRaisesRegex(OSError, "second output failure"):
                        restore(source, output, other)
                for destination, present in zip((output, other), existing):
                    self.assertEqual(destination.exists(), present)
                    if present:
                        self.assertEqual(destination.read_text(), "previous " + destination.name)
                self.assertEqual(source.read_bytes(), original)
                self.assertEqual(set(Path(directory).iterdir()),
                                 {source, *(p for p, present in zip((output, other), existing) if present)})

    def test_directory_destination_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output, other = [Path(directory) / name for name in
                                     ("export.tsv", "native.tsv", "unavailable.tsv")]
            self.write_supported_input(source)
            output.write_text("previous supported rows\n")
            other.mkdir()
            with self.assertRaises(IsADirectoryError):
                restore(source, output, other)
            self.assertEqual(output.read_text(), "previous supported rows\n")
            self.assertTrue(other.is_dir())
            self.assertEqual(set(Path(directory).iterdir()), {source, output, other})

    def test_restore_native_coordinates_and_partition_unavailable_rows(self):
        supported = dict(sequence_id="full", sequence_alignment="..AC-GT", germline_alignment="..ACAGT",
                         germline_reconstruction_status="available", v_germline_start="12", v_germline_end="15",
                         d_germline_start="2", d_germline_end="4", d_cigar="1I2=",
                         changeo_v_germline_start="1", changeo_v_germline_end="20",
                         changeo_d_germline_start="2", changeo_d_germline_end="4", changeo_d_cigar="1D2=")
        unavailable = dict(supported, sequence_id="reverse", germline_reconstruction_status="unavailable_reverse_d")
        with tempfile.TemporaryDirectory() as directory:
            source, output, other = [Path(directory) / name for name in ("export.tsv", "native.tsv", "unavailable.tsv")]
            with source.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=supported, delimiter="\t")
                writer.writeheader(); writer.writerows((supported, unavailable))
            original = source.read_bytes()
            self.assertEqual(restore(source, output, other), (1, 1))
            with output.open(newline="") as handle: row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["v_germline_start"], "1")
            self.assertEqual(row["v_germline_end"], "20")
            self.assertEqual(row["d_cigar"], "1D2=")
            self.assertEqual(row["sequence_alignment"], supported["sequence_alignment"])
            self.assertEqual(row["germline_alignment"], supported["germline_alignment"])
            with other.open(newline="") as handle: self.assertEqual(next(csv.DictReader(handle, delimiter="\t")), unavailable)
            self.assertEqual(source.read_bytes(), original)
            expected = (output.read_bytes(), other.read_bytes())
            self.assertEqual(restore(source, output, other), (1, 1))
            self.assertEqual((output.read_bytes(), other.read_bytes()), expected)
            self.assertEqual(set(Path(directory).iterdir()), {source, output, other})

    def test_invalid_input_preserves_existing_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output, other = [Path(directory) / name for name in ("export.tsv", "native.tsv", "unavailable.tsv")]
            source.write_text("sequence_id\ngermline_missing\n")
            output.write_text("previous supported rows\n")
            other.write_text("previous unavailable rows\n")
            with self.assertRaisesRegex(ValueError, "native reconstruction"):
                restore(source, output, other)
            self.assertEqual(output.read_text(), "previous supported rows\n")
            self.assertEqual(other.read_text(), "previous unavailable rows\n")
