"""Offline bundle tests using tiny cached family files; no downloads or jobs."""

import gzip
import json
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mipyrna import utility as mu
from mipyrna.rfamily import FilterRNAfamlies
from scripts.build_rfam_reference import add_family, build_reference


class RfamBundleTests(unittest.TestCase):
    def test_bundle_excludes_mirnas_preserves_memberships_and_is_reproducible(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            cache = root / "cache"
            (cache / "database_files").mkdir(parents=True)
            (cache / "fasta_files").mkdir()
            rows = []
            for accession, name, kind, seed, full in (
                ("RF00001", "familyA", "Gene; rRNA;", 1, 1),
                ("RF00002", "familyB", "Gene; tRNA;", 0, 2),
                ("RF00003", "mir-1", "Gene; miRNA;", 1, 0),
            ):
                fields = [""] * 19
                fields[0], fields[1], fields[18] = accession, name, kind
                fields[14], fields[15] = str(seed), str(full)
                rows.append("\t".join(fields))
            with gzip.open(cache / "database_files/family.txt.gz", "wt", encoding="latin1") as handle:
                handle.write("\n".join(rows) + "\n")
            (cache / "database_files/family.sql").write_text("fixture schema\n")
            (cache / "README").write_text("fixture Rfam 15.1\n")
            (cache / "COPYING").write_text("CC0\n")
            for accession, records in (
                ("RF00001", ">one\nacgu\n>two\nTTTT\n"),
                ("RF00002", ">three\nACGT\n"),
            ):
                with gzip.open(cache / f"fasta_files/{accession}.fa.gz", "wt") as handle:
                    handle.write(records)
            pins = {
                name: mu.file_sha256(cache / name)
                for name in ("database_files/family.txt.gz", "database_files/family.sql")
            }
            with (
                patch("scripts.build_rfam_reference.METADATA", pins),
                patch("scripts.build_rfam_reference.mu.download_file") as download,
            ):
                first = build_reference(cache, root / "first", jobs=2)
                build_reference(cache, root / "second", jobs=1)
                download.assert_not_called()
            self.assertEqual(first["included_families"], 2)
            self.assertEqual(first["excluded_miRNA_families"], 1)
            self.assertEqual(first["source_records"], 3)
            self.assertEqual(first["unique_sequences"], 2)
            self.assertEqual(first["multi_family_sequences"], 1)
            self.assertEqual(first["metadata_count_differences"], {"RF00002": {"metadata": 2, "fasta": 1}})
            self.assertEqual(mu.file_sha256(root / "first/ncRNA_rfam.tar.xz"), mu.file_sha256(root / "second/ncRNA_rfam.tar.xz"))
            records = mu.read_fasta(root / "first/ncRNA_rfam.fa")
            self.assertEqual(records[0], ("rfam_1|RF00001:familyA|RF00002:familyB", "ACGT"))
            self.assertNotIn("mir-1", (root / "first/ncRNA_rfam.fa").read_text())
            with tarfile.open(root / "first/ncRNA_rfam.tar.xz") as archive:
                manifest = archive.extractfile("rfam_manifest.json")
                assert manifest is not None
                with manifest:
                    self.assertEqual(json.load(manifest)["fasta_sha256"], mu.file_sha256(root / "first/ncRNA_rfam.fa"))
                self.assertIn("COPYING", archive.getnames())

    def test_truncated_gzip_is_not_silently_accepted(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "RF00001.fa.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(">one\nACGT\n")
            path.write_bytes(path.read_bytes()[:-8])
            with self.assertRaises(EOFError):
                add_family({"accession": "RF00001", "name": "family", "expected_records": 2}, path, {})

    @unittest.skipUnless(shutil.which("bowtie") and shutil.which("bowtie-build"), "Requires Bowtie 1")
    def test_real_bowtie_retains_a_read_with_more_than_five_reference_hits(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            resources = root / "resources"
            (resources / "data").mkdir(parents=True)
            (resources / "param").mkdir()
            config = Path(__file__).resolve().parents[1] / "mipyrna/param/bowtie_rfam.ini"
            shutil.copyfile(config, resources / "param/bowtie_rfam.ini")
            mature = "TGAGGTAGTAGGTTGTATAGTT"
            reference = root / "ncRNA_rfam.fa"
            mu.write_fasta([
                (f"rfam_{i}|RF00005:tRNA", "A" * i + mature + "C" * i)
                for i in range(1, 9)
            ], str(reference))
            manifest = root / "rfam_manifest.json"
            manifest.write_text(json.dumps({"release": "15.1", "fasta_sha256": mu.file_sha256(reference)}))
            with tarfile.open(resources / "data/ncRNA_rfam.tar.xz", "w:xz") as archive:
                archive.add(reference, arcname=reference.name)
                archive.add(manifest, arcname=manifest.name)
            reads = str(root / "reads.fa")
            mu.write_fasta([("s_1_x7", mature), ("s_2_x3", "C" * 22)], reads)
            screen = FilterRNAfamlies(samples={"s": ["s", "group", reads]}, outdir=temp)
            with patch("mipyrna.rfamily.ir.files", return_value=resources):
                try:
                    hits, filtered = screen.align_rfam(cpu=1, jobs=1, filter_reads=True)
                except Exception as exc:
                    errors = "\n".join(path.read_text() for path in Path(screen.outdir).rglob("*.err"))
                    raise AssertionError(errors) from exc
            self.assertEqual(hits["s"][2]["read_name"].tolist(), ["s_1_x7"])
            self.assertEqual(hits["s"][2]["name"].tolist(), ["RF00005:tRNA"])
            self.assertEqual(mu.read_fasta(filtered["s"][2]), [("s_2_x3", "C" * 22)])


if __name__ == "__main__":
    unittest.main()
