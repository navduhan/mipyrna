"""Offline packaged Rfam bundle and Bowtie screening tests."""

import json
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mipyrna import utility as mu
from mipyrna.rfamily import FilterRNAfamlies


class RfamBundleTests(unittest.TestCase):
    def test_packaged_rfam_archive_checksum_and_manifest_are_valid(self):
        data_dir = Path(__file__).resolve().parents[1] / "mipyrna/data"
        archive_path = data_dir / "ncRNA_rfam.tar.xz"
        sha_path = data_dir / "ncRNA_rfam.tar.xz.sha256"
        self.assertTrue(archive_path.is_file())
        self.assertTrue(sha_path.is_file())
        expected_sha = sha_path.read_text().split()[0]
        self.assertEqual(mu.file_sha256(archive_path), expected_sha)
        with tarfile.open(archive_path) as archive:
            names = set(archive.getnames())
            self.assertIn("ncRNA_rfam.fa", names)
            self.assertIn("rfam_manifest.json", names)
            self.assertIn("COPYING", names)
            manifest_fh = archive.extractfile("rfam_manifest.json")
            assert manifest_fh is not None
            with manifest_fh:
                manifest = json.load(manifest_fh)
            self.assertEqual(manifest.get("release"), "15.1")
            self.assertGreater(manifest.get("included_families", 0), 0)
            self.assertGreater(manifest.get("excluded_miRNA_families", 0), 0)

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
