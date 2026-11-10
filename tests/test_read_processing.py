"""Focused read-processing tests; no external tools or workflow runs are launched."""

import io
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from mipyrna import utility as mu
from mipyrna.cli import build_parser
from mipyrna.clustering import Cluster_reads
from mipyrna.known_scan import KnownScan, merge_known_calls
from mipyrna.novel import Novel_miRNA
from mipyrna.reads import Read_process
from mipyrna.rfamily import FilterRNAfamlies


def aligned_frame():
    """Synthetic parsed reads, including tied abundances and two chromosomes."""
    rows = [
        ("r0_x10", "ACGT" * 5, "2", "+", 50, 70, 20, "0", 10, 1),
        ("r1_x10", "ACGT" * 5, "1", "+", 100, 120, 20, "1", 10, 2),
        ("r2_x10", "CGTA" * 5, "1", "+", 101, 121, 20, "2", 10, 3),
        ("r3_x5", "GTAC" * 5, "1", "+", 102, 122, 20, "3", 5, 4),
        ("r4_x5", "GTAC" * 5, "2", "+", 51, 71, 20, "4", 5, 5),
        ("r5_x5", "TACG" * 5, "2", "+", 52, 72, 20, "5", 5, 6),
    ]
    return pd.DataFrame(rows, columns=[
        "read_name", "seq", "name", "strand", "start", "end", "length", "read_id", "reads", "align_id",
    ])


class ClusterOrderingTests(unittest.TestCase):
    def test_no_retained_reads_gives_empty_cluster_tables(self):
        frame = aligned_frame().iloc[:0]
        _, result = Cluster_reads()._cluster_sample("s", ["s", "group", frame])
        self.assertTrue(result[2].empty)
        self.assertTrue(result[3].empty)
        self.assertIn("cluster", result[2])
        self.assertIn("clust_size", result[3])

    def test_shuffled_sam_order_gives_identical_clusters_and_anchors(self):
        frame = aligned_frame()
        expected = Cluster_reads()._cluster_sample("s", ["s", "group", frame])[1]
        predictor = Novel_miRNA(read_cutoff=0)
        for seed in range(12):
            shuffled = frame.sample(frac=1, random_state=seed).reset_index(drop=True)
            # The parser assigns alignment ids in SAM order, not biological order.
            shuffled["align_id"] = range(1, len(shuffled) + 1)
            actual = Cluster_reads()._cluster_sample("s", ["s", "group", shuffled])[1]
            pd.testing.assert_frame_equal(expected[2], actual[2])
            pd.testing.assert_frame_equal(expected[3], actual[3])
            with patch.object(predictor, "get_putative_precursor") as precursor:
                predictor._score_cluster((1, 22, actual[2].query("cluster == 1"), "s"))
                reads = precursor.call_args.kwargs["filtered_cluster"]
                self.assertEqual(reads.iloc[0]["read_name"], "r1_x10")


class RfamTests(unittest.TestCase):
    def test_cli_default_and_config(self):
        self.assertFalse(build_parser().parse_args(["workflow"]).rfam_filter)
        self.assertTrue(build_parser().parse_args(["workflow", "--rfam-filter"]).rfam_filter)
        self.assertTrue(build_parser({"rfam_filter": "true"}).parse_args(["workflow"]).rfam_filter)

    def test_exact_screening_does_not_suppress_multimappers(self):
        config = mu.parse_config_file(Path(__file__).resolve().parents[1] / "mipyrna/param/bowtie_rfam.ini")
        options = " ".join(config["alignment"][2:]).split()
        self.assertIn("-v", options)
        self.assertEqual(options[options.index("-v") + 1], "0")
        self.assertEqual(options[options.index("-k") + 1], "1")
        self.assertNotIn("-m", options)
        self.assertNotIn("--norc", options)

    def test_alignment_and_optional_fasta_filter(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "reads.fa"
            mu.write_fasta([("s_1_x7", "ACGT"), ("s_2_x3", "TGCA")], str(source))
            resources = root / "resources"
            (resources / "data").mkdir(parents=True)
            content = b">ref1\nACGT\n"
            with tarfile.open(resources / "data/ncRNA_rfam.tar.gz", "w:gz") as archive:
                entry = tarfile.TarInfo("ncRNA_rfam.fa")
                entry.size = len(content)
                archive.addfile(entry, io.BytesIO(content))
            samples = {"s": ["s", "group", str(source)]}
            frame = pd.DataFrame({"read_name": ["s_1_x7", "s_1_x7"], "name": ["ref1", "ref2"], "reads": [7, 7]})
            parsed = {"s": ["s", "group", frame]}
            screen = FilterRNAfamlies(samples=samples, outdir=temp)
            with (
                patch("mipyrna.rfamily.ir.files", return_value=resources),
                patch("mipyrna.rfamily.Bowtie_Aligner") as aligner,
                patch.object(Read_process, "aligned_reads", return_value=parsed),
            ):
                hits, kept = screen.align_rfam(cpu=4, jobs=2, filter_reads=False)
                self.assertIs(kept, samples)
                self.assertIs(hits, parsed)
                self.assertFalse((Path(screen.outdir) / "rfam_filtered_reads").exists())
                call = aligner.return_value.run_alignment.call_args
                assert call is not None
                options = call.kwargs
                self.assertEqual(options["outType"], "SAM")
                self.assertEqual((options["cpu"], options["jobs"]), (4, 2))
                _, filtered = screen.align_rfam(filter_reads=True)
            self.assertEqual(list(mu.read_fasta(filtered["s"][2])), [("s_2_x3", "TGCA")])
            self.assertEqual((Path(screen.outdir) / "ncRNA_rfam.fa").read_bytes(), content)

    def test_empty_rfam_hits_preserve_reads(self):
        with tempfile.TemporaryDirectory() as temp:
            source = str(Path(temp) / "reads.fa")
            mu.write_fasta([("s_1_x7", "ACGT")], source)
            samples = {"s": ["s", "group", source]}
            screen = FilterRNAfamlies(samples=samples, outdir=temp)
            filtered = screen._extract_rfam_unmapped(
                samples, {"s": ["s", "group", pd.DataFrame(columns=["read_name"])]}
            )
            self.assertEqual(list(mu.read_fasta(filtered["s"][2])), [("s_1_x7", "ACGT")])

    def test_discovery_alert_counts_each_supporting_read_once(self):
        frame = aligned_frame().query("name == '1'").copy()
        frame["cluster"] = 1
        frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
        entry = {
            "mature_check": "ok", "precursor": "ACGT" * 20, "precursor_length": 80,
            "start": 90, "end": 170, "star": "TGCA", "position": "5-prime",
            "struct": "." * 80, "mfe": -20.0,
        }
        predictor = Novel_miRNA()
        with (
            patch.object(Read_process, "get_sequence", return_value="ACGT" * 5),
            patch("mipyrna.novel.fold_hairpins", return_value=[entry]),
            patch("mipyrna.novel.Precursor_Features") as features,
        ):
            features.return_value.generate_features.return_value = {"feature": 1.0}
            result = predictor.get_putative_precursor(
                frame, rfam={"r1_x10": {"ref2", "ref1"}}
            )
            self.assertIsNotNone(result)
            assert result is not None
            self.assertEqual(result.iloc[0]["rfam_reads"], 10)
            self.assertEqual(result.iloc[0]["Rfam"], "ref1,ref2")


class KnownCountTests(unittest.TestCase):
    def test_scan_alert_counts_and_tied_modes_are_order_independent(self):
        table = pd.DataFrame([{
            "Known_ID": "ath-miR1", "hairpin": "ath-MIR1", "hairpin_source": "miRBase",
            "hairpin_start": 10, "Database": "miRBase", "precursor": "A" * 60,
            "precursor_length": 60, "mature": "A" * 20, "mature_length": 20,
            "chrom": "1", "start": 100, "end": 160, "mature_start": 110,
            "strand": "+", "struct": "." * 60, "mfe": -10.0, "mature_check": "ok",
            "position": "5-prime", "star": "T" * 20,
        }])
        frame = pd.DataFrame([
            ("s_1_x7", "ath-miR1", 11, 32, 7),
            ("s_2_x7", "ath-miR1", 10, 30, 7),
            ("s_1_x7", "ath-miR1", 11, 32, 7),
        ], columns=["read_name", "name", "start", "end", "reads"])
        scan = KnownScan(rfam={"s": {"s_1_x7": {"ref1", "ref2"}}})
        expected = scan._count(table, {"s": frame}, {})
        for seed in range(8):
            actual = scan._count(table, {"s": frame.sample(frac=1, random_state=seed)}, {})
            pd.testing.assert_frame_equal(expected, actual)
        row = expected.iloc[0]
        self.assertEqual(row["total_reads"], 14)
        self.assertEqual(row["rfam_reads"], 7)
        self.assertEqual(row["Rfam"], "ref1,ref2")
        self.assertEqual(row["modal_5p_offset"], 0)
        self.assertEqual(row["length_mode"], 20)
        self.assertTrue(row["detected"])
        # Annotation-only calls keep the alert after merging into the known table.
        expected["features"] = ""
        expected["score"] = 1.0
        expected["ml_score"] = 1.0
        merged, added, _, _ = merge_known_calls(pd.DataFrame(), expected, "ath")
        self.assertEqual(added, 1)
        self.assertEqual(merged.iloc[0]["rfam_reads"], 7)
        self.assertEqual(merged.iloc[0]["Rfam"], "ref1,ref2")


if __name__ == "__main__":
    unittest.main()
