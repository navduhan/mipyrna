"""Focused read-processing tests; no external tools or workflow runs are launched."""

import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from mipyrna import utility as mu
from mipyrna.cli import build_parser
from mipyrna.clustering import Cluster_reads
from mipyrna.differential import _split_contrast
from mipyrna.known_scan import KnownScan, drop_detected_novel, merge_known_calls
from mipyrna.novel import (
    Novel_miRNA,
    _dinucleotide_shuffle,
    check_mature_position,
    estimate_star_sequence,
)
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
        self.assertNotIn("--strata", options)
        self.assertNotIn("--norc", options)

    def test_alignment_and_optional_fasta_filter(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "reads.fa"
            mu.write_fasta([("s_1_x7", "ACGT"), ("s_2_x3", "TGCA")], str(source))
            resources = root / "resources"
            (resources / "data").mkdir(parents=True)
            content = b">rfam_1|RF00001:5S_rRNA|RF00002:5_8S_rRNA\nACGT\n"
            with tarfile.open(resources / "data/ncRNA_rfam.tar.xz", "w:xz") as archive:
                entry = tarfile.TarInfo("ncRNA_rfam.fa")
                entry.size = len(content)
                archive.addfile(entry, io.BytesIO(content))
                reference = root / "reference.fa"
                reference.write_bytes(content)
                manifest = json.dumps({"release": "15.1", "fasta_sha256": mu.file_sha256(reference)}).encode()
                entry = tarfile.TarInfo("rfam_manifest.json")
                entry.size = len(manifest)
                archive.addfile(entry, io.BytesIO(manifest))
            samples = {"s": ["s", "group", str(source)]}
            frame = pd.DataFrame({
                "read_name": ["s_1_x7"], "name": ["rfam_1|RF00001:5S_rRNA|RF00002:5_8S_rRNA"], "reads": [7],
            })
            screen = FilterRNAfamlies(samples=samples, outdir=temp)
            with (
                patch("mipyrna.rfamily.ir.files", return_value=resources),
                patch("mipyrna.rfamily.Bowtie_Aligner") as aligner,
                patch.object(Read_process, "aligned_reads", side_effect=[
                    {"s": ["s", "group", frame.copy()]}, {"s": ["s", "group", frame.copy()]},
                ]),
            ):
                hits, kept = screen.align_rfam(cpu=4, jobs=2, filter_reads=False)
                self.assertIs(kept, samples)
                self.assertEqual(set(hits["s"][2]["name"]), {"RF00001:5S_rRNA", "RF00002:5_8S_rRNA"})
                self.assertEqual(hits["s"][2].drop_duplicates("read_name")["reads"].sum(), 7)
                self.assertFalse((Path(screen.outdir) / "rfam_filtered_reads").exists())
                call = aligner.return_value.run_alignment.call_args
                assert call is not None
                options = call.kwargs
                self.assertEqual(options["outType"], "SAM")
                self.assertEqual((options["cpu"], options["jobs"]), (4, 2))
                _, filtered = screen.align_rfam(filter_reads=True)
                aligner.reset_mock()
                (Path(screen.outdir) / "ncRNA_rfam.fa").write_text("stale reference\n")
                with self.assertRaisesRegex(ValueError, "differs from the bundle"):
                    screen.align_rfam()
                aligner.assert_not_called()
            self.assertEqual(list(mu.read_fasta(filtered["s"][2])), [("s_2_x3", "TGCA")])
            self.assertEqual(screen.reference_info["release"], "15.1")

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

    def test_opposite_arm_star_reads_are_counted_and_deduplicated(self):
        table = pd.DataFrame([{
            "Known_ID": "ath-miR156a-5p", "hairpin": "ath-MIR156a", "hairpin_source": "miRBase",
            "hairpin_start": 10, "Database": "miRBase", "precursor": "A" * 120,
            "precursor_length": 120, "mature": "A" * 21, "mature_length": 21,
            "chrom": "2", "start": 1000, "end": 1120, "mature_start": 1010,
            "strand": "+", "struct": "(" * 50 + "." * 20 + ")" * 50, "mfe": -40.0,
            "mature_check": "ok", "position": "5-prime", "star": "T" * 21,
        }])
        frame = pd.DataFrame([
            ("s_1_x12", "ath-miR156a-5p", 10, 31, 12),
            ("s_2_x5", "ath-miR156a-5p", 85, 106, 5),
            ("s_2_x5", "ath-miR156a-5p", 85, 106, 5),
        ], columns=["read_name", "name", "start", "end", "reads"])
        counts = KnownScan()._count(table, {"s": frame}, {"ath-miR156a-5p": (85, 21)})
        self.assertEqual(counts.iloc[0]["total_reads"], 12.0)
        self.assertEqual(counts.iloc[0]["star_reads"], 5.0)

    def test_merge_keeps_opposite_arms_distinct_and_merges_same_arm_variants(self):
        known = pd.DataFrame([{
            "Novel_ID": "ath_miPyRNA_2_miR_1", "Known_ID": "ath-miR156a-5p",
            "Known_IDs": "ath-miR156a-5p", "chrom": "2", "mature": "CTGACAGAAGAGAGTGAGCACA",
            "mature_length": 22, "precursor": "HP_SEQ", "precursor_length": 120,
            "start": 100, "end": 220, "mature_start": 110, "strand": "+",
            "struct": "." * 120, "mfe": -40.0, "mature_check": "ok", "features": "",
            "mature_reads": 50, "position": "5-prime", "star_reads": 5,
            "inconsistent_fraction": 0.0, "star": "TGCTCACTGCTCTTTCTGTCAGAT",
            "randfold_pvalue": 0.01, "cluster": 1, "rfam_reads": 0.0, "Rfam": "-",
            "score": 0.9, "ml_score": 0.9, "heuristic_score": 2.0,
            "Database": "miRBase", "raw_ID": "r1",
        }])
        scan = pd.DataFrame([
            {
                "Known_ID": "Ath-miR156a", "Database": "PmiREN", "precursor": "HP_SEQ",
                "precursor_length": 120, "mature": "TGACAGAAGAGAGTGAGCACA", "mature_length": 21,
                "chrom": "2", "start": 100, "end": 220, "mature_start": 111, "strand": "+",
                "struct": "." * 120, "mfe": -40.0, "mature_check": "ok", "features": "",
                "position": "5-prime", "star": "TGCTCACTGCTCTTTCTGTCAGAT",
                "total_reads": 50.0, "unique_reads": 10.0, "star_reads": 5.0,
                "rfam_reads": 0.0, "Rfam": "-", "score": 0.9, "ml_score": 0.9, "detected": True,
            },
            {
                "Known_ID": "ath-miR156a-3p", "Database": "miRBase", "precursor": "HP_SEQ",
                "precursor_length": 120, "mature": "TGCTCACTGCTCTTTCTGTCAGAT", "mature_length": 24,
                "chrom": "2", "start": 100, "end": 220, "mature_start": 185, "strand": "+",
                "struct": "." * 120, "mfe": -40.0, "mature_check": "ok", "features": "",
                "position": "3-prime", "star": "CTGACAGAAGAGAGTGAGCACA",
                "total_reads": 20.0, "unique_reads": 20.0, "star_reads": 50.0,
                "rfam_reads": 0.0, "Rfam": "-", "score": 0.9, "ml_score": 0.9, "detected": True,
            },
        ])
        merged, added, _, siblings = merge_known_calls(known, scan, "ath")
        self.assertEqual(added, 1)
        self.assertEqual(siblings, 1)
        self.assertEqual(merged["Known_ID"].tolist(), ["ath-miR156a-5p", "ath-miR156a-3p"])
        self.assertEqual(merged.iloc[0]["Known_IDs"], "ath-miR156a-5p;Ath-miR156a")

    def test_drop_detected_novel_handles_unmapped_and_paralog_loci(self):
        scan = pd.DataFrame([
            {
                "Known_ID": "ath-miR999", "chrom": "", "strand": "", "start": "", "end": "",
                "mature_start": "", "precursor": "UNMAPPED_PRECURSOR", "detected": True,
            },
            {
                "Known_ID": "ath-miR156a-5p", "chrom": "2", "strand": "+", "start": 100, "end": 220,
                "mature_start": 110, "precursor": "HP_A", "detected": True,
            },
        ])
        scan.attrs["all_loci"] = {"ath-miR156a-5p": [("4", "+", 510, 500, 620)]}
        novel = pd.DataFrame([
            {"Novel_ID": "n1", "chrom": "4", "strand": "+", "start": 505, "end": 615, "mature_start": 510, "precursor": "HP_B"},
            {"Novel_ID": "n2", "chrom": "1", "strand": "-", "start": 900, "end": 990, "mature_start": 910, "precursor": "UNIQUE_NOVEL"},
        ])
        filtered, removed = drop_detected_novel(novel, scan)
        self.assertEqual(removed, 1)
        self.assertEqual(filtered["Novel_ID"].tolist(), ["n2"])


class UtilityAndContrastTests(unittest.TestCase):
    def test_split_contrast_supports_hyphenated_groups(self):
        c1, c2, label = _split_contrast("Col-0-ago1-27", {"Col-0", "ago1-27"})
        self.assertEqual((c1, c2, label), ("Col-0", "ago1-27", "Col-0-ago1-27"))

    def test_read_input_file_skips_blank_lines_and_deduplicates_pairs(self):
        with tempfile.TemporaryDirectory() as temp:
            sheet = Path(temp) / "samples.txt"
            sheet.write_text("\n1 s1 Col-0 s1.fq\n\n2 s2 ago1-27 s2.fq\n   \n")
            parsed = mu.read_input_file(str(sheet), temp, paired=False)
            self.assertEqual(list(parsed["samples"]), ["s1", "s2"])
            self.assertEqual(parsed["combinations"], ["Col-0-ago1-27"])

    def test_local_run_raises_on_nonzero_exit(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(RuntimeError):
                mu.local_run("false_cmd", "exit 2", outdir=temp)

    def test_allocated_cpus_honors_slurm_cpus_per_task(self):
        with patch.dict("os.environ", {"SLURM_CPUS_PER_TASK": "12"}):
            self.assertEqual(mu._allocated_cpus(), 12)
            self.assertEqual(mu.get_cpu(), 9)

    def test_reference_index_and_integer_chromosome_lookup(self):
        import pickle

        from mipyrna.reads import Read_process

        with tempfile.TemporaryDirectory() as temp:
            fa = Path(temp) / "ref.fa"
            fa.write_text(">MT\nAAAAAAAAAA\n>1\nACGTACGTAA\n")
            reader = Read_process(reference=str(fa))
            reader.ensure_index()
            self.assertTrue((Path(temp) / "ref.fa.fai").is_file())
            self.assertEqual(reader.get_sequence((1, 0, 4, "+")), "ACGT")
            restored = pickle.loads(pickle.dumps(reader))
            self.assertEqual(restored.get_sequence(("1", 0, 4, "-")), "ACGT")


class LiteratureBiogenesisTests(unittest.TestCase):
    def test_single_variant_collapsed_stack_and_strand_separation(self):
        rows = [
            ("s_1_x250", "TGACAGAAGAGAGTGAGCACA", "1", "+", 100, 121, 21, "1", 250, 1),
            ("s_2_x80", "TGTGCTCACTCTCTTCTGTCA", "1", "-", 105, 126, 21, "2", 80, 2),
        ]
        frame = pd.DataFrame(rows, columns=[
            "read_name", "seq", "name", "strand", "start", "end", "length", "read_id", "reads", "align_id",
        ])
        _, out = Cluster_reads()._cluster_sample("s", ["s", "WT", frame])
        clusts = out[3]
        self.assertEqual(len(clusts), 2)
        self.assertEqual(set(clusts["strand"]), {"+", "-"})
        self.assertEqual(clusts["clust_size"].tolist(), [21, 21])
        self.assertEqual(sorted(clusts["reads"].tolist()), [80, 250])

    def test_duplex_check_and_two_nt_overhang_star_projection(self):
        # Canonical uninterrupted stem-loop with 2-nt 3' overhangs on both arms:
        # 5' flank (4 nt) + 21-nt 5p mature (19 bp + 2-nt 3' overhang) + 8-nt loop + 21-nt 3p star
        stem_5p = "CTGACAGAAGAGAGTGAGC"  # 19 nt paired core
        over_5p = "AC"
        mature_5p = stem_5p + over_5p
        comp = stem_5p[::-1].translate(str.maketrans("ACGT", "TGCA"))
        star_3p = comp + "GT"  # 19 nt paired to stem_5p + 2-nt 3' overhang
        precursor = "AAAA" + mature_5p + "TTTTAAAA" + star_3p + "TT"
        self.assertEqual(check_mature_position(precursor, mature_5p), "ok")
        self.assertEqual(check_mature_position(precursor, star_3p), "ok")
        self.assertEqual(estimate_star_sequence(precursor, mature_5p), star_3p)
        self.assertEqual(estimate_star_sequence(precursor, star_3p), mature_5p)
        loop_span = precursor[18:39]
        self.assertEqual(check_mature_position(precursor, loop_span), "mature in loop")

    def test_dinucleotide_shuffle_preserves_exact_dinucleotides(self):
        from collections import Counter

        import numpy as np

        rng = np.random.default_rng(7)
        seq = "UGACAGAAGAGAGUGAGCACAUAAUGGCUAAAUUAGCUUGCTCACUGCUCUUUCUGUCAGA"
        expected = Counter(zip(seq, seq[1:]))
        for _ in range(50):
            shuffled = _dinucleotide_shuffle(seq, rng)
            self.assertEqual(Counter(zip(shuffled, shuffled[1:])), expected)
            self.assertEqual((shuffled[0], shuffled[-1]), (seq[0], seq[-1]))

    def test_hairpin_wide_star_and_precision_across_separated_arms(self):
        stem_5p = "CTGACAGAAGAGAGTGAGC"
        mature_5p = stem_5p + "AC"
        star_3p = stem_5p[::-1].translate(str.maketrans("ACGT", "TGCA")) + "GT"
        precursor = "AAAA" + mature_5p + "TTTTTTTTTTTT" + star_3p + "AAAA"
        entry = {
            "mature_check": "ok", "precursor": precursor, "precursor_length": len(precursor),
            "start": 100, "end": 100 + len(precursor), "hairpin_start": 4,
            "star": star_3p, "position": "5-prime",
            "struct": "." * 4 + "(" * 19 + "." * 14 + ")" * 19 + "." * 6, "mfe": -32.0,
        }
        star_start = 100 + precursor.index(star_3p)
        mature_cluster = pd.DataFrame([
            ("r_m1_x40", mature_5p, "1", "+", 104, 125, 21, "1", 40, 1, 1),
            ("r_m2_x10", mature_5p[:-1], "1", "+", 105, 125, 20, "2", 10, 2, 1),
        ], columns=[
            "read_name", "seq", "name", "strand", "start", "end", "length", "read_id", "reads", "align_id", "cluster",
        ])
        hairpin_reads = pd.DataFrame([
            ("r_m1_x40", mature_5p, "1", "+", 104, 125, 21, "1", 40, 1),
            ("r_m2_x10", mature_5p[:-1], "1", "+", 105, 125, 20, "2", 10, 2),
            ("r_s1_x8", star_3p, "1", "+", star_start, star_start + 21, 21, "3", 8, 3),
            ("r_deg_x2", "A" * 20, "1", "+", 126, 146, 20, "4", 2, 4),
        ], columns=[
            "read_name", "seq", "name", "strand", "start", "end", "length", "read_id", "reads", "align_id",
        ])
        predictor = Novel_miRNA()
        with (
            patch.object(Read_process, "get_sequence", return_value=precursor),
            patch("mipyrna.novel.fold_hairpins", return_value=[entry]),
            patch("mipyrna.novel.Precursor_Features") as features,
        ):
            features.return_value.generate_features.return_value = {"f": 1.0}
            out = predictor.get_putative_precursor(mature_cluster, hairpin_reads=hairpin_reads)
            assert out is not None
            row = out.iloc[0]
            self.assertEqual(row["star_reads"], 8)
            self.assertAlmostEqual(row["five_prime_consistency"], 40 / 50)
            self.assertAlmostEqual(row["inconsistent_fraction"], 2 / 60)

    def test_known_survives_strict_without_star_while_novel_enforces_literature_rules(self):
        predictor = Novel_miRNA(species="ath", species_type="plants", require_star=True)
        pooled = pd.DataFrame([
            {
                "chrom": "1", "mature": "CTGACAGAAGAGAGTGAGCAC", "mature_length": 21,
                "score": 0.92, "mature_check": "ok", "mature_reads": 30, "star_reads": 0,
                "randfold_pvalue": 0.01, "inconsistent_fraction": 0.0,
                "five_prime_consistency": 0.95, "heuristic_score": 2.5,
            },
            {
                "chrom": "1", "mature": "CTGACAGAAGAGAGTGAGCAC", "mature_length": 21,
                "score": 0.92, "mature_check": "ok", "mature_reads": 25, "star_reads": 6,
                "randfold_pvalue": 0.01, "inconsistent_fraction": 0.05,
                "five_prime_consistency": 0.90, "heuristic_score": 8.5,
            },
            {
                "chrom": "2", "mature": "TTGGACTGAAGGGAGCTCCCT", "mature_length": 21,
                "score": 0.88, "mature_check": "ok", "mature_reads": 40, "star_reads": 0,
                "randfold_pvalue": 0.02, "inconsistent_fraction": 0.0,
                "five_prime_consistency": 0.92, "heuristic_score": 2.0,
            },
        ])
        dedup = predictor.filter_miRNA(pooled)
        # Star-supported sample wins the tie on identical precursor score.
        self.assertEqual(
            int(dedup.loc[dedup["mature"] == "CTGACAGAAGAGAGTGAGCAC", "star_reads"].iloc[0]), 6
        )
        # Strict pre-BLAST filter keeps star-less candidates so known miRNAs are not lost.
        strict = predictor.apply_strict_filters(dedup)
        self.assertEqual(len(strict), 2)
        # Post-BLAST novel Dicer filter enforces star support (and 24-nt plant star rule).
        kept_novel = predictor._dicer_filter(strict)
        self.assertEqual(kept_novel["mature"].tolist(), ["CTGACAGAAGAGAGTGAGCAC"])

        relaxed = Novel_miRNA(species="ath", species_type="plants", require_star=False)
        plant_24 = pd.DataFrame([
            {
                "mature": "A" * 24, "mature_length": 24, "star_reads": 0,
                "inconsistent_fraction": 0.0, "five_prime_consistency": 0.9, "heuristic_score": 3.0,
            },
            {
                "mature": "C" * 21, "mature_length": 21, "star_reads": 0,
                "inconsistent_fraction": 0.0, "five_prime_consistency": 0.9, "heuristic_score": 3.0,
            },
        ])
        self.assertEqual(relaxed._dicer_filter(plant_24)["mature"].tolist(), ["C" * 21])

    def test_cli_literature_defaults_and_overrides(self):
        default_args = build_parser().parse_args(["workflow"])
        self.assertTrue(default_args.require_star)
        self.assertEqual(default_args.read_cutoff, 20)
        self.assertEqual(default_args.max_inconsistent, 0.25)
        self.assertEqual(default_args.min_5p_consistency, 0.5)
        no_star = build_parser().parse_args(["workflow", "--no-require-star"])
        self.assertFalse(no_star.require_star)
        ini_no_star = build_parser({"require_star": "false"}).parse_args(["workflow"])
        self.assertFalse(ini_no_star.require_star)
        ini_override = build_parser({"require_star": "false"}).parse_args(["workflow", "--require-star"])
        self.assertTrue(ini_override.require_star)


if __name__ == "__main__":
    unittest.main()
