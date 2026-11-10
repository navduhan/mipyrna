"""Workflow wiring tests with all analysis stages mocked; no jobs are launched."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from mipyrna import utility as mu
from mipyrna.novel import Novel_miRNA
from mipyrna.workflow import _check_read_policy, run_workflow


class ResumePolicyTests(unittest.TestCase):
    def test_unchanged_policy_is_reusable(self):
        with tempfile.TemporaryDirectory() as temp:
            _check_read_policy(temp, False, False)
            (Path(temp) / ".resume/align.json").write_text("{}")
            _check_read_policy(temp, True, False)

    def test_filter_toggle_is_rejected_in_both_directions(self):
        for original in (False, True):
            with tempfile.TemporaryDirectory() as temp:
                _check_read_policy(temp, False, original)
                (Path(temp) / ".resume/align.json").write_text("{}")
                previous = (Path(temp) / ".resume/read_policy.json").read_text()
                with self.assertRaisesRegex(ValueError, "new --outdir"):
                    _check_read_policy(temp, True, not original)
                self.assertEqual((Path(temp) / ".resume/read_policy.json").read_text(), previous)

    def test_legacy_scores_are_rejected_but_collapse_only_can_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            _check_read_policy(temp, True, True)
            (Path(temp) / ".resume/read_policy.json").unlink()
            scores = Path(temp) / "3.miRNA_Discovery/temp"
            scores.mkdir(parents=True)
            (scores / "s_scored.pkl").touch()
            with self.assertRaises(ValueError):
                _check_read_policy(temp, True, True)


class WorkflowRfamWiringTests(unittest.TestCase):
    def test_filter_setting_applies_to_all_downstream_read_consumers(self):
        for filter_reads in (False, True):
            with self.subTest(filter_reads=filter_reads), tempfile.TemporaryDirectory() as temp:
                samples = {"s": ["s", "group", str(Path(temp) / "collapsed.fa")]}
                filtered = {"s": ["s", "group", str(Path(temp) / "filtered.fa")]}
                reads = pd.DataFrame({
                    "read_name": ["s_1_x7", "s_1_x7"], "name": ["ref1", "ref2"], "reads": [7, 7],
                })
                known = pd.DataFrame([{
                    "Known_ID": "ath-miR1", "Known_IDs": "ath-miR1", "Novel_ID": "n1",
                    "precursor": "ACGT" * 15, "mature": "ACGT" * 5,
                    "rfam_reads": 0 if filter_reads else 7, "Rfam": "-" if filter_reads else "ref1,ref2",
                }])
                novel = pd.DataFrame(columns=["Novel_ID", "precursor", "mature", "rfam_reads", "Rfam"])
                predictor = Novel_miRNA()

                def discovery_results(known=known, novel=novel):
                    blast = Path(temp) / "out/3.miRNA_Discovery/blast_results"
                    blast.mkdir(parents=True)
                    mu.write_fasta([("ath-miR1", "ACGT" * 5)], str(blast / "ath_known_mature.fa"))
                    return known.copy(), novel.copy(), None, None

                with (
                    patch("mipyrna.workflow.mu.read_input_file", return_value={
                        "samples": samples, "combinations": [], "targets": pd.DataFrame(),
                    }),
                    patch("mipyrna.workflow.provenance", return_value={}),
                    patch("mipyrna.workflow.set_log_file"),
                    patch("mipyrna.workflow._write_collapsed_fasta", return_value=samples),
                    patch("mipyrna.workflow.FilterRNAfamlies") as screen,
                    patch("mipyrna.workflow.Bowtie_Aligner") as aligner,
                    patch("mipyrna.workflow.Read_process") as reader,
                    patch("mipyrna.workflow.Cluster_reads") as clusterer,
                    patch("mipyrna.workflow.Novel_miRNA", return_value=predictor),
                    patch.object(predictor, "get_novel_miRNA", side_effect=discovery_results),
                    patch("mipyrna.workflow.KnownScan") as scan,
                    patch("mipyrna.workflow.Quantifier") as quantifier,
                    patch("mipyrna.workflow.runDESeq2", return_value=pd.DataFrame()),
                    patch("mipyrna.workflow.generate_report", return_value=("report.pdf", "metadata.json")),
                ):
                    screen.return_value.align_rfam.return_value = ({"s": ["s", "group", reads]}, filtered)
                    screen.return_value.reference_info = {"release": "15.1", "fasta_sha256": "test-reference"}
                    aligner.return_value.run_alignment.return_value = {"s": ["s", "group", "reads.sam"]}
                    reader.return_value.aligned_reads.return_value = {}
                    clusterer.return_value.get_clusters.return_value = {}
                    scan.return_value.scan.return_value = pd.DataFrame([{
                        "Known_ID": "ath-miR1", "total_reads": 10, "detected": True,
                    }])
                    quantifier.return_value.quantify_expression.return_value = pd.DataFrame({"mature_id": ["ath-miR1"], "s": [10]})
                    meta = run_workflow(
                        "samples.txt", temp, "genome.fa", "ath", outdir=str(Path(temp) / "out"),
                        run_qc=False, run_trim=False, known_scan=True, rfam_filter=filter_reads,
                        cpu=4, jobs=1,
                    )
                    expected = filtered if filter_reads else samples
                    self.assertEqual(aligner.return_value.run_alignment.call_args.kwargs["samplesDict"], expected)
                    self.assertEqual(scan.call_args.kwargs["samples"], expected)
                    self.assertEqual(quantifier.call_args.kwargs["samples"], expected)
                    self.assertEqual(meta["rfam"]["samples"]["s"], {
                        "flagged_reads": 7.0, "flagged_unique": 1,
                        "removed_reads": 7.0 if filter_reads else 0.0,
                        "removed_unique": 1 if filter_reads else 0,
                    })
                    self.assertEqual(meta["parameters"]["rfam_filter"], filter_reads)
                    table = pd.read_csv(meta["known_file"], sep="\t")
                    self.assertEqual(table["Known_ID"].tolist(), ["ath-miR1"])
                    self.assertIn("rfam_reads", table)
                    self.assertIn("Rfam", table)
        mu.set_resume(False)


if __name__ == "__main__":
    unittest.main()
