"""miRNA Quantification Module

This module aligns reads to miRNA hairpins and builds the raw count matrix
used for differential expression.

Features:
    - Alignment of collapsed reads to the hairpins of known and novel miRNAs with Bowtie
    - Counting of read abundance per mature miRNA, for reads inside a window around the mature position
    - Unweighted counts, or counts shared across the miRNAs a read maps to
    - Merged count matrix across samples written to Excel

Configuration:
    Inputs are the miRNA list (ID, hairpin, mature), the collapsed read samples and
    the output directory. Alignment settings are in ``bowtie_mirna.ini``. Output is
    ``Raw_read_counts.xlsx`` and a ``quant_align`` directory.

Dependencies:
    - HTSeq: SAM parsing
    - pandas: Count tables
    - Bowtie via mipyrna.aligner

Classes / Functions / Exceptions:
    - Quantifier: Count reads per miRNA across samples.

:Author: Naveen Duhan
"""

import importlib.resources as ir
import os
from collections import defaultdict
from functools import reduce

import HTSeq
import pandas as pd

from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner


class Quantifier:
    """Count reads per miRNA across samples."""

    def __init__(self, miRNAs=None, samples=None, slurm=False, outdir=".", upstream=2, downstream=5, weighted=False):
        """Initialize the quantifier and write the hairpin reference.

        Each miRNA gets its own reference entry named by its ID and holding its hairpin,
        so miRNAs that share a hairpin are counted separately.

        Args:
            miRNAs (list): ``(ID, hairpin, mature)`` tuples. Matures not found in their hairpin are skipped.
            samples (dict): Sample name to ``[group, replicate, collapsed reads path]``.
            slurm (bool): Use SLURM for alignment. Defaults to False.
            outdir (str): Output directory. Defaults to the current directory.
            upstream (int): Nucleotides a read may start before the mature. Defaults to 2.
            downstream (int): Nucleotides a read may end after the mature. Defaults to 5.
            weighted (bool): Divide each read count by the number of miRNAs it maps to. Defaults to False.
        """

        self.samples = samples
        self.slurm = slurm
        self.outdir = outdir
        self.upstream = upstream
        self.downstream = downstream
        self.weighted = weighted
        self.windows = {}
        self.sRNA = os.path.join(outdir, "miRNA_hairpins.fa")
        with open(self.sRNA, "w") as fh:
            for name, hairpin, mature in miRNAs or []:
                hairpin = str(hairpin).upper().replace("U", "T")
                mature = str(mature).upper().replace("U", "T")
                pos = hairpin.find(mature)
                if pos < 0 or name in self.windows:
                    continue
                self.windows[name] = (pos, pos + len(mature))
                fh.write(f">{name}\n{hairpin}\n")

    def _align_reads_against_miRNAs(self, cpu, mem, output):
        """Align the sample reads to the miRNA hairpins.

        Args:
            cpu (int): Threads for alignment.
            mem (int): Memory for alignment.
            output (str): Output directory for the index and alignments.

        Returns:
            dict: Sample name to ``[group, replicate, SAM path]``.
        """

        config = str(ir.files("mipyrna").joinpath("param/bowtie_mirna.ini"))
        aln = Bowtie_Aligner(ref_genome=self.sRNA, configFile=config, outdir=output, slurm=self.slurm)
        aln.build_index()

        results = aln.run_alignment(samplesDict=self.samples, outType="SAM", cpu=cpu, mem=mem)

        return results

    def _get_counts(self, aligned_files):
        """Sum aligned read counts per miRNA for each sample.

        Read abundance is taken from the ``x<count>`` suffix of the collapsed read name. A
        read counts for a miRNA when it lies within the window around the mature position.

        Args:
            aligned_files (dict): Sample name to ``[group, replicate, SAM path]``.

        Returns:
            dict: Sample name to ``[group, replicate, counts DataFrame]``.
        """
        all_counts = {}
        for k, sample in aligned_files.items():
            hits = defaultdict(set)
            copies = {}
            for line in HTSeq.SAM_Reader(str(sample[2])):
                if not line.aligned:
                    continue
                mstart, mend = self.windows[line.iv.chrom]
                if line.iv.start >= mstart - self.upstream and line.iv.end <= mend + self.downstream:
                    hits[line.read.name].add(line.iv.chrom)
                    copies[line.read.name] = int(line.read.name.rsplit("_x", 1)[1])

            totals = defaultdict(float)
            for read, names in hits.items():
                share = copies[read] / len(names) if self.weighted else copies[read]
                for name in names:
                    totals[name] += share

            result_df = pd.DataFrame({"mature_id": list(totals), "counts": list(totals.values())})
            all_counts[k] = [sample[0], sample[1], result_df]

        return all_counts

    def quantify_expression(self, cpu=8, mem=20):
        """Align, count and write the merged raw count matrix to ``Raw_read_counts.xlsx``.

        Args:
            cpu (int): Threads for alignment. Defaults to 8.
            mem (int): Memory for alignment. Defaults to 20.

        Returns:
            pandas.DataFrame: Integer counts with a ``mature_id`` column and one column per sample.
        """

        aligdir = mu.make_directory(os.path.join(self.outdir, "quant_align"))

        aligned_results = self._align_reads_against_miRNAs(cpu=cpu, mem=mem, output=aligdir)

        sample_counts = self._get_counts(aligned_files=aligned_results)

        sample_counts_dfs = []

        for k, s in sample_counts.items():
            ucounts = s[2][["mature_id", "counts"]]
            ucounts.columns = ["mature_id", k]

            sample_counts_dfs.append(ucounts)

        final_counts = reduce(
            lambda left, right: pd.merge(left, right, on=["mature_id"], how="outer"),
            sample_counts_dfs,
        ).fillna(0)
        final_counts["mature_id"] = final_counts["mature_id"].astype(str)
        final_counts.iloc[:, 1:] = final_counts.iloc[:, 1:].round().astype(int)
        final_counts.to_excel(f"{self.outdir}/Raw_read_counts.xlsx", index=False)

        return final_counts
