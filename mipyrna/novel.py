"""
Novel miRNA Prediction Module

This module predicts novel miRNAs from clustered small RNA reads. It extracts
candidate precursors around read clusters, folds them with ViennaRNA, scores them
with the bundled machine learning model, estimates randfold significance, finds
star sequences, and separates known from novel miRNAs by BLAST against the
known mature miRNA reference.

Features:
    - Candidate 5-prime and 3-prime precursor extraction around read clusters
    - RNA folding with ViennaRNA and mature position checks via forgi bulge graphs
    - Star sequence estimation and read counting
    - Dicer-consistency check of the reads on each hairpin
    - Randfold-like significance by shuffled-sequence minimum free energies
    - Model-based precursor scoring and a read-based heuristic score
    - 'all' and 'strict' candidate filtering modes
    - Optional parallel scoring of clusters (processes) and samples (threads)
    - Known versus novel classification with BLAST and FASTA/Excel output

Configuration:
    Set through the `Novel_miRNA` constructor: samples dictionary, reference genome,
    read and cluster cutoffs, precursor step, species and species type, model score
    cutoff, filter criteria, randfold settings, star/mature read minimums, random
    seed and worker counts. Writes pooled tables, FASTA files, BLAST results and
    known/novel Excel tables to the output directory.

Dependencies:
    - ViennaRNA (RNA): Folding
    - forgi: Bulge graph analysis
    - NCBI BLAST+: `makeblastdb` and `blastn`
    - pandas / numpy: Data handling
    - mipyrna.model_store, mipyrna.features, mipyrna.reads, mipyrna.knowndb: Scoring, features, reads, reference

Classes / Functions / Exceptions:
    - Novel_miRNA: Predicts and classifies novel miRNAs.
    - _init_worker: Store the predictor in a worker process.
    - _run_task: Score one cluster task in a worker process.

:Author: Naveen Duhan
"""

import contextlib
import copy
import multiprocessing
import os
import re
import shutil
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import forgi.graph.bulge_graph as fgb
import numpy as np
import pandas as pd
import RNA

from mipyrna import model_store as ms
from mipyrna import utility as mu
from mipyrna.features import Precursor_Features
from mipyrna.knowndb import get_known_reference, summarize_sources
from mipyrna.logger import MiPyRNALogger
from mipyrna.reads import Read_process

log = MiPyRNALogger(mode="a", log="novel")

_worker = None

FLANKS = (5, 15, 30, 50)
LOCAL_EXTENSION = 100


def _enclosing_pair(pt, i, j):
    """
    Find the base pair that directly encloses a base pair.

    Args:
        pt (tuple): ViennaRNA pair table (1-based, pt[k] is the partner of k or 0).
        i (int): 5-prime position of the inner pair.
        j (int): 3-prime position of the inner pair.

    Returns:
        tuple or None: (5-prime, 3-prime) of the enclosing pair, or None if the pair is not enclosed.
    """
    k = i - 1
    while k >= 1:
        if pt[k] == 0:
            k -= 1
        elif pt[k] < k:
            k = pt[k] - 1
        else:
            return (k, pt[k])
    return None


def _stem_loop(pt, anchor, max_length):
    """
    Climb from a base pair to the outer pair of the stem-loop that contains it.

    The climb continues through bulges and interior loops and stops below a
    multi-branch loop, the exterior loop or the length limit.

    Args:
        pt (tuple): ViennaRNA pair table (1-based).
        anchor (int): Paired position inside the mature.
        max_length (int): Maximum span of the stem-loop.

    Returns:
        tuple: (5-prime, 3-prime) positions of the outer base pair, 1-based.
    """
    i, j = sorted((anchor, pt[anchor]))
    while True:
        outer = _enclosing_pair(pt, i, j)
        if outer is None or outer[1] - outer[0] + 1 > max_length:
            return i, j
        branches, k = 0, outer[0] + 1
        while k < outer[1]:
            if pt[k] > k:
                branches += 1
                k = pt[k] + 1
            else:
                k += 1
        if branches > 1:
            return i, j
        i, j = outer


def _init_worker(predictor):
    """
    Store the predictor in a module-level variable of a worker process.

    Args:
        predictor (Novel_miRNA): Predictor used by the worker to score clusters.
    """
    global _worker
    _worker = predictor


def _run_task(task):
    """
    Score one cluster task with the predictor of the worker process.

    Args:
        task (tuple): (cluster ID, cluster size, cluster reads DataFrame).

    Returns:
        pandas.DataFrame or None: Candidate precursors of the cluster, or None if none were found.
    """
    assert _worker is not None
    return _worker._score_cluster(task)


def _run_randfold(task):
    """
    Compute the randfold p-value of one candidate with the predictor of the worker process.

    Args:
        task (tuple): (random seed sequence, precursor sequence).

    Returns:
        float: Randfold p-value.
    """
    assert _worker is not None
    seed, precursor = task
    _worker.rng = np.random.default_rng(seed)
    return _worker.estimate_randfold_pvalue(precursor)


class Novel_miRNA:
    """
    Predict novel miRNAs from clustered small RNA reads.
    """

    def __init__(
        self,
        samples=None,
        genome=None,
        read_cutoff=50,
        cluster_cutoff=28,
        max_extension=260,
        species=None,
        score=0.5,
        filter_criteria="all",
        species_type="plants",
        outdir=".",
        randfold_iterations=100,
        randfold_pvalue_cutoff=0.05,
        min_mature_reads=10,
        min_star_reads=1,
        require_star=False,
        max_inconsistent=0.1,
        min_novel_score=1.0,
        random_seed=42,
        workers=1,
        jobs=1,
        resume=False,
        slurm=False,
        mem=None,
    ):
        """
        Initialize the novel miRNA predictor.

        Args:
            samples (dict, optional): Sample name mapped to sample information (name at index 0, reads at 2, clusters at 3).
            genome (str, optional): Reference genome FASTA path.
            read_cutoff (int): Minimum reads of a cluster or mature candidate. Defaults to 50.
            cluster_cutoff (int): Cluster size below which the whole cluster is the mature. Defaults to 28.
            max_extension (int): Nucleotides folded on each side of the mature to find its stem-loop. Defaults to 260.
            species (str, optional): miRBase species code.
            score (float): Minimum model score of a miRNA. Defaults to 0.5.
            filter_criteria (str): 'all' or 'strict'. Defaults to 'all'.
            species_type (str): 'plants' or 'animals'. Defaults to 'plants'.
            outdir (str): Output directory. Defaults to the current directory.
            randfold_iterations (int): Shuffles for the randfold p-value. Defaults to 100.
            randfold_pvalue_cutoff (float): Maximum randfold p-value. Defaults to 0.05.
            min_mature_reads (int): Minimum mature reads in strict mode. Defaults to 10.
            min_star_reads (int): Minimum star reads when a star is required. Defaults to 1.
            require_star (bool): Require star reads in strict mode. Defaults to False.
            max_inconsistent (float): Maximum fraction of the reads on a novel hairpin outside the mature and star positions. Defaults to 0.1.
            min_novel_score (float): Minimum heuristic score of a novel miRNA. Defaults to 1.0.
            random_seed (int): Seed for shuffling. Defaults to 42.
            workers (int): Processes for scoring clusters. Defaults to 1.
            jobs (int): Threads for scoring samples. Defaults to 1.
            resume (bool): Reuse the scored candidates saved for a sample by an earlier run. Defaults to False.
            slurm (bool): Score each sample as an element of a Slurm job array, each with ``workers`` cores. Defaults to False.
            mem (int, optional): Memory in GB of each Slurm job. Defaults to the cluster setting.
        """

        self.genome = genome
        self.readsCutoff = read_cutoff
        self.clusterCutoff = cluster_cutoff
        self.species = species
        self.species_type = species_type
        self.outdir = outdir
        self.samples = samples
        self.max_extension = max_extension
        self.filterCriteria = filter_criteria
        self.score = score
        self.randfold_iterations = randfold_iterations
        self.randfold_pvalue_cutoff = randfold_pvalue_cutoff
        self.min_mature_reads = min_mature_reads
        self.min_star_reads = min_star_reads
        self.require_star = require_star
        self.max_inconsistent = max_inconsistent
        self.min_novel_score = min_novel_score
        self.random_seed = random_seed
        self.workers = workers
        self.jobs = jobs
        self.resume = resume
        self.slurm = slurm
        self.mem = mem
        self.rng = np.random.default_rng(random_seed)

    def get_precursor_structure(self, seq):
        """
        Fold a sequence and return its structure and minimum free energy.

        Args:
            seq (str): RNA sequence.

        Returns:
            tuple: Dot-bracket structure and minimum free energy.

        Raises:
            TypeError, ValueError, RuntimeError: If folding fails.
        """

        try:
            x = RNA.fold(seq)
        except (TypeError, ValueError, RuntimeError) as exc:
            log.warning(f"RNA folding failed: {exc}")
            raise
        return x

    def get_bulge(self, seq):
        """
        Fold a sequence and return its forgi bulge graph.

        Args:
            seq (str): DNA or RNA sequence; T is converted to U.

        Returns:
            forgi.graph.bulge_graph.BulgeGraph: Graph of the folded structure.
        """

        seq = seq.replace("T", "U")
        pstructure, _pmfe = RNA.fold(seq)

        bulge_graph = fgb.BulgeGraph.from_dotbracket(pstructure, seq)

        return bulge_graph

    def close_match(self, preseq, seq):
        """
        Find where seq matches preseq with one mismatch, or -1.

        Args:
            preseq (str): Sequence to search in.
            seq (str): Sequence to find.

        Returns:
            int: Start index of the first one-mismatch match, or -1 if none.
        """

        string_to_match = re.compile(
            "|".join(seq[:i] + "." + seq[i + 1 :] for i in range(len(seq)))
        )
        m = string_to_match.findall(preseq)
        if len(m) == 0:
            return -1
        idx = preseq.find(m[0])
        return idx

    def find_subseq(self, seq, string):
        """
        Locate a mature sequence in a precursor, falling back to a one-mismatch match.

        Args:
            seq (str): Precursor sequence.
            string (str): Mature sequence.

        Returns:
            int: Start index of the mature in the precursor, or -1 if not found.
        """

        for i in range(16, 4, -4):
            c = seq.find(string[:i])
            if c != -1:
                return c
        return self.close_match(seq, string)

    def check_mature(self, precursor, mature, bg=None):
        """
        Check where the mature sequence sits in the precursor structure.

        Args:
            precursor (str): Precursor sequence.
            mature (str): Mature sequence.
            bg (forgi.graph.bulge_graph.BulgeGraph, optional): Bulge graph of the precursor fold; folded here if not given.

        Returns:
            str or None: "ok", or the reason the check fails (e.g. "mature in loop"); None if no bulge graph is available.
        """

        if bg is None:
            bg = self.get_bulge(precursor)
        start = self.find_subseq(precursor, mature) + 1
        end = start + len(mature)
        loops = list(bg.hloop_iterator())
        if len(loops) == 0:
            return "too many loops"
        l = next(iter(bg.define_range_iterator("h0")))
        if (start < l[0] and end > l[0]) or (start < l[1] and end > l[1]):
            return "mature in loop"
        p = next(iter(bg.stem_bp_iterator("s0")))

        if start == 0:
            return "mature not found"
        if end > p[1] + 2:
            return "3p mature outside"
        if p[0] > start + 2:
            return "mature unpaired"
        return "ok"

    def find_star_sequence(self, precursor, mature, bg=None):
        """
        Estimate the star sequence from a given mature and precursor.

        Args:
            precursor (str): Precursor sequence.
            mature (str): Mature sequence.
            bg (forgi.graph.bulge_graph.BulgeGraph, optional): Bulge graph of the precursor fold; folded here if not given.

        Returns:
            str or None: Estimated star sequence, or None if it cannot be determined.
        """

        if bg is None:
            bg = self.get_bulge(precursor)
        start = self.find_subseq(precursor, mature) + 1

        for p in bg.adjacent_stem_pairs_iterator():
            pp = [p[0], p[2]]
            for s in pp:
                xx = list(bg.stem_bp_iterator(s))
                for i in xx:
                    if start == i[0]:
                        star = precursor[i[1] - len(mature) + 2 : i[1] + 2]
                        return star
                    elif start == i[1]:
                        star = precursor[i[0] - len(mature) + 1 : i[0] + 2]
                        return star
        return

    def estimate_randfold_pvalue(self, precursor):
        """
        Approximate randfold significance via shuffled-sequence minimum free energies.

        Args:
            precursor (str): Precursor sequence.

        Returns:
            float: Smoothed p-value, the share of shuffles with an MFE as low as the observed one; 1.0 for an empty precursor.
        """
        if not precursor:
            return 1.0
        seq = precursor.replace("T", "U")
        _, observed_mfe = RNA.fold(seq)
        better_or_equal = 0
        for _ in range(self.randfold_iterations):
            shuffled = "".join(self.rng.permutation(list(seq)))
            _, shuffled_mfe = RNA.fold(shuffled)
            if shuffled_mfe <= observed_mfe:
                better_or_equal += 1
        # +1 smoothing avoids zero p-values on small permutations.
        pvalue = (better_or_equal + 1) / (self.randfold_iterations + 1)
        return float(pvalue)

    def heuristic_score(self, row):
        """
        Calculate a heuristic score from read, star, structure and randfold evidence.

        Without star reads the score is lowered and the mature read term is dropped, so a
        candidate needs the star evidence of Dicer processing to score high.

        Args:
            row (pandas.Series): Candidate with mature_reads, star_reads, mature_check, mfe, randfold_pvalue and star values.

        Returns:
            float: Score rounded to four decimals.
        """
        score = 0.0
        star_reads = float(row.get("star_reads", 0))
        if star_reads > 0:
            score += 3.9 + np.log2(float(row.get("mature_reads", 0)) + 1.0) + 0.5 * np.log2(star_reads + 1.0)
        else:
            score -= 1.3
        if row.get("mature_check") == "ok":
            score += 1.0
        mfe = float(row.get("mfe", 0.0))
        if mfe < 0:
            score += min(3.0, abs(mfe) / 20.0)
        if float(row.get("randfold_pvalue", 1.0)) <= self.randfold_pvalue_cutoff:
            score += 1.5
        if isinstance(row.get("star"), str) and len(row.get("star")) > 0:
            score += 0.5
        return round(score, 4)

    def predict_precursor(self, df):
        """
        Score precursor features with the bundled model.

        Args:
            df (pandas.DataFrame): Candidates with a "features" column.

        Returns:
            pandas.DataFrame: Columns "miRNA" and "Non-miRNA" with the probability of each class.
        """

        prob = ms.predict(df["features"].tolist(), self.species_type)

        return pd.DataFrame({"miRNA": prob, "Non-miRNA": 1.0 - prob})

    def get_putative_precursor(self, filtered_cluster=None, other=None):
        """
        Extract and score the candidate hairpins of one read cluster.

        The mature and its flanks are folded once. The stem-loop that holds the
        mature is cut out, refolded and checked, and the model features are only
        built for hairpins with a good mature position. Randfold p-values and the
        heuristic score are added later, once the model has scored the candidates.

        Args:
            filtered_cluster (pandas.DataFrame): Reads of the mature cluster.
            other (pandas.DataFrame, optional): Remaining cluster reads used to count star reads.

        Returns:
            pandas.DataFrame or None: Candidate hairpins, or None if no valid candidate was found.

        Raises:
            ValueError: If filtered_cluster is not provided.
        """

        if filtered_cluster is None:
            raise ValueError("filtered_cluster is required")

        x = filtered_cluster.iloc[0]
        chrom, start, end, strand = x["name"], x.start, x.end, x.strand
        reader = Read_process(reference=self.genome)
        refseq = reader.get_sequence((chrom, max(0, start - 10), end + 10, strand))
        mature = reader.get_consensus_read(
            ref_precursor=refseq, aligned_cluster=filtered_cluster
        )
        if mature is None:
            mature = x.seq
        seqlen = len(mature)
        maturecounts = filtered_cluster.reads.sum()
        putative_precursor = []
        seen = set()

        for extension in sorted(
            {min(LOCAL_EXTENSION, self.max_extension), self.max_extension}
        ):
            max_length = extension + seqlen
            wstart = max(0, start - extension)
            window = reader.get_sequence([chrom, wstart, end + extension, strand])
            wend = wstart + len(window)
            offset = start - wstart if strand == "+" else wend - end
            pos = window.find(mature, max(0, offset - 5))
            if pos < 0:
                pos = self.find_subseq(window, mature)
            if pos < 0:
                continue

            try:
                md = RNA.md()
                md.max_bp_span = max_length
                wide, _ = RNA.fold_compound(window.replace("T", "U"), md).mfe()
            except (TypeError, ValueError, RuntimeError) as exc:
                log.warning(f"Skipped cluster at {chrom}:{start}: {exc}")
                continue
            pt = RNA.ptable(wide)

            centre = pos + seqlen // 2 + 1
            paired = [k for k in range(pos + 1, pos + seqlen + 1) if pt[k] != 0]
            if not paired:
                continue
            i, j = _stem_loop(
                pt, min(paired, key=lambda k: abs(k - centre)), max_length
            )

            for flank in FLANKS:
                lo = max(0, min(i - 1, pos) - flank)
                hi = min(len(window), max(j, pos + seqlen) + flank)
                precursor_seq = window[lo:hi]
                if precursor_seq in seen or "N" in precursor_seq:
                    continue
                seen.add(precursor_seq)
                try:
                    pre = precursor_seq.replace("T", "U")
                    struct, sc = RNA.fold(pre)
                    bg = fgb.BulgeGraph.from_dotbracket(struct, pre)
                    if self.check_mature(precursor_seq, mature, bg) != "ok":
                        continue
                except (TypeError, ValueError, RuntimeError) as exc:
                    log.warning(
                        f"Skipped precursor candidate at {chrom}:{start}: {exc}"
                    )
                    continue

                m = Precursor_Features(precursor_seq).generate_features()
                if m is None:
                    continue

                if strand == "+":
                    pstart, pend = wstart + lo, wstart + hi
                else:
                    pstart, pend = wend - hi, wend - lo
                star = self.find_star_sequence(precursor_seq, mature, bg)
                starcount = 0
                inconsistent = 0.0
                if other is not None and star is not None:
                    s = self.find_subseq(precursor_seq, star)
                    if s >= 0:
                        ss = pstart + s if strand == "+" else pend - s - len(star)
                        sreads = other[
                            (other.start >= ss - 2) & (other.end <= ss + len(star) + 3)
                        ]
                        starcount = sreads.reads.sum()
                if other is not None:
                    on_hairpin = other[(other.start >= pstart) & (other.end <= pend)].reads.sum()
                    inconsistent = max(0, on_hairpin - starcount) / (maturecounts + on_hairpin)
                putative_precursor.append(
                    {
                        "chrom": chrom,
                        "mature": mature,
                        "mature_length": seqlen,
                        "precursor": precursor_seq,
                        "precursor_length": len(precursor_seq),
                        "start": pstart,
                        "end": pend,
                        "mature_start": start,
                        "strand": strand,
                        "struct": struct,
                        "mfe": sc,
                        "mature_check": "ok",
                        "features": m,
                        "mature_reads": maturecounts,
                        "position": "5-prime" if centre < (i + j) / 2 else "3-prime",
                        "star_reads": starcount,
                        "inconsistent_fraction": inconsistent,
                        "star": star,
                        "randfold_pvalue": 1.0,
                        "cluster": x.cluster,
                    }
                )

        if len(putative_precursor) == 0:
            return
        return pd.DataFrame(putative_precursor)

    def filter_miRNA(self, df):
        """
        Keep the best-scoring precursor per mature sequence and assign novel IDs.

        Args:
            df (pandas.DataFrame): Scored candidates with "mature" and "score" columns.

        Returns:
            pandas.DataFrame: One row per mature sequence, with a Novel_ID column.
        """
        # Keep the single best-scoring precursor per mature sequence.
        idx = df.groupby("mature")["score"].idxmax().tolist()
        final_miRNAs = df.loc[idx].reset_index(drop=True)

        with contextlib.suppress(Exception):
            final_miRNAs["Novel_ID"] = [
                f"{self.species}_mipyrna_{final_miRNAs.iloc[i]['chrom']}_miR_{i + 1}"
                for i in range(len(final_miRNAs))
            ]

        return final_miRNAs

    def apply_strict_filters(self, df):
        """
        Apply the structure, score, read count and randfold filters.

        Args:
            df (pandas.DataFrame): Candidate miRNAs.

        Returns:
            pandas.DataFrame: Candidates passing all filters (and the star read minimum if required).
        """
        mask = (
            (df["mature_check"] == "ok")
            & (df["score"] >= self.score)
            & (df["mature_reads"] >= self.min_mature_reads)
            & (df["randfold_pvalue"] <= self.randfold_pvalue_cutoff)
        )
        if self.require_star:
            mask = mask & (df["star_reads"] >= self.min_star_reads)
        return df[mask].copy()

    @staticmethod
    def _overlaps(row, loci):
        """Check whether a hairpin overlaps any locus on the same chromosome and strand.

        Args:
            row (pandas.Series): Candidate with chrom, strand, start and end.
            loci (dict): ``(chrom, strand)`` to a list of ``(start, end)`` tuples.

        Returns:
            bool: True if the hairpin overlaps a locus.
        """
        return any(
            a <= row["end"] and b >= row["start"]
            for a, b in loci.get((row["chrom"], row["strand"]), [])
        )

    def _drop_known_loci(self, novel_df, known_df):
        """Remove novel candidates whose hairpin overlaps a known miRNA hairpin.

        Their reads are counted under the known miRNA, so they are not novel.

        Args:
            novel_df (pandas.DataFrame): Novel candidates.
            known_df (pandas.DataFrame): Known miRNAs.

        Returns:
            pandas.DataFrame: Novel candidates away from known miRNA loci.
        """
        loci = defaultdict(list)
        for r in known_df.itertuples():
            loci[(r.chrom, r.strand)].append((r.start, r.end))
        keep = [not self._overlaps(r, loci) for _i, r in novel_df.iterrows()]
        log.info(
            f"Removed {len(keep) - sum(keep)} novel candidates overlapping known miRNA hairpins"
        )
        return novel_df[keep]

    def _dicer_filter(self, df):
        """Remove novel candidates whose reads do not look like Dicer products.

        A candidate is dropped when more than ``max_inconsistent`` of the reads on its
        hairpin fall outside the mature and star positions, or when its heuristic score
        is below ``min_novel_score``.

        Args:
            df (pandas.DataFrame): Novel candidates.

        Returns:
            pandas.DataFrame: Candidates with Dicer-like reads.
        """
        keep = (df["inconsistent_fraction"] <= self.max_inconsistent) & (df["heuristic_score"] >= self.min_novel_score)
        log.info(f"Removed {int((~keep).sum())} novel candidates with inconsistent reads or a low heuristic score")
        return df[keep]

    def _one_per_locus(self, df):
        """Keep one candidate for each set of overlapping hairpins.

        Candidates with star reads come first, then those with the most mature reads.

        Args:
            df (pandas.DataFrame): Novel candidates.

        Returns:
            pandas.DataFrame: Candidates with non-overlapping hairpins.
        """
        ranked = df.assign(_has_star=df["star_reads"] > 0).sort_values(
            ["_has_star", "mature_reads", "score"], ascending=False
        )
        loci = defaultdict(list)
        keep = []
        for i, r in ranked.iterrows():
            if self._overlaps(r, loci):
                continue
            loci[(r["chrom"], r["strand"])].append((r["start"], r["end"]))
            keep.append(i)
        log.info(f"Kept {len(keep)} of {len(df)} novel candidates, one per locus")
        return df.loc[sorted(keep)]

    def _write_fasta(self, df, output, type="known", seq="mature"):
        """
        Write mature or hairpin sequences of known or novel miRNAs to a FASTA file.

        Args:
            df (pandas.DataFrame): miRNA table.
            output (str): Output FASTA path.
            type (str): "known" uses Known_ID (written once each), anything else uses Novel_ID. Defaults to "known".
            seq (str): "mature" or "hairpin". Defaults to "mature".

        Returns:
            str: The output path.
        """
        mature_seq = []
        with open(output, "w") as fp:
            for _i, r in df.iterrows():
                if seq == "mature":
                    if type == "known":
                        if r["Known_ID"] in mature_seq:
                            continue
                        else:
                            fp.write(f">{r['Known_ID']}\n{r['mature']}\n")
                            mature_seq.append(r["Known_ID"])
                    else:
                        fp.write(f">{r['Novel_ID']}\n{r['mature']}\n")
                if seq == "hairpin":
                    if type == "known":
                        if r["Known_ID"] in mature_seq:
                            continue
                        else:
                            fp.write(f">{r['Known_ID']}\n{r['precursor']}\n")
                            mature_seq.append(r["Known_ID"])
                    else:
                        fp.write(f">{r['Novel_ID']}\n{r['precursor']}\n")
        return output

    def predict_known_miRNAs(self, novel_df):
        """
        Split predictions into known and novel miRNAs by BLAST against the known mature miRNA reference.

        Writes known and novel Excel tables, FASTA files and a database summary. Each
        known miRNA keeps its best match in ``Known_ID`` and all its matches, best
        first, in ``Known_IDs``.
        Exits the program if BLAST+ is not in the path.

        Args:
            novel_df (pandas.DataFrame): Filtered candidate miRNAs with a Novel_ID column.

        Returns:
            tuple: Known miRNAs table, novel miRNAs table, known FASTA path and novel FASTA path.
        """

        blast = mu.make_directory(os.path.join(self.outdir, "blast_results"))
        novel_mature = os.path.join(self.outdir, "all_sample_pooled_mature.fa")

        mirbase_data = get_known_reference(self.species, self.species_type, blast)

        execPATH = shutil.which("makeblastdb")

        if execPATH is None:
            log.error("ncbi-blast not found in path")

            sys.exit(1)
        else:
            blastdb_command = f"makeblastdb -in {mirbase_data['mature']} -dbtype 'nucl'"

        mu.local_run(
            program="makeblastdb",
            command=blastdb_command,
            outdir=blast,
            message=f"building database for {mirbase_data['mature']}",
        )

        execPATH = shutil.which("blastn")

        if execPATH is None:
            log.error("ncbi-blast+ not found in path")

            sys.exit(1)
        else:
            blast_out = os.path.join(blast, "blast_mature_vs_known.txt")
            blastn_command = f"blastn -db {mirbase_data['mature']} -query {novel_mature} -evalue 100 -out {blast_out} -outfmt 6 -word_size 4"
            mu.local_run(
                program="blastn",
                command=blastn_command,
                outdir=blast,
                message="similarity search against known miRNAs",
            )

        blast_results = pd.read_csv(
            blast_out,
            sep="\t",
            names=[
                "Novel_ID",
                "Known_ID",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ],
        )

        # A mature may differ from its miRBase entry by one base at the 5-prime end,
        # so allow a one-base offset while keeping at least 18 aligned bases. When a
        # candidate matches several entries, keep the best alignment.
        hits = pd.DataFrame(
            blast_results[
                (blast_results.qstart < 3)
                & (blast_results.length > 17)
                & (blast_results.pident > 94)
            ]
        ).sort_values(by=["length", "pident", "bitscore"], ascending=False)

        filterd_results = hits.drop_duplicates(subset=["Novel_ID"])

        all_matches = (
            hits.drop_duplicates(subset=["Novel_ID", "Known_ID"])
            .groupby("Novel_ID", sort=False)["Known_ID"]
            .agg(";".join)
            .to_frame("Known_IDs")
            .reset_index()
        )

        known_miRNAs_dataframe = filterd_results[["Novel_ID", "Known_ID"]].merge(
            all_matches, on="Novel_ID"
        )

        known_miRNAs = known_miRNAs_dataframe.merge(
            novel_df, on=["Novel_ID"], indicator=True
        )

        known_miRNAs = known_miRNAs.drop("_merge", axis=1)
        known_miRNAs["Database"] = known_miRNAs["Known_ID"].map(mirbase_data["sources"])
        if len(known_miRNAs):
            summary = summarize_sources(known_miRNAs["Database"])
            summary.to_csv(
                f"{self.outdir}/known_database_summary.tsv", sep="\t", index=False
            )
            log.info(
                f"Known miRNAs by database: {dict(zip(summary['database'], summary['miRNAs']))}"
            )

        novel_miRNAs = (
            pd.merge(
                known_miRNAs["Novel_ID"],
                novel_df,
                on=["Novel_ID"],
                how="outer",
                indicator=True,
            )
            .query('_merge == "right_only"')
            .drop("_merge", axis=1)
        )

        novel_miRNAs = self._one_per_locus(
            self._drop_known_loci(self._dicer_filter(novel_miRNAs), known_miRNAs)
        )

        try:
            known_miRNAs["raw_ID"] = known_miRNAs["Novel_ID"]
            known_miRNAs["Novel_ID"] = [
                f"{self.species}_miPyRNA_{known_miRNAs.iloc[i]['chrom']}_miR_{i + 1}"
                for i in range(len(known_miRNAs))
            ]
        except Exception as exc:
            log.warning(f"Could not rename known miRNA IDs: {exc}")

        try:
            novel_miRNAs["raw_ID"] = novel_miRNAs["Novel_ID"]
            novel_miRNAs["Novel_ID"] = [
                f"{self.species}_miPyRNA_{novel_miRNAs.iloc[i]['chrom']}_miR_{i + 1}"
                for i in range(len(novel_miRNAs))
            ]
        except Exception as exc:
            log.warning(f"Could not rename novel miRNA IDs: {exc}")

        known_miRNAs.to_excel(f"{self.outdir}/known_miRNAs.xlsx", index=False)
        novel_miRNAs.to_excel(f"{self.outdir}/novel_miRNAs.xlsx", index=False)

        novel_fasta = self._write_fasta(
            df=novel_miRNAs,
            output=os.path.join(f"{self.outdir}/novel_mature.fa"),
            type="novel",
        )
        known_fasta = self._write_fasta(
            df=known_miRNAs,
            output=os.path.join(f"{self.outdir}/known_mature.fa"),
            type="known",
        )

        return known_miRNAs, novel_miRNAs, known_fasta, novel_fasta

    def _score_cluster(self, task):
        """
        Extract and score the candidate precursors of one read cluster.

        Args:
            task (tuple): (cluster ID, cluster size, cluster reads DataFrame).

        Returns:
            pandas.DataFrame or None: Candidate precursors, or None if the cluster has too few reads.
        """
        cluster_id, clust_size, df = task
        self.rng = np.random.default_rng([self.random_seed, int(cluster_id)])
        df = df.sort_values("reads", ascending=False)
        if clust_size < self.clusterCutoff:
            df["mature"] = True
            return self.get_putative_precursor(filtered_cluster=df)
        anchor = df.iloc[0]
        mm = df.loc[
            (abs(df.start - anchor.start) <= 3) & (abs(df.end - anchor.end) <= 5)
        ].copy()
        if mm.reads.sum() <= self.readsCutoff:
            return None
        mm["mature"] = True
        other = df.loc[~df.index.isin(mm.index)].copy()
        other["mature"] = False
        return self.get_putative_precursor(filtered_cluster=mm, other=other)

    def _score_sample(self, sample):
        """
        Score the candidate precursors of all read clusters of one sample.

        Args:
            sample (list): Sample information (name at index 0, reads at 2, clusters at 3).

        Returns:
            list: Non-empty candidate precursor DataFrames, one per scored cluster.
        """
        clusts = sample[3][sample[3].reads >= self.readsCutoff]
        groups = dict(tuple(sample[2].groupby("cluster")))
        tasks = [
            (c.cluster, c.clust_size, groups[c.cluster])
            for _i, c in clusts.iterrows()
            if c.cluster in groups
        ]
        log.info(
            f"Scoring {len(tasks)} clusters of sample {sample[0]} with {self.workers} worker(s)"
        )
        return [
            c
            for c in self._score_clusters(tasks, sample[0])
            if c is not None and len(c) > 0
        ]

    def _add_randfold(self, df, name):
        """
        Compute randfold p-values for the candidates that the model and basic checks keep.

        Folding the shuffles dominates the run time, so it is limited to candidates with
        a good mature position, enough mature reads and a model score of at least the
        cutoff. The others keep a p-value of 1.0. Each candidate has its own random
        stream, so its p-value does not depend on which other candidates were folded.

        Args:
            df (pandas.DataFrame): Scored candidates; the ``randfold_pvalue`` column is updated in place.
            name (str): Sample name used in the log messages.
        """
        keep = (
            (df["mature_check"] == "ok")
            & (df["mature_reads"] >= self.min_mature_reads)
            & (df["score"] >= self.score)
        )
        log.info(
            f"Sample {name}: {len(df)} candidates, {int((df['mature_check'] == 'ok').sum())} with a good mature position, "
            f"{int(keep.sum())} also pass the model cutoff of {self.score}"
        )
        df["randfold_pvalue"] = 1.0
        idx = df.index[keep]
        tasks = [
            (
                [self.random_seed, int(df.at[i, "cluster"]), int(i)],
                df.at[i, "precursor"],
            )
            for i in idx
        ]
        if tasks:
            if self.workers * self.jobs <= 1 or len(tasks) < 2:
                light = copy.copy(self)
                light.samples = None
                _init_worker(light)
                pvalues = [_run_randfold(t) for t in tasks]
            else:
                light = copy.copy(self)
                light.samples = None
                with multiprocessing.Pool(
                    self.workers * self.jobs, _init_worker, (light,)
                ) as pool:
                    pvalues = pool.map(
                        _run_randfold,
                        tasks,
                        chunksize=max(
                            1, min(10, len(tasks) // (self.workers * self.jobs * 4))
                        ),
                    )
            df.loc[idx, "randfold_pvalue"] = pvalues
        log.info(
            f"Sample {name}: {int((df['randfold_pvalue'] <= self.randfold_pvalue_cutoff).sum())} of {len(tasks)} "
            f"candidates pass randfold (p <= {self.randfold_pvalue_cutoff})"
        )

    def _score_clusters(self, tasks, name="Clusters"):
        """
        Score clusters serially or in a process pool; results keep the input order.

        Progress is logged about every 5% of the clusters.

        Args:
            tasks (list): Cluster tasks of (cluster ID, cluster size, reads DataFrame).
            name (str): Sample name used in the progress messages. Defaults to "Clusters".

        Returns:
            list: One result per task, a DataFrame or None.
        """
        total = len(tasks)
        step = max(1, total // 20)
        results = []

        def collect(iterator):
            for result in iterator:
                results.append(result)
                if len(results) % step == 0 or len(results) == total:
                    log.info(f"Scored {len(results)}/{total} clusters of sample {name}")

        if self.workers <= 1 or total < 2:
            collect(copy.copy(self)._score_cluster(t) for t in tasks)
            return results
        light = copy.copy(self)
        light.samples = None
        chunk = max(1, min(10, total // (self.workers * 8)))
        with multiprocessing.Pool(self.workers, _init_worker, (light,)) as pool:
            collect(pool.imap(_run_task, tasks, chunksize=chunk))
        return results

    def _cache_path(self, name):
        """Path of the scored candidates saved for a sample.

        Args:
            name (str): Sample name.

        Returns:
            str: Pickle path in the ``temp`` directory.
        """
        temp = os.path.join(self.outdir, "temp")
        os.makedirs(temp, exist_ok=True)
        return os.path.join(temp, f"{name}_scored.pkl")

    def _finish_sample(self, name, candidates):
        """Score the candidates of a sample with the model and randfold, and save them.

        Args:
            name (str): Sample name.
            candidates (list): Candidate precursor DataFrames of the sample.
        """
        final_data = pd.concat(candidates, ignore_index=True)
        final_data["score"] = self.predict_precursor(final_data)["miRNA"]
        final_data["ml_score"] = final_data["score"]
        self._add_randfold(final_data, name)
        final_data["heuristic_score"] = final_data.apply(self.heuristic_score, axis=1)
        cache = self._cache_path(name)
        final_data.to_csv(f"{os.path.dirname(cache)}/{name}_raw_miRNA.txt", sep="\t", index=False)
        final_data.to_pickle(cache + ".tmp")
        os.replace(cache + ".tmp", cache)

    def _score_one(self, sample):
        """Score one sample from start to finish; run by one element of a Slurm job array.

        Args:
            sample (list): Sample information (name at index 0, reads at 2, clusters at 3).
        """
        candidates = self._score_sample(sample)
        if candidates:
            self._finish_sample(sample[0], candidates)

    def _score_on_slurm(self, samples):
        """Score samples as one Slurm job array, one element per sample.

        Args:
            samples (list): Samples to score.
        """
        light = copy.copy(self)
        light.samples = None
        light.jobs = 1
        light.slurm = False
        mu.run_array(
            [(light, "_score_one", (sample,)) for sample in samples],
            outdir=self.outdir,
            name="novel_scoring",
            cpu=self.workers,
            mem=self.mem,
            limit=self.jobs,
        )

    def get_novel_miRNA(self):
        """
        Run novel miRNA prediction for all samples and return the final tables.

        Returns:
            tuple: Known miRNAs table, novel miRNAs table, known FASTA path and novel FASTA path; empty tables and None paths if nothing was found.

        Raises:
            ValueError: If samples is not provided.
        """

        if self.samples is None:
            raise ValueError("samples is required")

        outputs = []

        mu.make_directory(os.path.join(self.outdir, "temp"))

        samples = list(self.samples.values())
        cache = {s[0]: self._cache_path(s[0]) for s in samples}
        todo = [s for s in samples if not (self.resume and os.path.exists(cache[s[0]]))]
        if self.slurm and todo:
            self._score_on_slurm(todo)
            scored = {}
        else:
            with ThreadPoolExecutor(max_workers=max(1, self.jobs)) as pool:
                scored = dict(zip((s[0] for s in todo), pool.map(self._score_sample, todo), strict=True))

        for sample in samples:
            if sample[0] in scored:
                N = scored[sample[0]]
                if len(N) == 0:
                    log.warning(f"No precursor candidates produced for sample {sample[0]}")
                    continue
                self._finish_sample(sample[0], N)
            elif not os.path.exists(cache[sample[0]]):
                log.warning(f"No precursor candidates produced for sample {sample[0]}")
                continue
            else:
                log.info(f"Loading scored candidates of sample {sample[0]}")
            outputs.append(pd.read_pickle(cache[sample[0]]))

        if len(outputs) == 0:
            log.warning("No novel miRNA candidates were detected")
            return pd.DataFrame(), pd.DataFrame(), None, None

        final = pd.concat(outputs, ignore_index=True)

        final_miRNA = self.filter_miRNA(final)

        final_miRNA.to_csv(
            os.path.join(self.outdir, "all_sample_pooled_miRNAs.txt"),
            sep="\t",
            index=False,
        )

        if self.filterCriteria == "all":
            final_miRNA = final_miRNA[final_miRNA["score"] >= self.score]
            final_miRNA.to_csv(
                os.path.join(self.outdir, "all_sample_filtered_miRNAs.txt"),
                sep="\t",
                index=False,
            )
            mature_out = os.path.join(self.outdir, "all_sample_pooled_mature.fa")
            hairpin_out = os.path.join(self.outdir, "all_sample_pooled_hairpin.fa")
            self._write_fasta(final_miRNA, mature_out, type="raw", seq="mature")
            self._write_fasta(final_miRNA, hairpin_out, type="raw", seq="hairpin")

        if self.filterCriteria == "strict":
            final_miRNA = self.apply_strict_filters(final_miRNA)
            final_miRNA.to_csv(
                os.path.join(self.outdir, "all_sample_filtered_miRNAs.txt"),
                sep="\t",
                index=False,
            )
            mature_out = os.path.join(self.outdir, "all_sample_pooled_mature.fa")
            hairpin_out = os.path.join(self.outdir, "all_sample_pooled_hairpin.fa")
            self._write_fasta(final_miRNA, mature_out, type="raw", seq="mature")
            self._write_fasta(final_miRNA, hairpin_out, type="raw", seq="hairpin")

        if final_miRNA.empty:
            log.warning("No candidates left after filtering")
            return pd.DataFrame(), pd.DataFrame(), None, None

        known_miRNAs, novel_miRNAs, known_fasta, novel_fasta = (
            self.predict_known_miRNAs(novel_df=final_miRNA)
        )

        return known_miRNAs, novel_miRNAs, known_fasta, novel_fasta
