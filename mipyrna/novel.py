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
    - Randfold-like significance by shuffled-sequence minimum free energies
    - Model-based precursor scoring and a miRDeep2-style heuristic score
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
        precursor_step=5,
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
        random_seed=42,
        workers=1,
        jobs=1,
        resume=False,
    ):
        """
        Initialize the novel miRNA predictor.

        Args:
            samples (dict, optional): Sample name mapped to sample information (name at index 0, reads at 2, clusters at 3).
            genome (str, optional): Reference genome FASTA path.
            read_cutoff (int): Minimum reads of a cluster or mature candidate. Defaults to 50.
            cluster_cutoff (int): Cluster size below which the whole cluster is the mature. Defaults to 28.
            precursor_step (int): Step for extending candidate precursor length. Defaults to 5.
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
            random_seed (int): Seed for shuffling. Defaults to 42.
            workers (int): Processes for scoring clusters. Defaults to 1.
            jobs (int): Threads for scoring samples. Defaults to 1.
            resume (bool): Reuse the scored candidates saved for a sample by an earlier run. Defaults to False.
        """

        self.genome = genome
        self.readsCutoff = read_cutoff
        self.clusterCutoff = cluster_cutoff
        self.species = species
        self.species_type = species_type
        self.outdir = outdir
        self.samples = samples
        self.precursorStep = precursor_step
        self.filterCriteria = filter_criteria
        self.score = score
        self.randfold_iterations = randfold_iterations
        self.randfold_pvalue_cutoff = randfold_pvalue_cutoff
        self.min_mature_reads = min_mature_reads
        self.min_star_reads = min_star_reads
        self.require_star = require_star
        self.random_seed = random_seed
        self.workers = workers
        self.jobs = jobs
        self.resume = resume
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

    def check_mature(self, precursor, mature):
        """
        Check where the mature sequence sits in the precursor structure.

        Args:
            precursor (str): Precursor sequence.
            mature (str): Mature sequence.

        Returns:
            str or None: "ok", or the reason the check fails (e.g. "mature in loop"); None if no bulge graph is available.
        """

        bg = self.get_bulge(precursor)
        if bg is None:
            return
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

    def find_star_sequence(self, precursor, mature):
        """
        Estimate the star sequence from a given mature and precursor.

        Args:
            precursor (str): Precursor sequence.
            mature (str): Mature sequence.

        Returns:
            str or None: Estimated star sequence, or None if it cannot be determined.
        """

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

    def mirdeep2_like_score(self, row):
        """
        Calculate a heuristic miRDeep2-style score from read, structure and randfold evidence.

        Args:
            row (pandas.Series): Candidate with mature_reads, star_reads, mature_check, mfe, randfold_pvalue and star values.

        Returns:
            float: Score rounded to four decimals.
        """
        score = 0.0
        score += np.log2(float(row.get("mature_reads", 0)) + 1.0)
        score += 0.5 * np.log2(float(row.get("star_reads", 0)) + 1.0)
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

    def get_putative_precursor(
        self,
        filtered_cluster=None,
        other=None,
        seqlength=20,
        loop=20,
        precursor_step=10,
    ):
        """
        Extract and score candidate precursors around each read cluster.

        Candidates extend 5-prime and 3-prime of the mature, are folded and checked,
        and receive star reads, randfold p-values and a miRDeep2-style score.

        Args:
            filtered_cluster (pandas.DataFrame): Reads of the mature cluster.
            other (pandas.DataFrame, optional): Remaining cluster reads used to count star reads.
            seqlength (int): Unused. Defaults to 20.
            loop (int): Assumed loop length. Defaults to 20.
            precursor_step (int): Step for extending candidate length. Defaults to 10.

        Returns:
            pandas.DataFrame or None: Candidate precursors, or None if no valid candidate was found.

        Raises:
            ValueError: If filtered_cluster is not provided.
        """

        if filtered_cluster is None:
            raise ValueError("filtered_cluster is required")

        putative_precursor = []

        x = filtered_cluster.iloc[0]
        ref_coords = (x["name"], max(0, x.start - 10), x.end + 10, x.strand)
        # first get a sequence from reference genome based on coordinates
        refseq = Read_process(reference=self.genome).get_sequence(ref_coords)
        # check whether any sequencing reads is in sequence above
        mature_mirna = Read_process(reference=self.genome).get_consensus_read(
            ref_precursor=refseq, aligned_cluster=filtered_cluster
        )
        coords = (x["name"], x.start, x.end, x.strand)
        chrom, start, end, strand = coords
        if mature_mirna is None:
            mature_mirna = x.seq

        seqlen = len(mature_mirna)

        # get 5' precursor

        for i in range(10, 60, precursor_step):
            start5 = max(0, start - 1)
            end5 = start + 2 * seqlen - 1 + loop + i
            if end5 <= start5:
                continue
            coords = [chrom, start5, end5, strand]

            precursor_seq = Read_process(reference=self.genome).get_sequence(coords)

            if precursor_seq == None or "N" in precursor_seq:
                continue
            try:
                struct, sc = self.get_precursor_structure(precursor_seq)

                mature_status = self.check_mature(precursor_seq, mature_mirna)

                if mature_status is None:
                    continue
            except (TypeError, ValueError, RuntimeError) as exc:
                log.warning(f"Skipped precursor candidate at {chrom}:{start}: {exc}")
                continue

            k = Precursor_Features(precursor_seq)

            m = k.generate_features()

            if m == None:
                continue

            maturecounts = filtered_cluster.reads.sum()
            putative_precursor.append(
                {
                    "chrom": chrom,
                    "mature": mature_mirna,
                    "mature_length": len(mature_mirna),
                    "precursor": precursor_seq,
                    "precursor_length": len(precursor_seq),
                    "start": start5,
                    "end": end5,
                    "mature_start": start,
                    "strand": strand,
                    "struct": struct,
                    "mfe": sc,
                    "mature_check": mature_status,
                    "features": m,
                    "mature_reads": maturecounts,
                    "position": "5-prime",
                }
            )

        for i in range(10, 60, precursor_step):
            start3 = max(0, start - (loop + seqlen + i))
            end3 = end + seqlen + 1
            if end3 <= start3:
                continue
            coords = [chrom, start3, end3, strand]

            precursor_seq = Read_process(reference=self.genome).get_sequence(coords)

            if precursor_seq == None or "N" in precursor_seq:
                continue
            try:
                struct, sc = self.get_precursor_structure(precursor_seq)

                mature_status = self.check_mature(precursor_seq, mature_mirna)

                if mature_status is None:
                    continue
            except (TypeError, ValueError, RuntimeError) as exc:
                log.warning(f"Skipped precursor candidate at {chrom}:{start}: {exc}")
                continue

            k = Precursor_Features(sequence=precursor_seq)

            m = k.generate_features()

            if m == None:
                continue
            maturecounts = filtered_cluster.reads.sum()

            putative_precursor.append(
                {
                    "chrom": chrom,
                    "mature": mature_mirna,
                    "mature_length": len(mature_mirna),
                    "precursor": precursor_seq,
                    "precursor_length": len(precursor_seq),
                    "start": start3,
                    "end": end3,
                    "mature_start": start,
                    "strand": strand,
                    "struct": struct,
                    "mfe": sc,
                    "mature_check": mature_status,
                    "features": m,
                    "mature_reads": maturecounts,
                    "position": "3-prime",
                }
            )
        putative_precursor = pd.DataFrame(putative_precursor)

        if len(putative_precursor) == 0:
            return

        stars = []
        star_counts = []
        pvalues = []
        for _, row in putative_precursor.iterrows():
            star = self.find_star_sequence(row.precursor, row.mature)
            stars.append(star)
            starcount = 0
            if other is not None and star is not None:
                s = self.find_subseq(row.precursor, star)
                if s >= 0:
                    ss = row.start + s
                    se = ss + len(star)
                    sreads = other[(other.start >= ss - 2) & (other.end <= se + 3)]
                    starcount = sreads.reads.sum()
            star_counts.append(starcount)
            pvalues.append(self.estimate_randfold_pvalue(row.precursor))

        putative_precursor["star_reads"] = star_counts
        putative_precursor["star"] = stars
        putative_precursor["randfold_pvalue"] = pvalues
        putative_precursor["cluster"] = x.cluster
        putative_precursor["mirdeep2_score"] = putative_precursor.apply(
            self.mirdeep2_like_score, axis=1
        )

        return putative_precursor

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

        Writes known and novel Excel tables, FASTA files and a database summary.
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

        filterd_results = pd.DataFrame(
            blast_results[
                (blast_results.qstart < 2)
                & (blast_results.length > 19)
                & (blast_results.pident > 94)
            ]
        ).drop_duplicates(subset=["Novel_ID"])

        known_miRNAs_dataframe = filterd_results[["Novel_ID", "Known_ID"]]

        known_miRNAs = known_miRNAs_dataframe.merge(
            novel_df, on=["Novel_ID"], indicator=True
        )

        known_miRNAs = known_miRNAs.drop("_merge", axis=1)
        known_miRNAs["Database"] = known_miRNAs["Known_ID"].map(mirbase_data["sources"])
        if len(known_miRNAs):
            summary = summarize_sources(known_miRNAs["Database"])
            summary.to_csv(f"{self.outdir}/known_database_summary.tsv", sep="\t", index=False)
            log.info(f"Known miRNAs by database: {dict(zip(summary['database'], summary['miRNAs']))}")

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
            return self.get_putative_precursor(
                filtered_cluster=df, precursor_step=self.precursorStep
            )
        anchor = df.iloc[0]
        mm = df.loc[
            (abs(df.start - anchor.start) <= 3) & (abs(df.end - anchor.end) <= 5)
        ].copy()
        if mm.reads.sum() <= self.readsCutoff:
            return None
        mm["mature"] = True
        other = df.loc[~df.index.isin(mm.index)].copy()
        other["mature"] = False
        return self.get_putative_precursor(
            filtered_cluster=mm, other=other, precursor_step=self.precursorStep
        )

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
        log.info(f"Scoring {len(tasks)} clusters of sample {sample[0]} with {self.workers} worker(s)")
        return [c for c in self._score_clusters(tasks, sample[0]) if c is not None and len(c) > 0]

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

        mirna_temp = mu.make_directory(os.path.join(self.outdir, "temp"))

        samples = list(self.samples.values())
        cache = {s[0]: os.path.join(mirna_temp, f"{s[0]}_scored.pkl") for s in samples}
        todo = [s for s in samples if not (self.resume and os.path.exists(cache[s[0]]))]
        with ThreadPoolExecutor(max_workers=max(1, self.jobs)) as pool:
            scored = dict(zip((s[0] for s in todo), pool.map(self._score_sample, todo), strict=True))

        for sample in samples:
            if sample[0] not in scored:
                log.info(f"Resuming: loading scored candidates of sample {sample[0]}")
                outputs.append(pd.read_pickle(cache[sample[0]]))
                continue

            N = scored[sample[0]]
            if len(N) == 0:
                log.warning(f"No precursor candidates produced for sample {sample[0]}")
                continue

            final_data = pd.concat(N, ignore_index=True)

            pre_prob = self.predict_precursor(final_data)

            final_data["score"] = pre_prob["miRNA"]
            final_data["ml_score"] = final_data["score"]
            final_data["mirdeep2_score"] = final_data.apply(
                self.mirdeep2_like_score, axis=1
            )

            outputs.append(final_data)

            final_data.to_csv(
                f"{mirna_temp}/{sample[0]}_raw_miRNA.txt", sep="\t", index=False
            )
            final_data.to_pickle(cache[sample[0]] + ".tmp")
            os.replace(cache[sample[0]] + ".tmp", cache[sample[0]])

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
