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
import shlex
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
MAX_DUPLEX_UNPAIRED = 5
MAX_ASYMMETRIC_BULGE = 4


def enclosing_pair(pt, i, j):
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


def stem_loop(pt, anchor, max_length):
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
        outer = enclosing_pair(pt, i, j)
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


def fold_precursor(seq):
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


def bulge_graph(seq):
    """
    Fold a sequence and return its forgi bulge graph.

    Args:
        seq (str): DNA or RNA sequence; T is converted to U.

    Returns:
        forgi.graph.bulge_graph.BulgeGraph: Graph of the folded structure.
    """

    seq = seq.replace("T", "U")
    pstructure, _pmfe = RNA.fold(seq)

    return fgb.BulgeGraph.from_dotbracket(pstructure, seq)


def find_close_match(preseq, seq):
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


def find_subseq(seq, string):
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
    return find_close_match(seq, string)


def check_mature_position(precursor, mature, bg=None, struct=None, start=None):
    """
    Check whether the mature sequence forms a valid miRNA/miRNA* duplex in the precursor fold.

    Enforces single-stem duplex criteria from published miRNA annotation guidelines:
    the mature cannot fold onto itself or span a loop/multi-branch junction, at most
    ``MAX_DUPLEX_UNPAIRED`` bases in the paired core (excluding the 2-nt 3' overhang) may
    be unpaired, opposite-arm pairing must be monotonic without nested stems, asymmetric
    bulges cannot exceed ``MAX_ASYMMETRIC_BULGE`` nt, and the terminal loop between arms
    must be at least 3 nt.

    Args:
        precursor (str): Precursor sequence.
        mature (str): Mature sequence.
        bg (forgi.graph.bulge_graph.BulgeGraph, optional): Bulge graph of the precursor fold.
        struct (str, optional): Dot-bracket structure of the precursor; folded when omitted.
        start (int, optional): 0-based start of the mature in ``precursor``.

    Returns:
        str: ``"ok"`` or the reason the check fails.
    """

    pos = start
    if pos is None:
        pos = precursor.find(mature)
        if pos < 0:
            pos = find_subseq(precursor, mature)
    if pos < 0:
        return "mature not found"
    if pos + len(mature) > len(precursor):
        return "3p mature outside"
    if struct is None:
        if bg is not None:
            struct = bg.to_dotbracket_string()
        else:
            struct, _ = fold_precursor(precursor.replace("T", "U"))
    pt = RNA.ptable(struct)
    m5 = pos + 1
    m3 = pos + len(mature)
    if any(m5 <= pt[k] <= m3 for k in range(m5, m3 + 1) if pt[k] != 0):
        return "mature in loop"
    down = [(k, pt[k]) for k in range(m5, m3 + 1) if pt[k] > m3]
    up = [(k, pt[k]) for k in range(m5, m3 + 1) if 0 < pt[k] < m5]
    if down and up:
        return "mature in loop"
    side = down or up
    if not side:
        return "mature unpaired"
    m3_core = max(m5, m3 - 2)
    core_len = m3_core - m5 + 1
    core_paired = sum(1 for k, _ in side if k <= m3_core)
    if core_len - core_paired > MAX_DUPLEX_UNPAIRED:
        return "mature unpaired"
    for (k1, p1), (k2, p2) in zip(side, side[1:]):
        if p1 <= p2:
            return "too many loops"
        if abs((k2 - k1) - (p1 - p2)) > MAX_ASYMMETRIC_BULGE:
            return "asymmetric bulge"
        if any(p2 < q < pt[q] < p1 for q in range(p2 + 1, p1) if pt[q] != 0):
            return "too many loops"
    if down:
        if side[-1][1] - m3 < 3:
            return "mature in loop"
        s3 = side[0][1] + (side[0][0] - m5) + 2
        if s3 > len(precursor) + 2:
            return "3p mature outside"
    else:
        if m5 - side[0][1] < 3:
            return "mature in loop"
    return "ok"


def estimate_star_sequence(precursor, mature, bg=None, struct=None, start=None):
    """
    Estimate the miRNA* (star) sequence on the opposite arm with a 2-nt 3' overhang.

    Args:
        precursor (str): Precursor sequence.
        mature (str): Mature sequence.
        bg (forgi.graph.bulge_graph.BulgeGraph, optional): Bulge graph of the precursor fold.
        struct (str, optional): Dot-bracket structure of the precursor; folded when omitted.
        start (int, optional): 0-based start of the mature in ``precursor``.

    Returns:
        str or None: Estimated star sequence, or None if it cannot be determined.
    """

    pos = start
    if pos is None:
        pos = precursor.find(mature)
        if pos < 0:
            pos = find_subseq(precursor, mature)
    if pos < 0 or pos + len(mature) > len(precursor):
        return None
    if struct is None:
        if bg is not None:
            struct = bg.to_dotbracket_string()
        else:
            struct, _ = fold_precursor(precursor.replace("T", "U"))
    pt = RNA.ptable(struct)
    m5 = pos + 1
    m3 = pos + len(mature)
    pairs = [(k, pt[k]) for k in range(m5, m3 + 1) if pt[k] != 0 and not (m5 <= pt[k] <= m3)]
    if not pairs:
        return None
    down = [(k, p) for k, p in pairs if p > m3]
    up = [(k, p) for k, p in pairs if p < m5]
    side = down if len(down) >= len(up) else up
    if not side:
        return None
    k5, p5 = side[0]
    m3_core = max(m5, m3 - 2)
    k3, p3 = min(side, key=lambda item: abs(item[0] - m3_core))
    s3 = p5 + (k5 - m5) + 2
    s5 = p3 - (m3_core - k3)
    s5 = max(1, s5)
    s3 = min(len(precursor), s3)
    if s3 - s5 + 1 < 15 or s3 - s5 + 1 > 30:
        return None
    if not (s3 <= pos or s5 - 1 >= pos + len(mature)):
        return None
    return precursor[s5 - 1 : s3]


def safe_star_sequence(precursor, mature, bg=None, struct=None, start=None):
    """
    Estimate a star sequence and return None when it cannot be determined.

    Args:
        precursor (str): Precursor sequence.
        mature (str): Mature sequence.
        bg (forgi.graph.bulge_graph.BulgeGraph, optional): Bulge graph of the precursor fold.
        struct (str, optional): Dot-bracket structure of the precursor.
        start (int, optional): 0-based start of the mature in ``precursor``.

    Returns:
        str or None: Estimated star sequence, or None if the fold does not allow one.
    """

    try:
        return estimate_star_sequence(precursor, mature, bg=bg, struct=struct, start=start)
    except (TypeError, ValueError, RuntimeError, StopIteration) as exc:
        log.debug(f"Could not estimate a star sequence: {exc}")
        return None


def _dinucleotide_shuffle(seq, rng):
    """
    Shuffle an RNA sequence while preserving exact dinucleotide counts and endpoints.

    Uses the Altschul-Erickson Eulerian trail construction with a uniform spanning
    arborescence sampled by loop-erased random walk (Bonnet et al., 2004).

    Args:
        seq (str): Input sequence.
        rng (numpy.random.Generator): Random number generator.

    Returns:
        str: Shuffled sequence with identical dinucleotide counts.
    """

    if len(seq) <= 2:
        return seq
    alpha = list(dict.fromkeys(seq))
    if len(alpha) == 1:
        return seq
    last = seq[-1]
    adj = {u: [] for u in alpha}
    for a, b in zip(seq, seq[1:]):
        adj[a].append(b)
    non_root = [u for u in alpha if u != last and adj[u]]
    in_tree = {last}
    last_edges = {}
    ok = True
    for u in non_root:
        if u in in_tree:
            continue
        path = {}
        cur = u
        for _ in range(512):
            nxt = adj[cur][int(rng.integers(len(adj[cur])))]
            path[cur] = nxt
            cur = nxt
            if cur in in_tree:
                break
        else:
            ok = False
            break
        cur = u
        while cur not in in_tree:
            in_tree.add(cur)
            last_edges[cur] = path[cur]
            cur = path[cur]
    if not ok:
        last_edges = {}
        for a, b in zip(seq, seq[1:]):
            if a != last:
                last_edges[a] = b

    out_lists = {}
    for u in alpha:
        items = list(adj[u])
        if u in last_edges:
            idx = items.index(last_edges[u])
            chosen = items.pop(idx)
            rng.shuffle(items)
            items.append(chosen)
        else:
            rng.shuffle(items)
        out_lists[u] = items

    pos_idx = dict.fromkeys(alpha, 0)
    res = [seq[0]]
    cur = seq[0]
    for _ in range(len(seq) - 1):
        nxt = out_lists[cur][pos_idx[cur]]
        pos_idx[cur] += 1
        res.append(nxt)
        cur = nxt
    return "".join(res)


def describe_precursor(precursor, mature, start=None, max_length=None):
    """
    Fold a precursor and describe where a mature sits in it.

    Args:
        precursor (str): Precursor sequence.
        mature (str): Mature sequence.
        start (int, optional): 0-based start of the mature in the precursor; located when not given.
        max_length (int, optional): Maximum span of the stem-loop used for the arm label. Defaults to the precursor length.

    Returns:
        dict: Dot-bracket structure, minimum free energy, mature-position check, arm label and star sequence.

    Raises:
        TypeError, ValueError, RuntimeError: If folding or the position check fails.
    """

    rna = precursor.replace("T", "U")
    struct, mfe = fold_precursor(rna)
    if start is None:
        pos = precursor.find(mature)
        if pos < 0:
            pos = find_subseq(precursor, mature)
    else:
        pos = start
    check = check_mature_position(precursor, mature, struct=struct, start=pos)
    span = max_length or len(precursor)
    centre = pos + len(mature) // 2 + 1
    pt = RNA.ptable(struct)
    paired = [k for k in range(pos + 1, pos + len(mature) + 1) if 1 <= k < len(pt) and pt[k] != 0]
    if paired:
        i, j = stem_loop(pt, min(paired, key=lambda k: abs(k - centre)), span)
    else:
        i, j = 1, span

    return {
        "struct": struct,
        "mfe": mfe,
        "mature_check": check,
        "position": "5-prime" if centre < (i + j) / 2 else "3-prime",
        "star": safe_star_sequence(precursor, mature, struct=struct, start=pos),
    }


def fold_hairpins(reader, chrom, start, end, strand, mature, max_extension=LOCAL_EXTENSION):
    """
    Fold genomic windows around a mature and yield the hairpins that hold it.

    The window is folded once per extension, the stem-loop around the mature is cut out and refolded
    at several flank lengths. Candidates without a base pair inside the mature are skipped.

    Args:
        reader (Read_process): Reference reader used to fetch the window sequence.
        chrom (str): Chromosome.
        start (int): Leftmost genomic coordinate of the mature.
        end (int): Rightmost genomic coordinate of the mature (exclusive).
        strand (str): "+" or "-".
        mature (str): Mature sequence in the orientation of the strand.
        max_extension (int): Nucleotides folded on each side of the mature. Defaults to 100.

    Yields:
        dict: Candidate hairpin with sequence, structure, free energy, mature-position check, arm label, star and coordinates.
    """

    seqlen = len(mature)
    seen = set()
    for extension in sorted({min(LOCAL_EXTENSION, max_extension), max_extension}):
        max_length = extension + seqlen
        wstart = max(0, start - extension)
        window = reader.get_sequence([chrom, wstart, end + extension, strand])
        if not window:
            continue
        wend = wstart + len(window)
        offset = start - wstart if strand == "+" else wend - end
        pos = window.find(mature, max(0, offset - 5))
        if pos < 0:
            pos = find_subseq(window, mature)
        if pos < 0:
            continue

        try:
            md = RNA.md()
            md.max_bp_span = max_length
            wide, _ = RNA.fold_compound(window.replace("T", "U"), md).mfe()
        except (TypeError, ValueError, RuntimeError) as exc:
            log.warning(f"Could not fold {chrom}:{start}: {exc}")
            continue
        pt = RNA.ptable(wide)
        centre = pos + seqlen // 2 + 1
        paired = [k for k in range(pos + 1, pos + seqlen + 1) if pt[k] != 0]
        if not paired:
            continue
        i, j = stem_loop(pt, min(paired, key=lambda k: abs(k - centre)), max_length)

        for flank in FLANKS:
            lo = max(0, min(i - 1, pos) - flank)
            hi = min(len(window), max(j, pos + seqlen) + flank)
            precursor = window[lo:hi]
            if "N" in precursor or precursor in seen:
                continue
            seen.add(precursor)
            try:
                entry = describe_precursor(precursor, mature)
            except (TypeError, ValueError, RuntimeError, StopIteration) as exc:
                log.warning(f"Skipped precursor at {chrom}:{start}: {exc}")
                continue
            # Keep the arm label of the wide-window fold, as the discovery path always has.
            entry["position"] = "5-prime" if centre < (i + j) / 2 else "3-prime"
            if strand == "+":
                pstart, pend = wstart + lo, wstart + hi
            else:
                pstart, pend = wend - hi, wend - lo
            yield {
                "chrom": chrom,
                "mature": mature,
                "mature_length": seqlen,
                "precursor": precursor,
                "precursor_length": len(precursor),
                "start": pstart,
                "end": pend,
                "mature_start": start,
                "strand": strand,
                "hairpin_start": pos - lo,
                **entry,
            }


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
        read_cutoff=20,
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
        require_star=True,
        max_inconsistent=0.25,
        min_5p_consistency=0.5,
        min_novel_score=1.0,
        random_seed=42,
        workers=1,
        jobs=1,
        resume=False,
        slurm=False,
        mem=None,
        rfam=None,
        aligned_reads=None,
    ):
        """
        Initialize the novel miRNA predictor.

        Args:
            samples (dict, optional): Sample name mapped to sample information (name at index 0, reads at 2, clusters at 3).
            genome (str, optional): Reference genome FASTA path.
            read_cutoff (int): Minimum reads of a cluster or mature candidate. Defaults to 20.
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
            require_star (bool): Require star reads on novel miRNA candidates. Defaults to True.
            max_inconsistent (float): Maximum fraction of the reads on a novel hairpin outside the mature and star positions. Defaults to 0.25.
            min_5p_consistency (float): Minimum share of mature-locus reads starting at the modal 5-prime position for novel miRNAs. Defaults to 0.5.
            min_novel_score (float): Minimum heuristic score of a novel miRNA. Defaults to 1.0.
            random_seed (int): Seed for shuffling. Defaults to 42.
            workers (int): Processes for scoring clusters. Defaults to 1.
            jobs (int): Threads for scoring samples. Defaults to 1.
            resume (bool): Reuse the scored candidates saved for a sample by an earlier run. Defaults to False.
            slurm (bool): Score each sample as an element of a Slurm job array, each with ``workers`` cores. Defaults to False.
            mem (int, optional): Memory in GB of each Slurm job. Defaults to the cluster setting.
            rfam (dict, optional): Sample name to a dictionary of read name to the set of Rfam
                sequences the read matched. Defaults to None.
            aligned_reads (dict, optional): Sample name to ``[group, replicate, aligned reads DataFrame]``
                used for hairpin-wide star and Dicer-precision counting. Defaults to None.
        """

        self.genome = genome
        self.readsCutoff = read_cutoff
        self.clusterCutoff = cluster_cutoff
        self.species = species
        self.species_type = species_type
        self.outdir = outdir
        self.samples = samples
        self.aligned_reads = aligned_reads
        self.max_extension = max_extension
        self.filterCriteria = filter_criteria
        self.score = score
        self.randfold_iterations = randfold_iterations
        self.randfold_pvalue_cutoff = randfold_pvalue_cutoff
        self.min_mature_reads = min_mature_reads
        self.min_star_reads = min_star_reads
        self.require_star = require_star
        self.max_inconsistent = max_inconsistent
        self.min_5p_consistency = min_5p_consistency
        self.min_novel_score = min_novel_score
        self.random_seed = random_seed
        self.workers = workers
        self.jobs = jobs
        self.resume = resume
        self.slurm = slurm
        self.mem = mem
        self.rfam = rfam or {}
        self.rng = np.random.default_rng(random_seed)
        self._reader = Read_process(reference=self.genome) if self.genome else None
        if self._reader is not None:
            self._reader.ensure_index()

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

        return fold_precursor(seq)

    def get_bulge(self, seq):
        """
        Fold a sequence and return its forgi bulge graph.

        Args:
            seq (str): DNA or RNA sequence; T is converted to U.

        Returns:
            forgi.graph.bulge_graph.BulgeGraph: Graph of the folded structure.
        """

        return bulge_graph(seq)

    def close_match(self, preseq, seq):
        """
        Find where seq matches preseq with one mismatch, or -1.

        Args:
            preseq (str): Sequence to search in.
            seq (str): Sequence to find.

        Returns:
            int: Start index of the first one-mismatch match, or -1 if none.
        """

        return find_close_match(preseq, seq)

    def find_subseq(self, seq, string):
        """
        Locate a mature sequence in a precursor, falling back to a one-mismatch match.

        Args:
            seq (str): Precursor sequence.
            string (str): Mature sequence.

        Returns:
            int: Start index of the mature in the precursor, or -1 if not found.
        """

        return find_subseq(seq, string)

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

        return check_mature_position(precursor, mature, bg)

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

        return estimate_star_sequence(precursor, mature, bg)

    def estimate_randfold_pvalue(self, precursor):
        """
        Estimate randfold significance via dinucleotide-shuffled minimum free energies.

        Args:
            precursor (str): Precursor sequence.

        Returns:
            float: Smoothed p-value, the share of dinucleotide shuffles with an MFE as low as
            the observed one; 1.0 for an empty precursor.
        """
        if not precursor:
            return 1.0
        seq = precursor.replace("T", "U")
        _, observed_mfe = RNA.fold(seq)
        better_or_equal = 0
        for _ in range(self.randfold_iterations):
            shuffled = _dinucleotide_shuffle(seq, self.rng)
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

    def get_putative_precursor(self, filtered_cluster=None, other=None, rfam=None, hairpin_reads=None):
        """
        Extract and score the candidate hairpins of one read cluster.

        The mature and its flanks are folded once. The stem-loop that holds the
        mature is cut out, refolded and checked, and the model features are only
        built for hairpins with a good mature position. Randfold p-values and the
        heuristic score are added later, once the model has scored the candidates.

        Args:
            filtered_cluster (pandas.DataFrame): Reads of the mature cluster.
            other (pandas.DataFrame, optional): Remaining cluster reads used when
                ``hairpin_reads`` is not provided.
            rfam (dict, optional): Read name to the set of Rfam sequences the read matched.
                Defaults to None.
            hairpin_reads (pandas.DataFrame, optional): All aligned reads near the locus on
                the same chromosome and strand, used for hairpin-wide star and Dicer-precision
                counting across separated 5p and 3p arms. Defaults to None.

        Returns:
            pandas.DataFrame or None: Candidate hairpins, or None if no valid candidate was found.

        Raises:
            ValueError: If filtered_cluster is not provided.
        """

        if filtered_cluster is None:
            raise ValueError("filtered_cluster is required")

        x = filtered_cluster.iloc[0]
        chrom, start, end, strand = str(x["name"]), x.start, x.end, x.strand
        reader = self._reader or Read_process(reference=self.genome)
        refseq = reader.get_sequence((chrom, max(0, start - 10), end + 10, strand))
        mature = reader.get_consensus_read(
            ref_precursor=refseq, aligned_cluster=filtered_cluster
        )
        if mature is None:
            mature = x.seq
        seqlen = len(mature)
        dedup_cols = [c for c in ("read_name", "start", "end") if c in filtered_cluster.columns]
        clean_cluster = filtered_cluster.drop_duplicates(dedup_cols) if dedup_cols else filtered_cluster
        maturecounts = clean_cluster.reads.sum()
        five_prime_col = "start" if strand == "+" else "end"
        cluster_consistency = (
            float(clean_cluster.groupby(five_prime_col)["reads"].sum().max() / maturecounts)
            if maturecounts > 0
            else 0.0
        )
        rfam = rfam or {}
        if rfam:
            hits = filtered_cluster[filtered_cluster["read_name"].isin(rfam)].drop_duplicates("read_name")
            rfam_reads = float(hits["reads"].sum())
            accessions = sorted({name for read in hits["read_name"] for name in rfam.get(read, ())})
        else:
            rfam_reads = 0.0
            accessions = []
        rfam_label = ",".join(accessions) if accessions else "-"
        putative_precursor = []

        for entry in fold_hairpins(
            reader, chrom, start, end, strand, mature, self.max_extension
        ):
            if entry["mature_check"] != "ok":
                continue

            m = Precursor_Features(entry["precursor"]).generate_features()
            if m is None:
                continue

            pstart, pend = entry["start"], entry["end"]
            hstart = entry.get("hairpin_start")
            if hstart is None:
                found_pos = entry["precursor"].find(mature)
                if found_pos < 0:
                    found_pos = self.find_subseq(entry["precursor"], mature)
                hstart = found_pos if found_pos >= 0 else None
            if hstart is not None:
                m_gstart = pstart + hstart if strand == "+" else pend - hstart - seqlen
            else:
                m_gstart = start
            m_gend = m_gstart + seqlen

            star = entry["star"]
            s_pos = -1
            if isinstance(star, str) and star:
                s_pos = entry["precursor"].find(star)
                if s_pos < 0:
                    s_pos = self.find_subseq(entry["precursor"], star)

            starcount = 0
            inconsistent = 0.0
            consistency = cluster_consistency
            if hairpin_reads is not None and len(hairpin_reads):
                hp = hairpin_reads[(hairpin_reads.start >= pstart) & (hairpin_reads.end <= pend)]
                if len(hp):
                    hp_cols = [c for c in ("read_name", "start", "end") if c in hp.columns]
                    if hp_cols:
                        hp = hp.drop_duplicates(hp_cols)
                    if strand == "+":
                        mature_mask = (abs(hp.start - m_gstart) <= 2) & (abs(hp.end - m_gend) <= 4)
                    else:
                        mature_mask = (abs(hp.end - m_gend) <= 2) & (abs(hp.start - m_gstart) <= 4)
                    m_sub = hp.loc[mature_mask]
                    m_locus_reads = m_sub.reads.sum()
                    if m_locus_reads > 0:
                        consistency = float(
                            m_sub.groupby(five_prime_col)["reads"].sum().max() / m_locus_reads
                        )
                    if s_pos >= 0 and star is not None:
                        ss = pstart + s_pos if strand == "+" else pend - s_pos - len(star)
                        se = ss + len(star)
                        if strand == "+":
                            raw_star = (abs(hp.start - ss) <= 3) & (abs(hp.end - se) <= 4)
                        else:
                            raw_star = (abs(hp.end - se) <= 3) & (abs(hp.start - ss) <= 4)
                        star_mask = raw_star & (~mature_mask)
                        starcount = hp.loc[star_mask, "reads"].sum()
                    else:
                        star_mask = pd.Series(False, index=hp.index)
                    total_hp = hp.reads.sum()
                    if total_hp > 0:
                        inconsistent = float(
                            hp.loc[~(mature_mask | star_mask), "reads"].sum() / total_hp
                        )
            elif other is not None:
                if s_pos >= 0 and star is not None:
                    ss = pstart + s_pos if strand == "+" else pend - s_pos - len(star)
                    se = ss + len(star)
                    if strand == "+":
                        sreads = other[(abs(other.start - ss) <= 3) & (abs(other.end - se) <= 4)]
                    else:
                        sreads = other[(abs(other.end - se) <= 3) & (abs(other.start - ss) <= 4)]
                    starcount = sreads.reads.sum()
                on_hairpin = other[(other.start >= pstart) & (other.end <= pend)].reads.sum()
                denom = maturecounts + on_hairpin
                inconsistent = max(0.0, float(on_hairpin - starcount)) / denom if denom > 0 else 0.0

            putative_precursor.append(
                {
                    "chrom": chrom,
                    "mature": mature,
                    "mature_length": seqlen,
                    "precursor": entry["precursor"],
                    "precursor_length": entry["precursor_length"],
                    "start": pstart,
                    "end": pend,
                    "mature_start": m_gstart,
                    "strand": strand,
                    "struct": entry["struct"],
                    "mfe": entry["mfe"],
                    "mature_check": "ok",
                    "features": m,
                    "mature_reads": maturecounts,
                    "position": entry["position"],
                    "star_reads": starcount,
                    "inconsistent_fraction": inconsistent,
                    "five_prime_consistency": consistency,
                    "star": star,
                    "randfold_pvalue": 1.0,
                    "cluster": x.cluster,
                    "rfam_reads": rfam_reads,
                    "Rfam": rfam_label,
                }
            )

        if len(putative_precursor) == 0:
            return
        return pd.DataFrame(putative_precursor)

    def filter_miRNA(self, df):
        """
        Keep the best-supported precursor per mature sequence and assign novel IDs.

        Across samples and candidate flanks, valid precursors with star-strand support
        take precedence over star-less observations, followed by model score, heuristic
        score and mature read depth.

        Args:
            df (pandas.DataFrame): Scored candidates with "mature" and "score" columns.

        Returns:
            pandas.DataFrame: One row per mature sequence, with a Novel_ID column.
        """
        valid = (df["score"] >= self.score)
        if "mature_check" in df.columns:
            valid = valid & (df["mature_check"] == "ok")
        has_star = (
            df["star_reads"] >= self.min_star_reads
            if "star_reads" in df.columns
            else pd.Series(False, index=df.index)
        )
        sort_cols = ["_valid", "_has_star", "score"]
        for extra in ("heuristic_score", "mature_reads"):
            if extra in df.columns:
                sort_cols.append(extra)
        ranked = df.assign(_valid=valid, _has_star=has_star).sort_values(
            sort_cols, ascending=False, kind="mergesort"
        )
        idx = (
            ranked.drop_duplicates(subset=["mature"], keep="first")
            .sort_values("mature", kind="mergesort")
            .index.tolist()
        )
        final_miRNAs = df.loc[idx].reset_index(drop=True)

        with contextlib.suppress(Exception):
            final_miRNAs["Novel_ID"] = [
                f"{self.species}_mipyrna_{final_miRNAs.iloc[i]['chrom']}_miR_{i + 1}"
                for i in range(len(final_miRNAs))
            ]

        return final_miRNAs

    def apply_strict_filters(self, df):
        """
        Apply the structure, score, read count and randfold filters prior to known/novel split.

        Star-strand and Dicer-precision rules are enforced downstream in :meth:`_dicer_filter`
        on novel candidates only, so BLAST-confirmed known miRNAs without star reads are
        retained.

        Args:
            df (pandas.DataFrame): Candidate miRNAs.

        Returns:
            pandas.DataFrame: Candidates passing structure, score, mature read and randfold checks.
        """
        mask = (
            (df["mature_check"] == "ok")
            & (df["score"] >= self.score)
            & (df["mature_reads"] >= self.min_mature_reads)
            & (df["randfold_pvalue"] <= self.randfold_pvalue_cutoff)
        )
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
            chrom = getattr(r, "chrom", "")
            rstart = getattr(r, "start", None)
            rend = getattr(r, "end", None)
            if (
                chrom != ""
                and rstart is not None
                and rend is not None
                and str(rstart) not in ("", "nan")
                and str(rend) not in ("", "nan")
            ):
                loci[(chrom, r.strand)].append((int(rstart), int(rend)))
        keep = [not self._overlaps(r, loci) for _i, r in novel_df.iterrows()]
        log.info(
            f"Removed {len(keep) - sum(keep)} novel candidates overlapping known miRNA hairpins"
        )
        return novel_df[keep]

    def _dicer_filter(self, df):
        """Remove novel candidates whose reads do not meet published Dicer biogenesis criteria.

        Enforces hairpin-wide Dicer precision (``inconsistent_fraction <= max_inconsistent``,
        Axtell & Meyers 2018), 5-prime cleavage homogeneity (``five_prime_consistency >= min_5p_consistency``),
        heuristic score, kingdom mature length bounds, star-strand read support when
        ``require_star`` is enabled, and mandatory star reads for 23-24 nt plant candidates
        to exclude heterochromatic siRNAs.

        Args:
            df (pandas.DataFrame): Novel candidates.

        Returns:
            pandas.DataFrame: Candidates passing Dicer biogenesis checks.
        """
        if df.empty:
            return df
        keep = (df["inconsistent_fraction"] <= self.max_inconsistent) & (
            df["heuristic_score"] >= self.min_novel_score
        )
        if "five_prime_consistency" in df.columns:
            keep = keep & (df["five_prime_consistency"] >= self.min_5p_consistency)
        if self.require_star:
            keep = keep & (df["star_reads"] >= self.min_star_reads)
        if "mature_length" in df.columns:
            if self.species_type == "plants":
                keep = keep & df["mature_length"].between(20, 24)
                hc_len = df["mature_length"].isin((23, 24))
                keep = keep & (~hc_len | (df["star_reads"] >= max(1, self.min_star_reads)))
            elif self.species_type == "animals":
                keep = keep & df["mature_length"].between(18, 25)
        log.info(
            f"Removed {int((~keep).sum())} novel candidates failing Dicer precision, 5' consistency, star or score checks"
        )
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
        records = []
        seen = set()
        for _i, r in df.iterrows():
            if seq == "mature":
                sequence = r["mature"]
            elif seq == "hairpin":
                sequence = r["precursor"]
            else:
                continue
            if type == "known":
                if r["Known_ID"] in seen:
                    continue
                seen.add(r["Known_ID"])
                name = r["Known_ID"]
            else:
                name = r["Novel_ID"]
            records.append((name, sequence))
        return mu.write_fasta(records, output)

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
            blastdb_command = f"makeblastdb -in {shlex.quote(mirbase_data['mature'])} -dbtype 'nucl'"

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
            blastn_command = (
                f"blastn -db {shlex.quote(mirbase_data['mature'])} "
                f"-query {shlex.quote(novel_mature)} -strand plus "
                f"-evalue 100 -out {shlex.quote(blast_out)} -outfmt 6 -word_size 4"
            )
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
        # so allow a one-base offset while keeping at least 18 forward-strand aligned bases.
        # When a candidate matches several entries, keep the best alignment.
        hits = pd.DataFrame(
            blast_results[
                (blast_results.qstart < 3)
                & (blast_results.sstart < blast_results.send)
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
            task (tuple): ``(cluster ID, cluster size, cluster reads DataFrame, sample name)``
                optionally followed by ``local_reads`` DataFrame for hairpin-wide counting.

        Returns:
            pandas.DataFrame or None: Candidate precursors, or None if the cluster has too few reads.
        """
        if len(task) == 5:
            cluster_id, clust_size, df, name, local_reads = task
        else:
            cluster_id, clust_size, df, name = task
            local_reads = None
        self.rng = np.random.default_rng([self.random_seed, int(cluster_id)])
        rfam = self.rfam.get(name, {})
        df = df.sort_values("reads", ascending=False, kind="mergesort")
        if clust_size < self.clusterCutoff:
            df["mature"] = True
            if local_reads is not None:
                return self.get_putative_precursor(
                    filtered_cluster=df, rfam=rfam, hairpin_reads=local_reads
                )
            return self.get_putative_precursor(filtered_cluster=df, rfam=rfam)
        anchor = df.iloc[0]
        mm = df.loc[
            (abs(df.start - anchor.start) <= 3) & (abs(df.end - anchor.end) <= 5)
        ].copy()
        if mm.reads.sum() < self.readsCutoff:
            return None
        mm["mature"] = True
        other = df.loc[~df.index.isin(mm.index)].copy()
        other["mature"] = False
        if local_reads is not None:
            return self.get_putative_precursor(
                filtered_cluster=mm, other=other, rfam=rfam, hairpin_reads=local_reads
            )
        return self.get_putative_precursor(filtered_cluster=mm, other=other, rfam=rfam)

    def _score_sample(self, sample, all_reads=None):
        """
        Score the candidate precursors of all read clusters of one sample.

        Args:
            sample (list): Sample information (name at index 0, reads at 2, clusters at 3).
            all_reads (pandas.DataFrame, optional): Unfiltered genome-aligned reads of the
                sample for hairpin-wide star and Dicer-precision counting.

        Returns:
            list: Non-empty candidate precursor DataFrames, one per scored cluster.
        """
        clusts = sample[3][sample[3].reads >= self.readsCutoff]
        keep_clusters = set(clusts["cluster"])
        kept_reads = (
            sample[2][sample[2]["cluster"].isin(keep_clusters)]
            if keep_clusters
            else sample[2].iloc[:0]
        )
        groups = dict(tuple(kept_reads.groupby("cluster"))) if len(kept_reads) else {}
        if all_reads is None and self.aligned_reads and sample[0] in self.aligned_reads:
            all_reads = self.aligned_reads[sample[0]][2]
        if all_reads is None:
            all_reads = sample[2]
        by_locus = {}
        if all_reads is not None and len(all_reads):
            for (chrom, strand), sub in all_reads.groupby(["name", "strand"], sort=False):
                ordered = sub.sort_values("start", kind="mergesort")
                by_locus[(str(chrom), str(strand))] = (ordered, ordered["start"].to_numpy())

        pad = self.max_extension + 10
        tasks = []
        for _i, c in clusts.iterrows():
            if c.cluster not in groups:
                continue
            entry = by_locus.get((str(c["name"]), str(c.strand)))
            if entry is not None:
                ordered, starts = entry
                win_start = max(0, int(c.start) - pad)
                win_end = int(c.end) + pad
                lo = int(np.searchsorted(starts, max(0, win_start - 40), side="left"))
                hi = int(np.searchsorted(starts, win_end, side="right"))
                sub = ordered.iloc[lo:hi]
                local_reads = sub[(sub.end > win_start) & (sub.start < win_end)]
            else:
                local_reads = groups[c.cluster]
            tasks.append((c.cluster, c.clust_size, groups[c.cluster], sample[0], local_reads))

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
                light.aligned_reads = None
                _init_worker(light)
                pvalues = [_run_randfold(t) for t in tasks]
            else:
                light = copy.copy(self)
                light.samples = None
                light.aligned_reads = None
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
            tasks (list): Cluster tasks of (cluster ID, cluster size, reads DataFrame, sample name[, local_reads]).
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
        light.aligned_reads = None
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

    def _score_one(self, sample, all_reads=None):
        """Score one sample from start to finish; run by one element of a Slurm job array.

        Args:
            sample (list): Sample information (name at index 0, reads at 2, clusters at 3).
            all_reads (pandas.DataFrame, optional): Unfiltered genome-aligned reads of the sample.
        """
        candidates = self._score_sample(sample, all_reads=all_reads)
        if candidates:
            self._finish_sample(sample[0], candidates)

    def _score_on_slurm(self, samples):
        """Score samples as one Slurm job array, one element per sample.

        Args:
            samples (list): Samples to score.
        """
        light = copy.copy(self)
        light.samples = None
        light.aligned_reads = None
        light.jobs = 1
        light.slurm = False
        aligned = self.aligned_reads or {}
        mu.run_array(
            [
                (
                    light,
                    "_score_one",
                    (sample, aligned[sample[0]][2] if sample[0] in aligned else None),
                )
                for sample in samples
            ],
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
