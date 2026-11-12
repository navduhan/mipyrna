"""Read Clustering Module

This module groups aligned small RNA reads into candidate loci using a
bx-python ClusterTree, then filters the clusters by read length and size.

Features:
    - Read length filtering before clustering
    - Per-chromosome, per-strand interval clustering with a configurable distance
    - Per-cluster summary of read counts, maximum length and size
    - Samples clustered in parallel

Configuration:
    Distance, minimum cluster size and read length range are set in the
    constructor. Input is a dictionary of samples holding aligned reads.

Dependencies:
    - bx-python: ClusterTree interval clustering
    - pandas: Table handling
    - tqdm: Progress bar

Classes / Functions / Exceptions:
    - Cluster_reads: Cluster aligned reads into candidate miRNA loci.

:Author: Naveen Duhan
"""

import collections
import multiprocessing

import pandas as pd
from bx.intervals.cluster import ClusterTree
from tqdm import tqdm


class Cluster_reads:
    """Cluster aligned reads into candidate miRNA loci."""

    def __init__(self, cluster_distance=10, min_size=1, min_length=20, max_length=24, key='align_id'):
        """Initialize the clustering parameters.

        Args:
            cluster_distance (int): Maximum distance between reads in a cluster. Defaults to 10.
            min_size (int): Minimum number of collapsed read alignments in a cluster. Defaults to 1.
            min_length (int): Minimum read length kept. Defaults to 20.
            max_length (int): Maximum read length kept. Defaults to 24.
            key (str): Column used as the read index when present. Defaults to ``'align_id'``.
        """


        self.cluster_distance = cluster_distance
        self.min_size = min_size
        self.key = key
        self.min = min_length
        self.max = max_length

    def _cluster_sample(self, key, sample):
        """Cluster the reads of one sample.

        Args:
            key (str): Sample name.
            sample (list): ``[group, replicate, aligned reads DataFrame]``.

        Returns:
            tuple: ``(key, [group, replicate, reads DataFrame, clusters DataFrame])``.
        """
        filtered_reads = sample[2][(sample[2].length<=self.max) & (sample[2].length>=self.min)]
        # Canonical ids also keep ClusterTree's returned member order independent
        # of the SAM line numbers assigned by the alignment parser.
        input_file = filtered_reads.sort_values(
            ["name", "start", "end", "strand", "read_name"], kind="mergesort"
        ).copy()
        input_file = input_file.reset_index(drop=True)
        if self.key in input_file.columns:
            input_file[self.key] = range(1, len(input_file) + 1)
            input_file = input_file.set_index(self.key)

        cluster_trees = collections.defaultdict(
            lambda: ClusterTree(self.cluster_distance, self.min_size)
        )
        for chrom, strand, start, end, name in zip(
            input_file["name"],
            input_file["strand"],
            input_file.start,
            input_file.end,
            input_file.index,
            strict=True,
        ):
            cluster_trees[(chrom, strand)].insert(start, end, name)
        member_ids = []
        cl_starts = []
        cl_ends = []
        cl_nums = []
        icount = 1
        for chrom, strand in sorted(cluster_trees):
            cltree = cluster_trees[(chrom, strand)]

            for start, end, ids in cltree.getregions():
                count = len(ids)
                member_ids.extend(ids)
                cl_starts.extend([start] * count)
                cl_ends.extend([end] * count)
                cl_nums.extend([icount] * count)
                icount += 1
        if member_ids:
            df = input_file.loc[member_ids].copy()
            df["cl_start"] = cl_starts
            df["cl_end"] = cl_ends
            df["cluster"] = cl_nums
        else:
            df = input_file.iloc[:0].copy()
            for column in ("cl_start", "cl_end", "cluster"):
                df[column] = pd.Series(index=df.index, dtype=int)

        clusts = df.groupby(['name','cluster','cl_start','cl_end','strand'])\
                    .agg({'reads':'sum','length':'max'})\
                    .reset_index()\
                    .rename(columns={'cl_start':'start','cl_end':'end'})

        clusts['clust_size'] = clusts.end-clusts.start
        return key, [sample[0], sample[1], df, clusts]

    def get_clusters(self, samples, jobs=1):
        """Cluster the reads of each sample.

        Args:
            samples (dict): Sample name to ``[group, replicate, aligned reads DataFrame]``.
            jobs (int): Number of samples clustered at once. Defaults to 1.

        Returns:
            dict: Sample name to ``[group, replicate, reads DataFrame, clusters DataFrame]``.
            The clusters DataFrame has ``name``, ``cluster``, ``start``, ``end``,
            ``strand``, ``reads``, ``length`` and ``clust_size`` columns.
        """
        tasks = list(samples.items())
        if jobs <= 1 or len(tasks) < 2:
            results = [self._cluster_sample(k, v) for k, v in tqdm(tasks, desc="Processing samples")]
        else:
            with multiprocessing.Pool(min(jobs, len(tasks))) as pool:
                results = list(tqdm(pool.starmap(self._cluster_sample, tasks), total=len(tasks), desc="Processing samples"))
        return dict(results)
