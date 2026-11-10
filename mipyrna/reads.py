"""Read Processing Module

This module collapses small RNA reads, parses alignments into tables and
extracts sequences from a reference.

Features:
    - Parallel collapsing of identical reads into FASTA with read counts
    - Parsing of SAM or Bowtie alignments into DataFrames
    - Consensus read selection for a cluster
    - Reference sequence extraction with strand handling

Configuration:
    An optional reference FASTA is given to the constructor. Collapsed read
    names have the form ``<sample>_<n>_x<count>``.

Dependencies:
    - HTSeq: SAM parsing
    - Biopython: FASTQ parsing and reverse complement
    - pyfaidx: Indexed reference access
    - pandas, tqdm

Classes / Functions / Exceptions:
    - Read_process: Collapse, parse and query small RNA reads.

:Author: Naveen Duhan
"""

import gzip
import multiprocessing
import os
import re

import HTSeq
import pandas as pd
import pyfaidx as pf
from Bio import Seq, SeqIO
from tqdm import tqdm

from mipyrna.logger import MiPyRNALogger

log = MiPyRNALogger(mode='a', log='reads')


class Read_process:
    """Collapse, parse and query small RNA reads."""

    def __init__(self,  reference=None):
        """Initialize the read processor.

        Args:
            reference (str, optional): Reference FASTA used by :meth:`get_sequence`.
        """

        self.ref = reference


    def _collapse_process_sample(self, key, sample, outdir):
        """Collapse identical reads of one sample into a FASTA with read counts.

        Args:
            key (str): Sample name.
            sample (list): ``[group, replicate, FASTQ path]``; gzip files are supported.
            outdir (str): Output directory.

        Returns:
            tuple: ``(key, [group, replicate, collapsed FASTA path])``.
        """
        collapsed_reads = {}
        
        if sample[2].endswith(".gz"):
            file = gzip.open(sample[2], "rt")
        else:
            file = open(sample[2])
            
        fastqs = SeqIO.parse(file, 'fastq')

        for fq in fastqs:
            seq = str(fq.seq)
            
            if seq in collapsed_reads:
                collapsed_reads[seq]['reads'] += 1
            else:
                collapsed_reads[seq] = {'name': fq.name, 'reads': 1}

        count = 1  
        output = os.path.join(outdir, f"{key}_collapsed_reads.txt")

        with open(output, 'w') as fp:
            for d in collapsed_reads:
                fp.write(f">{key}_{count}_x{collapsed_reads[d]['reads']}\n{d}\n")
                count += 1

        return key, [sample[0], sample[1], output]
    
    def collapse_reads(self, samples, outdir, jobs=4):
        """Collapse the reads of all samples in parallel.

        Args:
            samples (dict): Sample name to ``[group, replicate, FASTQ path]``.
            outdir (str): Output directory.
            jobs (int): Number of worker processes. Defaults to 4.

        Returns:
            dict: Sample name to ``[group, replicate, collapsed FASTA path]``.
        """
        outreads = {}

        with multiprocessing.Pool(processes=max(1, jobs)) as pool:
            
            results = list(tqdm(pool.starmap(self._collapse_process_sample, [(key, sample, outdir) for key, sample in samples.items()]), total=len(samples), desc='Processing samples'))

        for key, result in results:
            outreads[key] = result

        return outreads

    def process_sample(self, key, sample, alignType):
        """Parse the alignments of one sample into a DataFrame.

        Args:
            key (str): Sample name.
            sample (list): ``[group, replicate, alignment path]``.
            alignType (str): ``'SAM'`` or ``'BWT'`` (Bowtie output).

        Returns:
            tuple: ``(key, [group, replicate, DataFrame])``. The DataFrame has
            ``read_name``, ``seq``, ``name``, ``strand``, ``start``, ``end``,
            ``length``, ``read_id``, ``reads`` and ``align_id`` columns.
        """
        areads = []
        count = 1

        if alignType == 'SAM':
            for r in HTSeq.SAM_Reader(str(sample[2])):
                try:
                    if r.aligned:
                        seq = r.read.seq.decode()
                        name = r.read.name
                        read_id = name.rsplit("_x", 1)[0].rsplit("_", 1)[-1]
                        reads = int(name.rsplit("_x", 1)[1]) if "_x" in name else 1
                        areads.append([name, seq, r.iv.chrom, r.iv.strand, r.iv.start, r.iv.end, len(str(seq)), read_id, reads, count])
                        count += 1
                        
                except Exception as exc:
                    log.debug(f"Skipping unparsable alignment {r.read.name}: {exc}")

        elif alignType == 'BWT':
            df = pd.read_csv(sample[2], sep="\t", header=None, dtype={2: str})
            df.columns = ['read', 'strand', 'chrom', 'start', 'mature', 'distance', 'num_align', 'info']

            for _i, row in df.iterrows():
                mature = str(row['mature'])
                if str(row['strand']) == "-":
                    mature = mature[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
                gseq = list(mature)
                edit = ["m"] * len(mature)
                mm = 0
                data = row['info']
                if isinstance(data, str):
                    changes = data.split(",")
                    for change in changes:
                        if re.match(r'(\d+):(\w+)>\w+', change):
                            mm += 1
                            pos, ont = change.split(':')
                            gseq[int(pos)] = ont.split(">")[1]  # Convert to lowercase
                            edit[int(pos)] = "M"
                read_name = str(row.read).split()[0]
                if "_x" in read_name:
                    _, raw_reads = read_name.rsplit("_x", 1)
                    try:
                        read_count = int(raw_reads)
                    except ValueError:
                        read_count = 1
                else:
                    read_count = 1
                start = int(row.start) + 1
                end = start + len(mature)
                areads.append([
                    read_name,
                    mature,
                    row.chrom,
                    row.strand,
                    start,
                    end,
                    len(mature),
                    count,
                    read_count,
                    count
                ])
                count += 1

        return key, [sample[0], sample[1], pd.DataFrame(areads, columns=['read_name', 'seq', 'name', 'strand', 'start', 'end', 'length', 'read_id', 'reads', 'align_id'])]

    def aligned_reads(self, samples, alignType='SAM', jobs=4):
        """Parse the alignments of all samples in parallel.

        Args:
            samples (dict): Sample name to ``[group, replicate, alignment path]``.
            alignType (str): ``'SAM'`` or ``'BWT'``. Defaults to ``'SAM'``.
            jobs (int): Number of worker processes. Defaults to 4.

        Returns:
            dict: Sample name to ``[group, replicate, aligned reads DataFrame]``.
        """
        outaligned = {}

        with multiprocessing.Pool(processes=max(1, jobs)) as pool:
            results = list(tqdm(pool.starmap(self.process_sample, [(key, sample, alignType) for key, sample in samples.items()])))
        for key, result in results:
            outaligned[key] = result

        return outaligned
    
    def get_consensus_read(self, ref_precursor, aligned_cluster):
        """Get the most abundant read of a cluster, preferring one found in the precursor.

        Args:
            ref_precursor (str): Precursor sequence.
            aligned_cluster (pandas.DataFrame): Aligned reads of the cluster.

        Returns:
            str: Sequence of the consensus read.
        """

        mature_read = None
        seq_start = 'start' if aligned_cluster.iloc[0].strand == '+' else 'end'
        # group = aligned_cluster.groupby([seq_start, 'length']).agg({'reads':np.sum}).sort_values(by='reads', ascending=False).reset_index()
        group = aligned_cluster.groupby([seq_start, 'length'])['reads'].agg('sum').sort_values(ascending=False).reset_index()
        consensus = aligned_cluster[(aligned_cluster[seq_start]==group.iloc[0][seq_start]) &(aligned_cluster.length == group.iloc[0].length)].sort_values(by='reads', ascending=False)

        for _k, v in consensus.iterrows():

            if v.seq in ref_precursor:

                mature_read = v.seq

                break

            if mature_read is None:

                mature_read = consensus.iloc[0].seq

        return mature_read
    
    def get_sequence(self, coordinates):
        """Get the reference sequence of a region.

        Args:
            coordinates (list): ``[chrom, start, end, strand]``. The reverse
                complement is returned for strand ``'-'``.

        Returns:
            str: The region sequence.
        """
        
        sequences = pf.Fasta(self.ref)

        seq = str(sequences[coordinates[0]][coordinates[1]:coordinates[2]])

        if coordinates[3] == '-':

            seq = str(Seq.Seq(seq).reverse_complement())

        return seq


