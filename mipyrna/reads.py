"""Read processing helpers.

Defines Read_process, which collapses, filters and converts small RNA reads
and extracts sequences from a reference.

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
        """Set the optional reference FASTA."""

        self.ref = reference


    def _collapse_process_sample(self, key, sample, outdir):
        """Collapse identical reads of one sample into a FASTA with read counts."""
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
    
    def collapse_reads(self, samples, outdir):
        """Collapse the reads of all samples."""
        outreads = {}

        with multiprocessing.Pool(processes=4) as pool:
            
            results = list(tqdm(pool.starmap(self._collapse_process_sample, [(key, sample, outdir) for key, sample in samples.items()]), total=len(samples), desc='Processing samples'))

        for key, result in results:
            outreads[key] = result

        return outreads

    def process_sample(self, key, sample, alignType):
        """Parse the alignments of one sample (SAM or Bowtie output) into a dataframe."""
        areads = []
        count = 1

        if alignType == 'SAM':
            for r in HTSeq.SAM_Reader(str(sample[2])):
                try:
                    if r.aligned:
                        seq = r.read.seq.decode()
                        readcounts = r.read.name.split("_")
                        areads.append([r.read.name, seq, r.iv.chrom, r.iv.strand, r.iv.start, r.iv.end, len(str(seq)), readcounts[1], int(readcounts[2].replace('x', '')), count])
                        count += 1
                        
                except Exception as exc:
                    log.debug(f"Skipping unparsable alignment {r.read.name}: {exc}")

        elif alignType == 'BWT':
            df = pd.read_csv(sample[2], sep="\t", header=None)
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

    def aligned_reads(self, samples, alignType='SAM'):
        """Parse the alignments of all samples."""
        outaligned = {}

        with multiprocessing.Pool(processes=4) as pool:
            results = list(tqdm(pool.starmap(self.process_sample, [(key, sample, alignType) for key, sample in samples.items()])))
        for key, result in results:
            outaligned[key] = result

        return outaligned
    
    def get_consensus_read(self, ref_precursor, aligned_cluster):
        """Return the most abundant read of a cluster, preferring one found in the precursor."""

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
        """Return the reference sequence for [chrom, start, end]."""
        
        sequences = pf.Fasta(self.ref)

        seq = str(sequences[coordinates[0]][coordinates[1]:coordinates[2]])

        if coordinates[3] == '-':

            seq = str(Seq.Seq(seq).reverse_complement())

        return seq


