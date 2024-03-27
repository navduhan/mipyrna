from Bio import SeqIO, Seq
import HTSeq
import gzip
import os
import re
import numpy as np
import pandas as pd
import pyfaidx as pf
from tqdm import tqdm
import multiprocessing



class Read_process:

    def __init__(self,  reference=None):

        self.ref = reference

        return

    def _collapse_process_sample(self, key, sample, outdir):
        collapsed_reads = {}
        
        if sample[2].endswith(".gz"):
            file = gzip.open(sample[2], "rt")
        else:
            file = open(sample[2], "rt")
            
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
        outreads = {}

        with multiprocessing.Pool(processes=4) as pool:
            
            results = list(tqdm(pool.starmap(self._collapse_process_sample, [(key, sample, outdir) for key, sample in samples.items()]), total=len(samples), desc='Processing samples'))

        for key, result in results:
            outreads[key] = result

        return outreads

    def process_sample(slef, key, sample, alignType):
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
                        
                except Exception:
                    continue

        elif alignType == 'BWT':
            df = pd.read_csv(sample[2], sep="\t", header=None)
            df.columns = ['read', 'strand', 'chrom', 'start', 'mature', 'distance', 'num_align', 'info']

            for i, row in df.iterrows():
                if row.strand == "-":
                    row.mature = row.mature[::-1]
                    row.mature = row.mature.translate(str.maketrans("ACGTN", "TGCAN")) 
                gseq = list(row.mature)
                edit = ["m"] * len(row.mature)
                mm = 0
                if not pd.isna(row['info']):
                    data = row['info']
                    changes = row['info'].split(",")
                    for change in changes:
                        if re.match(r'(\d+):(\w+)>\w+', change):
                            mm += 1
                            pos, ont = change.split(':')
                            gseq[int(pos)] = ont.split(">")[1]  # Convert to lowercase
                            edit[int(pos)] = "M"
                id_parts = row.read.split()
                db_parts = row.read.split()
                areads.append([id_parts[0], len(row.mature), row.mature, row.chrom, int(row.start)+1, int(row.start)+len(row.mature), row.strand, mm, ''.join(edit), count, id_parts[0].split("_x")[0], int(id_parts[0].split("_x")[1])])
                count += 1

        return key, [sample[0], sample[1], pd.DataFrame(areads, columns=['read_name', 'seq', 'name', 'strand', 'start', 'end', 'length', 'read_id', 'reads', 'align_id'])]

    def aligned_reads(self, samples, alignType='SAM'):
        outaligned = {}

        with multiprocessing.Pool(processes=4) as pool:
            results = tqdm(pool.starmap(self.process_sample, [(key, sample, alignType) for key, sample in samples.items()]))
            print(results)
        for key, result in results:
            outaligned[key] = result

        return outaligned
    
    def get_consensus_read(self, ref_precursor, aligned_cluster):

        mature_read = None
        if aligned_cluster.iloc[0].strand == '+':
            seq_start = 'start'
        else:
            seq_start = 'end'
        # group = aligned_cluster.groupby([seq_start, 'length']).agg({'reads':np.sum}).sort_values(by='reads', ascending=False).reset_index()
        group = aligned_cluster.groupby([seq_start, 'length'])['reads'].agg('sum').sort_values(ascending=False).reset_index()
        consensus = aligned_cluster[(aligned_cluster[seq_start]==group.iloc[0][seq_start]) &(aligned_cluster.length == group.iloc[0].length)].sort_values(by='reads', ascending=False)

        for k, v in consensus.iterrows():

            if v.seq in ref_precursor:

                mature_read = v.seq

                break

            if mature_read is None:

                mature_read = consensus.iloc[0].seq

        return mature_read
    
    def get_sequence(self, coordinates):
        
        sequences = pf.Fasta(self.ref)

        seq = str(sequences[coordinates[0]][coordinates[1]:coordinates[2]])

        if coordinates[3] == '-':

            seq = str(Seq.Seq(seq).reverse_complement())

        return seq



