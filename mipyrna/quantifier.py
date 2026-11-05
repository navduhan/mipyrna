"""Quantification of miRNAs for differential expression.

Counts aligned reads per miRNA and merges them into a count matrix.

:Author: Naveen Duhan
"""

import os
from functools import reduce

import HTSeq
import pandas as pd

from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner


class Quantifier:
    """Count reads per miRNA across samples."""
    
    def __init__(self, sRNA_file=None, samples=None, slurm=False, outdir="."):
        """Set the miRNA FASTA, samples and output directory."""
        
        self.sRNA = sRNA_file
        self.samples = samples
        self.slurm = slurm
        self.outdir = outdir
        
    
    def _align_reads_against_miRNAs(self, cpu, mem, output):
        """Align the sample reads to the miRNA sequences."""
        
        aln = Bowtie_Aligner(ref_genome=self.sRNA, outdir=output, slurm=self.slurm)
        aln.build_index()
        
        results= aln.run_alignment(samplesDict=self.samples, outType='SAM', cpu=cpu,mem=mem)
            
        return results
    
    def _get_counts(self, aligned_files):
        """Sum aligned read counts per miRNA for each sample."""
        all_counts = {}
        for k, sample in aligned_files.items():
            
            sam_file = HTSeq.SAM_Reader(str(sample[2]))
            
            aln_data = []
            
            for line in sam_file:
                if line.aligned == True:
                    aln_seq = line.read.seq.decode()
                    readcount = line.read.name.split("x")[1]
                    aln_data.append([aln_seq, line.read.name, int(readcount), line.iv.chrom,line.iv.start, line.iv.end, line.iv.strand, len(str(aln_seq))])
                    
            counts = pd.DataFrame(aln_data, columns=['sequence', 'read', 'counts', 'mature_id', 'start', 'end', 'strand', 'length' ])
            # check if all start are <=1
            filter_counts = counts[counts['start']<=1]
            
            agg_funcs = {'counts': 'sum', 'sequence': 'first', 'start': 'first', 'end': 'first', 'strand': 'first', 'length': 'first'}
            result_df = filter_counts.groupby('mature_id').agg(agg_funcs).reset_index()
    
            all_counts[k]=[sample[0], sample[1], result_df]
                    
        return all_counts
    

    def quantify_expression(self, cpu=8, mem=20):
        """Align, count and write the merged raw count matrix to Raw_read_counts.xlsx."""
        
        aligdir = mu.make_directory(os.path.join(self.outdir, "quant_align"))
        
        aligned_results = self._align_reads_against_miRNAs(cpu=cpu, mem=mem, output=aligdir)
        
        sample_counts = self._get_counts(aligned_files=aligned_results)
        
        sample_counts_dfs =[]
        

        for k, s in sample_counts.items():
            
            ucounts = s[2][['mature_id', 'counts']]
            ucounts.columns =['mature_id', k]
           
            sample_counts_dfs.append(ucounts)
        
        final_counts = reduce(lambda left,right:pd.merge(left,right, on=['mature_id'], how='outer'),sample_counts_dfs).fillna(0)
        final_counts['mature_id'] = final_counts['mature_id'].astype(str)
        final_counts.iloc[:, 1:] = final_counts.iloc[:, 1:].astype(int)
        final_counts.to_excel(f"{self.outdir}/Raw_read_counts.xlsx", index=False)
        
        return final_counts
