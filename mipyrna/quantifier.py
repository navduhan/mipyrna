#!/usr/bin/python

'''
Title: Quantifier module for differential expression of miRNAs

Author: Naveen Duhan
'''

import os
from waiting import wait
from Bio import SeqIO
import HTSeq
import pandas as pd
from functools import reduce
from mipyrna.aligner import Bowtie_Aligner
from mipyrna import utility as mu
from mipyrna.reads import Read_process


class Quantifier():
    
    def __init__(self, sRNA_file=None, samples=None, slurm=False, outdir="."):
        
        self.sRNA = sRNA_file
        self.samples = samples
        self.slurm = slurm
        self.outdir = outdir
        
        return
    
    def _align_reads_against_miRNAs(self, cpu, mem, output):
        
        aln = Bowtie_Aligner(ref_genome=self.sRNA, outdir=output, slurm=self.slurm)
        aln.build_index()
        
        results= aln.run_alignment(samplesDict=self.samples, outType='SAM', cpu=cpu,mem=mem)
            
        return results
    
    def _get_counts(self, aligned_files):
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
        
        aligdir = mu.make_directory(os.path.join(self.outdir, "quant_align"))
        
        aligned_results = self._align_reads_against_miRNAs(cpu=cpu, mem=mem, output=aligdir)
        
        sample_counts = self._get_counts(aligned_files=aligned_results)
        
        cc
        
        return final_counts