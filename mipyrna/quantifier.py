"""miRNA Quantification Module

This module aligns reads to miRNA sequences and builds the raw count matrix
used for differential expression.

Features:
    - Alignment of collapsed reads to known and novel miRNAs with Bowtie
    - Counting of read abundance per miRNA, restricted to reads starting at position 1 or earlier
    - Merged count matrix across samples written to Excel

Configuration:
    Inputs are the miRNA FASTA, the collapsed read samples and the output
    directory. Output is ``Raw_read_counts.xlsx`` and a ``quant_align`` directory.

Dependencies:
    - HTSeq: SAM parsing
    - pandas: Count tables
    - Bowtie via mipyrna.aligner

Classes / Functions / Exceptions:
    - Quantifier: Count reads per miRNA across samples.

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
        """Initialize the quantifier.

        Args:
            sRNA_file (str): FASTA of miRNA sequences to quantify.
            samples (dict): Sample name to ``[group, replicate, collapsed reads path]``.
            slurm (bool): Use SLURM for alignment. Defaults to False.
            outdir (str): Output directory. Defaults to the current directory.
        """
        
        self.sRNA = sRNA_file
        self.samples = samples
        self.slurm = slurm
        self.outdir = outdir
        
    
    def _align_reads_against_miRNAs(self, cpu, mem, output):
        """Align the sample reads to the miRNA sequences.

        Args:
            cpu (int): Threads for alignment.
            mem (int): Memory for alignment.
            output (str): Output directory for the index and alignments.

        Returns:
            dict: Sample name to ``[group, replicate, SAM path]``.
        """
        
        aln = Bowtie_Aligner(ref_genome=self.sRNA, outdir=output, slurm=self.slurm)
        aln.build_index()
        
        results= aln.run_alignment(samplesDict=self.samples, outType='SAM', cpu=cpu,mem=mem)
            
        return results
    
    def _get_counts(self, aligned_files):
        """Sum aligned read counts per miRNA for each sample.

        Read abundance is taken from the ``x<count>`` suffix of the collapsed read name.

        Args:
            aligned_files (dict): Sample name to ``[group, replicate, SAM path]``.

        Returns:
            dict: Sample name to ``[group, replicate, counts DataFrame]``.
        """
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
