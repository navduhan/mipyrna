"""
Known miRNA Identification Module

This module identifies known miRNAs in small RNA-seq data by aligning processed
reads to the combined mature miRNA reference of a species (miRBase plus
MirGeneDB or PmiREN) with Bowtie, and by writing the reads that do not match.

Features:
    - Builds the species mature miRNA reference through `get_known_reference`
    - Bowtie indexing and alignment of reads for every sample
    - Per-sample tables of reads aligned to known miRNAs
    - Per-sample FASTA files of reads not matching any known miRNA
    - Optional SLURM execution of the alignment step

Configuration:
    Set through the `Known_miRNA` constructor: samples dictionary, miRBase species
    code, species type (plants or animals), output directory and SLURM flag.
    Results are written to `known_miRNA_results` inside the output directory.

Dependencies:
    - Bowtie: Used through `mipyrna.aligner.Bowtie_Aligner`
    - mipyrna.knowndb: Reference construction
    - mipyrna.reads: Alignment parsing
    - mipyrna.utility: Directory and FASTA helpers

Classes / Functions / Exceptions:
    - Known_miRNA: Aligns reads to known miRNAs and extracts the remaining reads.

:Author: Naveen Duhan
"""

import os

from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner
from mipyrna.knowndb import get_known_reference
from mipyrna.reads import Read_process


class Known_miRNA:
    """
    Find known miRNAs by aligning reads to the known mature miRNA reference.
    """

    def __init__(self, samples=None, species=None,species_type='plants', outdir=".", slurm = False):
        """
        Initialize the known miRNA finder.

        Args:
            samples (dict, optional): Sample name mapped to a list of sample information
                (first two entries are carried through, the third is the reads FASTA).
            species (str, optional): miRBase species code, e.g. ath.
            species_type (str): 'plants' or 'animals'. Defaults to 'plants'.
            outdir (str): Base output directory. Defaults to the current directory.
            slurm (bool): Run the alignment through SLURM. Defaults to False.
        """
        
        self.samples = samples
        self.species = species
        self.slurm = slurm
        self.species_type = species_type
        self.outdir = mu.make_directory(os.path.join(outdir, "known_miRNA_results"))

    
    def _extract_known_unmapped(self, samples, rfam_samples):
        """
        Write the reads that did not match a known miRNA to a FASTA file per sample.

        Args:
            samples (dict): Sample name mapped to sample information, where the third
                entry is the input reads FASTA.
            rfam_samples (dict): Sample name mapped to aligned read information, where
                the third entry is a table with a 'read_name' column.

        Returns:
            dict: Sample name mapped to [first entry, second entry, filtered FASTA path].
        """
        
        outfiltered = {}
        
        output = mu.make_directory(os.path.join(self.outdir,'known_filtered_reads'))
        
        for k, sample in samples.items():
            
            rfam_mapped_ids = rfam_samples[k][2]['read_name']
            
            outfile = f"{output}/{k}_filtered.fasta"

            mu.exclude_fasta_records(sample[2], rfam_mapped_ids, outfile)

            outfiltered[k]=[sample[0],sample[1],outfile]
            
        return outfiltered
    
    def align_known(self, cpu=8,mem=20):
        """
        Align reads to the species mature miRNAs and return the remaining reads per sample.

        Also writes the aligned reads of each sample to `<sample>_known.txt`.

        Args:
            cpu (int): Number of threads. Defaults to 8.
            mem (int): Memory in GB for the alignment job. Defaults to 20.

        Returns:
            dict: Sample name mapped to [first entry, second entry, filtered FASTA path]
            for reads not matching known miRNAs.
        """

        
        mirbase_file = get_known_reference(self.species, self.species_type, self.outdir)['mature']

        aln = Bowtie_Aligner(ref_genome=mirbase_file, outdir=self.outdir, slurm=self.slurm)
        
        aln.build_index(cpu=cpu)
        
        results = aln.run_alignment(samplesDict=self.samples, cpu=cpu,mem=mem)
                
        ncRNA_reads = Read_process().aligned_reads(samples=results, alignType='SAM')
        
        for key, row in ncRNA_reads.items():
        
            row[2].to_csv(f"{self.outdir}/{key}_known.txt", sep="\t", index=False)

        return  self._extract_known_unmapped(samples=self.samples, rfam_samples=ncRNA_reads)  
