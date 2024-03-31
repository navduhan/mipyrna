#!/usr/bin/python
'''
Title: Class to remove non-coding RNA families from processed reads

Author: Naveen Duhan

'''
import os
import shutil
from waiting import wait
from Bio import SeqIO
from mipyrna.aligner import Bowtie_Aligner
from mipyrna import utility as mu
from mipyrna.reads import Read_process
import pkg_resources
from mipyrna.mirbase import MirBase

class Known_miRNA():

    def __init__(self, samples=None, species=None,species_type='plants', outdir=".", slurm = False):
        
        self.samples = samples
        self.species = species
        self.slurm = slurm
        self.species_type = species_type
        self.outdir = mu.make_directory(os.path.join(outdir, "known_miRNA_results"))

        return
    
    def _extract_known_unmapped(self, samples, rfam_samples):
        
        outfiltered = {}
        
        output = mu.make_directory(os.path.join(self.outdir,'known_filtered_reads'))
        
        for k, sample in samples.items():
            
            rfam_mapped_ids = rfam_samples[k][2]['read_name']
            
            rfam_mapped_ids.to_csv(f"{k}.txt", index=False)
    
            outfile=  f"{output}/{k}_filtered.fasta"
            
            execPATH = shutil.which('faSomeRecords')
            
            os.system(f"{execPATH} -exclude {sample[2]} {k}.txt {outfile}")
            
            os.remove(f"{k}.txt")
            
            outfiltered[k]=[sample[0],sample[1],outfile]
            
        return outfiltered
    
    def align_rfam(self, cpu=8,mem=20):

        
        mir = MirBase(species=self.species, speciesClass=self.species_type, outdir= self.outdir)

        mirbase_data = mir.get_species()
        
        aln = Bowtie_Aligner(ref_genome=mirbase_data['mature'], outdir=self.outdir, slurm=self.slurm)
        
        aln.build_index(cpu=cpu)
        
        results = aln.run_alignment(samplesDict=self.samples, cpu=cpu,mem=mem)
                
        ncRNA_reads = Read_process().aligned_reads(samples=results, alignType='SAM')
        
        for key, row in ncRNA_reads.items():
        
            row[2].to_csv(f"{self.outdir}/{key}_known.txt", sep="\t", index=False)

        return  self._extract_rfam_unmapped(samples=self.samples, rfam_samples=ncRNA_reads)  