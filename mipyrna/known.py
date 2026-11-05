"""Identification of known miRNAs.

Defines Known_miRNA, which aligns processed reads to miRBase and reports
reads matching known mature miRNAs and their precursors.

:Author: Naveen Duhan
"""

import os
import shutil
import subprocess

from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner
from mipyrna.mirbase import MirBase
from mipyrna.reads import Read_process


class Known_miRNA:
    """Find known miRNAs by aligning reads to miRBase."""

    def __init__(self, samples=None, species=None,species_type='plants', outdir=".", slurm = False):
        """Set the samples, species and output directory."""
        
        self.samples = samples
        self.species = species
        self.slurm = slurm
        self.species_type = species_type
        self.outdir = mu.make_directory(os.path.join(outdir, "known_miRNA_results"))

    
    def _extract_known_unmapped(self, samples, rfam_samples):
        """Write the reads that did not match a known miRNA to a FASTA file per sample."""
        
        outfiltered = {}
        
        output = mu.make_directory(os.path.join(self.outdir,'known_filtered_reads'))
        
        for k, sample in samples.items():
            
            rfam_mapped_ids = rfam_samples[k][2]['read_name']
            
            rfam_mapped_ids.to_csv(f"{k}.txt", index=False)
    
            outfile=  f"{output}/{k}_filtered.fasta"
            
            execPATH = shutil.which('faSomeRecords')

            if execPATH is None:
                raise FileNotFoundError("faSomeRecords (UCSC kent tools) was not found in PATH")

            subprocess.run([execPATH, "-exclude", sample[2], f"{k}.txt", outfile], check=True)
            
            os.remove(f"{k}.txt")
            
            outfiltered[k]=[sample[0],sample[1],outfile]
            
        return outfiltered
    
    def align_known(self, cpu=8,mem=20):
        """Align reads to the species mature miRNAs and return the remaining reads per sample."""

        
        mir = MirBase(species=self.species, speciesClass=self.species_type, outdir= self.outdir)

        mirbase_data = mir.get_species()
        
        mirbase_file = mu.replace_U_with_T(mirbase_data['mature'], os.path.join(self.outdir, "mature_d2r.fa"))
        
        aln = Bowtie_Aligner(ref_genome=mirbase_file, outdir=self.outdir, slurm=self.slurm)
        
        aln.build_index(cpu=cpu)
        
        results = aln.run_alignment(samplesDict=self.samples, cpu=cpu,mem=mem)
                
        ncRNA_reads = Read_process().aligned_reads(samples=results, alignType='SAM')
        
        for key, row in ncRNA_reads.items():
        
            row[2].to_csv(f"{self.outdir}/{key}_known.txt", sep="\t", index=False)

        return  self._extract_known_unmapped(samples=self.samples, rfam_samples=ncRNA_reads)  
