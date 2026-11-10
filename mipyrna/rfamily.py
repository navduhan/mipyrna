"""Rfam Filtering Module

This module screens reads against bundled Rfam ncRNA sequences and can write
the reads that did not match for optional filtering.

Features:
    - Extraction of the bundled Rfam ncRNA FASTA at runtime
    - Bowtie alignment of reads to Rfam sequences
    - FASTA output of reads that did not match Rfam

Configuration:
    Samples and the output directory are set in the constructor. Results go
    to ``ncRNA_results`` in the output directory.

Dependencies:
    - Bowtie via mipyrna.aligner
    - mipyrna.reads: Alignment parsing

Classes / Functions / Exceptions:
    - FilterRNAfamlies: Remove reads from known ncRNA families.

:Author: Naveen Duhan
"""

import importlib.resources as ir
import os
import shutil
import tarfile

from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner
from mipyrna.reads import Read_process


class FilterRNAfamlies:
    """Remove reads from known non-coding RNA families."""

    def __init__(self, samples=None,  outdir=".", slurm=False):
        """Initialize the filter and create the output directory.

        Args:
            samples (dict): Sample name to ``[group, replicate, reads FASTA path]``.
            outdir (str): Base output directory. Defaults to the current directory.
            slurm (bool): Use SLURM for alignment. Defaults to False.
        """
        
        self.samples = samples or {}
        
        self.slurm = slurm
        
        self.outdir = mu.make_directory(os.path.join(outdir, "ncRNA_results"))

    
    def _extract_rfam_unmapped(self, samples, rfam_samples):
        """Write the reads that did not match Rfam to a FASTA file per sample.

        Args:
            samples (dict): Sample name to ``[group, replicate, reads FASTA path]``.
            rfam_samples (dict): Sample name to ``[group, replicate, Rfam alignment DataFrame]``.

        Returns:
            dict: Sample name to ``[group, replicate, filtered FASTA path]``.
        """
        
        outfiltered = {}
        
        output = mu.make_directory(os.path.join(self.outdir,'rfam_filtered_reads'))
        
        for k, sample in samples.items():
            
            rfam_mapped_ids = rfam_samples[k][2]['read_name']
            
            outfile = f"{output}/{k}_filtered.fasta"

            mu.exclude_fasta_records(sample[2], rfam_mapped_ids, outfile)

            outfiltered[k]=[sample[0],sample[1],outfile]
            
        return outfiltered
    
    def align_rfam(self, cpu=8, mem=20, jobs=1, filter_reads=True):
        """Align reads to bundled Rfam sequences and optionally filter the matches.

        Args:
            cpu (int): Threads for alignment. Defaults to 8.
            mem (int): Memory for alignment. Defaults to 20.
            jobs (int): Samples aligned at once. Defaults to 1.
            filter_reads (bool): Write unmatched FASTA files. Otherwise return the
                original samples unchanged. Defaults to True for existing callers.

        Returns:
            tuple: ``(ncRNA_reads, filtered_samples)``: the parsed Rfam alignments
            and the reads that remain after removing Rfam matches, both keyed by sample.
        """

        ncRNA_file = os.path.join(self.outdir, "ncRNA_rfam.fa")

        if not os.path.exists(ncRNA_file):
            with tarfile.open(str(ir.files("mipyrna").joinpath("data/ncRNA_rfam.tar.gz")), "r:gz") as tar:
                source = tar.extractfile("ncRNA_rfam.fa")
                if source is None:
                    raise ValueError("The bundled Rfam archive has no ncRNA FASTA")
                with source, open(ncRNA_file, "wb") as target:
                    shutil.copyfileobj(source, target)
        
        config = str(ir.files("mipyrna").joinpath("param/bowtie_rfam.ini"))
        aln = Bowtie_Aligner(
            ref_genome=ncRNA_file, configFile=config, outdir=self.outdir, slurm=self.slurm
        )
        
        aln.build_index(cpu=cpu, mem=mem)
        
        results = aln.run_alignment(
            samplesDict=self.samples, fileType="fasta", outType="SAM", cpu=cpu, mem=mem, jobs=jobs
        )
        
        ncRNA_reads = Read_process().aligned_reads(samples=results, alignType='SAM', jobs=jobs)
        
        remaining = (
            self._extract_rfam_unmapped(samples=self.samples, rfam_samples=ncRNA_reads)
            if filter_reads else self.samples
        )
        return ncRNA_reads, remaining
