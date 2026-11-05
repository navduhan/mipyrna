"""miRBase Sequence Module

This module downloads miRBase hairpin and mature sequences and subsets them
for a species or for the other species of the same clade.

Features:
    - Validation of the species code against the bundled miRBase organisms table
    - Download of all mature and hairpin sequences from miRBase
    - Extraction of the sequences of one species
    - Extraction of the sequences of all or all other species of a clade

Configuration:
    The species code (for example ``ath``), the clade (``plants``, ``animals``,
    ``virus`` or ``fungi``) and the output directory are set in the constructor.

Dependencies:
    - Biopython: FASTA parsing
    - pandas: Organisms table
    - Internet access to www.mirbase.org

Classes / Functions / Exceptions:
    - MirBase: Download and subset miRBase sequences.

:Author: Naveen Duhan
"""

import importlib.resources as ir
import os

import pandas as pd
from Bio import SeqIO

from mipyrna import utility as mu


class MirBase:
    """Download and subset miRBase sequences for a species or clade."""
    def __init__(self, species=None, speciesClass= None, outdir='.'):
        """Initialize the miRBase helper and validate the species code.

        Args:
            species (str): miRBase species code, for example ``'ath'``.
            speciesClass (str): Clade: ``'plants'``, ``'animals'``, ``'virus'`` or ``'fungi'``.
            outdir (str): Output directory. Defaults to the current directory.

        Raises:
            ValueError: If ``speciesClass`` is None or the species code is not a known miRBase species.
        """

        if speciesClass is None:
            raise ValueError("speciesClass is required")

        self.species = species
        self.species_class = speciesClass
        self.outdir = outdir
        if not self.check_organim():
            raise ValueError(f"Invalid species code: {self.species}")
            


    def get_organism(self):
        """Get the path of the bundled miRBase organisms table.

        Returns:
            str: Path to ``data/organisms.txt`` inside the package.
        """

        return str(ir.files("mipyrna").joinpath("data/organisms.txt"))

    def check_organim(self):
        """Check that the species code is a known miRBase species.

        Returns:
            bool: True if the species code is in the organisms table.
        """
        
        file = str(ir.files("mipyrna").joinpath("data/organisms.txt"))
        df = pd.read_csv(file, sep="\t")

        df.columns = ['Species', 'Division', 'Name', 'Tree', 'TaxID']

        return self.species in df['Species'].values.tolist()


    def get_mature(self):
        """Download all mature miRBase sequences.

        An existing ``mature.fa`` in the output directory is replaced.

        Returns:
            str: Path to the downloaded ``mature.fa``.
        """

        out = os.path.join(self.outdir, 'mature.fa')


        if os.path.exists(out):
            os.remove(out)

        mu.download_file("https://www.mirbase.org/download/mature.fa", out)

        # subprocess.call(f"gunzip {os.path.join(self.outdir, 'mature.fa.gz')}", shell=True)

        return out

    def get_hairpin(self):
        """Download all miRBase hairpin sequences.

        An existing ``hairpin.fa`` in the output directory is replaced.

        Returns:
            str: Path to the downloaded ``hairpin.fa``.
        """

        out = os.path.join(self.outdir, 'hairpin.fa')

        if os.path.exists(out):
            os.remove(out)

        mu.download_file("https://www.mirbase.org/download/hairpin.fa", out)

        # subprocess.call(f"gunzip {os.path.join(self.outdir, 'hairpin.fa.gz')}", shell=True)

        return out

    def get_species(self):
        """Write the mature and hairpin sequences of the species.

        The full downloads are removed after filtering.

        Returns:
            dict: ``{'mature': path, 'hairpin': path}`` of the species FASTA files.
        """

        species_mature = os.path.join(self.outdir, f'{self.species}_mature.fa')
        species_precursor = os.path.join(self.outdir, f'{self.species}_hairpin.fa')

        mature = self.get_mature()
        mature_fasta = SeqIO.parse(mature, 'fasta')

        with open(species_mature, 'w') as fp:

            for fasta in mature_fasta:
                if self.species in fasta.id:
                    fp.write(f">{fasta.id}\n{fasta.seq}\n")

            fp.close()

        hairpin = self.get_hairpin()
        hairpin_fasta = SeqIO.parse(hairpin, 'fasta')

        with open(species_precursor, 'w') as fp:

            for fasta in hairpin_fasta:
                if self.species in fasta.id:
                    fp.write(f">{fasta.id}\n{fasta.seq}\n")

            fp.close()

        os.remove(mature)
        os.remove(hairpin)
        

        return {'mature': species_mature, 'hairpin': species_precursor}

    def get_other_species(self):
        """Write the sequences of the other species in the same clade.

        The selected species is excluded. The downloads and the extracted
        organisms table are removed afterwards.

        Returns:
            dict: ``{'other_mature': path, 'other_hairpin': path}`` of the FASTA files.
        """

        org_file = self.get_organism()

        df = pd.read_csv(org_file, sep="\t")

        df.columns = ['Species', 'Division', 'Name', 'Tree', 'TaxID']
        
        df_sp = df.values.tolist()

        other_species = []

        kingdom = {'plants': 'Viridiplantae', 'animals':'Metazoa', 'virus':'Viruses', 'fungi': 'Mycetozoa' }

        for d in df_sp:
            if d[3].split(";")[0] == kingdom[self.species_class.lower()] and  d[0] != self.species:
                other_species.append(d[0])

        species_mature = os.path.join(self.outdir, f'{self.species_class}_mature.fa')
        species_precursor = os.path.join(self.outdir, f'{self.species_class}_hairpin.fa')

        mature = self.get_mature()
        mature_fasta = SeqIO.parse(mature, 'fasta')

        with open(species_mature, 'w') as fp:

            for fasta in mature_fasta:
                for sp in other_species:
                    if sp in fasta.id:
                        fp.write(f">{fasta.id}\n{fasta.seq}\n")

            fp.close()

        hairpin = self.get_hairpin()
        hairpin_fasta = SeqIO.parse(hairpin, 'fasta')

        with open(species_precursor, 'w') as fp:

            for fasta in hairpin_fasta:
                for sp in other_species:
                    if sp in fasta.id:
                        fp.write(f">{fasta.id}\n{fasta.seq}\n")

            fp.close()

        os.remove(mature)
        os.remove(hairpin)
        os.remove(org_file)

        return {'other_mature': species_mature, 'other_hairpin': species_precursor}
    
    def get_all_species(self):
        """Write the sequences of every species in the clade.

        The downloads and the extracted organisms table are removed afterwards.

        Returns:
            dict: ``{'other_mature': path, 'other_hairpin': path}`` of the FASTA files.
        """

        org_file = self.get_organism()

        df = pd.read_csv(org_file, sep="\t")

        df.columns = ['Species', 'Division', 'Name', 'Tree', 'TaxID']
        
        df_sp = df.values.tolist()

        other_species = []

        kingdom = {'plants': 'Viridiplantae', 'animals':'Metazoa', 'virus':'Viruses', 'fungi': 'Mycetozoa' }

        for d in df_sp:
            if d[3].split(";")[0] == kingdom[self.species_class.lower()]:
                other_species.append(d[0])

        species_mature = os.path.join(self.outdir, f'{self.species_class}_mature.fa')
        species_precursor = os.path.join(self.outdir, f'{self.species_class}_hairpin.fa')

        mature = self.get_mature()
        mature_fasta = SeqIO.parse(mature, 'fasta')

        with open(species_mature, 'w') as fp:

            for fasta in mature_fasta:
                for sp in other_species:
                    if sp in fasta.id:
                        fp.write(f">{fasta.id}\n{fasta.seq}\n")

            fp.close()

        hairpin = self.get_hairpin()
        hairpin_fasta = SeqIO.parse(hairpin, 'fasta')

        with open(species_precursor, 'w') as fp:

            for fasta in hairpin_fasta:
                for sp in other_species:
                    if sp in fasta.id:
                        fp.write(f">{fasta.id}\n{fasta.seq}\n")

            fp.close()

        os.remove(mature)
        os.remove(hairpin)
        os.remove(org_file)

        return {'other_mature': species_mature, 'other_hairpin': species_precursor}
