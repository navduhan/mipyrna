#!/usr/bin/python3

"""
title: mipyrna miRbase utility
Author: Naveen Duhan
"""

import wget
import os 
import subprocess
import pandas as pd
from Bio import SeqIO
import pkg_resources
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

class MirBase:
    def __init__(self, species=None, speciesClass= None, outdir='.'):

        self.species = species
        self.species_class = speciesClass
        self.outdir = outdir
        if not self.check_organim():

            print("Please provide a valid species")
            


    def get_organism(self):

        organism_file = os.path.join(self.outdir, 'organisms.txt.gz')

        out = os.path.join(self.outdir, 'organisms.txt')

        if os.path.exists(out):
            os.remove(out)

        wget.download("https://www.mirbase.org/ftp/CURRENT/organisms.txt.gz", organism_file)

        subprocess.call(f"gunzip {os.path.join(self.outdir, 'organisms.txt.gz')}", shell=True)
            
        return out

    def check_organim(self):
        
        file = pkg_resources.resource_filename('mipyrna', "data/organisms.txt")
        df = pd.read_csv(file, sep="\t")

        df.columns = ['Species', 'Division', 'Name', 'Tree', 'TaxID']

        if self.species in df['Species'].values.tolist():
            return True

        return True


    def get_mature(self):

        outfile = os.path.join(self.outdir, 'mature.fa')
        out = os.path.join(self.outdir, 'mature.fa')


        if os.path.exists(out):
            os.remove(out)

        wget.download("https://www.mirbase.org/download/mature.fa", out )

        # subprocess.call(f"gunzip {os.path.join(self.outdir, 'mature.fa.gz')}", shell=True)

        return out

    def get_hairpin(self):

        outfile = os.path.join(self.outdir, 'hairpin.fa')

        out = os.path.join(self.outdir, 'hairpin.fa')

        if os.path.exists(out):
            os.remove(out)

        wget.download("https://www.mirbase.org/download/hairpin.fa", out)

        # subprocess.call(f"gunzip {os.path.join(self.outdir, 'hairpin.fa.gz')}", shell=True)

        return out

    def get_species(self):

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
