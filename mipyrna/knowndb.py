"""
Known miRNA Reference Database Module

This module builds the known mature miRNA reference of a species by combining
miRBase with MirGeneDB (animals) or PmiREN (plants) into one non-redundant FASTA,
and records which databases contain each sequence.

Features:
    - Mature miRNA retrieval from miRBase
    - Additional retrieval from PmiREN (plants) or MirGeneDB (animals)
    - Sequence-level de-duplication with per-miRNA database provenance
    - Summary counts of miRNAs per database and shared between databases
    - Non-fatal handling of failed extra-database downloads

Configuration:
    Inputs are a miRBase species code, species class (plants or animals) and an
    output directory. Writes `<species>_known_mature.fa`, the downloaded database
    FASTA files and `known_reference_sources.tsv` to the output directory.

Dependencies:
    - pandas: Species lookup and summaries
    - Biopython: FASTA parsing
    - mipyrna.mirbase: miRBase download
    - mipyrna.utility: File download helper

Classes / Functions / Exceptions:
    - _species_name: Latin name of a species code.
    - _read_mature: Read mature sequences from a FASTA file.
    - get_known_reference: Build the combined mature miRNA FASTA.
    - summarize_sources: Count miRNAs per source database.

:Author: Naveen Duhan
"""

import importlib.resources as ir
import os

import pandas as pd

from mipyrna import utility as mu
from mipyrna.logger import MiPyRNALogger
from mipyrna.mirbase import MirBase

log = MiPyRNALogger(mode='a', log='knowndb')

MIRGENEDB_URL = "https://mirgenedb.org/fasta/ALL?mat=1"
PMIREN_URL = "https://www.pmiren.com/ftp-download/v3.2/miRNA_info/{name}/{name}_mature.fa"


def _species_name(species):
    """
    Return the Latin name of a miRBase species code.

    Args:
        species (str): miRBase species code, e.g. ath.

    Returns:
        str: Latin species name from the bundled organisms table.

    Raises:
        IndexError: If the species code is not in the organisms table.
    """
    df = pd.read_csv(str(ir.files("mipyrna").joinpath("data/organisms.txt")), sep="\t")
    df.columns = ['Species', 'Division', 'Name', 'Tree', 'TaxID']
    return df.loc[df['Species'] == species, 'Name'].iloc[0]


def _read_mature(fasta, prefix=None):
    """
    Read mature miRNA sequences from a FASTA file.

    Args:
        fasta (str): Path to the FASTA file.
        prefix (str, optional): Keep only records whose ID starts with this prefix.

    Returns:
        list: (id, sequence) tuples with sequences upper-cased and U converted to T.
    """
    return mu.read_fasta(fasta, prefix=prefix)


def get_known_reference(species, species_class, outdir='.'):
    """
    Write the combined mature miRNA FASTA of a species and return it with the databases of each ID.

    Sequences shared between databases are written once under the first ID seen.
    A failure to get the PmiREN or MirGeneDB data is logged and skipped.

    Args:
        species (str): miRBase species code, e.g. ath.
        species_class (str): 'plants' or 'animals'.
        outdir (str): Output directory. Defaults to the current directory.

    Returns:
        dict: 'mature' holds the FASTA path and 'sources' maps each miRNA ID to a
        comma-separated string of databases containing it.
    """
    records = {'miRBase': _read_mature(MirBase(species=species, speciesClass=species_class, outdir=outdir).get_species()['mature'])}

    try:
        if species_class.lower() == 'plants':
            name = _species_name(species).replace(' ', '_')
            path = mu.download_file(PMIREN_URL.format(name=name), os.path.join(outdir, 'pmiren_mature.fa'))
            records['PmiREN'] = _read_mature(path)
        elif species_class.lower() == 'animals':
            path = mu.download_file(MIRGENEDB_URL, os.path.join(outdir, 'mirgenedb_mature.fa'))
            records['MirGeneDB'] = _read_mature(path, prefix=f"{species.capitalize()}-")
    except Exception as exc:
        log.warning(f"Could not get the extra database for {species}: {exc}")

    by_seq = {}
    for database, items in records.items():
        log.info(f"{database}: {len(items)} mature sequences")
        for mid, seq in items:
            entry = by_seq.setdefault(seq, [mid, []])
            if database not in entry[1]:
                entry[1].append(database)

    mature = os.path.join(outdir, f'{species}_known_mature.fa')
    sources = {}
    records = []
    for seq, (mid, databases) in by_seq.items():
        if mid in sources:
            continue
        records.append((mid, seq))
        sources[mid] = ','.join(databases)
    mu.write_fasta(records, mature)

    with open(os.path.join(outdir, 'known_reference_sources.tsv'), 'w') as fp:
        fp.write('Known_ID\tDatabase\n')
        for mid, databases in sorted(sources.items()):
            fp.write(f'{mid}\t{databases}\n')

    return {'mature': mature, 'sources': sources}


def summarize_sources(databases):
    """
    Count miRNAs per database, with a Shared row for those found in more than one.

    Args:
        databases (iterable): Comma-separated database names for each miRNA.

    Returns:
        pandas.DataFrame: Table with 'database' and 'miRNAs' columns.
    """
    counts = pd.Series(databases).map(lambda d: 'Shared' if ',' in d else f"{d} only").value_counts()
    return counts.rename_axis('database').reset_index(name='miRNAs')
