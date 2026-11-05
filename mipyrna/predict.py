#!/usr/bin/env python3

"""Precursor Prediction Module

This module scores a FASTA file of miRNA precursor sequences with the bundled
machine-learning model and reports the probability and call for each sequence.

Features:
    - Feature extraction for every FASTA record
    - Model scoring for plant or animal precursors
    - Skipped sequences are reported with the reason instead of being dropped

Configuration:
    Input is a FASTA path, the species type (``plants`` or ``animals``) and a
    probability threshold. Output is a pandas DataFrame.

Dependencies:
    - Biopython: FASTA parsing
    - pandas: Result table
    - mipyrna.features, mipyrna.model_store: Feature generation and model scoring

Classes / Functions / Exceptions:
    - predict_fasta: Score every sequence in a FASTA file.

:Author: Naveen Duhan
"""

import pandas as pd
from Bio import SeqIO

from mipyrna import model_store as ms
from mipyrna.features import Precursor_Features


def predict_fasta(fasta, species_type='plants', threshold=0.5):
    """Score the precursor sequences of a FASTA file.

    Args:
        fasta (str): Path to the FASTA file of precursor sequences.
        species_type (str): Model to use, ``'plants'`` or ``'animals'``. Defaults to ``'plants'``.
        threshold (float): Minimum P(miRNA) to call a sequence ``'miRNA'``. Defaults to 0.5.

    Returns:
        pandas.DataFrame: One row per sequence with columns ``id``, ``length``,
        ``P_miRNA``, ``call`` and ``skipped`` (the reason features could not be
        generated, otherwise None).
    """
    rows, feats, kept = [], [], []

    for record in SeqIO.parse(fasta, 'fasta'):
        generator = Precursor_Features(str(record.seq))
        features = generator.generate_features()
        row = {'id': record.id, 'length': len(record.seq), 'P_miRNA': None,
               'call': None, 'skipped': None}
        if features is None:
            row['skipped'] = generator.drop_reason or 'no features generated'
        else:
            feats.append(features)
            kept.append(len(rows))
        rows.append(row)

    for i, p in zip(kept, ms.predict(feats, species_type)):
        rows[i]['P_miRNA'] = float(p)
        rows[i]['call'] = 'miRNA' if p >= threshold else 'non-miRNA'

    return pd.DataFrame(rows)
