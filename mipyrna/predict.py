#!/usr/bin/env python3

"""Score a FASTA of miRNA precursor sequences with the bundled model."""

import pandas as pd
from Bio import SeqIO

from mipyrna import model_store as ms
from mipyrna.features import Precursor_Features


def predict_fasta(fasta, species_type='plants', threshold=0.5):
    """Return one row per sequence: P(miRNA), the call, or why it was skipped."""
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
