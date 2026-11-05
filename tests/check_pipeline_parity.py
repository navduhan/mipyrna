# Compares package features and predictions with the training code.
# Run from the repository root: python mipyrna/tests/check_pipeline_parity.py [n]
import os
import sys

import numpy as np
from Bio import SeqIO

root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(root, "mipyrna"))
sys.path.insert(0, os.path.join(root, "scripts-linux"))

from features import Precursor_Features as TrainFeatures  # pyright: ignore[reportMissingImports]

from mipyrna import model_store as ms
from mipyrna.features import Precursor_Features as PkgFeatures

n = int(sys.argv[1]) if len(sys.argv) > 1 else 60
data = os.path.join(root, "data", "plants_clean")

features = {}
mismatches = 0
for name in ("positive", "negative"):
    fasta = os.path.join(data, f"independent_test_{name}.fasta")
    seqs = [str(r.seq) for r in SeqIO.parse(fasta, "fasta")][:n]
    features[name] = []
    for s in seqs:
        a, b = PkgFeatures(s).generate_features(), TrainFeatures(s).generate_features()
        if a is None or b is None:
            mismatches += (a is None) != (b is None)
        elif list(a) != list(b) or not np.allclose(list(a.values()), list(b.values())):
            mismatches += 1
        else:
            features[name].append(a)

print(f"feature mismatches: {mismatches}")
pos = ms.predict(features["positive"], "plants") >= 0.5
neg = ms.predict(features["negative"], "plants") >= 0.5
print(f"positives called miRNA: {pos.sum()}/{len(pos)}")
print(f"negatives called miRNA: {neg.sum()}/{len(neg)}")
assert mismatches == 0 and pos.mean() > 0.9 and neg.mean() < 0.1
