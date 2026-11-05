#!/usr/bin/env python3

"""Bundled miPyRNA precursor models.

Chosen by validation MCC (train_*_results.txt); the held-out test was not used.
plants:  plants_clean fold 3
animals: animals_clean fold 5
"""

import hashlib
import importlib.resources as ir
from functools import cache
from typing import Any

import numpy as np

NUM_FEATURES = 129
MIRNA_CLASS = 0


def _models_dir():
    return ir.files("mipyrna").joinpath("data", "models")


def _checksums():
    lines = _models_dir().joinpath("SHA256SUMS").read_text().splitlines()
    return {l.split()[1]: l.split()[0] for l in lines if l.strip()}


@cache
def load_model_for(species_type) -> Any:
    if species_type not in ("plants", "animals"):
        raise ValueError(f"species_type must be 'plants' or 'animals', got {species_type!r}")
    from keras.models import load_model

    rel = f"{species_type}/model.keras"
    path = str(_models_dir().joinpath(rel))
    with open(path, "rb") as fh:
        if hashlib.sha256(fh.read()).hexdigest() != _checksums().get(rel):
            raise RuntimeError(f"checksum mismatch for {rel}")
    model = load_model(path, compile=False)
    assert model is not None
    return model


def model_info(species_type):
    rel = f"{species_type}/model.keras"
    return {rel: _checksums()[rel]}


def predict(feature_dicts, species_type):
    """P(miRNA) for each precursor feature dictionary."""
    if len(feature_dicts) == 0:
        return np.array([])
    x = np.array([list(d.values()) for d in feature_dicts], dtype=float)
    if x.shape[1] != NUM_FEATURES:
        raise ValueError(f"expected {NUM_FEATURES} features, got {x.shape[1]}")
    x = x.reshape(-1, 1, NUM_FEATURES, 1)
    model = load_model_for(species_type)
    return np.asarray(model.predict(x, verbose=0))[:, MIRNA_CLASS]
