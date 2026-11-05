#!/usr/bin/env python3

"""
Bundled miPyRNA Precursor Models

This module loads and runs the pre-trained Keras models shipped with miPyRNA
for scoring miRNA precursor candidates. Models were chosen by validation MCC
(train_*_results.txt); the held-out test set was not used.

Features:
    - Separate bundled models for plants (plants_clean fold 3) and animals (animals_clean fold 5)
    - SHA-256 integrity check of the model file before loading
    - Cached model loading (each model is loaded once per process)
    - Probability of the miRNA class for each precursor feature dictionary

Configuration:
    Models and the SHA256SUMS file are read from the packaged ``mipyrna/data/models``
    directory (``<species_type>/model.keras``). Input features are dictionaries with
    exactly NUM_FEATURES (129) values; the output is an array of P(miRNA).

Dependencies:
    - keras: Model loading and prediction
    - numpy: Feature array handling

Classes / Functions / Exceptions:
    - _models_dir: Packaged models directory
    - _checksums: Expected SHA-256 checksums of the packaged models
    - load_model_for: Load and verify the bundled model for plants or animals
    - model_info: Model file and checksum for the run report
    - predict: P(miRNA) for each precursor feature dictionary

:Author: Naveen Duhan
"""

import hashlib
import importlib.resources as ir
from functools import cache
from typing import Any

import numpy as np

NUM_FEATURES = 129
MIRNA_CLASS = 0


def _models_dir():
    """Return the packaged models directory.

    Returns:
        Traversable pointing to ``mipyrna/data/models``.
    """
    return ir.files("mipyrna").joinpath("data", "models")


def _checksums():
    """Return the SHA-256 checksums of the packaged models.

    Returns:
        Dictionary mapping relative model path to its expected SHA-256 hex digest,
        read from the SHA256SUMS file.
    """
    lines = _models_dir().joinpath("SHA256SUMS").read_text().splitlines()
    return {l.split()[1]: l.split()[0] for l in lines if l.strip()}


@cache
def load_model_for(species_type) -> Any:
    """Load the bundled model for plants or animals after verifying its checksum.

    The result is cached, so each species type is loaded once per process.

    Args:
        species_type: Either ``'plants'`` or ``'animals'``.

    Returns:
        The loaded Keras model (uncompiled).

    Raises:
        ValueError: If species_type is not ``'plants'`` or ``'animals'``.
        RuntimeError: If the model file checksum does not match SHA256SUMS.
    """
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
    """Return the model file and its checksum for the run report.

    Args:
        species_type: Either ``'plants'`` or ``'animals'``.

    Returns:
        Dictionary mapping the relative model path to its SHA-256 checksum.

    Raises:
        KeyError: If no checksum is listed for the model.
    """
    rel = f"{species_type}/model.keras"
    return {rel: _checksums()[rel]}


def predict(feature_dicts, species_type):
    """Score precursors with the bundled model.

    Args:
        feature_dicts: Sequence of feature dictionaries, each with NUM_FEATURES values.
        species_type: Either ``'plants'`` or ``'animals'``.

    Returns:
        NumPy array with P(miRNA) for each precursor (empty if no input).

    Raises:
        ValueError: If the number of features per precursor is not NUM_FEATURES,
            or species_type is invalid.
    """
    if len(feature_dicts) == 0:
        return np.array([])
    x = np.array([list(d.values()) for d in feature_dicts], dtype=float)
    if x.shape[1] != NUM_FEATURES:
        raise ValueError(f"expected {NUM_FEATURES} features, got {x.shape[1]}")
    x = x.reshape(-1, 1, NUM_FEATURES, 1)
    model = load_model_for(species_type)
    return np.asarray(model.predict(x, verbose=0))[:, MIRNA_CLASS]
