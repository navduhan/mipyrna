#!/usr/bin/env python3

"""
Functional annotation and enrichment helper utilities.
"""

import os
import requests
import pandas as pd

from mipyrna.logger import MiPyRNALogger
from mipyrna import utility as mu

log = MiPyRNALogger(mode='a', log='annotation')


def run_enrichment(targets_file=None, organism=None, outdir=".", sources=None, user_threshold=0.05):
    if targets_file is None:
        raise ValueError("targets_file is required")
    if organism is None:
        raise ValueError("organism is required (for example: athaliana, mmusculus, hsapiens)")

    if sources is None:
        sources = ["GO:BP", "GO:MF", "GO:CC", "KEGG"]

    annotation_dir = mu.make_directory(os.path.join(outdir, "functional_annotation"))
    df = pd.read_csv(targets_file, sep="\t")
    if "target_id" not in df.columns:
        raise ValueError("targets_file must contain a 'target_id' column")

    genes = sorted({str(x) for x in df["target_id"].dropna().tolist() if str(x).strip()})
    if len(genes) == 0:
        empty = pd.DataFrame(columns=["source", "native", "name", "p_value", "term_size", "query_size", "intersection_size"])
        out_path = os.path.join(annotation_dir, "enrichment_results.tsv")
        empty.to_csv(out_path, sep="\t", index=False)
        return {"enrichment_results": out_path, "count_terms": 0}

    payload = {
        "organism": organism,
        "query": genes,
        "sources": sources,
        "user_threshold": user_threshold,
    }
    response = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", json=payload, timeout=120)
    response.raise_for_status()
    data = response.json()
    results = data.get("result", [])

    out_rows = []
    for row in results:
        out_rows.append({
            "source": row.get("source"),
            "native": row.get("native"),
            "name": row.get("name"),
            "p_value": row.get("p_value"),
            "term_size": row.get("term_size"),
            "query_size": row.get("query_size"),
            "intersection_size": row.get("intersection_size")
        })
    out_df = pd.DataFrame(out_rows)
    out_path = os.path.join(annotation_dir, "enrichment_results.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)

    return {"enrichment_results": out_path, "count_terms": int(len(out_df))}
