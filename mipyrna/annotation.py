"""Functional Annotation Module

This module runs functional enrichment of miRNA target genes against the
g:Profiler web service.

Features:
    - Enrichment of GO:BP, GO:MF, GO:CC and KEGG terms by default
    - Configurable annotation sources and p-value threshold
    - Writes an empty table when there are no target genes

Configuration:
    Input is a tab-separated targets file with a ``target_id`` column and a
    g:Profiler organism code. Output is ``functional_annotation/enrichment_results.tsv``.

Dependencies:
    - requests: Calls the g:Profiler REST API (internet access required)
    - pandas: Table handling

Classes / Functions / Exceptions:
    - run_enrichment: Run enrichment and write the results table.

:Author: Naveen Duhan
"""

import os

import pandas as pd
import requests

from mipyrna import utility as mu
from mipyrna.logger import MiPyRNALogger

log = MiPyRNALogger(mode='a', log='annotation')


def run_enrichment(targets_file=None, organism=None, outdir=".", sources=None, user_threshold=0.05):
    """Run g:Profiler enrichment on the target genes.

    Writes ``functional_annotation/enrichment_results.tsv`` in ``outdir``.

    Args:
        targets_file (str): Tab-separated targets file with a ``target_id`` column.
        organism (str): g:Profiler organism code, for example ``athaliana``.
        outdir (str): Output directory. Defaults to the current directory.
        sources (list, optional): Annotation sources. Defaults to GO:BP, GO:MF, GO:CC and KEGG.
        user_threshold (float): P-value threshold. Defaults to 0.05.

    Returns:
        dict: ``{'enrichment_results': path, 'count_terms': number of terms}``.

    Raises:
        ValueError: If ``targets_file`` or ``organism`` is missing, or the file
            has no ``target_id`` column.
        requests.HTTPError: If the g:Profiler request fails.
    """
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
    log.info(f"Requesting g:Profiler enrichment for {len(genes)} genes ({organism})")
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

    log.info(f"Enrichment returned {len(out_df)} terms")
    return {"enrichment_results": out_path, "count_terms": len(out_df)}
