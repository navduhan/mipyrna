#!/usr/bin/env python3

"""Simple reproducible report generation for miPyRNA runs."""

import os
import json
import pandas as pd

from mipyrna import utility as mu


def _safe_count(path, sep="\t"):
    if path is None or not os.path.exists(path):
        return 0
    try:
        return len(pd.read_csv(path, sep=sep))
    except Exception:
        return 0


def generate_report(outdir=".", metadata=None):
    report_dir = mu.make_directory(os.path.join(outdir, "report"))
    report_path = os.path.join(report_dir, "analysis_report.md")
    meta_path = os.path.join(report_dir, "workflow_metadata.json")

    if metadata is None:
        metadata = {}

    with open(meta_path, "w") as jh:
        json.dump(metadata, jh, indent=2)

    known = metadata.get("known_file")
    novel = metadata.get("novel_file")
    targets = metadata.get("targets_file")
    enrichment = metadata.get("enrichment_file")
    deg = metadata.get("deg_file")

    known_count = _safe_count(known, sep="\t" if str(known).endswith(".txt") else ",")
    novel_count = _safe_count(novel, sep="\t" if str(novel).endswith(".txt") else ",")
    target_count = _safe_count(targets)
    enrich_count = _safe_count(enrichment)
    deg_count = _safe_count(deg)

    lines = [
        "# miPyRNA Analysis Report",
        "",
        "## Overview",
        "This report summarizes the miPyRNA workflow execution and generated artifacts.",
        "",
        "## Workflow Summary",
        f"- Run directory: `{os.path.abspath(outdir)}`",
        f"- Known miRNA entries: {known_count}",
        f"- Novel miRNA entries: {novel_count}",
        f"- Predicted miRNA-target interactions: {target_count}",
        f"- Functional enrichment terms: {enrich_count}",
        f"- Differential expression rows: {deg_count}",
        "",
        "## Key Files",
        f"- Metadata: `{meta_path}`",
        f"- Known miRNAs: `{known}`",
        f"- Novel miRNAs: `{novel}`",
        f"- Targets: `{targets}`",
        f"- Enrichment: `{enrichment}`",
        f"- Differential expression: `{deg}`",
        "",
        "## Notes",
        "miRNA discovery stage is miRDeep2-inspired and should be validated with dataset-level benchmarking before equivalence claims.",
        "",
    ]

    with open(report_path, "w") as fh:
        fh.write("\n".join(lines))

    return report_path, meta_path
