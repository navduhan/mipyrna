#!/usr/bin/python
"""
title: miRNA target prediction utilities

Author: Naveen Duhan
"""
import os
import shutil
import subprocess
import pkg_resources
import pandas as pd

from mipyrna.logger import MiPyRNALogger
from mipyrna import utility as mu

log = MiPyRNALogger(mode='a', log='miRNA_target')


def _parse_miranda_output(raw_output, parsed_output):
    records = []
    current = {"miRNA": None, "target_id": None, "score": None, "energy": None}
    with open(raw_output, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">>"):
                tokens = line[2:].split()
                current = {
                    "miRNA": tokens[0] if len(tokens) > 0 else None,
                    "target_id": tokens[1] if len(tokens) > 1 else None,
                    "score": None,
                    "energy": None
                }
                if len(tokens) > 2:
                    for tok in tokens[2:]:
                        if tok.lower().startswith("score"):
                            try:
                                current["score"] = float(tok.split(":")[-1])
                            except Exception:
                                pass
                        if tok.lower().startswith("energy"):
                            try:
                                current["energy"] = float(tok.split(":")[-1])
                            except Exception:
                                pass
                records.append(current.copy())
            elif "Score" in line or "Energy" in line:
                if not records:
                    continue
                parts = line.replace("=", " ").replace(":", " ").split()
                for i, token in enumerate(parts):
                    tk = token.lower()
                    if tk == "score" and i + 1 < len(parts):
                        try:
                            records[-1]["score"] = float(parts[i + 1])
                        except Exception:
                            pass
                    if tk == "energy" and i + 1 < len(parts):
                        try:
                            records[-1]["energy"] = float(parts[i + 1])
                        except Exception:
                            pass

    df = pd.DataFrame(records).dropna(subset=["miRNA", "target_id"], how="any")
    if len(df) == 0:
        df = pd.DataFrame(columns=["miRNA", "target_id", "score", "energy"])
    df.to_csv(parsed_output, sep="\t", index=False)
    return df


def run_targets(miRNA_file=None, configFile=None, mRNA_file=None, slurm=False, outdir=".", dryrun=False):
    if miRNA_file is None or mRNA_file is None:
        raise ValueError("Both miRNA_file and mRNA_file are required")

    if configFile is not None:
        config = mu.parse_config_file(configFile)
    else:
        stream = pkg_resources.resource_stream('mipyrna', "param/miranda.ini")
        config = mu.parse_config_file(stream.name)
        log.info("Using default config file miranda.ini")

    miranda_config = config[list(config.keys())[0]]
    args = " ".join(miranda_config)

    out = "targets_raw"
    if os.path.exists(outdir):
        output1 = os.path.join(outdir, out)
        output = mu.make_directory(output1, dryrun=dryrun)
    else:
        output = mu.make_directory(out, dryrun=dryrun)

    execPATH = shutil.which("miranda")
    if execPATH is None:
        raise RuntimeError("miranda command not found in PATH")

    raw_out = os.path.join(output, "miranda_targets.txt")
    parsed_out = os.path.join(output, "miranda_targets.tsv")
    cmd = f"{execPATH} {miRNA_file} {mRNA_file} {args} -out {raw_out}"

    if dryrun:
        return {"command": cmd, "raw_output": raw_out, "parsed_output": parsed_out}

    with open(os.path.join(output, "miranda.out"), "w+") as fout:
        with open(os.path.join(output, "miranda.err"), "w+") as ferr:
            rc = subprocess.call(cmd, shell=True, stdout=fout, stderr=ferr)
            if rc != 0:
                raise RuntimeError(f"miRanda failed with exit code {rc}")

    parsed = _parse_miranda_output(raw_out, parsed_out)
    summary = parsed.groupby("miRNA", as_index=False).agg(target_count=("target_id", "nunique"))
    summary_out = os.path.join(output, "miranda_summary.tsv")
    summary.to_csv(summary_out, sep="\t", index=False)

    return {
        "raw_output": raw_out,
        "parsed_output": parsed_out,
        "summary_output": summary_out,
        "count_predictions": int(len(parsed))
    }
