#!/usr/bin/env python3

"""Command-line entry point for miPyRNA."""

import argparse
import json


def build_parser():
    parser = argparse.ArgumentParser(
        prog="mipyrna",
        description="miPyRNA: small RNA-seq analysis toolkit",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="miPyRNA 0.1",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run a lightweight CLI sanity check and exit.",
    )

    sub = parser.add_subparsers(dest="command")

    wf = sub.add_parser("workflow", help="Run the end-to-end workflow")
    wf.add_argument("--input-file", required=True)
    wf.add_argument("--input-path", required=True)
    wf.add_argument("--genome", required=True)
    wf.add_argument("--species", required=True)
    wf.add_argument("--species-type", default="plants")
    wf.add_argument("--outdir", default="mipyrna_results")
    wf.add_argument("--paired", action="store_true")
    wf.add_argument("--skip-qc", action="store_true")
    wf.add_argument("--skip-trim", action="store_true")
    wf.add_argument("--mrna-file")
    wf.add_argument("--enrichment-organism")
    wf.add_argument("--run-enrichment", action="store_true")
    wf.add_argument("--cpu", type=int, default=8)
    wf.add_argument("--mem", type=int, default=20)

    tgt = sub.add_parser("targets", help="Run miRanda target prediction")
    tgt.add_argument("--mirna-file", required=True)
    tgt.add_argument("--mrna-file", required=True)
    tgt.add_argument("--config")
    tgt.add_argument("--outdir", default=".")

    enr = sub.add_parser("enrich", help="Run functional enrichment using g:Profiler")
    enr.add_argument("--targets-file", required=True)
    enr.add_argument("--organism", required=True)
    enr.add_argument("--outdir", default=".")

    rep = sub.add_parser("report", help="Generate a markdown analysis report")
    rep.add_argument("--outdir", default=".")
    rep.add_argument("--metadata-json")

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.check:
        print("miPyRNA CLI is available.")
        return 0

    if args.command == "workflow":
        from mipyrna.workflow import run_workflow
        meta = run_workflow(
            input_file=args.input_file,
            input_path=args.input_path,
            genome=args.genome,
            species=args.species,
            species_type=args.species_type,
            outdir=args.outdir,
            paired=args.paired,
            run_qc=not args.skip_qc,
            run_trim=not args.skip_trim,
            run_enrichment_stage=args.run_enrichment,
            mrna_file=args.mrna_file,
            enrichment_organism=args.enrichment_organism,
            cpu=args.cpu,
            mem=args.mem
        )
        print(json.dumps(meta, indent=2))
        return 0

    if args.command == "targets":
        from mipyrna.targets import run_targets
        result = run_targets(
            miRNA_file=args.mirna_file,
            mRNA_file=args.mrna_file,
            configFile=args.config,
            outdir=args.outdir
        )
        print(json.dumps(result, indent=2))
        return 0

    if args.command == "enrich":
        from mipyrna.annotation import run_enrichment
        result = run_enrichment(
            targets_file=args.targets_file,
            organism=args.organism,
            outdir=args.outdir
        )
        print(json.dumps(result, indent=2))
        return 0

    if args.command == "report":
        from mipyrna.report import generate_report
        metadata = None
        if args.metadata_json:
            with open(args.metadata_json, "r") as fh:
                metadata = json.load(fh)
        report, meta = generate_report(outdir=args.outdir, metadata=metadata)
        print(json.dumps({"report": report, "metadata": meta}, indent=2))
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
