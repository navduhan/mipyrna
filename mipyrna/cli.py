"""miPyRNA Command-Line Interface

This module defines the ``mipyrna`` command and dispatches its sub-commands to
the workflow, prediction, target, enrichment and report modules.

Features:
    - ``workflow``: end-to-end small RNA-seq analysis
    - ``predict``: score precursor sequences from a FASTA file
    - ``targets``: miRanda target prediction
    - ``enrich``: g:Profiler functional enrichment
    - ``report``: PDF report generation
    - ``--check`` and ``--version`` global options

Configuration:
    Options are given on the command line. Results of the sub-commands are
    printed as JSON, or written to the file given by ``--output`` for ``predict``.
    The workflow options can also be read from a ``[workflow]`` section of an INI
    file with ``mipyrna -c config.ini`` (or ``mipyrna workflow -c config.ini``);
    flags given on the command line override the file. See ``config.ini.example``.

Dependencies:
    - argparse: Argument parsing
    - mipyrna.workflow, predict, targets, annotation, report, check: Sub-command implementations

Classes / Functions / Exceptions:
    - build_parser: Build the argument parser.
    - main: Run the command-line interface.

:Author: Naveen Duhan
"""

import argparse
import configparser
import json
import sys
from importlib.metadata import PackageNotFoundError, version


def _version():
    """Get the installed miPyRNA version.

    Returns:
        str: The version string, or ``"unknown"`` if the package is not installed.
    """
    try:
        return version("miPyRNA")
    except PackageNotFoundError:
        return "unknown"


def load_config(path):
    """Read the ``[workflow]`` section of a config file.

    Args:
        path (str): Path to an INI file.

    Returns:
        dict: Option names (underscored) mapped to their raw string values.

    Raises:
        SystemExit: If the file cannot be read or has no ``[workflow]`` section.
    """
    cp = configparser.ConfigParser()
    if not cp.read(path):
        raise SystemExit(f"Cannot read config file: {path}")
    if not cp.has_section("workflow"):
        raise SystemExit(f"No [workflow] section in {path}")
    return {k.replace("-", "_"): v for k, v in cp.items("workflow")}


def build_parser(defaults=None):
    """Build the argument parser with its sub-commands.

    Args:
        defaults (dict, optional): ``workflow`` option values read from a config file.

    Returns:
        argparse.ArgumentParser: Parser with the ``workflow``, ``predict``,
        ``targets``, ``enrich`` and ``report`` sub-commands.
    """
    parser = argparse.ArgumentParser(
        prog="mipyrna",
        description="miPyRNA: small RNA-seq analysis toolkit",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"miPyRNA {_version()}",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check external tools, library versions and bundled models, then exit.",
    )

    parser.add_argument("-c", "--config-file", help="INI file with a [workflow] section; runs the workflow")

    sub = parser.add_subparsers(dest="command")

    wf = sub.add_parser("workflow", help="Run the end-to-end workflow")
    wf.add_argument("--input-file")
    wf.add_argument("--input-path")
    wf.add_argument("--genome")
    wf.add_argument("--species")
    wf.add_argument("--species-type", default="plants")
    wf.add_argument("--outdir", default="mipyrna_results")
    wf.add_argument("--paired", action="store_true")
    wf.add_argument("--skip-qc", action="store_true")
    wf.add_argument("--skip-trim", action="store_true")
    wf.add_argument("--mrna-file")
    wf.add_argument("--enrichment-organism")
    wf.add_argument("--run-enrichment", action="store_true")
    wf.add_argument("--cpu", type=int, default=8, help="Total cores for the run")
    wf.add_argument("--jobs", type=int, help="Samples to run at once; each gets cpu // jobs cores (default: one per 4 cores)")
    wf.add_argument("--mem", type=int, default=20)
    wf.add_argument("--adapter", help="3' adapter: kit name or sequence (default: auto-detect)")
    wf.add_argument("--dryrun", action="store_true", help="Show the QC, trimming and alignment commands without running them")
    wf.add_argument("--resume", action="store_true", help="Reuse --outdir and skip the stages and samples an earlier run finished")
    wf.add_argument("--require-star", action="store_true", help="Keep only novel miRNAs with star reads")
    wf.add_argument("--min-star-reads", type=int, default=1, help="Star reads needed with --require-star (default: 1)")
    wf.add_argument("--trimmer", choices=["trimmomatic", "flexbar", "trim_galore"], default="trimmomatic")

    wf.add_argument("-c", "--config-file", default=argparse.SUPPRESS, help="INI file with a [workflow] section")
    if defaults:
        known = {a.dest: a for a in wf._actions}
        values = {}
        for key, raw in defaults.items():
            act = known.get(key)
            if act is None or key in ("help", "config_file"):
                raise SystemExit(f"Unknown option in config file: {key}")
            if act.nargs == 0:
                values[key] = raw.strip().lower() in ("1", "true", "yes", "on")
            else:
                values[key] = act.type(raw) if callable(act.type) else raw
        wf.set_defaults(**values)

    pred = sub.add_parser("predict", help="Score a FASTA of precursor sequences")
    pred.add_argument("--fasta", required=True)
    pred.add_argument("--species-type", default="plants", choices=("plants", "animals"))
    pred.add_argument("--threshold", type=float, default=0.5)
    pred.add_argument("--output", default="predictions.tsv")

    tgt = sub.add_parser("targets", help="Run miRanda target prediction")
    tgt.add_argument("--mirna-file", required=True)
    tgt.add_argument("--mrna-file", required=True)
    tgt.add_argument("--config")
    tgt.add_argument("--outdir", default=".")

    enr = sub.add_parser("enrich", help="Run functional enrichment using g:Profiler")
    enr.add_argument("--targets-file", required=True)
    enr.add_argument("--organism", required=True)
    enr.add_argument("--outdir", default=".")

    rep = sub.add_parser("report", help="Generate a PDF analysis report")
    rep.add_argument("--outdir", default=".")
    rep.add_argument("--metadata-json")

    return parser


def main(argv=None):
    """Run the command-line interface.

    Args:
        argv (list, optional): Command-line arguments. Defaults to ``sys.argv[1:]``.

    Returns:
        int: Process exit code (0 on success, otherwise the code from ``--check``).
    """
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("-c", "--config-file")
    config_file = pre.parse_known_args(argv)[0].config_file
    parser = build_parser(load_config(config_file) if config_file else None)
    args = parser.parse_args(argv)
    if args.command is None and args.config_file:
        args = parser.parse_args(["workflow", *(argv if argv is not None else sys.argv[1:])])
    if args.check:
        from mipyrna.check import run_check
        return run_check()

    if args.command == "workflow":
        missing = [o for o in ("input_file", "input_path", "genome", "species") if not getattr(args, o)]
        if missing:
            parser.error("workflow requires " + ", ".join("--" + o.replace("_", "-") for o in missing))
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
            jobs=args.jobs,
            mem=args.mem,
            adapter=args.adapter,
            require_star=args.require_star,
            min_star_reads=args.min_star_reads,
            trimmer=args.trimmer,
            dryrun=args.dryrun,
            resume=args.resume
        )
        print(json.dumps(meta, indent=2))
        return 0

    if args.command == "predict":
        from mipyrna.predict import predict_fasta
        result = predict_fasta(args.fasta, args.species_type, args.threshold)
        result.to_csv(args.output, sep="\t", index=False)
        print(f"{len(result)} sequences written to {args.output}")
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
            with open(args.metadata_json) as fh:
                metadata = json.load(fh)
        report, meta = generate_report(outdir=args.outdir, metadata=metadata)
        print(json.dumps({"report": report, "metadata": meta}, indent=2))
        return 0

    parser.print_help()
    return 0

