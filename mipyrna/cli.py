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
    - ``--slurm``: submit the workflow as a Slurm job that scores samples as a job array

Configuration:
    Options are given on the command line. Results of the sub-commands are
    printed as JSON, or written to the file given by ``--output`` for ``predict``.
    The workflow options can also be read from a ``[workflow]`` section of an INI
    file with ``mipyrna -c config.ini`` (or ``mipyrna workflow -c config.ini``);
    flags given on the command line override the file. A ``[slurm]`` section holds
    ``partition``, ``account``, ``time``, ``qos``, ``constraint``, ``reservation``, ``gres`` and ``extra``. See ``config.ini.example``.

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
import os
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


def load_slurm(path):
    """Read the optional ``[slurm]`` section of a config file.

    Args:
        path (str): Path to an INI file.

    Returns:
        dict: Scheduler option names mapped to their string values; empty without the section.
    """
    cp = configparser.ConfigParser()
    cp.read(path)
    return dict(cp.items("slurm")) if cp.has_section("slurm") else {}


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
    wf.add_argument("--cpu", type=int, help="Total cores for the run (default: 8, or 20 per sample with --slurm)")
    wf.add_argument("--jobs", type=int, help="Samples to run at once; each gets cpu // jobs cores (default: one per 4 cores)")
    wf.add_argument("--mem", type=int, default=20)
    wf.add_argument("--adapter", help="3' adapter: kit name or sequence (default: auto-detect)")
    wf.add_argument("--dryrun", action="store_true", help="Show the QC, trimming and alignment commands without running them")
    wf.add_argument("--resume", action="store_true", help="Reuse --outdir and skip the stages and samples an earlier run finished")
    wf.add_argument("--require-star", action="store_true", help="Keep only novel miRNAs with star reads")
    wf.add_argument("--weighted", action="store_true", help="Share each read across the miRNAs it maps to when counting")
    wf.add_argument("--read-cutoff", type=int, default=50, help="Minimum reads of a cluster to build candidates (default: 50)")
    wf.add_argument("--min-mature-reads", type=int, default=10, help="Minimum mature reads of a candidate (default: 10)")
    wf.add_argument("--min-star-reads", type=int, default=1, help="Star reads needed with --require-star (default: 1)")
    wf.add_argument("--max-inconsistent", type=float, default=0.1, help="Maximum fraction of hairpin reads outside the mature and star positions (default: 0.1)")
    wf.add_argument("--min-novel-score", type=float, default=1.0, help="Minimum heuristic score of a novel miRNA (default: 1.0)")
    wf.add_argument("--known-scan", action="store_true", help="Detect known miRNAs from reads aligned to known precursors and merge them into the known table")
    wf.add_argument("--known-min-reads", type=int, default=10, help="Minimum reads on a known precursor arm to call it detected (default: 10)")
    wf.add_argument("--known-5p-tolerance", type=int, default=2, help="Allowed offset of a read 5' end from the annotated known mature start (default: 2)")
    wf.add_argument("--known-unique", action="store_true", help="Count a read for the known scan only when it maps to one precursor")
    wf.add_argument("--rfam-filter", action="store_true", help="Remove reads that align to the bundled Rfam ncRNA set before genome alignment and clustering (the Rfam match is always reported)")
    wf.add_argument("--slurm", action="store_true", help="Submit the run as a Slurm job; samples are scored as a job array")
    wf.add_argument("--partition", help="Slurm partition")
    wf.add_argument("--account", help="Slurm account")
    wf.add_argument("--time", help="Slurm time limit of every job, for example 3-00:00:00")
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


def submit_controller(args, argv):
    """Submit the workflow itself as a Slurm job and print its job id.

    The job gets ``cpu // jobs`` cores, 20 by default, and runs the same command inside Slurm.

    Args:
        args (argparse.Namespace): Parsed workflow arguments.
        argv (list, optional): The original command-line arguments.

    Returns:
        int: Process exit code.
    """
    import shlex

    from mipyrna import utility as mu

    n = len(mu.read_input_file(infile=args.input_file, inpath=args.input_path, paired=args.paired)["samples"])
    cpu = mu.resolve_cpu(args.cpu, args.jobs, n, True)
    per_sample = max(1, cpu // max(1, min(args.jobs or n, n)))
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    command = " ".join(shlex.quote(a) for a in [sys.executable, "-m", "mipyrna", *(argv if argv is not None else sys.argv[1:])])
    job_id = mu.clusterRun(
        job_name="mipyrna",
        sout=os.path.join(outdir, "mipyrna_slurm"),
        serror=os.path.join(outdir, "mipyrna_slurm"),
        command=command,
        cpu=per_sample,
        mem=args.mem,
    )
    print(f"Submitted mipyrna as Slurm job {job_id}; output in {outdir}/mipyrna_slurm.out")
    return 0


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
        if args.slurm:
            from mipyrna import utility as mu
            mu.slurm_preflight()
            options = load_slurm(config_file) if config_file else {}
            options.update({k: getattr(args, k) for k in ("partition", "account", "time") if getattr(args, k)})
            mu.set_slurm(**options)
            if not os.environ.get("SLURM_JOB_ID") and not args.dryrun:
                return submit_controller(args, argv)
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
            read_cutoff=args.read_cutoff,
            min_mature_reads=args.min_mature_reads,
            min_star_reads=args.min_star_reads,
            weighted=args.weighted,
            max_inconsistent=args.max_inconsistent,
            min_novel_score=args.min_novel_score,
            known_scan=args.known_scan,
            known_min_reads=args.known_min_reads,
            known_5p_tolerance=args.known_5p_tolerance,
            known_unique=args.known_unique,
            rfam_filter=args.rfam_filter,
            slurm=args.slurm,
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

