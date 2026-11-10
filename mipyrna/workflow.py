#!/usr/bin/env python3

"""miPyRNA Workflow Module

This module orchestrates the end-to-end small RNA-seq analysis, from raw reads
to miRNA discovery, quantification, differential expression, target prediction
and reporting.

Features:
    - Optional FastQC before and after trimming
    - Trimming with Trimmomatic, Flexbar or Trim Galore
    - Read collapsing, Bowtie alignment and read clustering
    - Known and novel miRNA identification
    - Quantification and DESeq2 differential expression
    - Optional miRanda target prediction and g:Profiler enrichment
    - PDF report and metadata JSON

Configuration:
    Inputs are a sample sheet, the reads directory, the genome FASTA and the
    miRBase species code. Outputs are written under ``outdir``. Optional stages
    that fail are logged and skipped.

Dependencies:
    - Bowtie, FastQC, a read trimmer and miRanda (external tools)
    - R with DESeq2 (through rpy2)
    - mipyrna submodules for each pipeline stage

Classes / Functions / Exceptions:
    - run_workflow: Run the full pipeline.

:Author: Naveen Duhan
"""

import json
import os
import time

import pandas as pd

from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner
from mipyrna.annotation import run_enrichment
from mipyrna.check import provenance
from mipyrna.clustering import Cluster_reads
from mipyrna.differential import runDESeq2
from mipyrna.logger import MiPyRNALogger, set_log_file
from mipyrna.novel import Novel_miRNA
from mipyrna.quality import fastqcRun
from mipyrna.quantifier import Quantifier
from mipyrna.reads import Read_process
from mipyrna.report import generate_report
from mipyrna.targets import run_targets
from mipyrna.trimming import flexbarRun, trim_galoreRun, trimmomaticRun

log = MiPyRNALogger(mode='a', log='workflow')


def _stage_dir(output, name, dryrun=False):
    """Create (or reuse when resuming) a numbered stage directory in the run directory.

    Args:
        output (str): Run output directory.
        name (str): Stage directory name, for example ``'2.Alignment'``.
        dryrun (bool): If True, return the path without creating it.

    Returns:
        str: Stage directory path.
    """
    return mu.make_directory(os.path.join(output, name), dryrun=dryrun)


def _save_stage(output, stage, results):
    """Record a finished stage and its results so a resumed run can skip it.

    Args:
        output (str): Run output directory.
        stage (str): Stage name.
        results (dict): Sample name to ``[group, replicate, file(s)...]`` produced by the stage.
    """
    path = os.path.join(output, ".resume", f"{stage}.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(results, fh)


def _load_stage(resume, output, stage):
    """Get the saved results of a finished stage when resuming.

    Args:
        resume (bool): Whether the run is a resume.
        output (str): Run output directory.
        stage (str): Stage name.

    Returns:
        dict or None: The saved results, or None if the stage must run (not resuming,
        never finished, or an output file is missing).
    """
    path = os.path.join(output, ".resume", f"{stage}.json")
    if not resume or not os.path.exists(path):
        return None
    with open(path) as fh:
        results = json.load(fh)
    if all(os.path.exists(f) for v in results.values() for f in v[2:]):
        log.info(f"Resuming: {stage} already finished")
        return results
    return None


def _write_collapsed_fasta(samples, outdir, jobs=4):
    """Collapse the reads of every sample into FASTA files.

    Args:
        samples (dict): Sample name to ``[group, replicate, read files...]``.
        outdir (str): Output directory.
        jobs (int): Number of parallel processes. Defaults to 4.

    Returns:
        dict: Sample name to ``[group, replicate, collapsed FASTA path]``.
    """
    collapsed = Read_process().collapse_reads(samples=samples, outdir=outdir, jobs=jobs)
    return collapsed


def run_workflow(
    input_file,
    input_path,
    genome,
    species,
    species_type='plants',
    outdir=".",
    paired=False,
    run_qc=True,
    run_trim=True,
    run_enrichment_stage=False,
    mrna_file=None,
    enrichment_organism=None,
    cpu=None,
    mem=20,
    jobs=None,
    adapter=None,
    trimmer="trimmomatic",
    require_star=False,
    read_cutoff=50,
    min_mature_reads=10,
    min_star_reads=1,
    weighted=False,
    max_inconsistent=0.1,
    min_novel_score=1.0,
    slurm=False,
    dryrun=False,
    resume=False
):
    """Run the full pipeline from raw reads to differential expression, targets and report.

    Args:
        input_file (str): Sample sheet describing samples and comparisons.
        input_path (str): Directory containing the read files.
        genome (str): Reference genome FASTA.
        species (str): miRBase species code.
        species_type (str): ``'plants'`` or ``'animals'``. Defaults to ``'plants'``.
        outdir (str): Output directory. Defaults to the current directory.
        paired (bool): True for paired-end reads. Defaults to False.
        run_qc (bool): Run FastQC before and after trimming. Defaults to True.
        run_trim (bool): Trim reads. Defaults to True.
        run_enrichment_stage (bool): Run functional enrichment of targets. Defaults to False.
        mrna_file (str, optional): mRNA FASTA enabling target prediction.
        enrichment_organism (str, optional): g:Profiler organism code for enrichment.
        cpu (int, optional): Total cores. Defaults to 8, or 20 per sample with ``slurm``.
        mem (int): Memory in GB. Defaults to 20.
        jobs (int, optional): Samples to run at once. Derived from ``cpu`` when None.
        adapter (str, optional): 3' adapter name or sequence. Auto-detected when None.
        trimmer (str): ``'trimmomatic'``, ``'flexbar'`` or ``'trim_galore'``. Defaults to ``'trimmomatic'``.
        require_star (bool): Keep only novel miRNAs with star reads. Defaults to False.
        read_cutoff (int): Minimum reads of a cluster to build candidates. Defaults to 50.
        min_mature_reads (int): Minimum mature reads of a candidate. Defaults to 10.
        min_star_reads (int): Star reads needed when ``require_star`` is set. Defaults to 1.
        weighted (bool): Share each read across the miRNAs it maps to when counting. Defaults to False.
        max_inconsistent (float): Maximum fraction of the reads on a novel hairpin outside the mature and star positions. Defaults to 0.1.
        min_novel_score (float): Minimum heuristic score of a novel miRNA. Defaults to 1.0.
        slurm (bool): Run each sample of the heavy steps as an element of a Slurm job array with ``cpu // jobs`` cores.
            Trimming, FastQC and the alignments are submitted as arrays too; collapsing and clustering run in
            this process, which should itself be a Slurm job. Defaults to False.
        dryrun (bool): Log the QC, trimming and alignment commands without running them or
            creating any files, then stop. Defaults to False.
        resume (bool): Reuse ``outdir`` and skip the QC, trimming, collapsing and alignment
            stages and the scored samples that an earlier run finished. Defaults to False.

    Returns:
        dict: Run metadata with provenance and the paths of the generated files
        (counts, known, novel, DEG, targets, enrichment and report files). With
        ``dryrun`` only the output directory and the planned steps.
    """
    start = time.time()
    mu.set_resume(resume)
    output = mu.make_directory(outdir, dryrun=dryrun)
    if not dryrun:
        set_log_file(os.path.join(output, "mipyrna.log"))
    parsed = mu.read_input_file(infile=input_file, inpath=input_path, paired=paired)
    samples = parsed["samples"]
    combinations = parsed["combinations"]
    targets = parsed["targets"]
    metadata = {
        "outdir": output,
        "samples": {k: v[1:] for k, v in samples.items()},
        "comparisons": combinations,
        "parameters": {
            "input_file": os.path.abspath(input_file),
            "input_path": os.path.abspath(input_path),
            "genome": os.path.abspath(genome),
            "species": species,
            "species_type": species_type,
            "paired": paired,
            "trimmer": trimmer if run_trim else None,
            "adapter": adapter,
        },
        "stage_dirs": {},
        **provenance(species_type),
    }
    cpu = mu.resolve_cpu(cpu, jobs, len(samples), slurm)
    score_jobs, score_threads = mu.split_cpu(cpu, jobs, len(samples))
    if slurm:
        jobs, _ = mu.split_cpu(score_threads, None, len(samples))
        threads = score_threads
        log.info(f"Using Slurm: {score_jobs} scoring job(s) at once with {score_threads} core(s) each")
    else:
        jobs, threads = score_jobs, score_threads
        log.info(f"Using {cpu} cores: {jobs} sample(s) at once, {threads} thread(s) each")
    metadata["parameters"].update(cpu=cpu, samples_at_once=score_jobs, threads_per_sample=score_threads, slurm=slurm)

    mu.use_array(score_threads if slurm and not dryrun else None, mem, score_jobs)
    qc_dir = _stage_dir(output, "1.Quality_and_trimming", dryrun=dryrun) if run_qc or run_trim else output
    if run_qc:
        try:
            raw_qc = _load_stage(resume, output, "qc_raw")
            if raw_qc is None:
                _, raw_qc = fastqcRun(sampleDict=samples, outDir=qc_dir, pairedEND=paired, cpu=threads, jobs=jobs, dryrun=dryrun)
                if not dryrun:
                    _save_stage(output, "qc_raw", raw_qc)
            metadata["raw_fastqc"] = str(raw_qc)
        except Exception as exc:
            log.warning(f"Raw FastQC step skipped due to error: {exc}")

    metadata["stage_dirs"]["quality"] = qc_dir
    working_samples = samples
    if run_trim:
        trimmed = _load_stage(resume, output, "trim")
        if trimmed is None:
            log.info("Trimming reads")
            trim_run = {"trimmomatic": trimmomaticRun, "flexbar": flexbarRun, "trim_galore": trim_galoreRun}[trimmer]
            trimmed, _ = trim_run(sampleDict=samples, outDir=qc_dir, paired=paired, adapter=adapter, cpu=threads, jobs=jobs, dryrun=dryrun)
            if not dryrun:
                _save_stage(output, "trim", trimmed)
        working_samples = trimmed
        if run_qc:
            try:
                trim_qc = _load_stage(resume, output, "qc_trim")
                if trim_qc is None:
                    _, trim_qc = fastqcRun(sampleDict=working_samples, outDir=qc_dir, pairedEND=paired, afterTrim=True, cpu=threads, jobs=jobs, dryrun=dryrun)
                    if not dryrun:
                        _save_stage(output, "qc_trim", trim_qc)
                metadata["trim_fastqc"] = str(trim_qc)
            except Exception as exc:
                log.warning(f"Trimmed FastQC step skipped due to error: {exc}")

    align_dir = _stage_dir(output, "2.Alignment", dryrun=dryrun)
    metadata["stage_dirs"]["alignment"] = align_dir
    if dryrun:
        log.info("[dryrun] Collapsing reads")
        collapsed = {k: [v[0], v[1], os.path.join(align_dir, "collapsed_reads", f"{k}_collapsed_reads.txt")] for k, v in working_samples.items()}
    else:
        collapsed = _load_stage(resume, output, "collapse")
        if collapsed is None:
            log.info("Collapsing reads")
            collapse_dir = mu.make_directory(os.path.join(align_dir, "collapsed_reads"))
            collapsed = _write_collapsed_fasta(working_samples, collapse_dir, jobs)
            _save_stage(output, "collapse", collapsed)

    log.info("Aligning reads to the genome")
    aligned_files = None if dryrun else _load_stage(resume, output, "align")
    if aligned_files is None:
        aln = Bowtie_Aligner(ref_genome=genome, outdir=align_dir, slurm=False, dryrun=dryrun)
        aln.build_index(cpu=cpu, mem=mem)
        aligned_files = aln.run_alignment(samplesDict=collapsed, fileType='fasta', outType='SAM', cpu=threads, mem=mem, jobs=jobs)
        if not dryrun:
            _save_stage(output, "align", aligned_files)
    if dryrun:
        log.info("[dryrun] Next: cluster reads, predict novel miRNAs, quantify, DESeq2, targets and report")
        return metadata

    log.info("Reading aligned reads")
    aligned_reads = Read_process().aligned_reads(samples=aligned_files, alignType='SAM', jobs=jobs)
    log.info("Clustering reads")
    clusters = Cluster_reads().get_clusters(aligned_reads, jobs=jobs)

    mirna_dir = _stage_dir(output, "3.miRNA_Discovery")
    metadata["stage_dirs"]["discovery"] = mirna_dir
    log.info("Predicting novel miRNAs")
    novel = Novel_miRNA(
        samples=clusters,
        genome=genome,
        species=species,
        species_type=species_type,
        outdir=mirna_dir,
        filter_criteria='strict',
        require_star=require_star,
        read_cutoff=read_cutoff,
        min_mature_reads=min_mature_reads,
        min_star_reads=min_star_reads,
        max_inconsistent=max_inconsistent,
        min_novel_score=min_novel_score,
        workers=score_threads,
        jobs=score_jobs,
        resume=resume,
        slurm=slurm and not dryrun,
        mem=mem,
    )
    metadata["parameters"].update(
        filter_criteria=novel.filterCriteria,
        read_cutoff=novel.readsCutoff,
        cluster_cutoff=novel.clusterCutoff,
        max_extension=novel.max_extension,
        score_cutoff=novel.score,
        randfold_iterations=novel.randfold_iterations,
        randfold_pvalue_cutoff=novel.randfold_pvalue_cutoff,
        min_mature_reads=novel.min_mature_reads,
        require_star=novel.require_star,
        min_star_reads=novel.min_star_reads,
        random_seed=novel.random_seed,
    )
    known_df, novel_df, known_fasta, novel_fasta = novel.get_novel_miRNA()
    db_summary = os.path.join(mirna_dir, "known_database_summary.tsv")
    if os.path.exists(db_summary):
        metadata["database_summary"] = db_summary
    if isinstance(known_df, pd.DataFrame):
        known_txt = os.path.join(mirna_dir, "known_miRNAs.tsv")
        known_df.to_csv(known_txt, sep="\t", index=False)
        metadata["known_file"] = known_txt
    if isinstance(novel_df, pd.DataFrame):
        novel_txt = os.path.join(mirna_dir, "novel_miRNAs.tsv")
        novel_df.to_csv(novel_txt, sep="\t", index=False)
        metadata["novel_file"] = novel_txt

    log.info("Quantifying miRNAs")
    quant_dir = _stage_dir(output, "4.Quantification")
    metadata["stage_dirs"]["quantification"] = quant_dir
    mirnas = []
    if isinstance(known_df, pd.DataFrame):
        known_unique = known_df.drop_duplicates(subset=["Known_ID"])
        mirnas += list(zip(known_unique["Known_ID"], known_unique["precursor"], known_unique["mature"]))
    if isinstance(novel_df, pd.DataFrame):
        mirnas += list(zip(novel_df["Novel_ID"], novel_df["precursor"], novel_df["mature"]))
    quant_ref = os.path.join(quant_dir, "miRNA_mature.fa")
    with open(quant_ref, "w") as fh:
        for name, _hairpin, mature in mirnas:
            fh.write(f">{name}\n{mature}\n")
    quant = Quantifier(miRNAs=mirnas, samples=collapsed, outdir=quant_dir, slurm=False, weighted=weighted)
    counts = quant.quantify_expression(cpu=cpu, mem=mem)
    counts_file = os.path.join(quant_dir, "Raw_read_counts.tsv")
    counts.to_csv(counts_file, sep="\t", index=False)
    metadata["counts_file"] = counts_file

    log.info("Running differential expression")
    try:
        deg_dir = _stage_dir(output, "5.Differential_Expression")
        metadata["stage_dirs"]["differential"] = deg_dir
        deg = runDESeq2(countDF=counts.copy(), targetFile=targets, combination=combinations, mirna_column='mature_id')
        deg_file = os.path.join(deg_dir, "DESeq2_results.tsv")
        deg.to_csv(deg_file, sep="\t", index=False)
        metadata["deg_file"] = deg_file
    except Exception as exc:
        log.warning(f"Differential expression step skipped due to error: {exc}")

    if mrna_file is not None and os.path.exists(mrna_file):
        log.info("Predicting targets")
        try:
            target_dir = _stage_dir(output, "6.Target_Prediction")
            metadata["stage_dirs"]["targets"] = target_dir
            target_results = run_targets(miRNA_file=quant_ref, mRNA_file=mrna_file, outdir=target_dir, dryrun=False)
            metadata["targets_file"] = target_results.get("parsed_output")
            if run_enrichment_stage and enrichment_organism:
                annotation_dir = _stage_dir(output, "7.Functional_Annotation")
                metadata["stage_dirs"]["annotation"] = annotation_dir
                enrich = run_enrichment(
                    targets_file=target_results.get("parsed_output"),
                    organism=enrichment_organism,
                    outdir=annotation_dir
                )
                metadata["enrichment_file"] = enrich.get("enrichment_results")
        except Exception as exc:
            log.warning(f"Target/enrichment stage skipped due to error: {exc}")

    log.info("Writing report")
    report_dir = _stage_dir(output, "8.Report")
    metadata["stage_dirs"]["report"] = report_dir
    metadata["elapsed_seconds"] = round(time.time() - start)
    report_path, meta_path = generate_report(outdir=report_dir, metadata=metadata)
    metadata["report_file"] = report_path
    metadata["metadata_file"] = meta_path
    return metadata
