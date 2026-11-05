#!/usr/bin/env python3

"""End-to-end workflow orchestration for miPyRNA."""

import os

import pandas as pd

from mipyrna import utility as mu
from mipyrna.aligner import Bowtie_Aligner
from mipyrna.annotation import run_enrichment
from mipyrna.check import provenance
from mipyrna.clustering import Cluster_reads
from mipyrna.differential import runDESeq2
from mipyrna.logger import MiPyRNALogger
from mipyrna.novel import Novel_miRNA
from mipyrna.quality import fastqcRun
from mipyrna.quantifier import Quantifier
from mipyrna.reads import Read_process
from mipyrna.report import generate_report
from mipyrna.targets import run_targets
from mipyrna.trimming import trimmomaticRun

log = MiPyRNALogger(mode='a', log='workflow')


def _write_collapsed_fasta(samples, outdir):
    """Collapse the reads of every sample into FASTA files."""
    collapsed = Read_process().collapse_reads(samples=samples, outdir=outdir)
    return collapsed


def _combine_fastas(files, output):
    """Concatenate the existing FASTA files into one."""
    with open(output, "w") as out:
        for fpath in files:
            if fpath is None or not os.path.exists(fpath):
                continue
            with open(fpath) as fh:
                out.write(fh.read().rstrip() + "\n")
    return output


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
    cpu=8,
    mem=20
):
    """Run the full pipeline from raw reads to differential expression, targets and report."""
    output = mu.make_directory(outdir)
    parsed = mu.read_input_file(infile=input_file, inpath=input_path, paired=paired)
    samples = parsed["samples"]
    combinations = parsed["combinations"]
    targets = parsed["targets"]
    metadata = {"outdir": output, **provenance(species_type)}

    if run_qc:
        try:
            _, raw_qc = fastqcRun(sampleDict=samples, outDir=output, pairedEND=paired)
            metadata["raw_fastqc"] = str(raw_qc)
        except Exception as exc:
            log.warning(f"Raw FastQC step skipped due to error: {exc}")

    working_samples = samples
    if run_trim:
        log.info("Trimming reads")
        trimmed, _ = trimmomaticRun(sampleDict=samples, outDir=output, paired=paired)
        working_samples = trimmed
        if run_qc:
            try:
                _, trim_qc = fastqcRun(sampleDict=working_samples, outDir=output, pairedEND=paired, afterTrim=True)
                metadata["trim_fastqc"] = str(trim_qc)
            except Exception as exc:
                log.warning(f"Trimmed FastQC step skipped due to error: {exc}")

    log.info("Collapsing reads")
    collapse_dir = mu.make_directory(os.path.join(output, "collapsed_reads"))
    collapsed = _write_collapsed_fasta(working_samples, collapse_dir)

    log.info("Aligning reads to the genome")
    aln = Bowtie_Aligner(ref_genome=genome, outdir=output, slurm=False)
    aln.build_index(cpu=cpu, mem=mem)
    aligned_files = aln.run_alignment(samplesDict=collapsed, fileType='fasta', outType='SAM', cpu=cpu, mem=mem)
    aligned_reads = Read_process().aligned_reads(samples=aligned_files, alignType='SAM')
    clusters = Cluster_reads().get_clusters(aligned_reads)

    log.info("Predicting novel miRNAs")
    novel = Novel_miRNA(
        samples=clusters,
        genome=genome,
        species=species,
        species_type=species_type,
        outdir=output,
        filter_criteria='strict'
    )
    known_df, novel_df, known_fasta, novel_fasta = novel.get_novel_miRNA()
    if isinstance(known_df, pd.DataFrame):
        known_txt = os.path.join(output, "known_miRNAs.tsv")
        known_df.to_csv(known_txt, sep="\t", index=False)
        metadata["known_file"] = known_txt
    if isinstance(novel_df, pd.DataFrame):
        novel_txt = os.path.join(output, "novel_miRNAs.tsv")
        novel_df.to_csv(novel_txt, sep="\t", index=False)
        metadata["novel_file"] = novel_txt

    log.info("Quantifying miRNAs")
    quant_ref = _combine_fastas([known_fasta, novel_fasta], os.path.join(output, "combined_miRNA.fa"))
    quant = Quantifier(sRNA_file=quant_ref, samples=collapsed, outdir=output, slurm=False)
    counts = quant.quantify_expression(cpu=cpu, mem=mem)
    counts_file = os.path.join(output, "Raw_read_counts.tsv")
    counts.to_csv(counts_file, sep="\t", index=False)
    metadata["counts_file"] = counts_file

    log.info("Running differential expression")
    try:
        deg = runDESeq2(countDF=counts.copy(), targetFile=targets, combination=combinations, mirna_column='mature_id')
        deg_file = os.path.join(output, "DESeq2_results.tsv")
        deg.to_csv(deg_file, sep="\t", index=False)
        metadata["deg_file"] = deg_file
    except Exception as exc:
        log.warning(f"Differential expression step skipped due to error: {exc}")

    if mrna_file is not None and os.path.exists(mrna_file):
        log.info("Predicting targets")
        try:
            target_results = run_targets(miRNA_file=quant_ref, mRNA_file=mrna_file, outdir=output, dryrun=False)
            metadata["targets_file"] = target_results.get("parsed_output")
            if run_enrichment_stage and enrichment_organism:
                enrich = run_enrichment(
                    targets_file=target_results.get("parsed_output"),
                    organism=enrichment_organism,
                    outdir=output
                )
                metadata["enrichment_file"] = enrich.get("enrichment_results")
        except Exception as exc:
            log.warning(f"Target/enrichment stage skipped due to error: {exc}")

    log.info("Writing report")
    report_path, meta_path = generate_report(outdir=output, metadata=metadata)
    metadata["report_file"] = report_path
    metadata["metadata_file"] = meta_path
    return metadata
