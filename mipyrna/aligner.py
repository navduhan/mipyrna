"""
Bowtie Aligner Module

This module provides Bowtie index building and read alignment wrappers for
small RNA data. It builds Bowtie indexes from a reference genome and aligns
reads to a reference or to the miRBase hairpin and mature sets.

Features:
    - Bowtie index building from plain or gzipped FASTA references
    - Alignment of FASTA/FASTQ reads, single-end or paired-end, per sample or pooled
    - Local multi-job execution and SLURM cluster submission
    - Dry-run mode that returns commands or output names without executing
    - Thread count adjustment from the configuration file

Configuration:
    Configured through ``bowtie1.ini`` (sections ``[index]`` and ``[alignment]``), either
    the packaged default or a user-supplied file. Inputs are a reference FASTA and a
    samples dictionary from ``utility.read_input_file``; outputs are the index directory
    and SAM files in ``bowtie_results``.

Dependencies:
    - Bowtie 1 (``bowtie`` and ``bowtie-build``): Must be on PATH
    - waiting: Waiting for SLURM jobs to finish
    - SLURM (optional): Cluster execution

Classes / Functions / Exceptions:
    - Bowtie_Aligner: Builds Bowtie indexes and aligns reads.

:Author: Naveen Duhan
"""

import gzip
import importlib.resources as ir
import os
import shutil
import subprocess
import sys

from waiting import wait

from mipyrna import utility as mu
from mipyrna.logger import MiPyRNALogger

# Intialize the logger

log = MiPyRNALogger(mode="a", log="aligner")


class Bowtie_Aligner:
    """Build Bowtie indexes and align reads with Bowtie.

    Attributes:
        genome: Reference genome FASTA file.
        outdir: Output directory.
        slurm: Whether jobs are submitted through SLURM.
        dryrun: Whether commands are only generated, not executed.
        config: Parsed Bowtie configuration (sections ``index`` and ``alignment``).
    """

    def __init__(
        self, ref_genome=None, configFile=None, outdir=".", slurm=False, dryrun=False
    ):
        """Set the reference, output directory and execution mode (local or SLURM).

        Args:
            ref_genome: Reference genome FASTA file.
            configFile: User-defined Bowtie config file. Defaults to the packaged ``bowtie1.ini``.
            outdir: Output directory.
            slurm: Submit jobs with the SLURM scheduler.
            dryrun: Generate commands without executing them.

        Raises:
            ValueError: If ref_genome is not given.
        """
        if ref_genome is None:
            raise ValueError("ref_genome is required")

        self.genome = ref_genome
        self.outdir = outdir
        self.slurm = slurm
        self.dryrun = dryrun

        if configFile != None:
            try:
                self.config = mu.parse_config_file(configFile)

                log.info(f"Using  config file {configFile}")

            except Exception:
                log.error("Please provide a valid config file")
        else:
            self.config = mu.parse_config_file(
                ir.files("mipyrna").joinpath("param/bowtie1.ini")
            )

            log.info("Using default config file bowtie.ini")

    def build_index(self, mem=20, tasks=1, cpu=8, dep=""):
        """Build a Bowtie index for the reference genome.

        The genome (decompressed if gzipped) is copied to the index directory and
        ``bowtie-build`` is run locally or through SLURM.

        Args:
            mem: Memory in GB for the SLURM job.
            tasks: Number of tasks for the SLURM job.
            cpu: Threads to use.
            dep: SLURM job id this job depends on.

        Returns:
            The ``bowtie-build`` command string in dry-run mode, otherwise None.

        System Exit:
            Exits with status 1 if the genome extension is not a FASTA extension
            or ``bowtie-build`` is not on PATH.
        """

        const = ["--threads", "-p"]

        if not self.slurm:
            config = mu.replace_cpu(self.config["index"], const, cpu=cpu)

        else:
            config = self.config["index"]

        directory = str(config[0]).split(" ")[0]

        indexName = str(config[1]).split(" ")[0]

        if self.dryrun:
            if os.path.exists(self.outdir):
                output = os.path.join(self.outdir, directory)

            else:
                output = os.path.join(self.outdir, directory)

        else:
            if os.path.exists(self.outdir):
                output1 = os.path.join(self.outdir, directory)

                output = mu.make_directory(output1)

            else:
                output = mu.make_directory(directory)

        genome_extensions = (".fa", ".fasta", ".fna", ".fa.gz", "fna.gz", ".fasta.gz")

        if self.genome.endswith(genome_extensions):
            if self.dryrun:
                pass

            elif self.genome.endswith(".gz"):
                target = os.path.join(output, os.path.basename(self.genome)[:-3])

                with (
                    gzip.GzipFile(self.genome, "rb") as fin,
                    open(target, "wb") as fout,
                ):
                    shutil.copyfileobj(fin, fout)

            else:
                shutil.copy(self.genome, output)

            if not self.dryrun:
                log.info(
                    f"{mu.get_basename(self.genome)} copied successfully in {directory}"
                )

        else:
            log.error(
                f"Please provide a valid genome fasta file with these {genome_extensions} extensions"
            )

            sys.exit(1)

        if indexName != "NA":
            basename = os.path.join(output, indexName)

        else:
            basename = os.path.join(output, mu.get_basename(self.genome))

        GenomeFasta = os.path.join(
            output, mu.get_basename(self.genome).removesuffix(".gz")
        )

        arg = " ".join(config[2:])

        execPATH = shutil.which("bowtie-build")

        if execPATH is None:
            log.error("bowtie-build not found in path")

            sys.exit(1)
        else:
            bowtie_command = f"{execPATH} -f {GenomeFasta} {arg} {basename}"

        if self.dryrun:
            log.info(f"[dryrun] {bowtie_command}")

            return bowtie_command
        else:
            if self.slurm:
                try:
                    job_id = mu.clusterRun(
                        job_name="bowtie-build",
                        sout=os.path.join(output, "bowtie_build.out"),
                        serror=os.path.join(output, "bowtie_build.err"),
                        command=bowtie_command,
                        mem=mem,
                        cpu=cpu,
                        tasks=tasks,
                        dep=dep,
                    )

                    log.info(
                        f"Job successfully submited for {GenomeFasta} with {job_id} for indexing"
                    )

                    wait(
                        lambda: mu.check_status(job_id),
                        waiting_for="alignment to finish",
                    )

                    log.info("Genome indexing job finished")
                except Exception:
                    log.error("Slurm job sumission failed")

            else:
                try:
                    logdir = os.path.join(output, "logs")
                    os.makedirs(logdir, exist_ok=True)

                    with open(os.path.join(logdir, "bowtie_build.out"), "w+") as fout:
                        with open(
                            os.path.join(logdir, "bowtie_build.err"), "w+"
                        ) as ferr:
                            job_id = subprocess.call(
                                bowtie_command, shell=True, stdout=fout, stderr=ferr
                            )
                            if job_id != 0:
                                log.warning(
                                    f"command exited with code {job_id}: {bowtie_command}"
                                )

                            log.info(
                                f"Job successfully submited for {GenomeFasta} for indexing"
                            )

                except Exception:
                    log.error("Job sumission failed")

        return

    def check_index(self):
        """Check whether a Bowtie index exists in ``<outdir>/bowtie_index``.

        Returns:
            Generator of paths of the ``.1.ebwt`` index files if the directory exists,
            otherwise False.
        """

        output = os.path.join(self.outdir, "bowtie_index")

        files = [r".1.ebwt", r".2.ebwt", r".3.ebwt", r".4.ebwt"]

        if os.path.exists(output) and os.path.isdir(output):
            for f in files:
                return (
                    os.path.join(output, _) for _ in os.listdir(output) if _.endswith(f)
                )

        return False

    def run_alignment(
        self,
        samplesDict=None,
        samples=None,
        fileType="fasta",
        outType="BWT",
        pairedEND=False,
        mem=20,
        cpu=8,
        tasks=1,
        dep="",
        jobs=1,
    ):
        """Align reads against the indexed reference genome.

        With ``samplesDict`` each sample is aligned separately; otherwise the single
        read file in ``samples`` is aligned as one pooled set.

        Args:
            samplesDict: Samples dictionary from ``utility.read_input_file``.
            samples: Sequence whose third element is the reads file, used when samplesDict is not given.
            fileType: Read format, ``'fasta'`` or ``'fastq'``.
            outType: Output type; ``'SAM'`` adds the ``--sam`` option.
            pairedEND: True if samples are paired-end.
            mem: Memory in GB for SLURM jobs.
            cpu: Threads per alignment job.
            tasks: Number of tasks for SLURM jobs.
            dep: SLURM job id this job depends on.
            jobs: Number of samples to align at once in local mode.

        Returns:
            Dictionary mapping sample key to ``[name, group, sam_file]`` when samplesDict
            is given, otherwise the path of the pooled SAM file.

        Raises:
            ValueError: If fileType is not ``'fastq'`` or ``'fasta'``, or neither
                samplesDict nor samples is given.

        System Exit:
            Exits if ``bowtie`` is not on PATH.
        """

        consta = ["--threads", "-p"]

        if self.slurm:
            config = self.config["alignment"]

        else:
            config = mu.replace_cpu(self.config["alignment"], consta, cpu=cpu)

        directory = str(config[0])

        reference = str(config[1])

        genomeIndex = os.path.join(self.outdir, directory, reference)

        if self.dryrun:
            if os.path.exists(self.outdir):
                output = os.path.join(self.outdir, "bowtie_results")

            else:
                output = os.path.join(self.outdir, "bowtie_results")

        else:
            if os.path.exists(self.outdir):
                output1 = os.path.join(self.outdir, "bowtie_results")

                output = mu.make_directory(output1)

            else:
                output = mu.make_directory("bowtie_results")

        arg = " ".join(config[2:])

        if outType == "SAM":
            arg += " --sam"

        outbowtie = {}
        local_cmds = {}

        job_id = []

        if samplesDict:
            for key, sample in samplesDict.items():
                if pairedEND:
                    outPrefix = os.path.join(output, sample[0])

                    outBAM = outPrefix + "_bowtie.sam"

                    outbowtie[key] = [sample[0], sample[1], outBAM]

                else:
                    outPrefix = os.path.join(output, sample[0])

                    outBAM = outPrefix + "_bowtie.sam"

                    outbowtie[key] = [sample[0], sample[1], outBAM]

                execPATH = shutil.which("bowtie")

                if execPATH is None:
                    log.error("bowtie aligner not found in path")
                    sys.exit(1)

                else:
                    if pairedEND:
                        if fileType == "fastq":
                            bowtie_command = f"{execPATH} -q -x {genomeIndex} {arg} -1 {sample[2]} -2 {sample[3]} > {outBAM}"

                        elif fileType == "fasta":
                            bowtie_command = f"{execPATH} -f -x {genomeIndex} {arg} -1 {sample[2]} -2 {sample[3]} > {outBAM}"

                        else:
                            raise ValueError(
                                f"fileType must be 'fastq' or 'fasta', got {fileType!r}"
                            )

                    else:
                        if fileType == "fastq":
                            bowtie_command = f"{execPATH} -q -x {genomeIndex} {arg} {sample[2]} > {outBAM}"

                        elif fileType == "fasta":
                            bowtie_command = f"{execPATH} -f -x {genomeIndex} {arg} {sample[2]} > {outBAM}"

                        else:
                            raise ValueError(
                                f"fileType must be 'fastq' or 'fasta', got {fileType!r}"
                            )

                    if self.dryrun:
                        log.info(f"[dryrun] {bowtie_command}")
                    else:
                        if self.slurm:
                            try:
                                job = mu.clusterRun(
                                    job_name="bowtie_align",
                                    sout=os.path.join(output, "bowtie.out"),
                                    serror=os.path.join(output, "bowtie.err"),
                                    command=bowtie_command,
                                    mem=mem,
                                    cpu=cpu,
                                    tasks=tasks,
                                    dep=dep,
                                )
                                job_id.append(job)
                                log.info(
                                    f"Job successfully submited for {outPrefix} with {job} for alignment"
                                )

                            except Exception:
                                log.error("Slurm job sumission failed")

                        else:
                            local_cmds[key] = bowtie_command

            if self.dryrun:
                return outbowtie

            if local_cmds:
                mu.local_run(
                    program="bowtie", command=local_cmds, outdir=output, jobs=jobs
                )

            if self.slurm:
                for job in job_id:
                    wait(
                        lambda job=job: mu.check_status(job),
                        waiting_for="alignment to finish",
                    )
                    log.info(f"Alignment completed for job {job}")

            return outbowtie

        else:
            if samples is None:
                raise ValueError("samples or samplesDict is required")

            outBAM = os.path.join(output, "all_reads_bowtie.sam")

            execPATH = shutil.which("bowtie")

            if execPATH is None:
                log.error("bowtie aligner not found in path")
                sys.exit(1)

            else:
                if fileType == "fastq":
                    bowtie_command = (
                        f"{execPATH} -q -x {genomeIndex} {arg} {samples[2]} > {outBAM}"
                    )

                elif fileType == "fasta":
                    bowtie_command = f"{execPATH} -f -x {genomeIndex} -p 40 -n 0 -e 80 -l 18 -a -m 5 --best --strata {samples[2]} > {outBAM}"

                else:
                    raise ValueError(
                        f"fileType must be 'fastq' or 'fasta', got {fileType!r}"
                    )

                if self.slurm:
                    try:
                        job = mu.clusterRun(
                            job_name="bowtie_align",
                            sout=os.path.join(output, "bowtie.out"),
                            serror=os.path.join(output, "bowtie.err"),
                            command=bowtie_command,
                            mem=mem,
                            cpu=cpu,
                            tasks=tasks,
                            dep=dep,
                        )
                        job_id.append(job)
                        log.info(
                            f"Job successfully submited for {outBAM} with {job} for alignment"
                        )

                    except Exception:
                        log.error("Slurm job sumission failed")

                else:
                    try:
                        logdir = os.path.join(output, "logs")
                        os.makedirs(logdir, exist_ok=True)

                        with open(os.path.join(logdir, "bowtie.out"), "w+") as fout:
                            with open(os.path.join(logdir, "bowtie.err"), "w+") as ferr:
                                job = subprocess.call(
                                    bowtie_command, shell=True, stdout=fout, stderr=ferr
                                )
                                if job != 0:
                                    log.warning(
                                        f"command exited with code {job}: {bowtie_command}"
                                    )
                                job_id.append(job)
                                log.info(
                                    f"Job successfully completed for {outBAM} for alignment"
                                )

                    except Exception:
                        log.exception("Job sumission failed")

            return outBAM
