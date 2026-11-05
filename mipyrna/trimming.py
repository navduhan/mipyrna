"""
Read Trimming Module

This module performs adapter and quality trimming of small RNA reads. It wraps
Trimmomatic, Flexbar and Trim Galore and detects the 3' adapter from common
small RNA library kits.

Features:
    - Trimming with Trimmomatic, Flexbar or Trim Galore
    - 3' adapter detection from the first reads of a FASTQ file
    - Single-end and paired-end support
    - Local multi-job execution and SLURM cluster submission
    - Dry-run mode that returns the expected output files

Configuration:
    Configured through ``flexbar.ini``, ``trimmomaticSE.ini`` / ``trimmomaticPE.ini`` or
    ``trim_galore.ini`` (packaged defaults or a user directory). Inputs are a samples
    dictionary from ``utility.read_input_file`` and an adapter name or sequence; outputs
    are trimmed FASTQ files per sample in the tool output directory.

Dependencies:
    - Trimmomatic, Flexbar, Trim Galore: Must be on PATH
    - SLURM (optional): Cluster execution

Classes / Functions / Exceptions:
    - ADAPTERS: Known 3' adapters of small RNA kits
    - detect_adapter: Detect the adapter in a FASTQ file
    - resolve_adapter: Resolve an adapter name, sequence or auto-detection
    - flexbarRun: Trim reads with Flexbar
    - trimmomaticRun: Trim reads with Trimmomatic
    - trim_galoreRun: Trim reads with Trim Galore

:Author: Naveen Duhan
"""

import gzip
import importlib.resources as ir
import os
import shutil
import sys

from mipyrna import utility as mu
from mipyrna.logger import MiPyRNALogger

# Intialize the logger

log = MiPyRNALogger(mode='a', log='qt')

# 3' adapters of common small RNA library kits.
ADAPTERS = {
    "truseq_small_rna": "TGGAATTCTCGGGTGCCAAGG",
    "illumina_v1.5": "CTGTAGGCACCATCAAT",
    "nebnext": "AGATCGGAAGAGCACACGTCT",
    "bgi": "AGTCGGAGGCCAAGCGGTC",
}


def detect_adapter(fastq, reads=100000, min_fraction=0.1):
    """Pick the known 3' adapter found in most of the first reads of a FASTQ file.

    Args:
        fastq: FASTQ file, optionally gzipped.
        reads: Number of reads to inspect.
        min_fraction: Minimum fraction of reads that must carry the adapter.

    Returns:
        The adapter sequence.

    Raises:
        ValueError: If no known adapter is found in enough reads.
    """
    opener = gzip.open if str(fastq).endswith(".gz") else open
    seeds = {name: seq[:12] for name, seq in ADAPTERS.items()}
    hits = dict.fromkeys(ADAPTERS, 0)
    total = 0
    with opener(fastq, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            total += 1
            for name, seed in seeds.items():
                if seed in line:
                    hits[name] += 1
            if total >= reads:
                break
    best = max(hits, key=lambda k: hits[k])
    if total == 0 or hits[best] / total < min_fraction:
        raise ValueError(f"No known adapter found in {fastq}; provide one with adapter=<name or sequence>")
    log.info(f"Detected adapter {best} in {hits[best] / total:.0%} of {total} reads")
    return ADAPTERS[best]


def resolve_adapter(adapter, fastq):
    """Return an adapter sequence from a kit name, a sequence, or auto-detection.

    Args:
        adapter: Kit name, adapter sequence, or None/``'auto'`` to detect it.
        fastq: FASTQ file used for auto-detection.

    Returns:
        The adapter sequence.

    Raises:
        ValueError: If auto-detection finds no known adapter.
    """
    if adapter in (None, "auto"):
        return detect_adapter(fastq)
    return ADAPTERS.get(adapter, adapter.upper())


def flexbarRun(sampleDict,  configFile=None, adapter=None, slurm=False, mem=10, cpu=8, task=1, paired=False, outDir=".", dryrun=False, dep='', jobs=1):
    """Trim adapters and low-quality bases from reads with Flexbar.

    See https://github.com/seqan/flexbar.

    Args:
        sampleDict: Samples dictionary from ``utility.read_input_file``.
        configFile: Directory containing the config file; the packaged default is used if absent or invalid.
        adapter: Adapter kit name, sequence, or None/``'auto'`` to detect it from the first sample.
        slurm: Submit jobs with SLURM.
        mem: Memory in GB for SLURM jobs.
        cpu: Total threads to use.
        task: Number of tasks for SLURM jobs.
        paired: True if samples are paired-end.
        outDir: Output directory for results.
        dryrun: Build commands without running them.
        dep: SLURM job id this job depends on.
        jobs: Number of samples to trim at once in local mode.

    Returns:
        Tuple of a dictionary of sample key to ``[name, group, trimmed file(s)]`` and
        the list of SLURM job ids.

    Raises:
        ValueError: If sampleDict is None.
        Exception: Re-raised if a sample entry lacks the expected input paths.

    System Exit:
        Exits if ``flexbar`` is not on PATH.
    """

    if sampleDict is None:
        raise ValueError("sampleDict is required")

    if configFile != None:

        try:
            log.info(f"Using  config file {os.path.join(configFile, 'flexbar.ini')}")

            config = mu.parse_config_file(os.path.join(configFile, 'flexbar.ini'))

        except Exception:

            log.error("Please provide a valid config file. Using default config file flexbar.ini")

            config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/flexbar.ini"))

            log.error("Please provide a valid config file. Using default config file flexbar.ini")

    else:

        log.info("Using default config file flexbar.ini")

        config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/flexbar.ini"))

        
    flexbar_config = config[next(iter(config.keys()))]

    output =  mu.get_basename(str(flexbar_config[0]).split(" ")[1])

    if os.path.exists(outDir) or dryrun:

        output1 = os.path.join(outDir, output)

        output = mu.make_directory(output1, dryrun=dryrun)

    else:
        _parent, base = os.path.split(output)
        output = mu.make_directory(base, dryrun=dryrun)

    const = ['threads', '-n']  # require when change number of CPUs

    if not slurm:

        flexbar_config = mu.replace_cpu(flexbar_config, const, cpu=cpu)

    adapter = resolve_adapter(adapter, next(iter(sampleDict.values()))[2])
    args = ' '.join(flexbar_config[1:]) + f" -as {adapter} -ao 3 -at RIGHT"

    outflex = {}  # initialize a dict to carry outfile names for next function

    job_id = []
    local_cmds = {}

    for key, value in sampleDict.items():  # Iterate thorough total number of samples

        if paired:
            try:
                input1 = value[2]
                input2 = value[3]

                outflex[key] = [value[0], value[1], os.path.join(
                    output, value[0]+"_1.fastq.gz"), os.path.join(output, value[0]+"_2.fastq.gz")]

            except Exception:

                log.error(
                    "please provide a paired END sample file or input Path is wrong")
                raise

            inputPair = mu.get_basename(input1)+" and " + mu.get_basename(input2)

            out = os.path.join(output, value[0])

        else:
            try:
                value[2]

            except Exception:

                log.error("please provide a valid input data path")
                raise

            outflex[key] = [value[0], value[1],
                            os.path.join(output, value[0]+".fastq.gz")]

            inputPair = mu.get_basename(value[2])

            out = os.path.join(output, value[0])

        execPATH = shutil.which("flexbar")

        if execPATH is None:  # get absolute path of flexbar
            
            log.error("flexbar command not found in path")

            sys.exit()

        else:
            
            if paired:

                flexbarCmd = f"{execPATH} -r {value[2]} -p {value[3]} -t {out} {args}"

            else:

                flexbarCmd = f"{execPATH} -r {value[2]} -t {out} {args}"
            if dryrun:
                log.info(f"[dryrun] {flexbarCmd}")
            else:

                if slurm == True:  # check if slurm job scheduling is enabled or not

                    try:
                        job = mu.clusterRun(job_name='flexbar', sout=os.path.join(output, "flexbar.out"), serror=os.path.join(
                            output, "flexbar.err"), command=flexbarCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)

                        job_id.append(job)

                        log.info(
                            f"Job submitted on slurm successfully for {inputPair} with {job}")

                    except Exception:

                        log.error("Slurm job sumission failed")

                else:

                    local_cmds[key] = flexbarCmd
    if local_cmds:
        mu.local_run(program="flexbar", command=local_cmds, outdir=output, jobs=jobs)
    return outflex, job_id


def trimmomaticRun(sampleDict=None, configFile=None, adapter=None, slurm=False, mem=10, cpu=8, task=1, paired=False, outDir=".", dryrun=False, dep='', jobs=1):
    """Trim adapters and low-quality bases from reads with Trimmomatic.

    The adapter is written to ``adapter.fa`` in the output directory and used for ILLUMINACLIP.

    Args:
        sampleDict: Samples dictionary from ``utility.read_input_file``.
        configFile: Directory containing the config file; the packaged default is used if absent or invalid.
        adapter: Adapter kit name, sequence, or None/``'auto'`` to detect it from the first sample.
        slurm: Submit jobs with SLURM.
        mem: Memory in GB for SLURM jobs.
        cpu: Total threads to use.
        task: Number of tasks for SLURM jobs.
        paired: True if samples are paired-end.
        outDir: Output directory for results.
        dryrun: Build commands without running them.
        dep: SLURM job id this job depends on.
        jobs: Number of samples to trim at once in local mode.

    Returns:
        Tuple of a dictionary of sample key to ``[name, group, trimmed file(s)]`` and
        the list of SLURM job ids.

    Raises:
        ValueError: If sampleDict is None.
        Exception: Re-raised if a sample entry lacks the expected input paths.

    System Exit:
        Exits if ``trimmomatic`` is not on PATH.
    """

    if sampleDict is None:
        raise ValueError("sampleDict is required")

    if configFile != None:

        try:

            if paired:
                
                log.info(f"Using  config file {os.path.join(configFile, 'trimmomaticPE.ini')}")

                config = mu.parse_config_file(os.path.join(configFile, 'trimmomaticPE.ini'))

            else:

                log.info(f"Using  config file {os.path.join(configFile, 'trimmomaticSE.ini')}")

                config = mu.parse_config_file(os.path.join(configFile, 'trimmomaticSE.ini'))

                
        except Exception:

            if paired:

                log.info("Please provide a valid config file. Using default config file trimmomaticPE.ini")
            
                config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/trimmomaticPE.ini"))
                

            else:

                log.info("Please provide a valid config file. Using default config file trimmomaticSE.ini")
                
                config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/trimmomaticSE.ini"))
                

    else:
        if paired:

            log.info("Using default config file trimmomaticPE.ini")

            config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/trimmomaticPE.ini"))
            

        else:

            log.info("Using default config file trimmomaticSE.ini")
            
            config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/trimmomaticSE.ini"))
            
           

    trimmomatic_config = config[next(iter(config.keys()))]

    output = mu.get_basename(str(trimmomatic_config[0]))

    if os.path.exists(outDir) or dryrun:

        output1 = os.path.join(outDir, output)

        output = mu.make_directory(output1, dryrun=dryrun)

    else:

        output = mu.make_directory(output, dryrun=dryrun)

    const = ['threads', '-n', 'cores']  # require when change number of CPUs

    if not slurm:

        trimmomatic_config = mu.replace_cpu(trimmomatic_config, const, cpu=cpu)

    adapter = resolve_adapter(adapter, next(iter(sampleDict.values()))[2])
    adapter_file = os.path.join(output, "adapter.fa")
    if not dryrun:
        with open(adapter_file, "w") as fh:
            fh.write(f">adapter\n{adapter}\n")
    options = [a for a in trimmomatic_config[1:] if a.startswith('-')]
    steps = [a for a in trimmomatic_config[1:] if not a.startswith('-')]
    clip = f"ILLUMINACLIP:{os.path.abspath(adapter_file)}:1:30:5"
    args = ' '.join(options + [clip] + steps)

    outtrimmomatic = {}  # initialize a dict to carry outfile names for next function

    job_id = []
    local_cmds = {}

    for key, value in sampleDict.items():  # Iterate thorough total number of samples

        if paired:
            try:
                inputFile1 = value[2]
                inputFile2 = value[3]

                outtrimmomatic[key] = [value[0], value[1], os.path.join(
                    output, value[0]+"_1.fastq.gz"), os.path.join(output, value[0]+"_2.fastq.gz")]

            except Exception:

                log.error(
                    "Please provide a paired END sample file or input Path is wrong")
                raise

            inputPair = mu.get_basename(
                inputFile1)+" and " + mu.get_basename(inputFile2)

        else:
            try:
                inputFile = value[2]

            except Exception:

                log.error("please provide a valid input data path")
                raise

            outtrimmomatic[key] = [value[0], value[1],
                                   os.path.join(output, value[0]+".fastq.gz")]

            inputPair = mu.get_basename(inputFile)

    

        execPATH = shutil.which('trimmomatic') # get absolute path of trimmomatic

        if execPATH is None:

            log.error("trimmomatic command not found in path")

            sys.exit()
        else:

            if paired:

                trimmomaticCmd = f"{execPATH} PE {value[2]} {value[3]} {outtrimmomatic[key][2]} {outtrimmomatic[key][3]} {args}"

            else:

                trimmomaticCmd = f"{execPATH} SE {value[2]} {outtrimmomatic[key][2]} {args}"

            if dryrun:
                log.info(f"[dryrun] {trimmomaticCmd}")
            else:
                if slurm == True:  # check if slurm job scheduling is enabled or not

                    try:
                        job = mu.clusterRun(job_name='trimmomatic', sout=os.path.join(output, "trimmomatic.out"), serror=os.path.join(
                            output, "trimmomatic.err"), command=trimmomaticCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)
                        job_id.append(job)
                        log.info(f"Job successfully submited for {inputPair} with {job}")

                    except Exception:

                        log.error("Slurm job sumission failed")

                else:

                    local_cmds[key] = trimmomaticCmd
    if local_cmds:
        mu.local_run(program="trimmomatic", command=local_cmds, outdir=output, jobs=jobs)
    try:
        shutil.move("LEADING:3", output+"/LEADING:3")
        shutil.move("TRAILING:3", output+"/TRAILING:3")
    except Exception:
        pass
    return outtrimmomatic, job_id


def trim_galoreRun(sampleDict=None, configFile=None, adapter=None, slurm=False, mem=10, cpu=8, task=1, paired=False, outDir=".", dryrun = False, dep='', jobs=1):
    """Trim adapters and low-quality bases from reads with Trim Galore.

    Args:
        sampleDict: Samples dictionary from ``utility.read_input_file``.
        configFile: Directory containing the config file; the packaged default is used if absent or invalid.
        adapter: Adapter kit name, sequence, or None/``'auto'`` to detect it from the first sample.
        slurm: Submit jobs with SLURM.
        mem: Memory in GB for SLURM jobs.
        cpu: Total threads to use.
        task: Number of tasks for SLURM jobs.
        paired: True if samples are paired-end.
        outDir: Output directory for results.
        dryrun: Build commands without running them.
        dep: SLURM job id this job depends on.
        jobs: Number of samples to trim at once in local mode.

    Returns:
        Tuple of a dictionary of sample key to ``[name, group, trimmed file(s)]`` and
        the list of SLURM job ids.

    Raises:
        ValueError: If sampleDict is None.
        Exception: Re-raised if a sample entry lacks the expected input paths.

    System Exit:
        Exits if ``trim_galore`` is not on PATH.
    """

    if sampleDict is None:
        raise ValueError("sampleDict is required")

    if configFile != None:

        try:

            log.info(f"Using config file {os.path.join(configFile, 'trim_galore.ini')}")

            config = mu.parse_config_file(os.path.join(configFile, 'trim_galore.ini'))
            

        except Exception:

            log.error("Please provide a valid config file. Using default config file trim_galore.ini")

            config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/trim_galore.ini"))

    else:
        log.info("Using default config file trim_galore.ini")

        config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/trim_galore.ini"))
        

    trim_galore_config = config[next(iter(config.keys()))]

    output =  mu.get_basename(str(trim_galore_config[0].split(" ")[1]))

    if os.path.exists(outDir) or dryrun:

        output1 = os.path.join(outDir, output)

        output = mu.make_directory(output1, dryrun=dryrun)

    else:

        output = mu.make_directory(output, dryrun=dryrun)

    # require when change number of CPUs
    const = ['threads', '-n', 'cores', '-j']

    if not slurm:

        trim_galore_config = mu.replace_cpu(trim_galore_config, const, cpu=cpu)

    adapter = resolve_adapter(adapter, next(iter(sampleDict.values()))[2])
    args = f"{' '.join(trim_galore_config[1:])} -a {adapter}"

    outtrim_galore = {}  # initialize a dict to carry outfile names for next function

    job_id = []
    local_cmds = {}

    for key, value in sampleDict.items():  # Iterate thorough total number of samples

        if paired:

            try:
                input1 = value[2]
                input2 = value[3]

                outtrim_galore[key] = [value[0], value[1], os.path.join(
                    output, value[0]+"_val_1.fq.gz"), os.path.join(output, value[0]+"_val_2.fq.gz")]
                
            except Exception:

                log.error(
                    "please provide a paired END sample file or input Path is wrong")
                raise

            out = value[0]
            inputPair = mu.get_basename(
                input1)+" and " + mu.get_basename(input2)

        else:
            try:
                value[2]

            except Exception:

                log.error("please provide a valid input data path")
                raise

            outtrim_galore[key] = [value[0], value[1],
                                   os.path.join(output, value[0]+"_trimmed.fq.gz")]

            out = value[0]

            inputPair = mu.get_basename(value[2])

        

        execPATH = shutil.which('trim_galore')  # get absolute path of trim_galore
        gzip = ''
        if not value[2].endswith(".gz"):
            gzip = '--gzip'
            

        if execPATH is None:

            log.error("trim_galore command not found in path")
            sys.exit()

        else:
            if paired:

                trim_galoreCmd = f"{execPATH} {args} {value[2]} {value[3]} --paired --basename {out} -o {output} {gzip}"

            else:

                trim_galoreCmd = f"{execPATH} {args} {value[2]} --basename {out} -o {output} {gzip}"

                # log.info(f"{trim_galoreCmd}")
                
            if dryrun:
                log.info(f"[dryrun] {trim_galoreCmd}")

            else:

                if slurm == True:  # check if slurm job scheduling is enabled or not

                    try:
                        job = mu.clusterRun(job_name='trim_galore', sout=os.path.join(output, "trim_galore.out"), serror=os.path.join(
                            output, "trim_galore.err"), command=trim_galoreCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)
                        job_id.append(job)
                        log.info(f"Job successfully submited for {inputPair} with {job}")

                    except Exception:

                        log.error("Slurm job sumission failed")

                else:

                    local_cmds[key] = trim_galoreCmd
    if local_cmds:
        mu.local_run(program="trim_galore", command=local_cmds, outdir=output, jobs=jobs)
    return outtrim_galore, job_id
