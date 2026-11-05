"""Read Quality Module

This module runs FastQC on raw or trimmed reads, locally or through SLURM.

Features:
    - Single-end and paired-end quality checks
    - Separate output directories for raw and trimmed reads
    - Local parallel execution or SLURM job submission
    - Dry-run mode that only returns the expected output paths

Configuration:
    FastQC parameters are read from ``param/fastqc.ini`` unless a config file
    is given. Results go to ``fastqc_results`` or ``trim_fastqc_results``.

Dependencies:
    - FastQC: Must be available on PATH

Classes / Functions / Exceptions:
    - fastqcRun: Run FastQC on a set of samples.

:Author: Naveen Duhan
"""

import importlib.resources as ir
import os
import shutil
import sys

from mipyrna import utility as mu
from mipyrna.logger import MiPyRNALogger

# Intialize the logger

log = MiPyRNALogger(mode='a', log='quality_check')

def fastqcRun(sampleDict=None, configFile=None,slurm=False, mem=10, cpu=8, task=1, pairedEND =False, afterTrim=False, dryrun = False, outDir=".", dep='', jobs=1):

    """Run FastQC on the samples (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

    Args:
        sampleDict (dict): Sample name to ``[group, replicate, read1, (read2)]``.
        configFile (str, optional): FastQC parameter file. Defaults to the bundled ``fastqc.ini``.
        slurm (bool): Submit jobs through SLURM. Defaults to False.
        mem (int): Memory per job when using SLURM. Defaults to 10.
        cpu (int): Threads per job. Defaults to 8.
        task (int): Tasks per CPU when using SLURM. Defaults to 1.
        pairedEND (bool): True if the samples are paired-end. Defaults to False.
        afterTrim (bool): True when checking trimmed reads. Defaults to False.
        dryrun (bool): Only build commands and return the output paths. Defaults to False.
        outDir (str): Output directory. Defaults to the current directory.
        dep (str): SLURM job dependency. Defaults to ``''``.
        jobs (int): Number of samples to run in parallel locally. Defaults to 1.

    Returns:
        tuple: ``(job_ids, fastqc_output)`` where ``fastqc_output`` maps sample
        names to ``[group, replicate, zip paths...]``.

    Raises:
        ValueError: If ``sampleDict`` is None.
        Exception: If read files are missing from a paired-end or single-end sample.
    """

    if sampleDict is None:
        raise ValueError("sampleDict is required")

    if configFile != None:

        config = mu.parse_config_file(configFile)

    else:
        config = mu.parse_config_file(ir.files("mipyrna").joinpath("param/fastqc.ini"))

        log.info("Using default config file fastqc.ini")

    fastqc_config = config[next(iter(config.keys()))]

    out = 'trim_fastqc_results' if afterTrim else "fastqc_results"

    if os.path.exists(outDir) or dryrun:

        output1 = os.path.join(outDir,out)

        output = mu.make_directory(output1, dryrun=dryrun)

    else:
        
        output = mu.make_directory(out, dryrun=dryrun)

    const = ['threads', '-t']  # require when change number of CPUs

    if not slurm:

        fastqc_config = mu.replace_cpu(fastqc_config, const, cpu=cpu)

    args = ' '.join(fastqc_config[0:])

    job_id = []
    local_cmds = {}

    fastqcOut = {}

    for key, value in sampleDict.items():  # Iterate thorough total number of samloges

        if pairedEND:
            try:
                input1 = value[2]
                input2 = value[3]

            except Exception:

                log.error(
                    "log please provide a paired END samloge file or input Path is wrong")
                raise

            inputPair = mu.get_basename(input1)+" and " + mu.get_basename(input2)

            
            
            if mu.get_file_extension(value[2]) == "gz":

                r1 = f"{output}/{mu.get_basename(mu.get_basename(value[2]))}.zip"

                r2 = f"{output}/{mu.get_basename(mu.get_basename(value[3]))}.zip"
            
                fastqcOut[key]=[value[0], value[1], r1, r2]

            else: 

                r1 = f"{output}/{mu.get_basename(value[2])}.zip"

                r2 = f"{output}/{mu.get_basename(value[3])}.zip"

                fastqcOut[key]=[value[0], value[1], r1, r2]
        else:
            try:
                value[2]

            except Exception:

                log.error("logease provide a valid input data path")
                raise
            
            inputPair = mu.get_basename(value[2])

            if mu.get_file_extension(inputPair) == "gz":
                 
                r1 = f"{output}/{mu.get_basename(inputPair)}.zip"

                fastqcOut[key]=[value[0], value[1], r1 ]

            else:

                r1 = f"{output}/{inputPair}.zip"

                fastqcOut[key]=[value[0], value[1], r1 ]


        execPATH = shutil.which('fastqc')
              # get absolute path of flexbar

        if execPATH is None:

            log.error("fastqc command not found in path")
            sys.exit()
        else:

            if pairedEND:

                fastqcCmd = f"{execPATH} -o {output} {args} {value[2]} {value[3]}"

            else:

                fastqcCmd = f"{execPATH} -o {output} {args} {value[2]}"
            
            if dryrun:
                log.info(f"[dryrun] {fastqcCmd}")
            else:
                if slurm == True:  # check if slurm job scheduling is enabled or not

                    try:
                        job = mu.clusterRun(job_name=
                            'fastqc', sout= os.path.join(output,"fastqc.out"), serror= os.path.join(output,"fastqc.err"), command=fastqcCmd, mem=mem, cpu=cpu, tasks=task, dep=dep)

                        job_id.append(job)

                        log.info(
                            f"Job submitted on slurm successfully for {inputPair} with {job}")

                    except Exception:

                        log.error("Slurm job sumission failed")

                else:

                    local_cmds[key] = fastqcCmd

    if local_cmds:
        mu.local_run(program="fastqc", command=local_cmds, outdir=output, jobs=jobs)
    return job_id, fastqcOut
