"""
miPyRNA Utility Module

This module provides shared helper functions for miPyRNA: sample file parsing,
configuration parsing, FASTA handling, job execution (local and SLURM), CPU
management and file and directory helpers.

Features:
    - Sample file parsing with differential expression combinations and targets
    - INI configuration parsing and thread option replacement
    - FASTA reading, subsetting, sequence extraction and U to T conversion
    - SLURM job submission and status checking
    - Parallel local command execution
    - File name, extension and directory helpers

Configuration:
    Tool parameters are read from ``<program>.ini`` files by ``parse_config_file``.
    The sample file is a whitespace-delimited table (columns: sample name, identifier,
    group, file names) parsed by ``read_input_file``.

Dependencies:
    - numpy, pandas: Array and table handling
    - biopython, pyfaidx: FASTA parsing and sequence extraction
    - psutil: CPU count
    - requests: File downloads
    - SLURM (optional): ``sbatch`` and ``squeue``

Classes / Functions / Exceptions:
    - read_input_file: Parse the sample file
    - fasta_to_df: FASTA to dataframe
    - get_sequence: Extract a sequence from a reference FASTA
    - parse_config_file: Parse a program INI file
    - clusterRun: Submit a SLURM job
    - local_run: Run commands locally
    - check_status: Check a SLURM job
    - get_cpu: CPU count to use
    - split_cpu: Split cores into jobs and threads
    - replace_cpu: Set the thread option of a tool config
    - get_basename, get_directory, get_parent, get_file_extension: Path helpers
    - make_directory: Create a directory
    - set_resume: Reuse existing directories
    - preprocess: Feature dictionaries to array
    - download_file: Download a URL
    - exclude_fasta_records: Remove FASTA records by name
    - replace_U_with_T: Convert U to T in a FASTA file

:Author: Naveen Duhan
"""

import configparser
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import psutil
import pyfaidx
import requests
from Bio import Seq, SeqIO

from mipyrna.logger import MiPyRNALogger

log = MiPyRNALogger(mode='a', log='utility')

def read_input_file(infile, inpath, paired = False):

    """Read the input sample file into a dictionary and build all comparison combinations.

    Args:
        infile: Input sample file describing the project.
        inpath: Directory containing the input read files.
        paired: True if reads are paired-end.

    Returns:
        Dictionary with ``samples`` (key to ``[name, group, file(s)]``), ``combinations``
        (group pairs for differential analysis) and ``targets`` (dataframe of samples and groups).

    Raises:
        Exception: Re-raised if the targets dataframe cannot be built.

    System Exit:
        Exits if the input file cannot be read.
    """
    samples = {} # dictionary to collect sample information from input sample file
    factors = [] # list for collecting the Identifier infromation from input sample file
    combinations = [] # list for genrating combinations for differential expression analysis

    try:
        with open(infile) as file:

            log.info("Reading input samples File ")

            for line in file:

                if not line.startswith("#") and not line.startswith("SampleName"):
                    line = line.strip()
                    lines = re.split(r'\s+', line.rstrip())

                    if paired:

                        samples[lines[1]] = [lines[1],lines[2],os.path.join(inpath,lines[3]), os.path.join(inpath,lines[4])]

                    else:

                        samples[lines[1]] = [lines[1], lines[2], os.path.join(inpath,lines[3])]

                    if lines[2] not in factors:

                        factors.append(lines[2])
    except Exception:

        log.error("Please provide a valid input file")

        sys.exit()

    finally:

        log.info("Input file %s read succesfully", infile)

    try:

        # create combinations from factors 
        for i in factors:
            for j in factors:

                if i != j and j+"-"+i not in combinations:

                    combinations.append(i+"-"+j)
    except Exception:

        log.error("Please provide a valid input file")
    
    finally:

        log.info("Combination created succesfully from %s", infile)

        samplename = []
        sample = []

    try:
        for k, s in samples.items():

            samplename.append(k)

            sample.append(s[1])

        targets = pd.DataFrame(samplename,index=list(samplename))

        targets = targets.assign(sample=sample)

    except Exception:

        log.error("Please provide a valid input file")
        raise
    
    finally:

        log.info("targets dataframe for differenatial created succesfully from %s", infile)

    return {'samples': samples, "combinations": combinations, "targets": targets}


def fasta_to_df(file):
    """Read a FASTA file into a dataframe with id, seq and description.

    Args:
        file: FASTA file.

    Returns:
        DataFrame with columns ``id``, ``seq`` and ``description``.
    """

    colums = ['id', 'seq', 'description']

    fastas = SeqIO.parse(file, 'fasta')

    sequences = [[f.id, str(f.seq), f.description] for f in fastas]

    seqdf = pd.DataFrame(sequences, columns=colums)

    return seqdf

def get_sequence(reference_fasta, coordinates):
        """Return the sequence for [chrom, start, end, strand] from a reference FASTA.

        Args:
            reference_fasta: Reference FASTA file.
            coordinates: Sequence of chromosome, start, end and strand; the
                reverse complement is returned for ``'-'``.

        Returns:
            The extracted sequence as a string.
        """
        
        sequences = pyfaidx.Fasta(reference_fasta)

        seq = str(sequences[coordinates[0]][coordinates[1]:coordinates[2]])

        if coordinates[3] == '-':
            seq = str(Seq.Seq(seq).reverse_complement())
        return seq


def parse_config_file(infile):
    """Parse a program configuration file.

    Options whose value contains ``NA`` are dropped.

    Args:
        infile: Path, file object or resource of a ``<program>.ini`` file.

    Returns:
        Dictionary mapping each section to the list of its option values
        (empty if the file could not be parsed).
    """
    sections_dict = {}

    config = configparser.ConfigParser()

    try:
        source_name = str(infile)
        if hasattr(infile, "read"):
            source_name = getattr(infile, "name", source_name)
            config.read_file(infile)
        elif hasattr(infile, "open") and not isinstance(infile, (str, os.PathLike)):
            with infile.open("r", encoding="utf-8") as handle:
                source_name = getattr(handle, "name", source_name)
                config.read_file(handle)
        else:
            with open(infile, encoding="utf-8") as handle:
                source_name = getattr(handle, "name", source_name)
                config.read_file(handle)

        sections = config.sections()

        for section in sections:

            options = config.options(section)

            temp_dict = {}
            voption = []

            for option in options:

                cc = config.get(section, option)

                temp_dict[option] = cc

            for value in temp_dict.values():

                if 'NA' not in value:

                    voption.append(value)

            sections_dict[section] = voption

        log.info("Config generated succesfully from %s", source_name)

    except Exception:

        log.error("Please provide a valid config file")
    
    return sections_dict

def clusterRun(job_name='pysirna',sout="pysirna", serror="pysirna", partition='compute', command='command', time=4, mem=10, cpu=8, tasks=1, dep=''):
    """Submit a job to a cluster with the SLURM scheduler.

    Args:
        job_name: SLURM job name.
        sout: Prefix of the standard output file (``.out`` is appended).
        serror: Prefix of the standard error file (``.err`` is appended).
        partition: SLURM partition.
        command: Command to execute.
        time: Time allotment in hours.
        mem: Memory in GB.
        cpu: CPUs per task.
        tasks: Number of tasks.
        dep: Job id this job depends on (``afterok``).

    Returns:
        The SLURM job id.

    Raises:
        Exception: Re-raised if job submission fails.
    """
    try:
        if dep != '':
            
            dep = f'--dependency=afterok:{dep} --kill-on-invalid-dep=yes '

        sbatch_command = f'sbatch -J {job_name} -o {sout}.out -e {serror}.err -t {time}:00:00 -p {partition} --mem={mem}000 --cpus-per-task={cpu} --ntasks={tasks} --wrap="{command}" {dep}'

        sbatch_response = subprocess.getoutput(sbatch_command)

        job_id = sbatch_response.split(' ')[-1].strip()

    except Exception:

        log.error("Job submission failed")
        raise

    return job_id


def local_run(program=None, command: str | dict = "", outdir=".", message=None, jobs=1):
    """Run shell commands locally, writing output to ``logs/<program>.out`` and ``logs/<program>.err``.

    Args:
        program: Program name used for output file names.
        command: A command, or a dict of sample key to command to run ``jobs`` at a time.
        outdir: Directory whose ``logs`` subdirectory receives the output and error files.
        message: Unused.
        jobs: Maximum number of commands running at once.

    Returns:
        Exit code, or a dict of sample key to exit code (None if a job could not run).
    """

    if message is None:
        message = {}

    logdir = os.path.join(outdir, "logs")
    os.makedirs(logdir, exist_ok=True)

    def run(name, cmd):
        """Run one command, logging its output to ``<name>.out`` and ``<name>.err``.

        Args:
            name: Output file prefix.
            cmd: Shell command.

        Returns:
            Exit code, or None if the job could not be run.
        """
        try:
            with open(os.path.join(logdir, f"{name}.out"), 'w+') as fout, \
                    open(os.path.join(logdir, f"{name}.err"), 'w+') as ferr:
                rc = subprocess.call(cmd, shell=True, stdout=fout, stderr=ferr)
        except Exception:
            log.error(f"Job submission failed for {name}")
            return None
        if rc != 0:
            log.warning(f"{name} exited with code {rc}; see logs/{name}.err")
        else:
            log.info(f"{name} finished")
        return rc

    if isinstance(command, dict):
        log.info(f"Running {program} on {len(command)} sample(s), {jobs} at a time")
        with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
            codes = pool.map(lambda k: run(f"{program}_{k}", command[k]), command)
            return dict(zip(command, codes))
    return run(program, command)


def check_status(job_id):
    """Check the status of a SLURM job.

    Args:
        job_id: SLURM job id.

    Returns:
        True if the job is no longer in the queue, otherwise False.
    """
    d = subprocess.check_output('squeue -j '+str(job_id), shell=True, universal_newlines=True)

    data = list(re.split(r"\s+ ", d))

    return len(data) == 6

def get_cpu():

    """Get the number of CPUs to use.

    Returns:
        80% of the CPU count of the system, rounded down.
    """

    return math.floor((psutil.cpu_count() or 1)*0.8)

def split_cpu(cpu, jobs=None, tasks=1):
    """Split a total core budget into parallel jobs and threads per job.

    Args:
        cpu: Total cores to use.
        jobs: Samples to run at once. Defaults to one job per 4 cores.
        tasks: Number of independent tasks (samples) available.

    Returns:
        Tuple of (jobs, threads per job).
    """
    cpu = max(1, int(cpu))
    if jobs is None:
        jobs = cpu // 4
    jobs = max(1, min(int(jobs), cpu, max(1, int(tasks))))
    return jobs, max(1, cpu // jobs)


def replace_cpu(args, args2, cpu=None):
    """Set the thread option in a tool config.

    Args:
        args: List of config option strings.
        args2: Option names that denote the thread count.
        cpu: Threads to use, capped at the machine core count. Defaults to the
            config value, capped at 80% of the machine.

    Returns:
        The config list with the thread option replaced.
    """

    mat = [i for i in args if any(j in i for j in args2)]

    opt , num = mat[0].split(" ")
    if cpu is not None:
        num = min(int(cpu), psutil.cpu_count() or 1)
    else:
        count = get_cpu()

        if int(num) > count:
            num = count
            log.warning("number of threads changed to available %s",count)

    mat2 = ' '.join([opt,str(num)])
    data = []
    for i in args:
        for j in args2:
            if j in i:
                i= mat2
        data.append(i)

    return data


def get_basename(filePATH):
    """Return the file name from a full path.

    Args:
        filePATH: Path to file.

    Returns:
        The base name of the file.
    """

    return os.path.basename(filePATH)


def get_directory(filePATH):
    """Return the directory of a file.

    Args:
        filePATH: Path to file.

    Returns:
        The directory part of the path.
    """

    return os.path.dirname(filePATH)


def get_parent(filePATH):
    """Return the path without the file extension.

    Args:
        filePATH: Path to file.

    Returns:
        The path with its extension removed.
    """

    return os.path.splitext(filePATH)[0]


def get_file_extension(filePATH):
    """Return the extension of a file.

    Args:
        filePATH: Path to file.

    Returns:
        The file extension including the dot.
    """

    return os.path.splitext(filePATH)[1]

_resume = False


def set_resume(flag=True):
    """Make :func:`make_directory` reuse existing directories instead of numbering new ones.

    Args:
        flag: True to reuse existing directories, False to restore the default.
    """
    global _resume
    _resume = flag


def make_directory(dir, dryrun=False):
    """Create a directory, adding a numeric suffix if it already exists.

    With :func:`set_resume` active an existing directory is reused as it is.

    Args:
        dir: Directory name.
        dryrun: If True, return the name without creating anything.

    Returns:
        Name of the created (or, in dry-run, requested) directory.
    """
    if dryrun:

        return dir
    
    else:
        outputdir = os.path.abspath(dir) 

        if _resume:
            os.makedirs(outputdir, exist_ok=True)
            return outputdir

        if os.path.exists(dir):
            parent, base = os.path.split(outputdir)
            counter = 0
            for sibdir in os.listdir(parent):
                if sibdir.startswith(base +'.'):
                    ext = sibdir[len(base)+1:]
                    if ext.isdecimal():
                        counter = max(counter, int(ext))
            outdir = os.path.join(parent, base+'.'+str(counter+1))

            os.mkdir(outdir)

            log.info(f"Succesfully created directory {outdir}")
            
        else:
            outdir = outputdir
            os.mkdir(outdir)
            log.info(f"Succesfully created directory {outdir}")
    

    return outdir

def preprocess(df, classLabel,size):
    """Convert a list of feature dictionaries to an array with a matching label list.

    Args:
        df: Iterable of feature dictionaries.
        classLabel: Class label assigned to every entry.
        size: Number of features per entry.

    Returns:
        Tuple of (array of shape ``(n, size)``, list of labels).
    """
    tmp = []
    _classLabels = []
    
    for dic in df:
        
        tmp.append(list(dic.values()))
        _classLabels.append(classLabel)
    tmp = np.array(tmp)



    aa = np.zeros((len(tmp),int(size)))

    for i,s in enumerate(tmp):

        aa[i]=list(s)

    return aa, _classLabels

def download_file(url, out):
    """Download a URL to a file.

    Args:
        url: URL to download.
        out: Output file path.

    Returns:
        The output file path.

    Raises:
        requests.HTTPError: If the request fails.
    """
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    with open(out, 'wb') as fp:
        fp.write(response.content)
    return out


def exclude_fasta_records(input_file, names, output_file):
    """Write the FASTA records whose first header word is not in names.

    Args:
        input_file: Input FASTA file.
        names: Record names to exclude.
        output_file: Output FASTA file.
    """
    names = set(names)
    keep = False
    with open(input_file) as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                keep = line[1:].split()[0] not in names
            if keep:
                f_out.write(line)

def replace_U_with_T(input_file, output_file):
    """Write a copy of a FASTA file with U replaced by T in the sequences.

    Args:
        input_file: Input FASTA file.
        output_file: Output FASTA file.

    Returns:
        The output file path.
    """
    with open(input_file) as f_in:
        fasta_lines = f_in.readlines()

    # Iterate through each line in the FASTA file
    modified_lines = []
    for line in fasta_lines:
        # Check if the line is a sequence line (starts with '>')
        if line.startswith('>'):
            modified_lines.append(line)  # If it's a sequence header, keep it unchanged
        else:
            # Replace 'U' with 'T' in the sequence
            modified_sequence = line.replace('U', 'T')
            modified_lines.append(modified_sequence)

    # Write the modified content to the output file
    with open(output_file, 'w') as f_out:
        f_out.write(''.join(modified_lines))
        
    return output_file
