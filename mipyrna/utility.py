#!/usr/bin/python3

"""
title: miPyRNA utility functions
Author: Naveen Duhan
"""

import os, sys
import re
import math
import pandas as pd
import numpy as np
import pyfaidx
import psutil
import subprocess
import configparser
import pandas as pd
from Bio import SeqIO, Seq

from mipyrna.logger import MiPyRNALogger

log = MiPyRNALogger(mode='a', log='utility')

def read_input_file(infile, inpath, paired = False):

    """This function reads input sample file and convert into a dictionary. It also make all possible combination for DEG analysis. Target dataframe for differential analysis.

    :param inputFile: input sample file containing the infromation about project

    :param inputPath: Path for input fastq files

    :param pairedEND: Check if reads are paired end]. Defaults to False

    :return:   samples, combinations and targets for differential expression
    :rtype: A dictionary
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
                    lines = re.split('\s+', line.rstrip())

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

                if i != j:

                    if j+"-"+i not in combinations:

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

        targets = pd.DataFrame(samplename,index=[i for i in samplename])

        targets = targets.assign(sample=sample)

    except Exception:

        log.error("Please provide a valid input file")
    
    finally:

        log.info("targets dataframe for differenatial created succesfully from %s", infile)

    return {'samples': samples, "combinations": combinations, "targets": targets}


def fasta_to_df(file):

    colums = ['id', 'seq', 'description']

    fastas = SeqIO.parse(file, 'fasta')

    sequences = [[f.id, str(f.seq), f.description] for f in fastas]

    seqdf = pd.DataFrame(sequences, columns=colums)

    return seqdf

def get_sequence(reference_fasta, coordinates):
        
        sequences = pyfaidx.Fasta(reference_fasta)

        seq = str(sequences[coordinates[0]][coordinates[1]:coordinates[2]])

        if coordinates[3] == '-':
            seq = str(Seq.Seq(seq).reverse_complement())
        return seq


def parse_config_file(infile):
    """
    This function parse the config file for all the programs used in pySeqRNA

    :param: configFile: <program>.ini config file containing arguments. 

    :retrun: Program specific arguments 

    :rtype: a dictionary
    """
    sections_dict = {}

    config = configparser.ConfigParser()

    try:

        config.read([infile])

        sections = config.sections()

        for section in sections:

            options = config.options(section)

            temp_dict = {}
            voption = []

            for option in options:

                cc = config.get(section, option)

                temp_dict[option] = cc

            for k, value in temp_dict.items():

                if 'NA' not in value:

                    voption.append(value)

            sections_dict[section] = voption

        log.info("Config generated succesfully from %s", infile)

    except Exception:

        log.error("Please provide a valid config file")
    
    return sections_dict

def clusterRun(job_name='pysirna',sout="pysirna", serror="pysirna", partition='compute', command='command', time=4, mem=10, cpu=8, tasks=1, dep=''):
    """
    This function is for submitting job on cluster with SLURM job Scheduling

    :param job_name: Slurm job name on HPC. Defaults to 'pyseqRNA'.

    :param command:  Command to excute on HPC.

    :param time: Slurm Job time allotment. 

    :param mem: Memory to use in GB.

    :param cpu: Number of CPUs to use for the job.

    :param tasks: Number of tasks to execute.

    :param dep: Slurm Job dependency. Defaults to ''.

    :returns:

        :rtype: Slurm sbatch ID
    """
    try:
        if dep != '':
            
            dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)

        sbatch_command = 'sbatch -J {} -o {}.out -e {}.err -t {}:00:00 -p {} --mem={}000 --cpus-per-task={} --ntasks={} --wrap="{}" {}'.format(
            job_name, sout, serror, time, partition, mem, cpu, tasks,  command, dep)

        sbatch_response = subprocess.getoutput(sbatch_command)

        job_id = sbatch_response.split(' ')[-1].strip()

    except Exception:

        log.error("Job submission failed")

    return job_id


def local_run(program=None, command=None, outdir=None, message={}):

    try:

        with open(os.path.join(outdir, f"{program}.out"), 'w+') as fout:
                with open(os.path.join(outdir, f"{program}.err"), 'w+') as ferr:
                    subprocess.call(command, shell=True,stdout=fout,stderr=ferr)
                    log.info(f"Job successfully submited for {message}")
    except Exception:
        
        log.error("Job sumission failed")
        
    return

def check_status(job_id):
    """
    This function is check status of slurm job

    :param job_id: slurm job id

    :returns: True/False

    :rtype: If job completed return True. Default False.
    """
    d = subprocess.check_output('squeue -j '+str(job_id), shell=True, universal_newlines=True)

    data = list(re.split("\s+ ",d))

    if len(data)==6:

        return True

    return False

def get_cpu():

    """
    This function get actual CPU count of the system 

    :returns: Integer 
    
    :rtype: int with 80 % of CPU count
    """

    return math.floor(psutil.cpu_count()*0.8)

def replace_cpu(args, args2):
    '''
    This function replace the actual CPU in config file.

    :returns: Change CPU count to 80% of available CPU 
    '''

    mat = [i for i in args if any(j in i for j in args2)]
 

    opt , num = mat[0].split(" ")
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
    """
    This function get the base name of the file from full path 

    :param  filePATH: Path to file.
    """

    return os.path.basename(filePATH)


def get_directory(filePATH):
    """
    This function retrun directory of a file 

    :param  filePATH: Path to file.
    """

    return os.path.dirname(filePATH)


def get_parent(filePATH):
    """
    This function return the file name without extension

    :param  filePATH: Path to file.
    """

    return os.path.splitext(filePATH)[0]


def get_file_extension(filePATH):
    """
    This function return the extension of file 

    :param  filePATH: Path to file.
    """

    return os.path.splitext(filePATH)[1]

def make_directory(dir, dryrun=False):
    """
    This function create a directory 

    :param  dir: Directory name.

    :returns: Name of created directory.
    """
    if dryrun:

        return dir
    
    else:
        outputdir = os.path.abspath(dir) 

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

def replace_U_with_T(input_file, output_file):
    with open(input_file, 'r') as f_in:
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