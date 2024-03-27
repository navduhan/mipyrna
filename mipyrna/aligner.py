#!/usr/bin/python
"""
title: mipyrna aligner module 

Author: Naveen Duhan

"""
import os
import sys
import shutil
import subprocess
import pkg_resources
from mipyrna.logger import MiPyRNALogger
from mipyrna import utility as mu
from waiting import wait

# Intialize the logger

log = MiPyRNALogger(mode='a', log='aligner')

class Bowtie_Aligner():
    """Class for Bowtie aligner

    :param configFile: User defined Bowtie config file  
    :param ref_genome: Reference genome file
    :param slurm: To execute jobs with slurm scheduler
    :param dryrun: print coomand without executing
    :param outdir: output directory. Default is from default parameters

    """

    def __init__(self, ref_genome=None, configFile=None, outdir=".", slurm=False, dryrun=False):
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
            stream = pkg_resources.resource_stream('mipyrna', "param/bowtie1.ini")

            self.config = mu.parse_config_file(stream.name)
           
            log.info("Using default config file bowtie.ini")

        return

    def build_index(self, mem=20, tasks = 1, cpu = 8, dep=""):

        const = ['--threads','-p']

        if not self.slurm:
            
            config = mu.replace_cpu(self.config['index'], const)

        else:

            config = self.config['index']

        directory = str(config[0]).split(" ")[0]

        indexName = str(config[1]).split(" ")[0]

        if self.dryrun:
            if os.path.exists(self.outdir):

                output =  os.path.join(self.outdir,directory)

            else:
                output =  os.path.join(self.outdir,directory)
        
        else:

            if os.path.exists(self.outdir):

                output1 =  os.path.join(self.outdir,directory)

                output = mu.make_directory(output1)

            else:
                
                output = mu.make_directory(directory)

        genome_extensions = tuple(['.fa', '.fasta', '.fna', '.fa.gz', 'fna.gz', '.fasta.gz'])
        
        if self.genome.endswith(genome_extensions):

            if self.genome.endswith(".gz"):

                os.system(' '.join(["gunzip", self.genome]))

            os.system(' '.join(["cp", self.genome, output]))

            log.info(f"{mu.get_basename(self.genome)} copied successfully in {directory}")
        
        else:

            log.error(f"Please provide a valid genome fasta file with these {genome_extensions} extensions")

            sys.exit(1)

        if indexName != 'NA':

            basename = os.path.join(output, indexName)

        else:

            basename = os.path.join(output,mu.get_basename(self.genome))

        GenomeFasta = os.path.join(output, mu.get_basename(self.genome))

        arg = ' '.join(config[2:])

        execPATH = shutil.which('bowtie-build')

        if execPATH is None:

            log.error("bowtie-build not found in path")

            sys.exit(1)
        else:

            bowtie_command = f"{execPATH} -f {GenomeFasta} {arg} {basename}"

        if self.dryrun:

                return bowtie_command           
        else:

            if self.slurm:
                try:

                    job_id = mu.clusterRun(job_name='bowtie-build',sout=os.path.join(output, "bowtie_build.out") , serror= os.path.join(output, "bowtie_build.err"), command= bowtie_command, mem=mem, cpu=cpu, tasks=tasks, dep=dep)

                    log.info("Job successfully submited for {} with {} for indexing".format(GenomeFasta, job_id))
                    
                    wait(lambda: mu.check_status(job_id), waiting_for="alignment to finish")

                    log.info("Genome indexing job finished")
                except Exception:

                    log.error("Slurm job sumission failed")
                    
            else:

                try:

                    with open(os.path.join(output, "bowtie_build.out"), 'w+') as fout:
                            with open(os.path.join(output, "bowtie_build.err"), 'w+') as ferr:
                                job_id = subprocess.call(bowtie_command, shell=True,stdout=fout,stderr=ferr)

                                log.info("Job successfully submited for {} for indexing".format(GenomeFasta))

                except Exception:
                    
                    log.error("Job sumission failed")

        return 


    def check_index(self):

        output = os.path.join(self.outdir, 'bowtie_index')

        files = [r'.1.ebwt', r'.2.ebwt', r'.3.ebwt', r'.4.ebwt']

        if (os.path.exists(output) and os.path.isdir(output)):
            for f in files:

                if not (os.path.join(output, _) for _ in os.listdir(output) if _.endswith(f)):

                    return False
                else:

                    return True

        return False

    def run_alignment(self, samplesDict=None, samples=None, fileType='fasta', outType='BWT' , pairedEND=False,  mem= 20, cpu=8, tasks=1, dep=''):


        """This function align reads against indexed reference genome.

            :param sampels: samples dictionary containing sample information.

            :param pairedEND: True if samples are paired.

            :param mem: Provide memory in GB to use. Default 20 GB.

            :param tasks: Number of cpu-tasks to run. Defaults to 1.

            :param cpu: Total number of threads to use. Default 8.

            :param dep: slurm job id on which this job depends.  Defaults to ''.
        """

        consta = ['--threads', '-p']

        if self.slurm:

            config = self.config['alignment']
            
        else:
            config = mu.replace_cpu(self.config['alignment'], consta)
            
        
        directory = str(config[0])

        reference = str(config[1])

        genomeIndex = os.path.join(self.outdir,directory,reference)

        if self.dryrun:

            if os.path.exists(self.outdir):

                output =  os.path.join(self.outdir,"bowtie_results")

            else:
                output =  os.path.join(self.outdir,"bowtie_results")
        
        else:

            if os.path.exists(self.outdir):

                output1 =  os.path.join(self.outdir,"bowtie_results")

                output = mu.make_directory(output1)

            else:
                
                output = mu.make_directory("bowtie_results")

        arg = ' '.join(config[2:])
        
        if outType =='SAM':
            arg += ' --sam'

        outbowtie = {}
        
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

                
                execPATH = shutil.which('bowtie')

                if execPATH is None:

                    log.error("bowtie aligner not found in path")
                    sys.exit(1)
                
                else:

                    if pairedEND:

                        if fileType == 'fastq':

                            bowtie_command = f"{execPATH} -q -x {genomeIndex} {arg} -1 {sample[2]} -2 {sample[3]} > {outBAM}"

                        if fileType == 'fasta':

                            bowtie_command = f"{execPATH} -f -x {genomeIndex} {arg} -1 {sample[2]} -2 {sample[3]} > {outBAM}"

                    else:

                        if fileType == 'fastq':

                            bowtie_command = f"{execPATH} -q -x {genomeIndex} {arg} {sample[2]} > {outBAM}"

                        if fileType == 'fasta':

                            bowtie_command = f"{execPATH} -f -x {genomeIndex} {arg} {sample[2]} > {outBAM}"

                            # print(bowtie_command)

                    if self.dryrun:
                            pass
                    else:
                
                        if self.slurm:
                            try:
                                job = mu.clusterRun(job_name='bowtie_align', sout=os.path.join(output, "bowtie.out") , serror=os.path.join(output, "bowtie.err") ,command= bowtie_command, mem=mem, cpu=cpu, tasks=tasks, dep=dep)
                                job_id.append(job)
                                log.info("Job successfully submited for {} with {} for alignment".format(outPrefix, job))

                            except Exception:

                                log.error("Slurm job sumission failed")

                        else:

                            try:
                                with open(os.path.join(output, "bowtie.out"), 'w+') as fout:
                                    with open(os.path.join(output, "bowtie.err"), 'w+') as ferr:
                                        job = subprocess.call(bowtie_command, shell=True,stdout=fout,stderr=ferr)
                                        job_id.append(job)
                                        log.info("Job successfully completed for {} for alignment".format(outPrefix))

                            except Exception:
                                    
                                    log.exception("Job sumission failed")
                if self.dryrun:

                    return outbowtie 
                
            if self.slurm:
                
                for job in job_id:     
                    
                    wait(lambda: mu.check_status(job), waiting_for="alignment to finish")
                    log.info(f"Alignment completed for job {job}")
                
            return outbowtie
        
        else:
            outBAM = os.path.join(output,"all_reads_bowtie.sam")

            execPATH = shutil.which('bowtie')

            if execPATH is None:

                log.error("bowtie aligner not found in path")
                sys.exit(1)
                
            else:

                if fileType == 'fastq':

                    bowtie_command = f"{execPATH} -q -x {genomeIndex} {arg} {samples[2]} > {outBAM}"

                if fileType == 'fasta':

                    bowtie_command = f"{execPATH} -f -x {genomeIndex} -p 40 -n 0 -e 80 -l 18 -a -m 5 --best --strata {samples[2]} > {outBAM}"

                if self.slurm:
                    try:
                        job = mu.clusterRun(job_name='bowtie_align', sout=os.path.join(output, "bowtie.out") , serror=os.path.join(output, "bowtie.err") ,command= bowtie_command, mem=mem, cpu=cpu, tasks=tasks, dep=dep)
                        job_id.append(job)
                        log.info("Job successfully submited for {} with {} for alignment".format(outPrefix, job))

                    except Exception:

                        log.error("Slurm job sumission failed")

                else:

                    try:
                        with open(os.path.join(output, "bowtie.out"), 'w+') as fout:
                            with open(os.path.join(output, "bowtie.err"), 'w+') as ferr:
                                job = subprocess.call(bowtie_command, shell=True,stdout=fout,stderr=ferr)
                                job_id.append(job)
                                log.info("Job successfully completed for {} for alignment".format(outPrefix))

                    except Exception:
                            
                            log.exception("Job sumission failed")


            return outBAM


