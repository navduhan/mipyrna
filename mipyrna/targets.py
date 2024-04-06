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

# Intialize the logger

log = MiPyRNALogger(mode='a', log='miRNA_target')

def run_targets(miRNA_file=None, configFile=None, mRNA_file=None, slurm=False, outdir=".", dryrun=False):
    
    if configFile != None:

        config = mu.parse_config_file(configFile)

    else:
        stream = pkg_resources.resource_stream('pysirna', "param/miranda.ini")

        config = mu.parse_config_file(stream.name)

        log.info("Using default config file miranda.ini")

    miranda_config = config[list(config.keys())[0]]

    out = "targets_raw"
    
    if os.path.exists(outdir):

        output1 = os.path.join(outdir,out)

        output = mu.make_directory(output1, dryrun=dryrun)

    else:
        
        output = mu.make_directory(out, dryrun=dryrun)
        
    
    
    return