# Use the official Miniconda base image
FROM ubuntu:24.04
FROM python:3.10
FROM continuumio/miniconda3:latest

# Author: Naveen Duhan
# Title: Dockerfile for miPyRNA

# Install system dependencies
RUN apt-get update && \
    apt-get install -y git bzip2 wget curl nano htop

# Clone the pySeqRNA repository
RUN git clone https://github.com/navduhan/mipyrna.git /mipyrna

# Set the working directory
WORKDIR /mipyrna

# Install dependencies from the YAML file
RUN conda env create -f mipyrna_environment.yml

# Activate the conda environment
# Note: The activation command is generally not needed in RUN instructions, so we directly use conda.
RUN echo "source activate mipyrna-0.2" > ~/.bashrc

# Install pySeqRNA package
RUN /opt/conda/bin/conda run -n mipyrna-0.2 pip install .

WORKDIR /home

# Create data and output directories
RUN mkdir -p /home/data /home/output

# Create a startup script to check architecture and append the STAR script at runtime


# Set the default command to run the startup script, remove it, and then start bash
CMD ["bash"]

##################################################################
##                                                              ##
##  Written by Naveen Duhan (naveen.duhan@outlook.com)          ##
##  Kaundal Bioinformatics Lab, Utah State University           ##
##  Released under the terms of GNU General Public Licence v3   ## 
##                                                              ##
##################################################################