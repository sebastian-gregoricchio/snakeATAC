#######################################
## SPACCa: Snakefile for CRAM to BAM ##
#######################################

import os
#conda_prefix = str(os.environ["CONDA_PREFIX"])

import sys
#sys.path.insert(1, conda_prefix+"/lib/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+"/site-packages")

from typing import List
import pathlib
import re
import numpy
import pandas as pd
import math
import glob
from itertools import combinations


# Define general variables
genome_fasta = str(config["genome_fasta"])


### working directory
home_dir = os.path.join(config["bam_out_directory"],"")
#shell('mkdir -p {home_dir}')
workdir: home_dir


# get the unique samples names and other variables
if not (os.path.exists(config["cram_directory"])):
    os.system("printf '\033[1;31m\\n!!! *cram_directory* does not exist !!!\\n\\n\033[0m'")


CRAMFILES = []
# Iterate directory
for file in os.listdir(config["cram_directory"]):
    # check only text files
    if file.endswith(".cram"):
        CRAMFILES.append(file)

CRAMNAMES = numpy.unique([re.sub("[.]cram$", "", i) for i in CRAMFILES])


### Generation of global wildcard_constraints
# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
            raise ValueError("constraint_to(): Expected a list, got str instead")
    return "({})".format("|".join(values))

wildcard_constraints:
    CRAMS = constraint_to(CRAMNAMES)


# ruleorder: fastQC_filtered_BAM > normalized_bigWig > raw_bigWig

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions

rule A_initialization:
    input:
        bam = expand(os.path.join(config["bam_out_directory"], "{crams}.bam"), crams = CRAMNAMES),
        bam_bai = expand(os.path.join(config["bam_out_directory"], "{crams}.bam.bai"), crams = CRAMNAMES)
    threads: 1
    shell:
        """
        printf '\033[1;36mCrams converted to bams!\\n\033[0m'
        """

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================


# samtools cram to bam -----------------------------------------------------------------------------
rule B_cram_to_bam:
    input:
        cram = os.path.join(config["cram_directory"], "{CRAMS}.cram")
    output:
        bam = os.path.join(config["bam_out_directory"], "{CRAMS}.bam"),
        bam_bai = os.path.join(config["bam_out_directory"], "{CRAMS}.bam.bai")
    params:
        sample = "{CRAMS}",
        genome_fasta = genome_fasta,
        log_dir = os.path.join(config["bam_out_directory"], "log_cramToBam/")
    threads:
        max(math.floor(workflow.cores/len(CRAMNAMES)), 1)
    log:
        out = os.path.join(config["bam_out_directory"], "log_cramToBam/{CRAMS}_cramToBam_samtools.log")
    shell:
        """
        printf '\033[1;36mProccessing: {params.sample}...\\n\033[0m'

        mkdir -p {params.log_dir}

        $CONDA_PREFIX/bin/samtools view -@ {threads} -b -T {params.genome_fasta} -o {output.bam} {input.cram} &> {log.out}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam} {output.bam_bai}
        """


# ------------------------------------------------------------------------------
#                                 END pipeline
# ------------------------------------------------------------------------------
