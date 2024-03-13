#######################################
## SPACCa: Snakefile for DNA mapping ##
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
from itertools import combinations


# Define general variables
genome_fasta = str(config["genome_fasta"])


### working directory
home_dir = os.path.join(config["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir



### Other variables
bwa_version = "bwa-mem2"
bwa_idx = ".bwt.2bit.64"


if (eval(str(config["remove_duplicates"])) == True):
    DUP = "dedup"
else:
    DUP = "mdup"



# get the unique samples names and other variables
if not (os.path.exists(config["fastq_directory"])):
    os.system("printf '\033[1;31m\\n!!! *fastq_directory* does not exist !!!\\n\\n\033[0m'")

FILENAMES = next(os.walk(config["fastq_directory"]))[2]
RUNNAMES = [re.sub(rf"{config['fastq_suffix']}$", "", i) for i in FILENAMES]
SAMPLENAMES = numpy.sort(numpy.unique([re.sub(rf"{config['read_suffix'][0]}|{config['read_suffix'][1]}.*$", "", i) for i in RUNNAMES]))


# Chromosome remove chr_remove_pattern
if (len(config["remove_other_chromosomes_pattern"]) > 0):
    chr_remove_pattern = '^chrM|^M|'+config["remove_other_chromosomes_pattern"]
else:
    chr_remove_pattern = '^chrM|^M'


### Optional analysis outputs
if (eval(str(config["run_fastq_qc"])) == True):
    multiqc_fastq = "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.html"
else:
    multiqc_fastq = []


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
    SAMPLE = constraint_to(SAMPLENAMES),
    RUNS = constraint_to(RUNNAMES)


# ruleorder: fastQC_filtered_BAM > normalized_bigWig > raw_bigWig

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        multiqc_bam_report = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.html",
        multiqc_fastq = multiqc_fastq
    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================

# Generate bwa index if required [OPTIONAL]  ---------------------------------------------
if not os.path.exists(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),bwa_idx])):
    rule generate_genome_index:
        input:
            genome = ancient(genome_fasta)
        output:
            genome_bwt_2bit_64 = ''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),bwa_idx])
        params:
            bwa = bwa_version
        threads:
            workflow.cores
        benchmark:
            "benchmarks/generate_genome_index/generate_genome_index---benchmark.txt"
        shell:
            """
            printf '\033[1;36mGenerating the genome index...\\n\033[0m'
            ${{CONDA_PREFIX}}/bin/{params.bwa} index {input.genome}
            samtools faidx {input.genome}
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """



# cutdapat -------------------------------------------------------------------------------
rule cutadapt_PE:
    input:
        R1 = os.path.join(config["fastq_directory"], "".join(["{SAMPLE}", config['read_suffix'][0], config['fastq_suffix']])),
        R2 = os.path.join(config["fastq_directory"], "".join(["{SAMPLE}", config['read_suffix'][1], config['fastq_suffix']]))
    output:
        R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"])),
        R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][1], "_trimmed.fastq.gz"]))
    params:
        sample = "{SAMPLE}",
        opts = str(config["cutadapt_trimm_options"]),
        fw_adapter_sequence = str(config["fw_adapter_sequence"]),
        rv_adapter_sequence = str(config["fw_adapter_sequence"])
    log:
        out = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.out",
        err = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.err"
    threads:
        workflow.cores
    benchmark:
        "benchmarks/cutadapt_PE/cutadapt_PE---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: reads trimming...\\n\033[0m'
        mkdir -p 01_trimmed_fastq/logs/

        ${{CONDA_PREFIX}}/bin/cutadapt \
        -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 \
        -a {params.fw_adapter_sequence} -A {params.rv_adapter_sequence} {params.opts} \
        -o {output.R1_trimm} -p {output.R2_trimm} {input.R1} {input.R2} > {log.out} 2> {log.err}
        """


# BWA mapping -----------------------------------------------------------------------------
if (eval(str(config["umi_present"])) == True):
    rule BWA_PE_UMI:
        input:
            R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"])),
            R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][1], "_trimmed.fastq.gz"])),
            genome_index = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),bwa_idx]))
        output:
            align_summary = "02_BAM/bwa_summary/{SAMPLE}.BWA_summary.txt",
            bam = temp("02_BAM/{SAMPLE}.sorted.bam")
        params:
            bwa_opts = str(config["bwa_options"]),
            sample = "{SAMPLE}",
            bwa = bwa_version,
            genome_fasta = genome_fasta
        log:
            out = "02_BAM/logs/{SAMPLE}.sort.log"
        threads:
            max((workflow.cores - 1), 1)
        benchmark:
            "benchmarks/BWA_PE_UMI/BWA_PE_UMI---{SAMPLE}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: bwa mapping...\\n\033[0m'
            mkdir -p 02_BAM/logs
            mkdir -p 02_BAM/bwa_summary

            ${{CONDA_PREFIX}}/bin/{params.bwa} mem \
            -t {threads} \
            -R '@RG\\tID:{params.sample}\\tDS:{params.sample}\\tPL:ILLUMINA\\tSM:{params.sample}' \
            -M {params.genome_fasta} {input.R1_trimm} {input.R2_trimm} |
            mawk -v OFS="\\t" '
            {{
        	    if(substr($0, 1, 1) != "@") {{
        	    # line below behaves like samtools view -bF 3840 -
        	    if ($2 >= 256) {{next}}
        	    # remove the umi from the readname and append it as RX tag
        	    n=length($1)
        	    for (p=0; substr($1, n-p, 1) != ":"; p++);
        	    u="\tRX:Z" substr($1, n-p, p+1)
        	    $1 = substr($1, 0, n-p-1)
        	    }}
              print $0 u
            }}' |
	          $CONDA_PREFIX/bin/samtools fixmate -m - - |
            $CONDA_PREFIX/bin/samtools sort -m 2G -T 02_BAM/{params.sample} -@ {threads} -O bam - > {output.bam} 2> {log.out};
            $CONDA_PREFIX/bin/samtools flagstat {output.bam} > {output.align_summary}
            """

else:
    rule BWA_PE:
        input:
            R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][0], "_trimmed.fastq.gz"])),
            R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", config['read_suffix'][1], "_trimmed.fastq.gz"])),
            genome_index = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),bwa_idx]))
        output:
            align_summary = "02_BAM/bwa_summary/{SAMPLE}.BWA_summary.txt",
            bam = temp("02_BAM/{SAMPLE}.sorted.bam")
        params:
            bwa_opts = str(config["bwa_options"]),
            sample = "{SAMPLE}",
            bwa = bwa_version,
            genome_fasta = genome_fasta
        threads:
            max((workflow.cores - 1), 1)
        log:
            out = "02_BAM/logs/{SAMPLE}.sort.log"
        benchmark:
            "benchmarks/BWA_PE/BWA_PE---{SAMPLE}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: bwa mapping...\\n\033[0m'
            mkdir -p 02_BAM/bwa_summary

            ${{CONDA_PREFIX}}/bin/{params.bwa} mem \
            -t {threads} \
            -M {params.genome_fasta} {input.R1_trimm} {input.R2_trimm} |
	          $CONDA_PREFIX/bin/samtools fixmate -m - - |
            $CONDA_PREFIX/bin/samtools sort -m 2G -T 02_BAM/{params.sample} -@ {threads} -O bam - > {output.bam} 2> {log.out};
            $CONDA_PREFIX/bin/samtools flagstat {output.bam} > {output.align_summary}
            """


# samtools mapq filter -----------------------------------------------------------------------------
rule MAPQ_MT_filter:
    input:
        source_bam = "02_BAM/{SAMPLE}.sorted.bam"
    output:
        bam_mapq_only = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), ".bam"]))),
        bam_mapq_only_sorted = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"]))),
        bam_mapq_only_sorted_index = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bai"]))),
        idxstats_file = "02_BAM/reads_per_chromosome/{SAMPLE}_idxstats_read_per_chromosome.txt",
        bam_mapq_only_sorted_woMT = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam"]))),
        bam_mapq_only_sorted_woMT_index = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam.bai"])))
    params:
        sample = "{SAMPLE}",
        MAPQ_threshold = config["MAPQ_threshold"],
        chr_remove_pattern = chr_remove_pattern
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/MAPQ_MT_filter/MAPQ_MT_filter---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: filtering MAPQ and remove mitochondrial chromosome...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools view -@ {threads} -h -q {params.MAPQ_threshold} {input.source_bam} -o {output.bam_mapq_only}

        $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted} {output.bam_mapq_only_sorted_index}

        $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} > {output.idxstats_file}

        printf '\033[1;36m{params.sample}: Removing MT from BAM...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} | cut -f 1 | grep -v -E '{params.chr_remove_pattern}' | xargs ${{CONDA_PREFIX}}/bin/samtools view -@ {threads} -b {output.bam_mapq_only_sorted} > {output.bam_mapq_only_sorted_woMT}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_woMT} {output.bam_mapq_only_sorted_woMT_index}
        """


if (eval(str(config["umi_present"])) == True):
    rule gatk4_markdups_umiAware:
        input:
            bam_mapq_only_sorted = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam"])),
            bam_mapq_only_sorted_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam.bai"]))
        output:
            bam_mdup = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, ".bam"])),
            bai_mdup = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, ".bai"])),
            umi_metrics = "02_BAM/umi_metrics/{SAMPLE}_UMI_metrics.txt",
            dup_metrics = "02_BAM/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
            flagstat_filtered = os.path.join("02_BAM/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_flagstat.txt"]))
        params:
            remove_duplicates = (str(config["remove_duplicates"])).lower(),
            sample = "{SAMPLE}"
        log:
            out = "02_BAM/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
            err = "02_BAM/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
        threads:
            workflow.cores
        benchmark:
            "benchmarks/gatk4_markdups_umiAware/gatk4_markdups_umiAware---{SAMPLE}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: UMI-aware gatk MarkDuplicates...\\n\033[0m'

            mkdir -p 02_BAM/umi_metrics
            mkdir -p 02_BAM/MarkDuplicates_metrics
            mkdir -p 02_BAM/MarkDuplicates_logs
            mkdir -p 02_BAM/flagstat

            $CONDA_PREFIX/bin/gatk UmiAwareMarkDuplicatesWithMateCigar \
            --INPUT {input.bam_mapq_only_sorted} \
            --OUTPUT {output.bam_mdup} \
            --REMOVE_DUPLICATES {params.remove_duplicates} \
            --MAX_EDIT_DISTANCE_TO_JOIN 1 \
            --UMI_METRICS_FILE {output.umi_metrics} \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --UMI_TAG_NAME RX \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY STRICT \
            --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """
else: # no-UMI dedup
    rule gatk4_markdups:
        input:
            bam_mapq_only_sorted = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam"])),
            bam_mapq_only_sorted_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT.bam.bai"]))
        output:
            bam_mdup = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, ".bam"])),
            bai_mdup = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, ".bai"])),
            dup_metrics = "02_BAM/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
            flagstat_filtered = os.path.join("02_BAM/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_flagstat.txt"]))
        params:
            remove_duplicates = (str(config["remove_duplicates"])).lower(),
            sample = "{SAMPLE}"
        log:
            out = "02_BAM/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
            err = "02_BAM/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
        threads:
            workflow.cores
        benchmark:
            "benchmarks/gatk4_markdups/gatk4_markdups---{SAMPLE}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: 'standard' gatk MarkDuplicates...\\n\033[0m'

            mkdir -p 02_BAM/MarkDuplicates_metrics
            mkdir -p 02_BAM/MarkDuplicates_logs
            mkdir -p 02_BAM/flagstat

            $CONDA_PREFIX/bin/gatk MarkDuplicatesWithMateCigar \
            --INPUT {input.bam_mapq_only_sorted} \
            --OUTPUT {output.bam_mdup} \
            --REMOVE_DUPLICATES {params.remove_duplicates} \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY LENIENT \
            --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """


# multiQC aligned bams
rule multiQC_bam:
    input:
        flagstat_filtered = expand(os.path.join("02_BAM/flagstat/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_flagstat.txt"])), sample = SAMPLENAMES)
    output:
        multiqc_bam_report = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.html",
    params:
        out_directory = "03_quality_controls/multiQC_bam_filtered/",
        multiqc_bam_report_name = "multiQC_bam_filtered.html"
    log:
        out = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.out",
        err = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.err"
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/multiQC_bam/multiQC_bam---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPerforming multiQC on filtered bam flagstat/picard...\\n\033[0m'

        mkdir -p {params.out_directory}

        $CONDA_PREFIX/bin/multiqc -f \
        -o {params.out_directory} \
        -n {params.multiqc_bam_report_name} \
        --dirs 02_BAM/MarkDuplicates_metrics 02_BAM/flagstat > {log.err} 2> {log.out}
        """



# fastQC on fastq raw data ----------------------------------------------------------------------------------------
rule fastQC_trimmed_fastq:
    input:
        fastq_trimm = expand(os.path.join("01_trimmed_fastq", "".join(["{runs}_trimmed.fastq.gz"])), runs = RUNNAMES)
    output:
        fastqc_html = expand(os.path.join("03_quality_controls/trimmed_fastq_fastqc","{runs}_trimmed_fastqc.html"), runs = RUNNAMES),
        fastqc_zip = expand(os.path.join("03_quality_controls/trimmed_fastq_fastqc","{runs}_trimmed_fastqc.zip"), runs = RUNNAMES)
    threads:
        workflow.cores
    benchmark:
        "benchmarks/fastQC_trimmed_fastq/fastQC_trimmed_fastq---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPerforming fastQC on trimmed fastq...\\n\033[0m'

        mkdir -p 03_quality_controls/trimmed_fastq_fastqc
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir 03_quality_controls/trimmed_fastq_fastqc 01_trimmed_fastq/*_trimmed.fastq.gz
        """



rule multiQC_trimmed_fastq:
    input:
        fastqc_zip = expand(os.path.join("03_quality_controls/trimmed_fastq_fastqc","{runs}_trimmed_fastqc.zip"), runs = RUNNAMES)
    output:
        multiqc_fastqc_report = "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.html"
    params:
        fastqc_zip_dir = os.path.join(home_dir, "03_quality_controls/trimmed_fastq_fastqc/"),
        out_directory = os.path.join(home_dir, "03_quality_controls/trimmed_fastq_multiQC/"),
        home_dir = home_dir,
        cutadapt_logs = os.path.join(home_dir, "01_trimmed_fastq/logs/"),
        multiqc_fastqc_report_name = "multiQC_report_trimmed_fastq.html"
    log:
        out = os.path.join(home_dir, "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.out"),
        err = os.path.join(home_dir, "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.err")
    threads: 1
    benchmark:
        "benchmarks/multiQC_trimmed_fastq/multiQC_trimmed_fastq---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPerforming multiQC of trimmed fastq...\\n\033[0m'
        printf '\033[1;36mGenerating multiQC report for trimmed fastq...\\n\033[0m'

        cd {params.fastqc_zip_dir}

        $CONDA_PREFIX/bin/multiqc -f \
        -o {params.out_directory} \
        -n {params.multiqc_fastqc_report_name} \
        --dirs-depth 2 \
        --dirs ./ {params.cutadapt_logs} > {log.err} 2> {log.out}

        cd {params.home_dir}
        """

# ------------------------------------------------------------------------------
#                                 END pipeline
# ------------------------------------------------------------------------------
