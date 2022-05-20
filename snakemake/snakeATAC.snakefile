from typing import List
import pathlib
import re
import numpy
import os

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

# working diirectory
workdir: config["output_directory"]


# get the unique samples names and other variables
FILENAMES = next(os.walk(config["runs_directory"]))[2]
RUNNAMES = [re.sub(rf"{config['fastq_extension']}$", "", i) for i in FILENAMES]
SAMPLENAMES = numpy.unique([re.sub(rf"{config['runs_suffix'][0]}|{config['runs_suffix'][1]}.*$", "", i) for i in RUNNAMES])

if (eval(str(config["perform_HMCan_correction"])) == True):
    BINS=config["smallBinLength"]
else:
    BINS=config["bigWig_binSize"]


if (eval(str(config["call_summits"])) == True):
    SUMMITS="--call-summits"
else:
    SUMMITS=""


if (config["perform_HMCan_correction"] == True):
    PEAKSDIR = "04_Normalization/HMCan_output/"
    SUMMARYDIR = "05_Overall_quality_and_info/"
    GATKDIR = "06_Variant_calling/"
    PEAKCALLER = "HMCAN"
else:
    PEAKSDIR = "05_Peaks_MACS3/"
    SUMMARYDIR = "06_Overall_quality_and_info/"
    GATKDIR = "07_Variant_calling/"
    PEAKCALLER = "MACS3"

# reference sample definition
#if config["reference_sample"] in SAMPLENAMES:
#    REFERENCE_SAMPLE = config["reference_sample"]
#else:
#    REFERENCE_SAMPLE = SAMPLENAMES[0]


# generation of global wildcard_constraints
wildcard_constraints:
    RUNS=constraint_to(RUNNAMES),
    SAMPLES=constraint_to(SAMPLENAMES)



# Generate optional inputs for HMCan (1) vs Sequencing depth normalization
norm_outputs = []
if (eval(str(config["perform_HMCan_correction"]))):
    norm_outputs.append("04_Normalization/HMCan_output/CONFIGURATION_file_HMCan.txt") #HMCan_config
    norm_outputs.append(expand(os.path.join("04_Normalization/HMCan_output/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_regions.bed"])), sample=SAMPLENAMES)) # HMCan_regions
    norm_outputs.append(expand(os.path.join("04_Normalization/HMCan_output/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_peaks.narrowPeak"])), sample=SAMPLENAMES)) # HMCan_peaks
    norm_outputs.append(expand(os.path.join("04_Normalization/HMCan_output/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_CNV_profile.txt"])), sample=SAMPLENAMES)) # HMCan_CNV_profile

else:
    norm_outputs.append(expand(os.path.join("05_Peaks_MACS3/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_FDR{fdr}_peaks.xls"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), fdr=str(config["FDR_cutoff"]))) # peaks_xls MACS3
    norm_outputs.append(expand(os.path.join("05_Peaks_MACS3/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_FDR{fdr}_peaks.narrowPeak"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), fdr=str(config["FDR_cutoff"]))) # narrowPeaks MACS3


# GATK outputs
if (eval(str(config["call_variants"]))):
    bsqr_table = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_bsqr.table"), sample=SAMPLENAMES))
    dedup_BAM_bsqr = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_bsqr.bam"), sample=SAMPLENAMES))
    dedup_BAM_bsqr_index = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_bsqr.bam.bai"), sample=SAMPLENAMES))
    vcf = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_gatk.vcf.gz"), sample=SAMPLENAMES))
else:
    dedup_BAM_bsqr = []
    dedup_BAM_bsqr_index = []
    vcf = []


if (eval(str(config["call_SNPs"]))):
    snp = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_gatk-snp.vcf"), sample=SAMPLENAMES))
else:
    snp = []


if (eval(str(config["call_indels"]))):
    indels = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_gatk.vcf.gz"), sample=SAMPLENAMES))
else:
    indels = []


correlation_outputs = []
if (len(SAMPLENAMES) > 1):
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "PCA_on_BigWigs_wholeGenome.pdf"))) # PCA
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"))) # heatamap_spearman
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"))) # hetamap_pearson
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"))) # scatterplot_spearman
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf"))) # scatterplot_pearson

# ========================================================================================
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ========================================================================================



# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        # Rule A
        fastQC_raw_zip = ancient(expand(os.path.join("01_fastQC_raw", "{run}_fastqc.zip"), run=RUNNAMES)),

        # Rule B
        multiQC_raw_html = ancient("01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"),

        # Rule C
        #SAM = ancient(expand(os.path.join("01b_SAM_tempFolder/", "{sample}.sam"), sample=SAMPLENAMES)),

        # Rule D
        filtBAM_sorted_woMT = ancient(expand(os.path.join("02_BAM/", "{sample}_mapQ{MAPQ}_sorted_woMT.bam"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        filtBAM_sorted_woMT_index = ancient(expand(os.path.join("02_BAM/", "{sample}_mapQ{MAPQ}_sorted_woMT.bam.bai"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        flagstat_unfiltered_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_UNfiltered_bam.txt"), sample=SAMPLENAMES)),
        flagstat_on_filtered_woMT_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_filtered_bam_woMT.txt"), sample=SAMPLENAMES)),

        # Rule E
        # dedup_BAM = ancient(expand(os.path.join("03_BAM_dedup/unshifted_bams/", "{sample}_mapQ{MAPQ}_sorted_woMT_dedup.bam"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        # dedup_BAM_index = ancient(expand(os.path.join("03_BAM_dedup/unshifted_bams/", "{sample}_mapQ{MAPQ}_sorted_woMT_dedup.bai"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        # dedup_BAM_metrics = ancient(expand(os.path.join("03_BAM_dedup/metrics", "{sample}_metrics_woMT_dedup_bam.txt"), sample=SAMPLENAMES)),
        # dedup_BAM_flagstat = ancient(expand(os.path.join("03_BAM_dedup/flagstat/", "{sample}_flagstat_filtered_bam_woMT_dedup.txt"), sample=SAMPLENAMES)),

        # Rule F
        dedup_BAM_shifted_sorted = ancient(expand(os.path.join("03_BAM_dedup/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_sorted.bam"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        dedup_BAM_shifted_sorted_index = ancient(expand(os.path.join("03_BAM_dedup/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_sorted.bam.bai"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        dedup_BAM_flagstat_shifted_sorted = ancient(expand(os.path.join("03_BAM_dedup/flagstat/", "{sample}_flagstat_woMT_dedup_shifted_sorted.txt"), sample=SAMPLENAMES)),

        # Rule G
        fastQC_zip_BAM = ancient(expand(os.path.join("03_BAM_dedup/fastQC/", "{sample}_mapQ{MAPQ}_sorted_woMT_dedup_fastqc.zip"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),

        # Rule H
        multiQC_BAM_html = ancient("03_BAM_dedup/fastQC/multiQC_dedup_bams/multiQC_report_BAMs_dedup.html"),

        # Rule I
        fragmentSizePlot = expand(os.path.join("03_BAM_dedup/fragmentSizeDistribution_plots/", "{sample}_fragment_size_distribution.pdf"), sample=SAMPLENAMES),
        report_pdf = "03_BAM_dedup/fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf",

        # Rule J
        #scalingFactors_txt_result = "04_Normalization/scalingFactor/scalingFactor_results.txt",

        # Rule_normalization
        norm_bw = ancient(expand(os.path.join("04_Normalization/normalized_bigWigs/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_normalized_bs{binSize}.bw"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), binSize=str(BINS))),

        # Rule peakCalling
        norm_outputs = norm_outputs,

        # Rule counts_summary
        #temp_file_counts = ancient(expand(os.path.join(SUMMARYDIR, "{sample}_counts_summary.temp"), sample=SAMPLENAMES)),
        summary_file = os.path.join(SUMMARYDIR, "counts_summary.txt"),

        # Rule PCA and plotCorrelation
        correlation_outputs = correlation_outputs,

        # Rule Heatmap zScores peaks
        rawScores_hetamap_MACS3 = ancient(os.path.join(SUMMARYDIR, "".join(["Heatmap_on_log1p.rawScores_for_", PEAKCALLER, ".peaks_union_population.pdf"]))),

        # Rules gatk variant calling
        dedup_BAM_bsqr = dedup_BAM_bsqr,
        dedup_BAM_bsqr_index = dedup_BAM_bsqr_index,
        vcf = vcf,
        snp = snp,
        indels = indels

    #params:
    #    summary_file = str(os.path.join(SUMMARYDIR, "counts_summary.txt"))

    shell:
        # """
        # mkdir -p 04_Normalization/HMCan_output/temp_GENOME
        # rm -r 04_Normalization/HMCan_output/temp_GENOME
        #
        # uniq -u {params.summary_file} > summary_file.temp
        # (head -n 1 summary_file.temp && tail -n +2 summary_file.temp | sort -k 1) > {params.summary_file}
        # rm summary_file.temp
        #
        # mkdir -p 01b_SAM_tempFolder
        # rm -R 01b_SAM_tempFolder
        #
        # pdfcombine 03_BAM_dedup/fragmentSizeDistribution_plots/*_fragment_size_distribution.pdf -o 03_BAM_dedup/fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf -sf
        #
        # printf '\033[1;36mPipeline ended!\\n\033[0m'
        # """
        """
        mkdir -p 04_Normalization/HMCan_output/temp_GENOME
        rm -r 04_Normalization/HMCan_output/temp_GENOME

        mkdir -p 01b_SAM_tempFolder
        rm -R 01b_SAM_tempFolder

        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """
# ========================================================================================


# ----------------------------------------------------------------------------------------
# Perform the FastQC on raw fastq.gz
rule A_fastQC_raw:
    input:
        fastq_gz = ancient(os.path.join(config["runs_directory"], "".join(["{RUNS}", config['fastq_extension']])))
    output:
        html = os.path.join("01_fastQC_raw","{RUNS}_fastqc.html"),
        zip =  os.path.join("01_fastQC_raw","{RUNS}_fastqc.zip")
    params:
        build_fastqcDir = os.path.dirname("01_fastQC_raw/multiQC_raw/"),
        fastQC_raw_outdir = os.path.join(config["output_directory"], "01_fastQC_raw"),
        run = "{RUNS}",
        CPUs = config["fastQC_threads"]
    threads:
        config["fastQC_threads"]
    shell:
        """
        printf '\033[1;36m{params.run}: Performing fastQC on raw fastq...\\n\033[0m'

        mkdir -p {params.build_fastqcDir}
        fastqc -t {params.CPUs} --outdir {params.fastQC_raw_outdir} {input.fastq_gz}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for raw fastq reports
rule B_multiQC_raw:
    input:
        fastqc_zip = ancient(expand(os.path.join("01_fastQC_raw", "{run}_fastqc.zip"), run=RUNNAMES))
    output:
        multiqcReportRaw = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"
    params:
        fastqc_Raw_reports = os.path.join("01_fastQC_raw", "*.zip"),
        multiQC_raw_outdir = os.path.join(config["output_directory"], "01_fastQC_raw/multiQC_raw/")
    shell:
        """
        printf '\033[1;36mGenerating multiQC report for fatsq quality test...\\n\033[0m'
        multiqc -f --outdir {params.multiQC_raw_outdir} -n multiQC_report_fastqRaw.html {params.fastqc_Raw_reports}
        """
# ----------------------------------------------------------------------------------------

if not os.path.exists("".join([config["genome_fasta"], ".fai"])):
    # ----------------------------------------------------------------------------------------
    # Reads alignement
    rule Cextra_generate_genome_index:
        input:
            multiqcReportRaw = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"
        output:
            genome_fai = "".join([config["genome_fasta"], ".fai"])
        params:
            genome = config["genome_fasta"],
        threads:
            config["bwa_threads"]
        shell:
            """
            printf '\033[1;36mGenerating the genome index...\\n\033[0m'
            bwa index {params.genome}
            samtools faidx {params.genome}
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
# Reads alignement
rule C_bwa_align:
    input:
        multiqcReportRaw = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html",
        R1 = ancient(os.path.join(config["runs_directory"], "".join(["{SAMPLES}", config['runs_suffix'][0], config['fastq_extension']]))),
        R2 = ancient(os.path.join(config["runs_directory"], "".join(["{SAMPLES}", config['runs_suffix'][1], config['fastq_extension']]))),
        genome_fai = "".join([config["genome_fasta"], ".fai"])
    output:
        SAM = os.path.join("01b_SAM_tempFolder/", "{SAMPLES}.sam")
    params:
        build_SAM = "01b_SAM_tempFolder/log/",
        genome = config["genome_fasta"],
        sample = "{SAMPLES}",
        CPUs = config["bwa_threads"]
    threads:
        config["bwa_threads"]
    log:
        out = os.path.join("01b_SAM_tempFolder/log/{SAMPLES}_bwa-mem_.out")
    shell:
        """
        mkdir -p {params.build_SAM}
        printf '\033[1;36m{params.sample}: alignment of the reads...\\n\033[0m'
        bwa mem -t {params.CPUs} {params.genome} {input.R1} {input.R2} > {output.SAM} 2> {log.out}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# SAM filtering for mapping quality and BAM generation | BAM MT-reads removal
rule D_sam_to_bam:
    input:
        SAM = ancient(os.path.join("01b_SAM_tempFolder/", "{SAMPLES}.sam"))
    output:
        filtBAM = temp(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), ".bam"]))),
        filtBAM_sorted = temp(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted.bam"]))),
        filtBAM_sorted_index = temp(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted.bam.bai"]))),
        flagstat_on_unfiltered_BAM = os.path.join("02_BAM/flagstat/", "{SAMPLES}_flagstat_UNfiltered_bam.txt"),
        filtBAM_sorted_woMT = os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT.bam"])),
        filtBAM_sorted_woMT_index = os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT.bam.bai"])),
        flagstat_on_filtered_woMT_BAM = os.path.join("02_BAM/flagstat/", "{SAMPLES}_flagstat_filtered_bam_woMT.txt")
    params:
        build_BAM = os.path.dirname("02_BAM/flagstat/"),
        mapq_cutoff = str(config["mapQ_cutoff"]),
        sample = "{SAMPLES}",
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    shell:
        """
        mkdir -p {params.build_BAM}

        printf '\033[1;36m{params.sample}: filtering SAM and generate BAM...\\n\033[0m'
        samtools view -@ {params.CPUs} -b -q {params.mapq_cutoff} {input.SAM} -o {output.filtBAM}

        printf '\033[1;36m{params.sample}: sorting BAM...\\n\033[0m'
        samtools sort -@ {params.CPUs} {output.filtBAM} -o {output.filtBAM_sorted}
        samtools index -@ {params.CPUs} -b {output.filtBAM_sorted} {output.filtBAM_sorted_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from unfiltered BAM...\\n\033[0m'
        samtools flagstat {output.filtBAM_sorted} -@ {params.CPUs} > {output.flagstat_on_unfiltered_BAM}

        printf '\033[1;36m{params.sample}: Removing MT reads from BAM...\\n\033[0m'
        samtools idxstats {output.filtBAM_sorted} | cut -f 1 | grep -v M | xargs samtools view -@ {params.CPUs} -b {output.filtBAM_sorted} > {output.filtBAM_sorted_woMT}
        samtools index -@ {params.CPUs} -b {output.filtBAM_sorted_woMT} {output.filtBAM_sorted_woMT_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from BAM without MT-DNA...\\n\033[0m'
        samtools flagstat {output.filtBAM_sorted_woMT} -@ {params.CPUs} > {output.flagstat_on_filtered_woMT_BAM}

        rm {input.SAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# BAM duplicates removal and relative flagstat | BAM reads shifting
rule E_bam_deduplication:
    input:
        BAM = ancient(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT.bam"])))
    output:
        dedup_BAM = os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bam"])),
        dedup_BAM_index = os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bai"])),
        dedup_BAM_metrics = os.path.join("03_BAM_dedup/metrics", "{SAMPLES}_metrics_woMT_dedup_bam.txt"),
        dedup_BAM_flagstat = os.path.join("03_BAM_dedup/flagstat/", "{SAMPLES}_flagstat_filtered_bam_woMT_dedup.txt")
    params:
        build_BAMdedup_metrics = os.path.dirname("03_BAM_dedup/metrics/"),
        build_BAMdedup_unshifted = os.path.dirname("03_BAM_dedup/unshifted_bams/"),
        build_BAMdedup_flagstat = os.path.dirname("03_BAM_dedup/flagstat/"),
        build_BAMdedup_multiQC = os.path.dirname("03_BAM_dedup/fastQC/multiQC_dedup_bams/"),
        build_BAMdedup_fragmentPlots = os.path.dirname("03_BAM_dedup/fragmentSizeDistribution_plots/"),
        sample = "{SAMPLES}",
        minFragmentLength = config["minFragmentLength"],
        maxFragmentLength = config["maxFragmentLength"],
        rm_dup = config["remove_duplicates"],
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    shell:
        """
        mkdir -p {params.build_BAMdedup_metrics}
        mkdir -p {params.build_BAMdedup_unshifted}
        mkdir -p {params.build_BAMdedup_flagstat}
        mkdir -p {params.build_BAMdedup_multiQC}
        mkdir -p {params.build_BAMdedup_fragmentPlots}

        printf '\033[1;36m{params.sample}: Removing BAM duplicates...\\n\033[0m'
        picard MarkDuplicates CREATE_INDEX=true INPUT={input.BAM} OUTPUT={output.dedup_BAM} METRICS_FILE={output.dedup_BAM_metrics} ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES={params.rm_dup} VALIDATION_STRINGENCY=STRICT
        samtools flagstat {output.dedup_BAM} -@ {params.CPUs} > {output.dedup_BAM_flagstat}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# BAM reads shifting
rule F_bam_shifting:
    input:
        dedup_BAM = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bam"]))),
        dedup_BAM_index = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bai"])))
    output:
        dedup_BAM_shifted_toSort = temp(os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup_shifted.ToSort.bam"]))),
        dedup_BAM_shifted_sorted = os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam"])),
        dedup_BAM_shifted_sorted_index = os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam.bai"])),
        dedup_BAM_shifted_sorted_flagstat = os.path.join("03_BAM_dedup/flagstat/", "{SAMPLES}_flagstat_woMT_dedup_shifted_sorted.txt")
    params:
        sample = "{SAMPLES}",
        minFragmentLength = str(config["minFragmentLength"]),
        maxFragmentLength = str(config["maxFragmentLength"]),
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    shell:
        """
        printf '\033[1;36m{params.sample}: Shifting reads in BAM...\\n\033[0m'
        alignmentSieve -p {params.CPUs} --ATACshift --bam {input.dedup_BAM} --outFile {output.dedup_BAM_shifted_toSort} --minFragmentLength {params.minFragmentLength} --maxFragmentLength {params.maxFragmentLength}

        printf '\033[1;36m{params.sample}: Sorting shifted BAM...\\n\033[0m'
        samtools sort -@ {params.CPUs} {output.dedup_BAM_shifted_toSort} -o {output.dedup_BAM_shifted_sorted}
        samtools index -@ {params.CPUs} -b {output.dedup_BAM_shifted_sorted} {output.dedup_BAM_shifted_sorted_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from shifted BAM...\\n\033[0m'
        samtools flagstat {output.dedup_BAM_shifted_sorted} -@ {params.CPUs} > {output.dedup_BAM_shifted_sorted_flagstat}

        echo '------------------------------------------------------------------------'
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# FastQC on BAMs
rule G_fastQC_BAMs:
    input:
        dedup_BAM = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bam"]))),
        dedup_BAM_index = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bai"])))
    output:
        html = os.path.join("03_BAM_dedup/fastQC/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup_fastqc.html"])),
        zip = os.path.join("03_BAM_dedup/fastQC/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join(config["output_directory"], "03_BAM_dedup/fastQC/"),
        sample = "{SAMPLES}",
        CPUs = config["fastQC_threads"]
    threads:
        config["fastQC_threads"]
    shell:
        """
        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        fastqc -t {params.CPUs} --outdir {params.fastQC_BAMs_outdir} {input.dedup_BAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for BAMs
rule H_multiQC_BAMs:
    input:
        BAM_fastqc_zip = ancient(expand(os.path.join("03_BAM_dedup/fastQC/", "{sample}_mapQ{MAPQ}_sorted_woMT_dedup_fastqc.zip"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"])))
    output:
        multiqcReportBAM = "03_BAM_dedup/fastQC/multiQC_dedup_bams/multiQC_report_BAMs_dedup.html"
    params:
        fastQC_BAM_reports = os.path.join("03_BAM_dedup/fastQC/", "*.zip"),
        multiQC_BAM_outdir = os.path.join(config["output_directory"], "03_BAM_dedup/fastQC/multiQC_dedup_bams/")
    shell:
        """
        printf '\033[1;36mGenerating multiQC report from deduplicated bam quality test...\\n\033[0m'
        multiqc -f --outdir {params.multiQC_BAM_outdir} -n multiQC_report_BAMs_dedup.html {params.fastQC_BAM_reports}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule I1_fragment_size_distribution:
    input:
        BAM = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bam"]))),
        BAM_index = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bai"]))),
        multiqcReportBAM = "03_BAM_dedup/fastQC/multiQC_dedup_bams/multiQC_report_BAMs_dedup.html"
    output:
        plot = os.path.join("03_BAM_dedup/fragmentSizeDistribution_plots/", "{SAMPLES}_fragment_size_distribution.pdf")
    params:
        sample = "{SAMPLES}",
        plotFormat = config["plot_format"],
        binSize = str(config["window_length"]),
        blacklist = config["blacklist_file"],
        CPUs = config["threads_bamPEFragmentSize"],
        max_fragment_length = config["max_fragment_length"]
    threads:
        config["threads_bamPEFragmentSize"]
    shell:
        """
        printf '\033[1;36m{params.sample}: Plotting the fragment size distribution...\\n\033[0m'

        bamPEFragmentSize \
        -p {params.CPUs} \
        -b {input.BAM} \
        --plotFileFormat {params.plotFormat} \
        --plotTitle {params.sample} \
        --samplesLabel {params.sample} \
        --binSize {params.binSize} \
        --maxFragmentLength {params.max_fragment_length} \
        --blackListFileName {params.blacklist} \
        -o {output.plot}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule I2_fragment_size_distribution_report:
    input:
        plots = ancient(expand(os.path.join("03_BAM_dedup/fragmentSizeDistribution_plots/", "{sample}_fragment_size_distribution.pdf"), sample = SAMPLENAMES))
    output:
        report_pdf ="03_BAM_dedup/fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"
    threads:
        1
    shell:
        """
        printf '\033[1;36mMerging fragmentSizeDistribution reports in a unique PDF...\\n\033[0m'
        pdfcombine 03_BAM_dedup/fragmentSizeDistribution_plots/*_fragment_size_distribution.pdf -o {output.report_pdf} -sf
        """


# ----------------------------------------------------------------------------------------

# ****************************************************************************************
# START OF CONDITONAL NORMALIZATION SECTION: HMCan (1) vs Sequencing deepth (2)
# ****************************************************************************************
# Run HMCan normalization if required --> perform_HMCan_correction: "False"
# ****************************************************************************************
if (eval(str(config["perform_HMCan_correction"])) == True):

    #  Generate config file for HMCcan
    rule J1_HMCan_config_file_and_genome_splitting:
        input:
            dir = ancient(os.path.dirname("04_Normalization/normalized_bigWigs/"))
        output:
            HMCan_config = "04_Normalization/HMCan_output/CONFIGURATION_file_HMCan.txt"
        params:
            # multiBamSummary bins
            labels = ' '.join(SAMPLENAMES),
            threads = str(config["multiBamSummary_threads"]),
            blacklist = config["blacklist_file"],
            minFragmentLength = str(config["minFragmentLength"]),
            maxFragmentLength = str(config["maxFragmentLength"]),
            binSize = str(config["window_length"]),
            # faidx
            genome = config["genome_fasta"],
            working_dir = config["output_directory"],
            # HMCan config
            format = config["format"],
            GCIndex = config["GCIndex"],
            smallBinLength = config["smallBinLength"],
            largeBinLength = config["largeBinLength"],
            genomePath = os.path.join(config["output_directory"], "04_Normalization/HMCan_output/temp_GENOME/"),
            pvalueThreshold = config["pvalueThreshold"],
            mergeDistance = config["mergeDistance"],
            iterationThreshold = config["iterationThreshold"],
            finalThreshold = config["finalThreshold"],
            maxIter = config["maxIter"],
            PrintWig = config["PrintWig"],
            blackListFile = config["blacklist_file"],
            PrintPosterior = config["PrintPosterior"],
            PrintBedGraph = config["PrintBedGraph"],
            CallPeaks = config["CallPeaks"],
            pairedEnds = config["pairedEnds"],
            Type = config["Type"],
            GCmergeDistance = config["GCmergeDistance"],
            RemoveDuplicates = config["RemoveDuplicates"],
            CNAnormalization = config["CNAnormalization"]
        shell:
            """
            printf '\033[1;36mSplitting the genome in single chromosomes files for HMCan normalization...\\n\033[0m'

            mkdir -p {params.genomePath}
            cd {params.genomePath}
            faidx {params.genome} --split-files
            cd {params.working_dir}



            printf '\033[1;36mReading the HMcan configuration information...\\n\033[0m'

            mkdir -p 04_Normalization/HMCan_output/

            echo format {params.format} > {output.HMCan_config}
            echo GCIndex {params.GCIndex} >> {output.HMCan_config}
            echo smallBinLength {params.smallBinLength} >> {output.HMCan_config}
            echo largeBinLength {params.largeBinLength} >> {output.HMCan_config}
            echo genomePath {params.genomePath} >> {output.HMCan_config}
            echo pvalueThreshold {params.pvalueThreshold} >> {output.HMCan_config}
            echo mergeDistance {params.mergeDistance} >> {output.HMCan_config}
            echo iterationThreshold {params.iterationThreshold} >> {output.HMCan_config}
            echo finalThreshold {params.finalThreshold} >> {output.HMCan_config}
            echo maxIter {params.maxIter} >> {output.HMCan_config}
            echo PrintWig {params.PrintWig} >> {output.HMCan_config}
            echo blackListFile {params.blackListFile} >> {output.HMCan_config}
            echo PrintPosterior {params.PrintPosterior} >> {output.HMCan_config}
            echo PrintBedGraph {params.PrintBedGraph} >> {output.HMCan_config}
            echo CallPeaks {params.CallPeaks} >> {output.HMCan_config}
            echo pairedEnds {params.pairedEnds} >> {output.HMCan_config}
            echo Type {params.Type} >> {output.HMCan_config}
            echo GCmergeDistance {params.GCmergeDistance} >> {output.HMCan_config}
            echo RemoveDuplicates {params.RemoveDuplicates} >> {output.HMCan_config}
            echo CNAnormalization {params.CNAnormalization} >> {output.HMCan_config}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # CNV correction by HMCan
    rule J2_signal_correction_for_CNVs_HMCan:
        input:
            dedup_BAM_shifted_sorted = ancient(os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam"]))),
            dedup_BAM_shifted_sorted_index = ancient(os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam.bai"]))),
            config_file = "04_Normalization/HMCan_output/CONFIGURATION_file_HMCan.txt",
        output:
            HMCan_regions = os.path.join("04_Normalization/HMCan_output/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_regions.bed"])),
            HMCan_peaks = os.path.join("04_Normalization/HMCan_output/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_peaks.narrowPeak"])),
            HMCan_CNV_profile = os.path.join("04_Normalization/HMCan_output/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_CNV_profile.txt"])),
            HMCan_bedGraph = temp(os.path.join("04_Normalization/HMCan_output/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted.bedGraph"])))
        params:
            HMCan_path = config["HMCan_path"],
            basename = os.path.join("04_Normalization/HMCan_output/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted"])),
            sample = "{SAMPLES}"
        threads:
            config["HMCan_threads"]
        resources:
            mem_mb = 50
        shell:
            """
            printf '\033[1;36m{params.sample}: HMcan signal correction...\\n\033[0m'
            {params.HMCan_path} {input.dedup_BAM_shifted_sorted} - {input.config_file} {params.basename}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Computing the scaling factor
    rule J3_computing_scaling_factor:
        input:
            dedup_BAM_shifted_sorted = ancient(expand(os.path.join("03_BAM_dedup/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_sorted.bam"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
            dedup_BAM_shifted_sorted_index = ancient(expand(os.path.join("03_BAM_dedup/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_sorted.bam.bai"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
            fragment_size_report_pdf ="03_BAM_dedup/fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"
        output:
            npz_results = temp("04_Normalization/scalingFactor/scalingFactor_results.npz"),
            txt_result_notNorm = temp("04_Normalization/scalingFactor/scalingFactor_results_to_normalize.txt")
        params:
            # multiBamSummary bins
            build_normalization = os.path.dirname("04_Normalization/scalingFactor/"),
            build_normBigWig = os.path.dirname("04_Normalization/normalized_bigWigs/"),
            labels = ' '.join(SAMPLENAMES),
            blacklist = config["blacklist_file"],
            minFragmentLength = str(config["minFragmentLength"]),
            maxFragmentLength = str(config["maxFragmentLength"]),
            binSize = str(config["window_length"]),
            CPUs = config["multiBamSummary_threads"]
        threads:
            config["multiBamSummary_threads"]
        shell:
            """
            mkdir -p {params.build_normalization}
            mkdir -p {params.build_normBigWig}

            printf '\033[1;36mComputing the scaling factors for the intersample normalization...\\n\033[0m'

            multiBamSummary bins \
            -p {params.CPUs} \
            --bamfiles {input.dedup_BAM_shifted_sorted} \
            -o {output.npz_results} \
            --blackListFileName {params.blacklist} \
            --scalingFactors {output.txt_result_notNorm} \
            --minFragmentLength {params.minFragmentLength} \
            --maxFragmentLength {params.maxFragmentLength} \
            --binSize {params.binSize} \
            --labels {params.labels}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # Normalization scaling factor to reference sample
    rule J4_computing_scaling_factor_normalizationToReferenceSample:
        input:
            dedup_BAM_shifted_sorted = ancient(expand(os.path.join("03_BAM_dedup/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_sorted.bam"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
            txt_result_notNorm = ancient("04_Normalization/scalingFactor/scalingFactor_results_to_normalize.txt")
        output:
            txt_result = "04_Normalization/scalingFactor/scalingFactor_results.txt"
        params:
            ref_sample = config["reference_sample"],
            samples_names = SAMPLENAMES
        run:
            if params.ref_sample in SAMPLENAMES:
                reference_sample = params.ref_sample
            else:
                reference_sample = params.samples_names[0]

            shell("printf '\033[1;36mThe scaling factors are normalized to {params.ref_sample}\\n\033[0m'")

            import pandas as pd
            tb = pd.read_table(input.txt_result_notNorm)
            reference_factor = tb.loc[tb["sample"] == reference_sample, "scalingFactor"].values[0]
            tb["scalingFactor"] = tb["scalingFactor"].div(reference_factor)
            tb.to_csv(output.txt_result, encoding="utf-8", index=False, header=True, sep="\t")
    # ----------------------------------------------------------------------------------------



    # ----------------------------------------------------------------------------------------
    # Normalized bigWig generation
    rule J5_bigWig_CNV_adjusted_signal:
        input:
            HMCan_bedGraph = ancient(os.path.join("04_Normalization/HMCan_output/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted.bedGraph"]))),
            scaling_factors = ancient("04_Normalization/scalingFactor/scalingFactor_results.txt")
        output:
            norm_bdg = temp(os.path.join("04_Normalization/HMCan_output/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_NORMALIZED.bedGraph"]))),
            norm_bw = os.path.join("04_Normalization/normalized_bigWigs/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_normalized_bs", str(config["smallBinLength"]), ".bw"]))
        params:
            sample = "{SAMPLES}",
            chrSizes = config["chromosome_sizes_file"]
        threads:
            config["HMCan_threads"]
        run:
            shell("printf '\033[1;36m{params.sample}: Normalization of the corrected signal and generation of the bigWig file...\\n\033[0m'")

            # Import bedGraph and scaling factor
            import pandas as pd
            bdg = pd.read_csv(str(input.HMCan_bedGraph),  sep='\s+', engine='python', header=None)
            factor = pd.read_csv(str(input.scaling_factors),  sep='\s+', engine='python')
            factor = factor[factor['sample']==params.sample].iloc[0,1]

            # Normalize the bedGraph scores
            bdg[3] = factor * bdg[3]

            # Export normalized_bedGraph
            bdg.to_csv(output.norm_bdg, encoding='utf-8', index=False, header=False, sep='\t')

            # Converting normalized_bedGraph to BigWig
            shell("bedGraphToBigWig {output.norm_bgd} {params.chrSizes} {output.norm_bw}")
    # ----------------------------------------------------------------------------------------



# ********************************************************************************************
else: # Skip HMCan and perform a classical "sequencing-depth" normalization
# ********************************************************************************************
    # bigWig generation from BAM (not corrected by HMCan)
    rule J1_bigWig_normalization_woCNVcorrection:
        input:
            dedup_BAM_shifted_sorted = ancient(os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam"]))),
            dedup_BAM_shifted_sorted_index = ancient(os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam.bai"])))
        output:
            norm_bw = os.path.join("04_Normalization/normalized_bigWigs/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_normalized_bs", str(config["bigWig_binSize"]), ".bw"]))
        params:
            sample = "{SAMPLES}",
            binSize = config["bigWig_binSize"],
            normalization_method = config["normalization_method"],
            effective_genomeSize = config["effective_genomeSize"],
            ignore_for_normalization = config["ignore_for_normalization"],
            blacklist = config["blacklist_file"],
            CPUs = config["bamCoverage_threads"]
        threads:
            config["bamCoverage_threads"]
        shell:
            """
            printf "\033[1;36m{params.sample}: generation of the normalized bigWig file...\\n\033[0m"

            bamCoverage -p {params.CPUs} \
            --bam {input.dedup_BAM_shifted_sorted} \
            --binSize {params.binSize} \
            --normalizeUsing {params.normalization_method} \
            --effectiveGenomeSize {params.effective_genomeSize} \
            --ignoreForNormalization {params.ignore_for_normalization} \
            --extendReads \
            --blackListFileName {params.blacklist} \
            -of "bigwig" \
            -o {output.norm_bw}
            """
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # MACS3 peakCalling on uncorrected bams
    rule J2_MACS3_peakCalling:
        input:
            dedup_BAM_shifted_sorted = ancient(os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam"]))),
            dedup_BAM_shifted_sorted_index = ancient(os.path.join("03_BAM_dedup/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_sorted.bam.bai"]))),
            norm_bw = ancient(expand(os.path.join("04_Normalization/normalized_bigWigs/", ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES))
        output:
            peaks_xls = os.path.join("05_Peaks_MACS3/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_FDR", str(config["FDR_cutoff"]), "_peaks.xls"])),
            narrowPeaks = os.path.join("05_Peaks_MACS3/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_FDR", str(config["FDR_cutoff"]), "_peaks.narrowPeak"])),
            summits = os.path.join("05_Peaks_MACS3/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_FDR", str(config["FDR_cutoff"]), "_summits.bed"]))
        params:
            genomeSize = str(config["genome_size_MACS"]),
            FDR = str(config["FDR_cutoff"]),
            basename = ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_FDR", str(config["FDR_cutoff"])]),
            summits = SUMMITS,
            sample = "{SAMPLES}"
        log:
            log = os.path.join("05_Peaks_MACS3/log/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_FDR", str(config["FDR_cutoff"]), ".log"]))
        shell:
            """
            printf '\033[1;36m{params.sample}: Calling peaks by MACS3...\\n\033[0m'

            mkdir -p 05_Peaks_MACS3/log/

            macs3 callpeak \
            -t {input.dedup_BAM_shifted_sorted} \
            -g {params.genomeSize} \
	        -n {params.basename} \
            -q {params.FDR} \
            -f BAMPE \
	        --outdir 05_Peaks_MACS3 \
            --keep-dup all \
            --nolambda \
            {params.summits} 2> {log.log}
            """
    # ----------------------------------------------------------------------------------------
# ****************************************************************************************
# END OF NORMALIZATION PART (if else condition)
# ****************************************************************************************


# ----------------------------------------------------------------------------------------
# Computation of the counts summary table
rule K_counts_summary:
    input:
        R1 = ancient(expand(os.path.join(config["runs_directory"], "".join(["{sample}", config['runs_suffix'][0], config['fastq_extension']])), sample = SAMPLENAMES)),
        R2 = ancient(expand(os.path.join(config["runs_directory"], "".join(["{sample}", config['runs_suffix'][1], config['fastq_extension']])), sample = SAMPLENAMES)),
        flagstat_on_unfiltered_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_UNfiltered_bam.txt"), sample = SAMPLENAMES)),
        flagstat_on_filtered_woMT_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_filtered_bam_woMT.txt"), sample = SAMPLENAMES)),
        dedup_BAM_flagstat = ancient(expand(os.path.join("03_BAM_dedup/flagstat/", "{sample}_flagstat_filtered_bam_woMT_dedup.txt"), sample = SAMPLENAMES)),
        dedup_BAM_shifted_sorted_flagstat = ancient(expand(os.path.join("03_BAM_dedup/flagstat/", "{sample}_flagstat_woMT_dedup_shifted_sorted.txt"), sample = SAMPLENAMES)),
        #scaling_factors = ancient("04_Normalization/scalingFactor/scalingFactor_results.txt"),
        peaks_file = ancient(expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapQ", str(config["mapQ_cutoff"]), "_woMT_dedup_shifted_FDR", str(config["FDR_cutoff"]), "_peaks.narrowPeak"])), sample = SAMPLENAMES))
    output:
        summary_file = os.path.join(SUMMARYDIR, "counts_summary.txt")
    params:
        build_summary_directory = os.path.dirname(SUMMARYDIR),
        R1_suffix = config['runs_suffix'][0],
        R2_suffix = config['runs_suffix'][1],
        sample_list = str(' '.join(SAMPLENAMES)),
        multiQC_report = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw_data/multiqc_general_stats.txt",
        peaks_dir = PEAKSDIR
    threads: 1
    shell:
        """
        mkdir -p {params.build_summary_directory}

        printf '\033[1;36mGeneration of a general counts summary table...\\n\033[0m'
        printf Sample'\\t'Reads_R1'\\t'Reads_R2'\\t'Reads_total'\\t'unfiltered_BAM'\\t'Percentage_MT'\\t'dedup_BAM'\\t'duplicated_reads'\\t'shifted_BAM'\\t'loss_post_shifting'\\t'n.peaks'\\n' > {output.summary_file}

        for NAME in {params.sample_list}
        do
            printf '\033[1;36m     - %s: adding stats to summary table...\\n\033[0m' $NAME

            R1=$(grep ${{NAME}}{params.R1_suffix} {params.multiQC_report} | cut -f 6 | sed 's/\.[0-9]\+//')
            R2=$(grep ${{NAME}}{params.R2_suffix} {params.multiQC_report} | cut -f 6 | sed 's/\.[0-9]\+//')
            TOTAL=$((R1 + R2))

            unfilteredBAM=$(grep mapped 02_BAM/flagstat/${{NAME}}_flagstat_UNfiltered_bam.txt | head -n 1 | cut -f 1 -d ' ')
            woMT_BAM=$(grep mapped 02_BAM/flagstat/${{NAME}}_flagstat_filtered_bam_woMT.txt | head -n 1 | cut -f 1 -d ' ')
            percMT=$(echo "scale=1; (100 - (($woMT_BAM/$unfilteredBAM) * 100))" | bc)

            dedupBAM=$(grep mapped 03_BAM_dedup/flagstat/${{NAME}}_flagstat_filtered_bam_woMT_dedup.txt | head -n 1 | cut -f 1 -d ' ')
            dedupREADS=$((woMT_BAM - dedupBAM))

            shiftedBAM=$(grep mapped 03_BAM_dedup/flagstat/${{NAME}}_flagstat_woMT_dedup_shifted_sorted.txt | head -n 1 | cut -f 1 -d ' ')
            lossReads=$((dedupBAM - shiftedBAM))

            peaks=$(wc -l {params.peaks_dir}${{NAME}}*.*Peak | cut -f 1 -d ' ')

            printf ${{NAME}}'\\t'$R1'\\t'$R2'\\t'$TOTAL'\\t'$unfilteredBAM'\\t'$percMT'\\t'$dedupBAM'\\t'$dedupREADS'\\t'$shiftedBAM'\\t'$lossReads'\\t'$peaks'\\n' >> {output.summary_file}
        done

        uniq -u {output.summary_file} > summary_file.temp
        (head -n 1 summary_file.temp && tail -n +2 summary_file.temp | sort -k 1) > {output.summary_file}
        rm summary_file.temp
        """


# ----------------------------------------------------------------------------------------
# Generation of samples PCA and Heatmap
rule L1_multiBigwigSummary:
    input:
        norm_bw = ancient(expand(os.path.join("04_Normalization/normalized_bigWigs/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_normalized_bs{binSize}.bw"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), binSize=str(BINS)))
    output:
        matrix = temp(os.path.join(SUMMARYDIR, "temp_multiBigWigSummary_matrix.npz"))
    params:
        labels = ' '.join(SAMPLENAMES),
        window = config["binning_window_size"],
        blacklist = config["blacklist_file"],
        CPUs = config["threads_multiBigwigSummary"]
    threads:
        config["threads_multiBigwigSummary"]
    shell:
        """
        printf '\033[1;36mComparing the whole signal among samples...\\n\033[0m'

        multiBigwigSummary bins -p {params.CPUs} -b {input.norm_bw} --labels {params.labels} --binSize {params.window} --blackListFileName {params.blacklist} -o {output.matrix}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of samples PCA and Heatmap
rule L2_PCA_and_samples_correlation:
    input:
        matrix = os.path.join(SUMMARYDIR, "temp_multiBigWigSummary_matrix.npz")
    output:
        PCA = os.path.join(SUMMARYDIR, "PCA_on_BigWigs_wholeGenome.pdf"),
        hetamap_spearman = os.path.join(SUMMARYDIR, "Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        hetamap_pearson = os.path.join(SUMMARYDIR, "Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"),
        scatterplot_spearman = os.path.join(SUMMARYDIR, "Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        scatterplot_pearson = os.path.join(SUMMARYDIR, "Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf")
    params:
        heatmap_color = config["heatmap_color"]
    shell:
        """
        printf '\033[1;36mPlotting the correlation and variability of the whole signal among samples...\\n\033[0m'

        printf '\033[1;36m    - plotting PCA...\\n\033[0m'
        plotPCA -in {input.matrix} -o {output.PCA} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf'


        printf '\033[1;36m    - plotting Spearman correlation heatmap...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.hetamap_spearman}

        printf '\033[1;36m    - plotting Pearson correlation heatmap...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.hetamap_pearson}



        printf '\033[1;36m    - plotting Spearman correlation scatterplot...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_spearman}

        printf '\033[1;36m    - plotting Pearson correlation scatterplot...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_pearson}
        """


# ----------------------------------------------------------------------------------------
# Absolute peaks file and relative matrix score generation for MACS3/HMCan peaks
rule M_all_peaks_file_and_score_matrix:
    input:
        norm_bw = ancient(expand(os.path.join("04_Normalization/normalized_bigWigs/", "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_normalized_bs{binSize}.bw"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), binSize=str(BINS))),
        peaks_file = (expand(os.path.join(str(PEAKSDIR), "{sample}_mapQ{MAPQ}_woMT_dedup_shifted_FDR{fdr}_peaks.narrowPeak"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), fdr=str(config["FDR_cutoff"])))
    output:
        concatenation_bed = temp(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation.bed")),
        concatenation_bed_sorted = temp(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_sorted.bed")),
        concatenation_bed_collapsed = temp(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed.bed")),
        concatenation_bed_collapsed_sorted = temp(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        score_matrix_peaks = temp(os.path.join(SUMMARYDIR, "temp_peaks_score_matrix_all_samples_MACS3.npz")),
        score_matrix_peaks_table = temp(os.path.join(SUMMARYDIR, ''.join(["temp_peaks_score_matrix_all_samples_table_", PEAKCALLER, ".tsv"])))
    params:
        peaks_dir = PEAKSDIR,
        labels = ' '.join(SAMPLENAMES),
        blacklist = config["blacklist_file"],
        CPUs = config["threads_multiBigwigSummary"]
    threads:
        config["threads_multiBigwigSummary"]
    shell:
        """
        printf '\033[1;36mGenerating a file result of the merge of all the MACS3 peaks...\\n\033[0m'
        cat {params.peaks_dir}*.*Peak >> {output.concatenation_bed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed} > {output.concatenation_bed_sorted}

        bedtools merge -i {output.concatenation_bed_sorted} | uniq > {output.concatenation_bed_collapsed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed_collapsed} > {output.concatenation_bed_collapsed_sorted}


        printf '\033[1;36mComputing the score matrix for all the MACS3 peaks per each sample...\\n\033[0m'

        multiBigwigSummary BED-file \
        -p {params.CPUs} \
        -b {input.norm_bw} \
        -o {output.score_matrix_peaks} \
        --BED {output.concatenation_bed_collapsed_sorted} \
        --blackListFileName {params.blacklist} \
        --outRawCounts {output.score_matrix_peaks_table} \
        --labels {params.labels}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Compute peaks z-scores and plot heatmap for MACS3/HMCan peaks
rule N_peaks_zScores_and_heatmap:
    input:
        score_matrix_peaks_table = ancient(os.path.join(SUMMARYDIR, "temp_peaks_score_matrix_all_samples_table_{PEAKCALLER}.tsv"))
    output:
        rawScores_hetamap = os.path.join(SUMMARYDIR, "Heatmap_on_log1p.rawScores_for_{PEAKCALLER}.peaks_union_population.pdf")
    params:
        heatmap_color = config["heatmap_color"],
        zScore_heatmap_color = config["zScore_heatmap_color"],
        heatmap_basename_rawScores = os.path.join(SUMMARYDIR, "Heatmap_on_log1p.rawScores_for_{PEAKCALLER}.peaks_union_population"),
        heatmap_basename_zScore = os.path.join(SUMMARYDIR, "Heatmap_on_zScores_for_{PEAKCALLER}.peaks_union_population"),
        n_samples = len(SAMPLENAMES),
        peak_caller = PEAKCALLER
    threads:
        config["threads_multiBigwigSummary"]
    run:
        # Messege
        shell("printf '\033[1;36mLoading {PEAKCALLER} peak raw score table...\\n\033[0m'")

        # Import multiBigWig summary table
        import pandas as pd
        matrix = pd.read_csv(str(input.score_matrix_peaks_table),  sep='\s+', engine='python')

        # Use peak coordinates as ID ofr each row
        matrix["peak_ID"] = matrix[matrix.columns[:3]].apply(lambda x: '_'.join(x.dropna().astype(str)),axis=1)
        matrix = matrix[matrix.columns[3:]]
        matrix = matrix.set_index('peak_ID')

        # Conversion in log1p
        import numpy as np
        matrix_log1p = np.log1p(matrix)

        # Required to avoid the use of X11
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        # Required to avoid errors in scipy iterations
        import sys
        sys.setrecursionlimit(100000)


        # Messege
        shell("printf '\033[1;36mPlotting the {PEAKCALLER} peak raw score hetamaps...\\n\033[0m'")

        # Generation of the rawScore heatmap and clustering
        from bioinfokit import analys, visuz
        visuz.gene_exp.hmap(df=matrix_log1p,
                            cmap=params.heatmap_color,
                            rowclus=True,
                            colclus=(params.n_samples > 1),
                            figtype="pdf",
                            ylabel=False,
                            figname=str(params.heatmap_basename_rawScores),
                            dim=(params.n_samples, 9), # W * H
                            tickfont=(6, 4))


        ## Generation of the zScore heatmap and clustering
        #from bioinfokit import analys, visuz
        #visuz.gene_exp.hmap(df=matrix,
                            #cmap=params.zScore_heatmap_color,
                            #rowclus=True,
                            #colclus=True,
                            #zscore=0,
                            #figtype ="pdf",
                            #ylabel=False,
                            #figname=str(params.heatmap_basename_zScore),
                            #dim=(6, 12),
                            #tickfont=(6, 4))

        # ---------------- Manual way to compute zScore heatmap Z=(rowScore - rowMean)/rowSD -----------------
        if params.n_samples > 1:
            # Messege
            shell("printf '\033[1;36mComputing the zScores for {PEAKCALLER} peaks...\\n\033[0m'")

            import pandas as pd
            matrix = pd.read_csv(str(input.score_matrix_peaks_table),  sep='\s+', engine='python')
            stat_tb = pd.DataFrame({'rowMeans': matrix[matrix.columns[3:]].mean(axis=1),
                                    'SD':  matrix[matrix.columns[3:]].std(axis=1)})

            scores = []
            for i in list(range(3,len(matrix.columns))):
                scores.append((matrix[matrix.columns[i]] - stat_tb["rowMeans"]) / stat_tb["SD"])


            zScores = pd.DataFrame(scores).transpose()
            zScores.columns = list(matrix.columns)[3:len(matrix.columns)]
            zScores = zScores.fillna(0)
            zScores['peak_ID'] = matrix[matrix.columns[:3]].apply(lambda x: '_'.join(x.dropna().astype(str)),axis=1)
            zScores = zScores.set_index('peak_ID')


            # Messege
            shell("printf '\033[1;36mPlotting the {PEAKCALLER} peak zScores hetamaps...\\n\033[0m'")

            visuz.gene_exp.hmap(df=zScores,
                                cmap=params.zScore_heatmap_color,
                                rowclus=True,
                                colclus=True,
                                figtype ="pdf",
                                ylabel=False,
                                figname=str(params.heatmap_basename_zScore),
                                dim=(params.n_samples, 9),
                                tickfont=(6, 4))
        #------- ------- ------- ------- ------- ------- ------- ------- ------- -------
# ----------------------------------------------------------------------------------------



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  VARIANT CALLING  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ----------------------------------------------------------------------------------------
# Create reference genome dictionary
rule O_Make_reference_genome_dictionary:
    input:
        genome = config["genome_fasta"]
    output:
        genome_dict = ''.join([re.sub("[.]([a-z]|[A-Z])*$", "",config["genome_fasta"]),'.dict'])
    shell:
        """
        picard CreateSequenceDictionary REFERENCE={input.genome} OUTPUT={output.genome_dict}
        """
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# run base score recalibration (BSQR) of the bams
rule P_GATK_bam_base_quality_score_recalibration:
    input:
        genome_dict = ancient(''.join([re.sub("[.]([a-z]|[A-Z])*$", "",config["genome_fasta"]),'.dict'])),
        dedup_BAM = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bam"]))),
        dedup_BAM_index = ancient(os.path.join("03_BAM_dedup/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_dedup.bai"])))
    params:
        sample = "{SAMPLES}",
        gatk_directory = GATKDIR,
        genome = config["genome_fasta"],
        dbsnp = config["dbsnp_file"],
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    output:
        bsqr_table = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.table"),
        dedup_BAM_bsqr = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam"),
        dedup_BAM_bsqr_index = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam.bai")
    shell:
        """
        printf '\033[1;36m{params.sample}: Base Quality Score Recalibration of the deduplicated unshifted bam...\\n\033[0m'
        mkdir -p {params.gatk_directory}{params.sample}

        gatk --java-options '-Xmx4G' BaseRecalibrator \
        --input {input.dedup_BAM} \
        --known-sites {params.dbsnp} \
        --output {output.bsqr_table} \
        --reference {params.genome}

        gatk --java-options '-Xmx4G' ApplyBQSR \
        -R {params.genome} \
        -I {input.dedup_BAM} \
        --bqsr-recal-file {output.bsqr_table} \
        -O {output.dedup_BAM_bsqr}

        printf '\033[1;36m{params.sample}: Indexing recalibrated bam...\\n\033[0m'
        samtools index -@ {params.CPUs} -b {output.dedup_BAM_bsqr} {output.dedup_BAM_bsqr_index}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# run gatk haplotype caller
rule Q_GATK_haplotype_calling:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        dedup_BAM_bsqr = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam")),
        dedup_BAM_bsqr_index = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam.bai"))
    params:
        genome = config["genome_fasta"],
        sample = "{SAMPLES}",
        to_copy_bed = os.path.join(GATKDIR, "all_samples_peaks_concatenation_collapsed_sorted.bed")
    output:
        gvcf = temp(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.g.vcf.gz"))
    shell:
        """
        cp {input.concatenation_bed_collapsed_sorted} {params.to_copy_bed}

        printf '\033[1;36m{params.sample}: GATK Haplotype calling...\\n\033[0m'

        gatk --java-options '-Xmx4g' HaplotypeCaller \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -I {input.dedup_BAM_bsqr} \
        -O {output.gvcf} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# correct the genotypes that come out of haplotype caller
rule R_GATK_haplotype_calling_correction:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        gvcf = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.g.vcf.gz"))
    params:
        sample = "{SAMPLES}",
        genome = config["genome_fasta"]
    output:
        vcf = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.vcf.gz")
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Haplotype call correction...\\n\033[0m'

        gatk --java-options '-Xmx4g' GenotypeGVCFs \
        --include-non-variant-sites \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.gvcf} \
        -O {output.vcf}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Select SNPs
rule S_GATK_call_SNPs:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        vcf = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.vcf.gz"))
    params:
        sample = "{SAMPLES}",
        genome = config["genome_fasta"]
    output:
        snp = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk-snp.vcf")
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK SNP calling...\\n\033[0m'

        gatk --java-options "-Xmx4g" SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type SNP \
        --select-type NO_VARIATION \
        --select-type-to-exclude INDEL \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.snp}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Call indels
rule T_GATK_call_indels:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        vcf = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.vcf.gz"))
    params:
        sample = "{SAMPLES}",
        genome = config["genome_fasta"]
    output:
        indels = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk-indel.vcf")
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Indel calling...\\n\033[0m'

        gatk --java-options "-Xmx4g" SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type INDEL \
        --select-type NO_VARIATION \
        --select-type-to-exclude SNP \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.indels}
        """
# ----------------------------------------------------------------------------------------





#    gatk-indel:
# cols:
#    info: [QD,FS,MQ,AC,ExcessHet]
#    format: [DP,GQ]


#snp_filter: ['gatk-snp~DP>10']
#indel_filter: ['gatk-indel~DP>10']
# ----------------------------------------------------------------------------------------
