import os
#conda_prefix = str(os.environ["CONDA_PREFIX"])

import sys
#sys.path.insert(1, conda_prefix+"/lib/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+"/site-packages")

from typing import List
import pathlib
import re
import math
import numpy
import pandas as pd

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
home_dir = os.path.join(config["workflow_configuration"]["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir


# get the unique samples names and other variables
if not (os.path.exists(config["workflow_configuration"]["runs_directory"])):
    os.system("printf '\033[1;31m\\n!!! *runs_directory* does not exist !!!\\n\\n\033[0m'")
else:
    BAMS = []
    for file in os.listdir(config["workflow_configuration"]["runs_directory"]):
        if file.endswith(config['bam_features']['bam_suffix']):
            BAMS.append(file)

    SAMPLENAMES = numpy.unique([re.sub(rf"{config['bam_features']['bam_suffix']}$", "", i) for i in BAMS])



#BINS = config["bigWig_binSize"]
PEAKCALLER = "MACS3"
PEAKSDIR = "04_"+PEAKCALLER+"_peaks/"
FOOTDIR = "05_Footprinting_TOBIAS/"
SUMMARYDIR = "06_Overall_quality_and_info/"
GATKDIR = "07a_Variant_calling/"
COPYWRITERDIR = "07b_CNV_detection/"


MAPQ = str(config["bam_features"]["MAPQ_threshold"])
genome_fasta = str(config["genomic_annotations"]["genome_fasta"])
genome_id = config["genomic_annotations"]["genome_id"]
BLACKLIST = str(config["genomic_annotations"]["blacklist"])

rm_duplicates = str(config["bam_features"]["remove_duplicates"]).lower()
if (eval(str(rm_duplicates.capitalize())) == True):
    DUP="dedup"
else:
    DUP="mdup"


if (eval(str(config["peak_calling"]["call_summits"])) == True):
    SUMMITS="--call-summits"
else:
    SUMMITS=""



# generation of global wildcard_constraints
wildcard_constraints:
    RUNS=constraint_to(BAMS),
    SAMPLES=constraint_to(SAMPLENAMES)

# Correlations and heatmaps outputs
correlation_outputs = []
if (len(SAMPLENAMES) > 1):
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA1.2_on_BigWigs_wholeGenome.pdf")) # PCA 1-2
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA2.3_on_BigWigs_wholeGenome.pdf")) # PCA 2-3
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf")) # heatamap_spearman
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf")) # hetamap_pearson
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf")) # scatterplot_spearman
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf")) # scatterplot_pearson

correlation_outputs_peaks = []
if (len(SAMPLENAMES) > 1):
    correlation_outputs_peaks.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/PCA/PCA1.2_on_BigWigs_peakUnion.pdf"))
    correlation_outputs_peaks.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/PCA/PCA2.3_on_BigWigs_peakUnion.pdf"))
    correlation_outputs_peaks.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_heatmap_on_BigWigs_peakUnion_spearmanMethod.pdf"))
    correlation_outputs_peaks.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_heatmap_on_BigWigs_peakUnion_pearsonMethod.pdf"))
    correlation_outputs_peaks.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_scatterplot_on_BigWigs_peakUnion_spearmanMethod.pdf"))
    correlation_outputs_peaks.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_scatterplot_on_BigWigs_peakUnion_pearsonMethod.pdf"))

if (eval(str(config["peak_calling"]["peak_score_heatmap"]["plot_heatmap"]))):
    hetamap_peaks = os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_log1p.rawScores_for_", PEAKCALLER, ".peaks_union_population.pdf"]))
else:
    hetamap_peaks = []


# CNA signal corrections
if (eval(str(config["copy_number_variation"]["call_CNV"]))):
    CNA_corrected_bw = expand(''.join(["03_Normalization/RPM_normalized_CNA.corrected/{sample}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized_CNA.corrected_bs", str(config["copy_number_variation"]["corrected_bigWig_binSize"]), ".bw"]), sample=SAMPLENAMES)
else:
    CNA_corrected_bw = []


# Somatic variants calling
if (eval(str(config["somatic_variants"]["call_SNPs"]))):
    snp = os.path.join(GATKDIR, "SV_count_plots/all.samples_InDel_counts_plot.pdf")
else:
    snp = []


if (eval(str(config["somatic_variants"]["call_indels"]))):
    indels = os.path.join(GATKDIR, "SV_count_plots/all.samples_SNP_counts_plot.pdf")
else:
    indels = []


# ========================================================================================
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ========================================================================================

# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        filtBAM_sorted_woMT = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])), sample=SAMPLENAMES),
        filtBAM_sorted_woMT_index = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])), sample=SAMPLENAMES),
        multiQC_BAM_html = os.path.join(SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/multiQC_report_BAMs_", DUP, ".html"])),
        report_pdf = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"),
        normalized_bigWig = expand(''.join(["03_Normalization/RPM_normalized/{sample}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized.bw"]), sample=SAMPLENAMES),
        narrowPeaks_peaks = expand(os.path.join(PEAKSDIR, "{sample}_mapq{mapq}_woMT_{dup}_qValue{qValue}_peaks.narrowPeak"), sample=SAMPLENAMES, mapq=MAPQ, dup=DUP, qValue=str(config["peak_calling"]["qValue_cutoff"])),
        narrowPeaks_peaks_chr = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks_chr.narrowPeak"])), sample=SAMPLENAMES),
        tobias_output = expand(os.path.join(FOOTDIR, ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_",DUP,"_{extension}"])), sample=SAMPLENAMES, extension = ["AtacBias.pickle", "atacorrect.pdf", "bias.bw", "corrected.bw", "expected.bw", "uncorrected.bw"]),
        summary_file = os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"),
        correlation_outputs = correlation_outputs,
        correlation_outputs_peaks = correlation_outputs_peaks,
        lorenz_plot = os.path.join(SUMMARYDIR, "Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf"),
        CNA_corrected_bw = CNA_corrected_bw,
        rawScores_hetamap_peaks = hetamap_peaks,
        snp = snp,
        indels = indels

    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """
# ========================================================================================


# -------------------------------------------------------------------------------------------

if not (os.path.exists(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))):
    # ----------------------------------------------------------------------------------------
    # Reads alignement
    rule generate_genome_index:
        input:
            genome = ancient(genome_fasta),
        output:
            genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
        threads: 1
        benchmark:
            "benchmarks/generate_genome_index/generate_genome_index---benchmark.txt"
        shell:
            """
            $CONDA_PREFIX/bin/samtools faidx {input.genome}
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """
# ----------------------------------------------------------------------------------------


if (eval(str(config["bam_features"]["skip_bam_filtering"])) == False):
    rule MAPQ_filter:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLES}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mapq_only = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, ".bam"]))),
            bam_mapq_only_sorted = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted.bam"]))),
            bam_mapq_only_sorted_index = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted.bai"]))),
            bam_mapq_only_sorted_woMT = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"]))),
            bam_mapq_only_sorted_woMT_index = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"])))
        params:
            sample = "{SAMPLES}",
            MAPQ_threshold = config["MAPQ_threshold"]
        threads:
            max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
        benchmark:
            "benchmarks/MAPQ_MT_filter/MAPQ_MT_filter---{SAMPLES}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: filtering MAPQ and reemove mitochondrial chromosome...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools view -@ {threads} -h -q {params.MAPQ_threshold} {input.source_bam} -o {output.bam_mapq_only}

            $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted} {output.bam_mapq_only_sorted_index}

            printf '\033[1;36m{params.sample}: Removing MT from BAM...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} | cut -f 1 | grep -v -E '^chrM|^M' | xargs ${{CONDA_PREFIX}}/bin/samtools view -@ {threads} -b {output.bam_mapq_only_sorted} > {output.bam_mapq_only_sorted_woMT}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_woMT} {output.bam_mapq_only_sorted_woMT_index}
            """


    if (eval(str(config["bam_features"]["umi_present"])) == True):
        rule gatk4_markdups_umiAware:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
                umi_metrics = "01_BAM_filtered/umi_metrics/{SAMPLES}_UMI_metrics.txt",
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLES}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, "_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["bam_features"]["remove_duplicates"])).lower(),
                sample = "{SAMPLES}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLES}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLES}_MarkDuplicates.err"
            threads:
                max((workflow.cores-1), 1)
            benchmark:
                "benchmarks/gatk4_markdups_umiAware/gatk4_markdups_umiAware---{SAMPLES}_benchmark.txt"
            shell:
                """
                printf '\033[1;36m{params.sample}: UMI-aware gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/umi_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

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
    else: # Single-end/no-UMI dedup
        rule gatk4_markdups:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLES}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, "_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["bam_features"]["remove_duplicates"])).lower(),
                sample = "{SAMPLES}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLES}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLES}_MarkDuplicates.err"
            threads:
                max((workflow.cores-1), 1)
            benchmark:
                "benchmarks/gatk4_markdups/gatk4_markdups---{SAMPLES}_benchmark.txt"
            shell:
                """
                printf '\033[1;36m{params.sample}: 'standard' gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

                $CONDA_PREFIX/bin/gatk MarkDuplicatesWithMateCigar \
                --INPUT {input.bam_mapq_only_sorted} \
                --OUTPUT {output.bam_mdup} \
                --REMOVE_DUPLICATES {params.remove_duplicates} \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --VALIDATION_STRINGENCY LENIENT \
                --CREATE_INDEX true \
                --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

                $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
                """
else:
    rule bam_link__skip_filtering:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLES}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
            bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, "_flagstat.txt"]))
        params:
            sample = "{SAMPLES}"
        threads:
            max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
        benchmark:
            "benchmarks/bam_link__skip_filtering/bam_link__skip_filtering---{SAMPLES}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample} (skip filtering): linking bam, indexing and computing flagstat...\\n\033[0m'

            mkdir -p 01_BAM_filtered/flagstat

            BAM_REAL=$(realpath {input.source_bam})
            ln -s $BAM_REAL {output.bam_mdup}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mdup} {output.bai_mdup}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """


# ----------------------------------------------------------------------------------------
# BAM reads shifting
rule bam_shifting_and_RPM_normalization:
    input:
        dedup_BAM = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_BAM_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
        genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
    output:
        #temp
        dedup_BAM_sortedByName = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByName.bam"]))),
        dedup_BEDPE_sortedByName = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByName.bedpe"]))),
        dedup_BED_sortedByName_shifted = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByName_shifted.bed"]))),
        dedup_BED_sortedByPos_shifted = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted.bed"]))),
        dedup_BED_sortedByPos_shifted_noBlack = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack.bed"]))),
        dedup_BED_sortedByPos_shifted_noBlack_noIgnoreChr = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack_noIgnoreChr.bed"]))),
        dedup_BedGraph_sortedByPos_shifted_noBlack_RPM = temp(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack_RPM.bedGraph"]))),
        chrSizes = temp("03_Normalization/RPM_normalized/temp_{SAMPLES}_chrSizes_from_genome.txt"),
        # keep
        norm_bw = ''.join(["03_Normalization/RPM_normalized/{SAMPLES}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized.bw"])
    params:
        sample = "{SAMPLES}",
        build_normalization = "03_Normalization/RPM_normalized/bamToBed_log",
        blacklist = BLACKLIST,
        mapq_cutoff = MAPQ,
        minFragmentLength = str(config["bam_features"]["minFragmentLength"]),
        maxFragmentLength = str(config["bam_features"]["maxFragmentLength"]),
        ignore_chr = '|'.join([re.sub('\..*$', '', i) for i in str(config["genomic_annotations"]["ignore_for_normalization"]).split(" ")])
    threads:
        max(math.floor(workflow.cores/len(BAMS)), 1)
    log:
        out = "03_Normalization/RPM_normalized/bamToBed_log/{SAMPLES}_bamToBed.log"
    benchmark:
        "benchmarks/bam_shifting_and_RPM_normalization/bam_shifting_and_RPM_normalization---{SAMPLES}_benchmark.txt"
    priority: 90
    shell:
        """
        printf '\033[1;36m{params.sample}: shifting reads and normalization of the data...\\n\033[0m'
        mkdir -p {params.build_normalization}

        printf '\033[1;36m{params.sample}: re-sorting by read name...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools sort -@ {threads} -n -o {output.dedup_BAM_sortedByName} {input.dedup_BAM}

        printf '\033[1;36m{params.sample}: Bam filtering and conversion to bedPE...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools view -@ {threads} -b -f 3 {output.dedup_BAM_sortedByName} | bedtools bamtobed -i stdin -cigar -bedpe > {output.dedup_BEDPE_sortedByName} 2> {log.out}

        printf '\033[1;36m{params.sample}: Shifting read fragments...\\n\033[0m'
        awk -v OFS='\\t' '{{if($9=="+"){{print $1,$2+4,$6-5,$7,1,$9}}else if($9=="-"){{print $1,$2-5,$6+4,$7,1,$9}}}}' {output.dedup_BEDPE_sortedByName} > {output.dedup_BED_sortedByName_shifted}
        sort -k1,1 -k2,2n {output.dedup_BED_sortedByName_shifted} | grep '+' > {output.dedup_BED_sortedByPos_shifted}

        printf '\033[1;36m{params.sample}: Filter blacklist and fragmentSize...\\n\033[0m'
        $CONDA_PREFIX/bin/bedtools intersect -a {output.dedup_BED_sortedByPos_shifted} -b {params.blacklist} -v | awk '(($3-$2 >= {params.minFragmentLength}) && ($3-$2 <= {params.maxFragmentLength}))' > {output.dedup_BED_sortedByPos_shifted_noBlack}

        printf '\033[1;36m{params.sample}: Compute and RPM-Normalize coverage...\\n\033[0m'
        cut -f 1,2,3 {output.dedup_BED_sortedByPos_shifted_noBlack} | grep -v -E '{params.ignore_chr}' > {output.dedup_BED_sortedByPos_shifted_noBlack_noIgnoreChr}
        TOTREADS=$(printf %d\\\\n $(wc -l < {output.dedup_BED_sortedByPos_shifted_noBlack_noIgnoreChr}))
        cut -f1,2 {input.genome_fai} > {output.chrSizes}
        $CONDA_PREFIX/bin/bedtools genomecov -bg -i {output.dedup_BED_sortedByPos_shifted_noBlack} -g {output.chrSizes} | awk -F "\\t" -v tot="$TOTREADS" '{{print $1,$2,$3,$4/((tot*2)/1000000)}}' OFMT="%.20f" > {output.dedup_BedGraph_sortedByPos_shifted_noBlack_RPM}

        printf '\033[1;36m{params.sample}: Computing RPM-normalized bigWig...\\n\033[0m'
        $CONDA_PREFIX/bin/bedGraphToBigWig {output.dedup_BedGraph_sortedByPos_shifted_noBlack_RPM} {output.chrSizes} {output.norm_bw}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# FastQC on BAMs
rule fastQC_BAMs:
    input:
        dedup_BAM = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_BAM_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"]))
    output:
        html = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, "_fastqc.html"])),
        zip = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, "_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join(config["workflow_configuration"]["output_directory"], "02_BAM_fastQC/"),
        sample = "{SAMPLES}"
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/fastQC_BAMs/fastQC_BAMs---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir {params.fastQC_BAMs_outdir} {input.dedup_BAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for BAMs
if (eval(str(config["bam_features"]["skip_bam_filtering"])) == False):
    rule multiQC_BAMs:
        input:
            BAM_fastqc_zip = expand(os.path.join("02_BAM_fastQC/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, "_fastqc.zip"])), sample=SAMPLENAMES),
            narrowPeaks = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])), sample=SAMPLENAMES)
        output:
            multiqcReportBAM = os.path.join(SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/multiQC_report_BAMs_", DUP, ".html"]))
        params:
            fastQC_BAM_reports_dir = "02_BAM_fastQC/",
            picard_metrics_dir = "01_BAM_filtered/metrics/",
            dedup_BAM_flagstat_dir = "01_BAM_filtered/flagstat/",
            macs_dir = PEAKSDIR,
            multiQC_BAM_outdir = os.path.join(config["workflow_configuration"]["output_directory"], SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/"]))
        benchmark:
            "benchmarks/fmultiQC_BAMs/multiQC_BAMs---benchmark.txt"
        shell:
            """
            printf '\033[1;36mGenerating multiQC report from deduplicated bam quality test...\\n\033[0m'

            $CONDA_PREFIX/bin/multiqc -f \
            --outdir {params.multiQC_BAM_outdir} \
            -n multiQC_report_BAMs_{DUP}.html \
            --dirs \
            {params.fastQC_BAM_reports_dir} \
            {params.picard_metrics_dir} \
            {params.dedup_BAM_flagstat_dir} \
            {params.macs_dir}
            """
else:
    rule multiQC_BAMs:
        input:
            BAM_fastqc_zip = expand(os.path.join("02_BAM_fastQC/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, "_fastqc.zip"])), sample=SAMPLENAMES),
            narrowPeaks = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])), sample=SAMPLENAMES)
        output:
            multiqcReportBAM = os.path.join(SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/multiQC_report_BAMs_", DUP, ".html"]))
        params:
            fastQC_BAM_reports_dir = "02_BAM_fastQC/",
            dedup_BAM_flagstat_dir = "01_BAM_filtered/flagstat/",
            macs_dir = PEAKSDIR,
            multiQC_BAM_outdir = os.path.join(config["workflow_configuration"]["output_directory"], SUMMARYDIR, ''.join(["multiQC_", DUP, "_bams/"]))
        benchmark:
            "benchmarks/multiQC_BAMs/multiQC_BAMs---benchmark.txt"
        shell:
            """
            printf '\033[1;36mGenerating multiQC report from deduplicated bam quality test...\\n\033[0m'

            $CONDA_PREFIX/bin/multiqc -f \
            --outdir {params.multiQC_BAM_outdir} \
            -n multiQC_report_BAMs_{DUP}.html \
            --dirs \
            {params.fastQC_BAM_reports_dir} \
            {params.dedup_BAM_flagstat_dir} \
            {params.macs_dir}
            """

# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule fragment_size_distribution:
    input:
        BAM = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        BAM_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
    output:
        plot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{SAMPLES}_fragment_size_distribution.pdf"),
        fragmentSize_RawFragmentLengths = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLES}_fragmentSize_RawFragmentLengths.txt"),
        fragmentSize_metrics = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLES}_fragmentSize_metrics.txt")
    params:
        build_summary_directory = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log"),
        build_summary_directory_table = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize"),
        sample = "{SAMPLES}",
        plotFormat = "pdf",
        binSize = str(config["quality_controls"]["fragmentSize_window_length"]),
        blacklist = BLACKLIST,
        maxFragmentLength = config["bam_features"]["maxFragmentLength"]
    log:
        out = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLES}_fragmentSize_log.out"),
        err = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLES}_fragmentSize_log.err")
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/fragment_size_distribution/fragment_size_distribution---{SAMPLES}_benchmark.txt"
    priority: -1
    shell:
        """
        printf '\033[1;36m{params.sample}: Plotting the fragment size distribution...\\n\033[0m'

        mkdir -p {params.build_summary_directory}
        mkdir -p {params.build_summary_directory_table}

        $CONDA_PREFIX/bin/bamPEFragmentSize \
        -p {threads} \
        -b {input.BAM} \
        --plotFileFormat {params.plotFormat} \
        --plotTitle {params.sample} \
        --samplesLabel {params.sample} \
        --binSize {params.binSize} \
        --maxFragmentLength {params.maxFragmentLength} \
        --blackListFileName {params.blacklist} \
        --outRawFragmentLengths {output.fragmentSize_RawFragmentLengths} \
        --table {output.fragmentSize_metrics} \
        --histogram {output.plot} > {log.out} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule fragment_size_distribution_report:
    input:
        plots = expand(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{sample}_fragment_size_distribution.pdf"), sample = SAMPLENAMES)
    output:
        report_pdf = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf")
    params:
        distribution_plots_pattern = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/*_fragment_size_distribution.pdf")
    threads: 1
    benchmark:
        "benchmarks/fragment_size_distribution_report/fragment_size_distribution_report---benchmark.txt"
    shell:
        """
        printf '\033[1;36mMerging fragmentSizeDistribution reports in a unique PDF...\\n\033[0m'
        $CONDA_PREFIX/bin/pdfcombine {params.distribution_plots_pattern} -o {output.report_pdf} -sf
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# PeakCalling on uncorrected bams (MACS3 callpeak)
rule peakCalling_MACS3:
    input:
        dedup_BAM_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_BAM_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
    output:
        peaks_xls = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.xls"])),
        narrowPeaks = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])),
        summits = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_summits.bed"]))
    params:
        genomeSize = str(config["genomic_annotations"]["effective_genomeSize"]),
        peak_caller = PEAKCALLER.lower(),
        peaks_dir = PEAKSDIR,
        qValue = str(config["peak_calling"]["qValue_cutoff"]),
        basename = ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"])]),
        summits = SUMMITS,
        sample = "{SAMPLES}"
    log:
        out = os.path.join(''.join([PEAKSDIR,"log/"]), ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), ".log"]))
    benchmark:
        "benchmarks/peakCalling_MACS3/peakCalling_MACS3---{SAMPLES}_benchmark.txt"
    priority: 100
    shell:
        """
        printf '\033[1;36m{params.sample}: Calling peaks by {params.peak_caller}...\\n\033[0m'

        $CONDA_PREFIX/bin/{params.peak_caller} callpeak \
        -t {input.dedup_BAM_sorted} \
        -g {params.genomeSize} \
        -n {params.basename} \
        -q {params.qValue} \
        -f BAMPE \
        --outdir {params.peaks_dir} \
        --keep-dup all \
        --nolambda \
        {params.summits} 2> {log.out}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Computation of the counts summary table
rule counts_summary:
    input:
        flagstat_filtered = expand(os.path.join("01_BAM_filtered/flagstat/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, "_flagstat.txt"])), sample = SAMPLENAMES),
        dedup_BED_sortedByPos_shifted_noBlack = expand(os.path.join("03_Normalization/RPM_normalized/", ''.join(["temp_{sample}_mapq", MAPQ, "_sorted_woMT_",DUP,"_sortedByPos_shifted_noBlack.bed"])), sample=SAMPLENAMES),
        peaks_file = expand(os.path.join(PEAKSDIR, ''.join(["{sample}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])), sample = SAMPLENAMES),
        norm_bw = expand("03_Normalization/RPM_normalized/{sample}_mapq{mapq}_woMT_{dup}_shifted_RPM.normalized.bw", sample=SAMPLENAMES, dup=DUP, mapq=MAPQ)
    output:
        summary_file = os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"),
        summary_file_temp = temp(os.path.join(SUMMARYDIR, "Counts/summary_file.temp"))
    params:
        build_summary_directory = os.path.dirname(SUMMARYDIR),
        sample_list = str(' '.join(SAMPLENAMES)),
        peaks_dir = PEAKSDIR,
        FRiP_threshold = config["peak_calling"]["FRiP_threshold"],
        MAPQ = MAPQ,
        DUP = DUP
    threads: 1
    benchmark:
        "benchmarks/counts_summary/counts_summary---benchmark.txt"
    priority: 80
    shell:
        """
        mkdir -p {params.build_summary_directory}/Counts/subread_featureCounts_output/

        printf '\033[1;36mGeneration of a general counts summary table...\\n\033[0m'
        printf Sample'\\t'dedup_BAM'\\t'shifted_readsM'\\t'loss_post_shifting'\\t'n.peaks'\\t'FRiP.perc'\\t'FRiP.quality'\\n' > {output.summary_file}

        for NAME in {params.sample_list}
        do
            printf '\033[1;36m     - %s: adding stats to summary table...\\n\033[0m' $NAME

            mkdir -p {params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/

            woMT_BAM=$(grep mapped 01_BAM_filtered/flagstat/${{NAME}}_mapq{params.MAPQ}_sorted_woMT_{params.DUP}_flagstat.txt | head -n 1 | cut -f 1 -d ' ')

            dedupBAM=$(grep mapped 01_BAM_filtered/flagstat/${{NAME}}_mapq{params.MAPQ}_sorted_woMT_{params.DUP}_flagstat.txt | head -n 1 | cut -f 1 -d ' ')

            halfShiftedBEDPE=$(printf %d\\\\n $(wc -l < 03_Normalization/RPM_normalized/temp_${{NAME}}_mapq{params.MAPQ}_sorted_woMT_{params.DUP}_sortedByPos_shifted_noBlack.bed))
            shiftedBEDPE=$((halfShiftedBEDPE + halfShiftedBEDPE))
            lossReads=$((dedupBAM - shiftedBEDPE))

            peaks=$(wc -l {params.peaks_dir}${{NAME}}*.*Peak | cut -f 1 -d ' ')

            awk 'BEGIN{{FS=OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' {params.peaks_dir}${{NAME}}*.*Peak > {params.peaks_dir}${{NAME}}.saf
            FEATURECOUNTSLOG={params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/${{NAME}}.readCountInPeaks.log
            featureCounts -p -a {params.peaks_dir}${{NAME}}.saf -F SAF -o {params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/${{NAME}}.readCountInPeaks 01_BAM_filtered/${{NAME}}_mapq*_sorted_woMT_{DUP}.bam 2> ${{FEATURECOUNTSLOG}}
            rm {params.peaks_dir}${{NAME}}.saf
            frip=$(grep 'Successfully assigned alignments' ${{FEATURECOUNTSLOG}} | sed -e 's/.*(//' | sed 's/%.*$//')
            fripScore=$(echo $frip | sed 's/\\..*$//')
            fripLabel=$(if [ $fripScore -ge {params.FRiP_threshold} ]; then echo 'good'; else echo 'bad'; fi)

            printf ${{NAME}}'\\t'$dedupBAM'\\t'$shiftedBEDPE'\\t'$lossReads'\\t'$peaks'\\t'$frip'\\t'$fripLabel'\\n' >> {output.summary_file}
        done

        uniq -u {output.summary_file} > {output.summary_file_temp}
        (head -n 1 {output.summary_file_temp} && tail -n +2 {output.summary_file_temp} | sort -k 1) > {output.summary_file}
        """



# Add chr to peak_suffix
rule add_chr_to_peaks:
    input:
        peaks = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])),
        summary_file = os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"),
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed")
    output:
        peaks_chr = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks_chr.narrowPeak"]))
    params:
        blacklist = BLACKLIST
    threads: 1
    benchmark:
        "benchmarks/add_chr_to_peaks/add_chr_to_peaks---{SAMPLES}_benchmark.txt"
    shell:
        """
        $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a {input.peaks} -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > {output.peaks_chr}
        """




# ----------------------------------------------------------------------------------------
# Generation of matrix scores
rule multiBigwigSummary:
    input:
        norm_bw = expand("03_Normalization/RPM_normalized/{sample}_mapq"+MAPQ+"_woMT_"+DUP+"_shifted_RPM.normalized.bw", sample=SAMPLENAMES)
    output:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.npz")
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "/Sample_comparisons/"])),
        labels = ' '.join(SAMPLENAMES),
        window = config["quality_controls"]["multiBigwigSummary_binning_window_size"],
        blacklist = BLACKLIST
    threads:
        max((workflow.cores-1), 1)
    benchmark:
        "benchmarks/multiBigwigSummary/multiBigwigSummary---benchmark.txt"
    shell:
        """
        printf '\033[1;36mComparing the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        $CONDA_PREFIX/bin/multiBigwigSummary bins -p {threads} -b {input.norm_bw} --labels {params.labels} --binSize {params.window} --blackListFileName {params.blacklist} -o {output.matrix}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of samples PCA and Heatmap
rule PCA_and_samples_correlation:
    input:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.npz")
    output:
        PCA_OneTwo = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA1.2_on_BigWigs_wholeGenome.pdf"),
        PCA_TwoThree = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA2.3_on_BigWigs_wholeGenome.pdf"),
        hetamap_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        hetamap_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"),
        scatterplot_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        scatterplot_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf")
    params:
        make_directory = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/"),
        heatmap_color = config["quality_controls"]["correlation_heatmap_color"]
    benchmark:
        "benchmarks/PCA_and_samples_correlation/PCA_and_samples_correlation---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting the correlation and variability of the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        printf '\033[1;36m    - plotting PCAs...\\n\033[0m'
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_OneTwo} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf'
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_TwoThree} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf' --PCs 2 3


        printf '\033[1;36m    - plotting Spearman correlation heatmap...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
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
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
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
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
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
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
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


# ----------------------------------------------------------------------------------------
# Generation of Lorenz curves
rule Lorenz_curve:
    input:
        dedup_bam_sorted = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])), sample=SAMPLENAMES),
        dedup_bam_sorted_index = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])), sample=SAMPLENAMES)
    output:
        lorenz_plot = os.path.join(SUMMARYDIR, "Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf")
    params:
        all_bams = ' '.join(expand(os.path.join(''.join(["01_BAM_filtered/{sample}_mapq", MAPQ ,"_sorted_woMT_", DUP, ".bam"])), sample=SAMPLENAMES)),
        labels = ' '.join(SAMPLENAMES),
        blacklist = BLACKLIST,
        binSize = config["quality_controls"]["plotFingerprint"]["binSize"],
        sampledRegions = config["quality_controls"]["plotFingerprint"]["sampledRegions"],
        extra_params = config["quality_controls"]["plotFingerprint"]["extra_parameters"]
    threads:
        max((workflow.cores-1), 1)
    benchmark:
        "benchmarks/Lorenz_curve/Lorenz_curve---benchmark.txt"
    priority: -10
    shell:
        """
        printf '\033[1;36mPlotting the Lorenz curves-Fingerprint for all samples...\\n\033[0m'

        $CONDA_PREFIX/bin/plotFingerprint \
        --bamfiles {params.all_bams} \
        --plotFile {output.lorenz_plot} \
        --labels {params.labels} \
        --blackListFileName {params.blacklist} \
        --binSize {params.binSize} \
        --numberOfSamples {params.sampledRegions}\
        -p {threads} {params.extra_params}
        """
# ----------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
###                           COPY NUMBER VARIATION                          ###
# ------------------------------------------------------------------------------
# generate copywriteR script
rule generate_copywriteR_Rscript:
    output:
        script = os.path.join(COPYWRITERDIR, "CopywriteR_script.R")
    threads: 1
    benchmark:
        "benchmarks/generate_copywriteR_Rscript/generate_copywriteR_Rscript---genomeMappability_benchmark.txt"
    run:
        shell("printf '\033[1;36mGenerating the Rscript used to compute and visualize the log2 fequencies (CopywriteR)...\\n\033[0m'")

        x = """# Load parameters
        args = commandArgs(trailingOnly = TRUE)

        sample_id = as.character(strsplit(grep('--SAMPLENAME*', args, value = TRUE), split = '=')[[1]][[2]])
        target_id = as.character(strsplit(grep('--BAMFILE*', args, value = TRUE), split = '=')[[1]][[2]])
        peaks = as.character(strsplit(grep('--PEAKS*', args, value = TRUE), split = '=')[[1]][[2]])
        cores = as.numeric(strsplit(grep('--CORES*', args, value = TRUE), split = '=')[[1]][[2]])
        data.folder = as.character(strsplit(grep('--COPYWRITERDIR*', args, value = TRUE), split = '=')[[1]][[2]])
        kb.resolution = as.numeric(strsplit(grep('--RESOLUTION*', args, value = TRUE), split = '=')[[1]][[2]])
        genome = as.character(strsplit(grep('--GENOME*', args, value = TRUE), split = '=')[[1]][[2]])
        CNA_threshold = as.numeric(strsplit(grep('--CNATHRESHOLD*', args, value = TRUE), split = '=')[[1]][[2]])
        CNA_color = as.character(strsplit(grep('--LINECOLOR*', args, value = TRUE), split = '=')[[1]][[2]])
        point_size = as.numeric(strsplit(grep('--POINTSIZE*', args, value = TRUE), split = '=')[[1]][[2]])
        point_alpha = as.numeric(strsplit(grep('--POINTALPHA*', args, value = TRUE), split = '=')[[1]][[2]])



        #########################
        # Load libs
        require(dplyr)
        require(ggplot2)


        ### Set-up copywriteR
        bp.param = BiocParallel::SnowParam(workers = cores, type = "SOCK")


        # Run CopyWriteR
        sample_dir = paste0(gsub("/$","",data.folder), "/", sample_id)

        CopywriteR::CopywriteR(sample.control = data.frame(target_id, target_id),
                               destination.folder = sample_dir,
                               reference.folder = paste0(data.folder,"/",genome,"_",kb.resolution,"kb"),
                               capture.regions.file = peaks,
                               bp.param = bp.param,
                               keep.intermediary.files = F)


        # Plot by copyWriteR
        CopywriteR::plotCNA(destination.folder = sample_dir)


        # Read segmentation
        loadRData = function(fileName){
          load(fileName)
          get(ls()[ls() != "fileName"])}

        segmentation = loadRData(paste0(sample_dir, "/CNAprofiles/segment.Rdata"))


        # Read log2_readCounts and plot it
        log2.tb =
          read.table(file = file.path(sample_dir, "CNAprofiles", "log2_read_counts.igv"), header = TRUE) %>%
          #dplyr::filter(Chromosome != "Y") %>%
          dplyr::mutate(mid.point = (Start+End)/2,
                        Chromosome = factor(Chromosome, levels = unique(Chromosome)))
        names(log2.tb)[5] = "log2.value"


        # Re-plot CNA
        ymax = ceiling(max(abs(boxplot(log2.tb$log2.value, na.rm = T)$stats[,1]))) + 1

        CNA.plot =
          ggplot(log2.tb,
                 aes(x = mid.point,
                     y = log2.value)) +
          geom_point(size = point_size,
                     stroke = NA,
                     alpha = point_alpha) +
          theme_classic() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(color = "black"),
                axis.text.y = element_text(color = "black"),
                panel.spacing = unit(0,'lines'),
                strip.background = element_blank(),
                strip.placement = "outside",
                panel.grid = element_blank(),
                panel.border = element_blank()) +
          geom_hline(yintercept = c(-1,1)*log2(CNA_threshold), linetype = 2, color = CNA_color) +
          geom_vline(xintercept = +Inf, linetype = 3, color = "gray60") +
          ylim(c(-1,1)*ymax) +
          xlab("Chromosome") +
          ylab("log2(copy number)") +
          ggtitle(sample_id) +
          geom_segment(data =
                         dplyr::filter(segmentation$output, abs(seg.mean) > 0) %>%
                         dplyr::mutate(chrom = gsub("23", "X", chrom)) %>%
                         dplyr::mutate(chrom = gsub("24", "Y", chrom)) %>%
                         dplyr::mutate(chrom = factor(chrom)) %>%
                         dplyr::rename(Chromosome = chrom,
                                       log2.value = seg.mean),
                       aes(x = loc.start,
                           xend = loc.end,
                           y = log2.value,
                           yend = log2.value),
                       show.legend = F,
                       inherit.aes = F,
                       color = CNA_color,
                       linewidth = 0.75) +
          facet_grid(~ Chromosome,
                     space = "free_x",
                     scales = "free_x",
                     switch = "x")


        pdf(paste0(sample_dir, "/CNA.plot_", sample_id, "_all.chr_",kb.resolution,"kb.pdf"), width = 16, height = 5)
        print(CNA.plot)
        invisible(dev.off())"""

        with open(str(output.script),"w+") as f:
          f.writelines(x)


# Calculate genome GC_mappability
rule generate_copywriteR_genome_map:
    output:
        genome_mappability = os.path.join(COPYWRITERDIR, ''.join([str(re.sub("_.*$", "", genome_id).lower()),"_",str(config["copy_number_variation"]["kb_bin_resolution"]),"kb/GC_mappability.rda"]))
    params:
        data_folder = os.path.join(home_dir, COPYWRITERDIR),
        resolution = str(config["copy_number_variation"]["kb_bin_resolution"]),
        genome = genome_id
    log:
        out = os.path.join(COPYWRITERDIR, "logs/genome_mappability_file_generation_copywriteR.log")
    benchmark:
        "benchmarks/generate_copywriteR_genome_map/generate_copywriteR_genome_map---genomeMappability_benchmark.txt"
    shell:
        """
        ### Generate genome index
        printf '\033[1;36mGenerating genome mappability files (CopywriteR)...\\n\033[0m'
        $CONDA_PREFIX/bin/Rscript -e 'CopywriteR::preCopywriteR(output.folder = "{params.data_folder}", bin.size = {params.resolution}*1000, ref.genome = "{params.genome}", prefix = "")' &> {log.out}
        """


# Running CopyWriteR
rule CopywriteR_CNA_detection:
    input:
        copywriteR_script = os.path.join(COPYWRITERDIR, "CopywriteR_script.R"),
        target_bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        target_bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
        peaks = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])),
        genome_mappability = os.path.join(COPYWRITERDIR, ''.join([str(re.sub("_.*$", "", genome_id).lower()),"_",str(config["copy_number_variation"]["kb_bin_resolution"]),"kb/GC_mappability.rda"])),
        norm_bw = ''.join(["03_Normalization/RPM_normalized/{SAMPLES}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized.bw"]),
        genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
    output:
        #CNA_plot = os.path.join(COPYWRITERDIR, ''.join(["{SAMPLES}/CNA.plot_{SAMPLES}_all.chr_", str(config["copy_number_variation"]["kb_bin_resolution"]), "kb.pdf"])),
        collapsed_peaks = temp(os.path.join(COPYWRITERDIR, "{SAMPLES}/{SAMPLES}_collapsed_peaks.bed")),
        log2_read_counts = os.path.join(COPYWRITERDIR, "{SAMPLES}/CNAprofiles/log2_read_counts.igv"),
        #collapsed_peaks = temp(os.path.join(COPYWRITERDIR, ''.join(["temp_{SAMPLES}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_collapsed.peaks.bed"]))),
        #CNA_correctionFactors_bedGraph_sorted = os.path.join(''.join([COPYWRITERDIR, "CNA_profiles/{SAMPLES}/{SAMPLES}_mapq",MAPQ,"_sorted_woMT_",DUP,"_linear.ratio_read_counts_filtered_sorted.bedGraph"])),
        #chrSizes = temp(''.join([COPYWRITERDIR, "temp_{SAMPLES}_chrSizes_from_genome.txt"])),
        #CNA_corrected_bw = os.path.join(''.join([COPYWRITERDIR, "CNA_corrected_signal/{SAMPLES}_mapq",MAPQ,"_sorted_woMT_",DUP,"_CNA.corrected.bw"]))
    params:
        outdir = os.path.join(home_dir, COPYWRITERDIR, "{SAMPLES}"),
        sample = "{SAMPLES}",
        target = os.path.join(home_dir, "01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        collapsed_peaks = os.path.join(home_dir, COPYWRITERDIR, "{SAMPLES}/{SAMPLES}_collapsed_peaks.bed"),
        data_folder = os.path.join(home_dir, COPYWRITERDIR),
        resolution = str(config["copy_number_variation"]["kb_bin_resolution"]),
        genome = re.sub("_.*$", "", genome_id).lower(),
        peak_suffix = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_",DUP,"_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks"])),
        CNA_threshold = abs(config["copy_number_variation"]["CNA_threshold"]),
        CNA_plot_line_colors = config["copy_number_variation"]["CNA_plot_line_colors"],
        CNA_plot_point_size = config["copy_number_variation"]["CNA_plot_point_size"],
        CNA_plot_point_alpha = config["copy_number_variation"]["CNA_plot_point_transparency"]
    threads:
        max((workflow.cores-1), 1)
    log:
        out = os.path.join(COPYWRITERDIR, "logs/{SAMPLES}_copywriteR.log")
    benchmark:
        "benchmarks/CopywriteR/CopywriteR---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Running CopywriteR...\\n\033[0m'

        mkdir -p {params.outdir}/CNAprofiles
        rm -r {params.outdir}/CNAprofiles
        mkdir -p {params.outdir}

        $CONDA_PREFIX/bin/bedtools merge -i {input.peaks} > {output.collapsed_peaks}

        $CONDA_PREFIX/bin/Rscript {input.copywriteR_script} \
        --SAMPLENAME={params.sample} \
        --BAMFILE={params.target} \
        --PEAKS={params.collapsed_peaks} \
        --CORES={threads} \
        --COPYWRITERDIR={params.data_folder} \
        --RESOLUTION={params.resolution} \
        --GENOME={params.genome} \
        --CNATHRESHOLD={params.CNA_threshold} \
        --LINECOLOR={params.CNA_plot_line_colors} \
        --POINTSIZE={params.CNA_plot_point_size} \
        --POINTALPHA={params.CNA_plot_point_alpha} &> {log.out}
        """


# Conversion of CNA counts (filtered to threshold) to bedGraph
rule convert_CNA_to_bedGraph:
    input:
        log2_read_counts = os.path.join(COPYWRITERDIR, "{SAMPLES}/CNAprofiles/log2_read_counts.igv")
    output:
        bedGraph_filtered = temp(os.path.join(COPYWRITERDIR, ''.join(["{SAMPLES}/CNAprofiles/{SAMPLES}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bedGraph"])))
    params:
        sample = "{SAMPLES}",
        CNA_threshold = int(config["copy_number_variation"]["CNA_threshold"])
    threads: 1
    benchmark:
        "benchmarks/convert_CNA_to_bedGraph/convert_CNA_to_bedGraph---{SAMPLES}_benchmark.txt"
    run:
        shell("printf '\033[1;36m{params.sample}: Filter and converting CNA counts to bedGraph...\\n\033[0m'")

        counts = pd.read_csv(input.log2_read_counts,  sep='\t+', engine='python', skiprows=1)
        counts.iloc[:,4] = 2**counts.iloc[:,4]

        bdg = counts.iloc[:,[0,1,2,4]]
        bdg_filt = bdg[abs(bdg.iloc[:,3]) >= params.CNA_threshold]
        bdg_filt.to_csv(output.bedGraph_filtered, header=False, index=False, sep='\t')



# Extract chromosome sizes fro converstion to bigWig
rule extract_chromosome_sizes:
    input:
        genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
    output:
        chrSizes = "03_Normalization/RPM_normalized/.genome_chromosome_sizes.txt"
    benchmark:
        "benchmarks/extract_chromosome_sizes/extract_chromosome_sizes---benchmark.txt"
    shell:
        """
        printf '\033[1;36mExtracting chromosome sizes from genome index (.fai)...\\n\033[0m'
        cut -f1,2 {input.genome_fai} > {output.chrSizes}
        """


# Correct CNV on RPM normalized bigWigs
rule correct_CNV_normalized_bigWigs:
    input:
        bedGraph_filtered = os.path.join(COPYWRITERDIR, ''.join(["{SAMPLES}/CNAprofiles/{SAMPLES}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bedGraph"])),
        chrSizes = "03_Normalization/RPM_normalized/.genome_chromosome_sizes.txt",
        normalized_bigWig = ''.join(["03_Normalization/RPM_normalized/{SAMPLES}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized.bw"])
    output:
        sorted_bedGraph = os.path.join(COPYWRITERDIR, ''.join(["{SAMPLES}/CNAprofiles/{SAMPLES}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts_sorted.bedGraph"])),
        CNA_corrected_bw = ''.join(["03_Normalization/RPM_normalized_CNA.corrected/{SAMPLES}_mapq", MAPQ, "_woMT_", DUP ,"_shifted_RPM.normalized_CNA.corrected_bs", str(config["copy_number_variation"]["corrected_bigWig_binSize"]), ".bw"])
    params:
        sample = "{SAMPLES}",
        bigWig_filtered = os.path.join(COPYWRITERDIR, ''.join(["{SAMPLES}/CNAprofiles/{SAMPLES}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bw"])),
        bw_binSize = config["copy_number_variation"]["corrected_bigWig_binSize"]
    threads:
        max((workflow.cores-1), 1)
    log:
        bedGraphToBigWig_out = os.path.join(COPYWRITERDIR, "logs/{SAMPLES}_bedGraphToBigWig.log"),
        bigwigCompare_out = os.path.join("03_Normalization/RPM_normalized_CNA.corrected/log/{SAMPLES}_CNA.correction_bigwig.out")
    benchmark:
        "benchmarks/correct_CNV_normalized_bigWigs/correct_CNV_normalized_bigWigs---{SAMPLES}_benchmark.txt"
    shell:
        """
        sort -k1,1 -k2,2n {input.bedGraph_filtered} > {output.sorted_bedGraph}
        NROWSBDG=$(wc -l {output.sorted_bedGraph} | head -n 1 | cut -f 1 -d ' ')

        if [[ $NROWSBDG -gt 0 ]]
        then
            printf '\033[1;36m{params.sample}: signal correction for CNA...\\n\033[0m'

            $CONDA_PREFIX/bin/bedGraphToBigWig {output.sorted_bedGraph} {input.chrSizes} {params.bigWig_filtered} &> {log.bedGraphToBigWig_out}

            $CONDA_PREFIX/bin/bigwigCompare \
            --bigwig1 {input.normalized_bigWig} \
            --bigwig2 {params.bigWig_filtered} \
            -o {output.CNA_corrected_bw} \
            --operation ratio \
            --binSize {params.bw_binSize} \
            -p {threads} \
            -of bigwig \
            --pseudocount 0 1 &> {log.bigwigCompare_out}
        else
            cp {input.normalized_bigWig} {output.CNA_corrected_bw}
        fi
        """



# ----------------------------------------------------------------------------------------
# Absolute peaks file and relative matrix score generation for called peaks
rule all_peaks_file_and_score_matrix:
    input:
        norm_bw = expand("03_Normalization/RPM_normalized/{sample}_mapq{mapq}_woMT_{dup}_shifted_RPM.normalized.bw", sample=SAMPLENAMES,  dup=DUP, mapq=MAPQ),
        peaks_file = expand(os.path.join(str(PEAKSDIR), "{sample}_mapq{mapq}_woMT_{dup}_qValue{qValue}_peaks.narrowPeak"), sample=SAMPLENAMES, mapq=MAPQ, dup=DUP, qValue=str(config["peak_calling"]["qValue_cutoff"]))
    output:
        concatenation_bed = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation.bed")),
        concatenation_bed_sorted = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation_sorted.bed")),
        concatenation_bed_collapsed = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation_collapsed.bed")),
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"),
        score_matrix_peaks = os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_", PEAKCALLER, ".npz"])),
        score_matrix_peaks_table = os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_table_", PEAKCALLER, ".tsv"]))
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Peak_comparison/"])),
        peak_caller = PEAKCALLER,
        peaks_dir = PEAKSDIR,
        labels = ' '.join(SAMPLENAMES),
        blacklist = BLACKLIST
    threads:
        max((workflow.cores-1), 1)
    benchmark:
        "benchmarks/all_peaks_file_and_score_matrix/all_peaks_file_and_score_matrix---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating a file result of the merge of all the {params.peak_caller} peaks...\\n\033[0m'

        mkdir -p {params.make_directory}

        cat {params.peaks_dir}*.*Peak >> {output.concatenation_bed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed} > {output.concatenation_bed_sorted}

        $CONDA_PREFIX/bin/bedtools merge -i {output.concatenation_bed_sorted} | uniq > {output.concatenation_bed_collapsed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed_collapsed} > {output.concatenation_bed_collapsed_sorted}


        printf '\033[1;36mComputing the score matrix for all the {params.peak_caller} peaks per each sample...\\n\033[0m'

        $CONDA_PREFIX/bin/multiBigwigSummary BED-file \
        -p {threads} \
        -b {input.norm_bw} \
        -o {output.score_matrix_peaks} \
        --BED {output.concatenation_bed_collapsed_sorted} \
        --blackListFileName {params.blacklist} \
        --outRawCounts {output.score_matrix_peaks_table} \
        --labels {params.labels}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Compute peaks z-scores and plot heatmap for called peaks
rule peaks_zScores_and_heatmap:
    input:
        score_matrix_peaks_table = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_table_"+PEAKCALLER+".tsv")
    output:
        rawScores_hetamap = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_log1p.rawScores_for_"+PEAKCALLER+".peaks_union_population.pdf")
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/"])),
        heatmap_color = config["peak_calling"]["peak_score_heatmap"]["raw_heatmap_color"],
        zScore_heatmap_color = config["peak_calling"]["peak_score_heatmap"]["zScore_heatmap_color"],
        heatmap_basename_rawScores = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_log1p.rawScores_for_"+PEAKCALLER+".peaks_union_population"),
        heatmap_basename_zScore = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_zScores_for_"+PEAKCALLER+".peaks_union_population"),
        n_samples = len(SAMPLENAMES),
        peak_caller = PEAKCALLER
    threads:
        max((workflow.cores-1), 1)
    benchmark:
        "benchmarks/peaks_zScores_and_heatmap/peaks_zScores_and_heatmap---benchmark.txt"
    run:
        # Messege
        shell("printf '\033[1;36mPeak-scores heatmaps computation:\\n\033[0m'")
        shell("printf '\033[1;36m    - Loading {params.peak_caller} peak raw score table...\\n\033[0m'")
        shell("mkdir -p {params.make_directory}")

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
        shell("printf '\033[1;36m    - Plotting the {params.peak_caller} peak raw score hetamaps...\\n\033[0m'")

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
            shell("printf '\033[1;36m    - Computing the zScores for {PEAKCALLER} peaks...\\n\033[0m'")

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
            shell("printf '\033[1;36m    - Plotting the {PEAKCALLER} peak zScores hetamaps...\\n\033[0m'")

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


# ----------------------------------------------------------------------------------------
# Compute peaks z-scores and plot heatmap for called peaks
rule peaks_correlation_and_PCA_at_peaks:
    input:
        matrix = os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_",PEAKCALLER,".npz"]))
    output:
        PCA_OneTwo = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/PCA/PCA1.2_on_BigWigs_peakUnion.pdf"),
        PCA_TwoThree = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/PCA/PCA2.3_on_BigWigs_peakUnion.pdf"),
        hetamap_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_heatmap_on_BigWigs_peakUnion_spearmanMethod.pdf"),
        hetamap_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_heatmap_on_BigWigs_peakUnion_pearsonMethod.pdf"),
        scatterplot_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_scatterplot_on_BigWigs_peakUnion_spearmanMethod.pdf"),
        scatterplot_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/Correlation_scatterplot_on_BigWigs_peakUnion_pearsonMethod.pdf")
    params:
        make_directory_corr = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Peak_comparison/Correlations/"])),
        make_directory_PCA = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Peak_comparison/PCA/"])),
        heatmap_color = config["quality_controls"]["correlation_heatmap_color"]
    benchmark:
        "benchmarks/peaks_correlation_and_PCA_at_peaks/peaks_correlation_and_PCA_at_peaks---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting the correlation and variability of the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory_corr}
        mkdir -p {params.make_directory_PCA}

        printf '\033[1;36m    - plotting PCAs...\\n\033[0m'
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_OneTwo} -T 'PCA on BigWigs (peak union)' --plotFileFormat 'pdf'
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_TwoThree} -T 'PCA on BigWigs (peak union)' --plotFileFormat 'pdf' --PCs 2 3


        printf '\033[1;36m    - plotting Spearman correlation heatmap...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs (peak union)" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.hetamap_spearman}

        printf '\033[1;36m    - plotting Pearson correlation heatmap...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs (peak union)" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.hetamap_pearson}



        printf '\033[1;36m    - plotting Spearman correlation scatterplot...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs - ln values (peak union)" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_spearman}

        printf '\033[1;36m    - plotting Pearson correlation scatterplot...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs - ln values (peak union)" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_pearson}
        """
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Correct signal for Tn5 bias
rule TOBIAS_ATACorrect:
    input:
        dedup_bam_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        dedup_bam_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
        peaks = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_", DUP, "_qValue", str(config["peak_calling"]["qValue_cutoff"]), "_peaks.narrowPeak"])),
        genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
    output:
        tobias_output = expand(os.path.join(FOOTDIR, ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_",DUP,"_{extension}"])), extension = ["AtacBias.pickle", "atacorrect.pdf", "bias.bw", "corrected.bw", "expected.bw", "uncorrected.bw"], allow_missing=True)
    params:
        sample = "{SAMPLES}",
        foot_folder = FOOTDIR,
        genome = genome_fasta,
        blacklist = BLACKLIST
    log:
        out = os.path.join(FOOTDIR, "log_tobias/{SAMPLES}_tobias.log")
    threads:
        max(math.floor(workflow.cores/4), 1)
    benchmark:
        "benchmarks/TOBIAS_ATACorrect/TOBIAS_ATACorrect---{SAMPLES}_benchmark.txt"
    priority: -5
    shell:
        """
        printf '\033[1;36m{params.sample}: generating footprint scores with TOBIAS...\\n\033[0m'

        mkdir -p {params.foot_folder}/log_tobias

        $CONDA_PREFIX/bin/TOBIAS ATACorrect \
        --bam {input.dedup_bam_sorted} \
        --genome {params.genome} \
        --peaks {input.peaks} \
        --blacklist {params.blacklist} \
        --outdir {params.foot_folder} \
        --cores {threads} &> {log.out}
        """
# ----------------------------------------------------------------------------------------


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  VARIANT CALLING  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Create reference genome dictionary
rule make_reference_genome_dictionary:
    input:
        genome = ancient(genome_fasta)
    output:
        genome_dict = ''.join([re.sub("[.]([a-z]|[A-Z])*$", "",genome_fasta),'.dict'])
    benchmark:
        "benchmarks/make_reference_genome_dictionary/make_reference_genome_dictionary---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating genome dictionary...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk CreateSequenceDictionary REFERENCE={input.genome} OUTPUT={output.genome_dict}
        """
# ----------------------------------------------------------------------------------------


# run base score recalibration (BSQR) of the bams
rule GATK_bam_base_quality_score_recalibration:
    input:
        genome_dict = ancient(''.join([re.sub("[.]([a-z]|[A-Z])*$", "",genome_fasta),'.dict'])),
        target_bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bam"])),
        target_bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_", DUP, ".bai"])),
    output:
        target_bam_withRG = temp(os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_withRG.bam"]))),
        bsqr_table = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["bsqr_tables/{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.table"])),
        target_bam_bsqr = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bam"])),
        target_bam_bsqr_index = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bai"]))
    params:
        sample = "{SAMPLES}",
        gatk_directory = GATKDIR,
        genome = genome_fasta,
        dbsnp = config["somatic_variants"]["dbsnp_file"]
    log:
        readGroup_log = os.path.join(GATKDIR, "recalibrated_bams/logs/{SAMPLES}_addReadGroup.log"),
        BaseRecalibrator_log = os.path.join(GATKDIR, "recalibrated_bams/logs/{SAMPLES}_BaseRecalibrator.log"),
        ApplyBQSR_log = os.path.join(GATKDIR, "recalibrated_bams/logs/{SAMPLES}_ApplyBQSR.log")
    threads:
        max(math.floor(workflow.cores/2), 1)
    benchmark:
        "benchmarks/GATK_bam_base_quality_score_recalibration/GATK_bam_base_quality_score_recalibration---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Adding Read Groups to filtered bams...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk AddOrReplaceReadGroups \
        -I {input.target_bam} \
        -O {output.target_bam_withRG} \
        -RGID 1 \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM {params.sample} &> {log.readGroup_log}


        printf '\033[1;36m{params.sample}: Base Quality Score Recalibration of the deduplicated bam...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk BaseRecalibrator \
        --input {output.target_bam_withRG} \
        --known-sites {params.dbsnp} \
        --output {output.bsqr_table} \
        --reference {params.genome} &> {log.BaseRecalibrator_log}

        $CONDA_PREFIX/bin/gatk ApplyBQSR \
        -R {params.genome} \
        -I {output.target_bam_withRG} \
        --bqsr-recal-file {output.bsqr_table} \
        -O {output.target_bam_bsqr} &> {log.ApplyBQSR_log}

        printf '\033[1;36m{params.sample}: Indexing recalibrated bam...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.target_bam_bsqr} {output.target_bam_bsqr_index}
        """

# ----------------------------------------------------------------------------------------


# compute the coverage
rule plotCoverage_merged_peaks:
    input:
        target_bam_bsqr = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bam"])),
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed")
    output:
        coverage_plot = os.path.join(GATKDIR, "coverage_plots/{SAMPLES}_plotCoverage.pdf")
    params:
        label = "{SAMPLES}",
        blacklist = BLACKLIST
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "coverage_plots/logs/{SAMPLES}_plotCoverage.log")
    benchmark:
        "benchmarks/plotCoverage_merged_peaks/plotCoverage_merged_peaks---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.label}: Plotting coverage at ChIP (merged) peaks...\\n\033[0m'

        $CONDA_PREFIX/bin/plotCoverage \
        --bamfiles {input.target_bam_bsqr} \
        --labels {params.label} \
        --BED {input.concatenation_bed_collapsed_sorted} \
        --blackListFileName {params.blacklist} \
        --plotFile {output.coverage_plot} \
        --plotTitle {params.label} \
        --numberOfProcessors {threads} &> {log.out}
        """

# ----------------------------------------------------------------------------------------

# run gatk haplotype caller
rule GATK_haplotype_calling:
    input:
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"),
        target_bam_bsqr = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bam"])),
        target_bam_bsqr_index = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{SAMPLES}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bai"]))
    params:
        genome = genome_fasta,
        sample = "{SAMPLES}",
        to_copy_bed = os.path.join(GATKDIR, "all_samples_peak_concatenation_collapsed_sorted.bed")
    output:
        gvcf = temp(os.path.join(GATKDIR, ''.join(["VCF/{SAMPLES}_", DUP, "_gatk.g.vcf.gz"]))),
        gvcf_idx = temp(os.path.join(GATKDIR, ''.join(["VCF/{SAMPLES}_", DUP, "_gatk.g.vcf.gz.tbi"])))
    threads:
        max(math.floor(workflow.cores/2), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/logs/{SAMPLES}_HaplotypeCaller.log")
    benchmark:
        "benchmarks/GATK_haplotype_calling/GATK_haplotype_calling---{SAMPLES}_benchmark.txt"
    shell:
        """
        cp {input.concatenation_bed_collapsed_sorted} {params.to_copy_bed}

        printf '\033[1;36m{params.sample}: GATK Haplotype calling...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk HaplotypeCaller \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -I {input.target_bam_bsqr} \
        -O {output.gvcf} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation &> {log.out}
        """


# correct the genotypes that come out of haplotype caller
rule GATK_haplotype_calling_correction:
    input:
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"),
        gvcf = os.path.join(GATKDIR, ''.join(["VCF/{SAMPLES}_", DUP, "_gatk.g.vcf.gz"])),
        gvcf_idx = os.path.join(GATKDIR, ''.join(["VCF/{SAMPLES}_", DUP, "_gatk.g.vcf.gz.tbi"]))
    params:
        sample = "{SAMPLES}",
        genome = genome_fasta
    output:
        vcf = os.path.join(GATKDIR, ''.join(["VCF/{SAMPLES}_", DUP, "_gatk.vcf.gz"]))
    threads:
        max(math.floor(workflow.cores/4), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/logs/{SAMPLES}_GenotypeGVCFs.log")
    benchmark:
        "benchmarks/GATK_haplotype_calling_correction/GATK_haplotype_calling_correction---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Haplotype call correction...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk GenotypeGVCFs \
        --include-non-variant-sites \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.gvcf} \
        -O {output.vcf} &> {log.out}
        """

# ----------------------------------------------------------------------------------------


# Call SNPs
rule GATK_call_SNPs:
    input:
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"),
        vcf = os.path.join(GATKDIR, ''.join(["VCF/{SAMPLES}_", DUP, "_gatk.vcf.gz"]))
    params:
        sample = "{SAMPLES}",
        genome = genome_fasta
    output:
        snp = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp.vcf"]))),
        snp_idx = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp.vcf.idx"])))
    threads:
        max(math.floor(workflow.cores/4), 1)
    log:
        out = os.path.join(GATKDIR, ''.join(["VCF/SNP/logs/{SAMPLES}_SNP_SelectVariants.log"]))
    benchmark:
        "benchmarks/GATK_call_SNPs/GATK_call_SNPs---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK SNP calling...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type SNP \
        --select-type NO_VARIATION \
        --select-type-to-exclude INDEL \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.snp} &> {log.out}
        """


# Filter SNPs
rule GATK_filter_SNPs:
    input:
        snp_vcf = os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp.vcf"]))
    output:
        filtered_snp_vcf = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), ".vcf"])))
    params:
        sample = "{SAMPLES}",
        DP_snp_threshold = config["somatic_variants"]["DP_snp_threshold"],
        QUAL_snp_threshold = config["somatic_variants"]["QUAL_snp_threshold"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/GATK_filter_SNPs/GATK_filter_SNPs---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift SNP filtering...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift filter '( DP > {params.DP_snp_threshold} & ( QUAL > {params.QUAL_snp_threshold} ))' {input.snp_vcf} > {output.filtered_snp_vcf}
        """


# Annotate SNP
rule SnpSift_annotate_SNPs:
    input:
        filtered_snp_vcf = os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), ".vcf"]))
    output:
        filtered_snp_vcf_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf"]))),
        filtered_snp_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf.gz"])),
        filtered_snp_vcf_gz_annotated_idx = os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf.gz.tbi"]))
    params:
        sample = "{SAMPLES}",
        dbsnp = config["somatic_variants"]["dbsnp_file"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/SNP/logs/{SAMPLES}_SnpSift_annotate_SNPs.log")
    benchmark:
        "benchmarks/SnpSift_annotate_SNPs/SnpSift_annotate_SNPs---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift SNP annotation...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift annotate -a {params.dbsnp} {input.filtered_snp_vcf} > {output.filtered_snp_vcf_annotated} 2> {log.out}

        printf '\033[1;36m{params.sample}: Filtered SNP vcf bgzipping...\\n\033[0m'
        $CONDA_PREFIX/bin/bgzip -@ {threads} -c {output.filtered_snp_vcf_annotated} > {output.filtered_snp_vcf_annotated_gz}
        $CONDA_PREFIX/bin/tabix -p vcf {output.filtered_snp_vcf_annotated_gz}
        """


# Export SNPs table
rule GATK_vcf2txt_SNPs:
    input:
        filtered_snp_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf.gz"]))
    output:
        filtered_snp_allGT_tb_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_allGT_annotated.txt"]))),
        filtered_snp_tb_annotated = os.path.join(GATKDIR, ''.join(["VCF/SNP/{SAMPLES}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"]))
    params:
        sample = "{SAMPLES}",
        fileds = ''.join(['"', '" "'.join(config["somatic_variants"]["SnpSift_vcf_fields_to_extract"]), '"'])
    threads:
        max(math.floor(workflow.cores/5), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/SNP/logs/{SAMPLES}_GATK_vcf2txt_SNPs.log")
    benchmark:
        "benchmarks/GATK_vcf2txt_SNPs/GATK_vcf2txt_SNPs---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Filtered SNP export to txt...\\n\033[0m'

        $CONDA_PREFIX/bin/SnpSift extractFields \
        -s "," \
        -e "." \
        {input.filtered_snp_vcf_annotated_gz} \
        {params.fileds} > {output.filtered_snp_allGT_tb_annotated} 2> {log.out}

        grep -v '0/0' {output.filtered_snp_allGT_tb_annotated} > {output.filtered_snp_tb_annotated}
        """


# Merge all SNPs tables
rule GATK_merge_SNP_tables:
    input:
        snp_txt = expand(os.path.join(GATKDIR, ''.join(["VCF/SNP/{sample}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"])), sample=SAMPLENAMES)
    params:
        sample = SAMPLENAMES,
        gatk_dir = GATKDIR,
        suffix_tb = ''.join(["_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"])
    output:
        merged_snps = os.path.join(GATKDIR, ''.join(["VCF/SNP/all.samples_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"]))
    threads: 1
    benchmark:
        "benchmarks/GATK_merge_SNP_tables/GATK_merge_SNP_tables---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mMerging all SNP txt tables in one...\\n\033[0m'")

        import pandas as pd

        merged_table = pd.DataFrame()

        for s in params.sample:
            if os.path.getsize(''.join([GATKDIR, "/VCF/SNP/", s, params.suffix_tb])) > 0:
                tb = pd.read_table(''.join([GATKDIR, "/VCF/SNP/", s, params.suffix_tb]))
                tb.insert(0, "sample_ID", s, allow_duplicates=True)
                merged_table = pd.concat([merged_table, tb], ignore_index = True, sort = False)

        merged_table.to_csv(output.merged_snps, encoding="utf-8", index=False, header=True, sep="\t")



# Plot SNP occurences
rule GATK_plot_SNPs:
    input:
        merged_SNPs = os.path.join(GATKDIR, ''.join(["VCF/SNP/all.samples_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"]))
    params:
        sample = SAMPLENAMES,
        plot_title = ''.join(["SNPs (on ", DUP, " bams): DP > ", str(config["somatic_variants"]["DP_snp_threshold"]), ", QUAL > ", str(config["somatic_variants"]["QUAL_snp_threshold"]), ", w\o 0|0"])
    output:
        snp_plot = os.path.join(GATKDIR, "SV_count_plots/all.samples_SNP_counts_plot.pdf")
    threads: 1
    benchmark:
        "benchmarks/GATK_plot_SNPs/GATK_plot_SNPs---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mPlotting SNP occurences per sample...\\n\033[0m'")

        import pandas as pd
        import numpy as np
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        tb = pd.read_table(input.merged_SNPs)

        occurences = []
        for s in np.flip(params.sample):
            occurences.append(len(tb[tb["sample_ID"] == s].index))

        occurences_tb = pd.DataFrame({'Sample':np.flip(params.sample), 'SNP.counts':occurences})
        plot = occurences_tb.plot.barh(x='Sample', y='SNP.counts', title = params.plot_title)
        plot.figure.savefig(output.snp_plot, bbox_inches='tight')


# ----------------------------------------------------------------------------------------


# Call InDels
rule GATK_call_InDels:
    input:
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"),
        vcf = os.path.join(GATKDIR, ''.join(["VCF/{SAMPLES}_", DUP, "_gatk.vcf.gz"]))
    output:
        indels = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel.vcf"]))),
        indels_idx = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel.vcf.idx"])))
    params:
        sample = "{SAMPLES}",
        genome = genome_fasta
    threads:
        max(math.floor(workflow.cores/2), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{SAMPLES}_InDel_SelectVariants.log")
    benchmark:
        "benchmarks/GATK_call_InDels/GATK_call_InDels---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Indel calling...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type INDEL \
        --select-type NO_VARIATION \
        --select-type-to-exclude SNP \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.indels} &> {log.out}
        """


# Filter InDels
rule GATK_filter_InDels:
    input:
        indel_vcf = os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel.vcf"]))
    output:
        filtered_indel_vcf = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), ".vcf"])))
    params:
        sample = "{SAMPLES}",
        DP_indel_threshold = config["somatic_variants"]["DP_indel_threshold"],
        QUAL_indel_threshold = config["somatic_variants"]["QUAL_indel_threshold"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{SAMPLES}_filter_InDels.log")
    benchmark:
        "benchmarks/GATK_filter_InDels/GATK_filter_InDels---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift InDel filtering...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift filter '( (DP > {params.DP_indel_threshold}) & ( QUAL > {params.QUAL_indel_threshold} ))' {input.indel_vcf} > {output.filtered_indel_vcf} 2> {log.out}
        """


# Annotate InDels
rule SnpSift_annotate_InDels:
    input:
        filtered_indel_vcf = os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), ".vcf"]))
    output:
        filtered_indel_vcf_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf"]))),
        filtered_indel_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf.gz"])),
        filtered_indel_vcf_gz_annotated_idx = os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf.gz.tbi"]))
    params:
        sample = "{SAMPLES}",
        dbsnp = config["somatic_variants"]["dbsnp_file"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{SAMPLES}_SnpSift_annotate_InDels.log")
    benchmark:
        "benchmarks/SnpSift_annotate_InDels/SnpSift_annotate_InDels---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift InDel annotation...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift annotate -a {params.dbsnp} {input.filtered_indel_vcf} > {output.filtered_indel_vcf_annotated} 2> {log.out}

        printf '\033[1;36m{params.sample}: Filtered InDel vcf bgzipping...\\n\033[0m'
        $CONDA_PREFIX/bin/bgzip -@ {threads} -c {output.filtered_indel_vcf_annotated} > {output.filtered_indel_vcf_annotated_gz}
        $CONDA_PREFIX/bin/tabix -p vcf {output.filtered_indel_vcf_annotated_gz}
        """


# Export InDels table
rule GATK_vcf2txt_InDels:
    input:
        filtered_indel_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf.gz"]))
    output:
        filtered_indel_allGT_tb_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_allGT_annotated.txt"]))),
        filtered_indel_tb_annotated = os.path.join(GATKDIR, ''.join(["VCF/InDel/{SAMPLES}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"]))
    params:
        sample = "{SAMPLES}",
        fileds = ''.join(['"', '" "'.join(config["somatic_variants"]["SnpSift_vcf_fields_to_extract"]), '"'])
    threads:
        max(math.floor(workflow.cores/5), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{SAMPLES}_GATK_vcf2txt_InDels.log")
    benchmark:
        "benchmarks/GATK_vcf2txt_InDels/GATK_vcf2txt_InDels---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Filtered InDel export to txt...\\n\033[0m'

        $CONDA_PREFIX/bin/SnpSift extractFields \
        -s "," \
        -e "." \
        {input.filtered_indel_vcf_annotated_gz} \
        {params.fileds} > {output.filtered_indel_allGT_tb_annotated} 2> {log.out}

        grep -v '0/0' {output.filtered_indel_allGT_tb_annotated} > {output.filtered_indel_tb_annotated}
        """



# Merge all InDels tables
rule GATK_merge_InDel_tables:
    input:
        indel_txt_annotated = expand(os.path.join(GATKDIR, ''.join(["VCF/InDel/{sample}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"])), sample=SAMPLENAMES)
    output:
        merged_indels_annotated = os.path.join(GATKDIR, ''.join(["VCF/InDel/all.samples_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"]))
    params:
        sample = SAMPLENAMES,
        gatk_dir = GATKDIR,
        suffix_tb = ''.join(["_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"])
    threads: 1
    benchmark:
        "benchmarks/GATK_merge_InDel_tables/GATK_merge_InDel_tables---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mMerging all InDel txt tables in one...\\n\033[0m'")

        import pandas as pd

        merged_table = pd.DataFrame()

        for s in params.sample:
            if os.path.getsize(''.join([GATKDIR, "/VCF/InDel/", s, params.suffix_tb])) > 0:
                tb = pd.read_table(''.join([GATKDIR, "/VCF/InDel/", s, params.suffix_tb]))
                tb.insert(0, "sample_ID", s, allow_duplicates=True)
                merged_table = pd.concat([merged_table, tb], ignore_index = True, sort = False)

        merged_table.to_csv(output.merged_indels_annotated, encoding="utf-8", index=False, header=True, sep="\t")



# Plot InDels occurences
rule GATK_plot_InDels:
    input:
        merged_indels_annotated = os.path.join(GATKDIR, ''.join(["VCF/InDel/all.samples_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"]))
    params:
        sample = SAMPLENAMES,
        plot_title = ''.join(["InDels (on ", DUP, " bams): DP > ", str(config["somatic_variants"]["DP_indel_threshold"]), ", QUAL > ", str(config["somatic_variants"]["QUAL_indel_threshold"]), ", w\o 0|0"])
    output:
        indel_plot = os.path.join(GATKDIR, "SV_count_plots/all.samples_InDel_counts_plot.pdf")
    threads: 1
    benchmark:
        "benchmarks/GATK_plot_InDels/GATK_plot_InDels---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mPlotting indel occurences per sample...\\n\033[0m'")

        import pandas as pd
        import numpy as np
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        tb = pd.read_table(input.merged_indels_annotated)

        occurences = []
        for s in np.flip(params.sample):
            occurences.append(len(tb[tb["sample_ID"] == s].index))

        occurences_tb = pd.DataFrame({'Sample':np.flip(params.sample), 'InDel.counts':occurences})
        plot = occurences_tb.plot.barh(x='Sample', y='InDel.counts', title = params.plot_title)
        plot.figure.savefig(output.indel_plot, bbox_inches='tight')



# ------------------------------------------------------------------------------
#                                 END PIPELINE                                 #
# ------------------------------------------------------------------------------
