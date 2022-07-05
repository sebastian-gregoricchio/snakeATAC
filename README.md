[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
![release](https://img.shields.io/github/v/release/sebastian-gregoricchio/snakeATAC)
![update](https://badges.pufler.dev/updated/sebastian-gregoricchio/snakeATAC)
![visits](https://badges.pufler.dev/visits/sebastian-gregoricchio/snakeATAC)
[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/snakeATAC/LICENSE.md/LICENSE)
[![forks](https://img.shields.io/github/forks/sebastian-gregoricchio/snakeATAC?style=social)](https://github.com/sebastian-gregoricchio/snakeATAC/fork)
<!---![downloads](https://img.shields.io/github/downloads/sebastian-gregoricchio/Rseb/total.svg)--->

# snakeATAC <img src="https://sebastian-gregoricchio.github.io/snakeATAC/resources/snakeATAC_logo.svg" align="right" height = 150/>
## Introduction
`SnakeATAC` is a snakemake based end-to-end pipeline to analyze ATAC-seq data. The input files required to run the pipeline are Paired-End fastq files. The pipeline include data quality check and normalization. It is included also a step of data reads shifting in order to take into account the Tn5 transposome insertion bias.

<br></br>

## Installation an dependencies
To install the pipeline it is required to download this repository and the installation of a conda environment is strongly suggested.
Follow the following steps for the installation:
* Place yourself in the directory where the repository should be downloaded with `cd </target/folder>`
* download the GitHub repository `git clone https://github.com/sebastian-gregoricchio/snakeATAC`
* install the conda environment from the env .yaml file contained in the repository:

`conda env create -f </target/folder>/snakeATAC/workflow/envs/snakeATAC_conda_env_stable.yamll`
* activate the environment: `conda activate snakeATAC` (if the env is not activated the pipeline won't work)

<br></br>

## How to run the pipeline
The snakemake pipeline requires only two files: a) the `.snakefile`, containing all the rules that will be run; b) the `configuration.yaml` file, in which the user can define and customize all the parameters for the different pipeline steps.
To partially avoid unexpected errors during the execution of the pipeline, a so called 'dry-run' is strongly recommended. Indeed, adding a `-n` at the end of the snakemake running command will allow snakemake to check that all links and file/parameters dependencies are satisfied before to run the "real" processes. This command will therefore help the debugging process.

```
snakemake -s </target/folder>/snakeATAC/workflow/snakeATAC.snakefile --configfile </target/folder>/snakeATAC/config/snakeATAC_config.yaml --cores 12 -n
```

If no errors occur, the pipeline can be run with the same command line without the terminal `-n`:

```
snakemake -s </target/folder>/snakeATAC/workflow/snakeATAC.snakefile --configfile </target/folder>/snakeATAC/config/snakeATAC_config.yaml --cores 12
```

Notice that the absence of errors does not mean that the pipeline will run without any issues; the "dry-run" is only checking whether all the resources are available.


### Snakefile
Explain the rules  + workflow



### Configuration file
The configuration file is a yaml-formatted file containing all the parameters that are passed to different steps of the pipelines such as the directory with the input files, reference genome, threads of the tools, etc.
The snakeATAC configuration file is divided in two sections: a) 'experiment-specific', with al the parameters that most likely are changing from experiment to experiment; b) 'common', with parameters that are quite stable independently of the experiments design. The latter should be changed only for very specific needs and is in turn compose by two sections depending on whether the copy number variation is performed or not.
Hereafter, the meaning of the different parameters is described.

<br></br>

#### Experiment-specific section
| Parameter   |   Description   |
|------------:|:----------------|
|*runs_directory*| The full path to the directory were the input fastq files are contained, e.g. `/home/user/ATAC/00_runs/`. Importantly, the name of the files, deprived of the read suffix (e.g., _R1/_R2) and file extension (e.g., .fastq.gz) will be used as sample name.|
|*output_directory*| The full path to the folder in which the results should be stored, e.g. `"/home/user/ATAC/"`. |
|*fastq_extension*| The extension of the input fastq files. This string will be removed from the input file names to obtain the sample names. Examples: `".fastq.gz"`, `".fq.gz"`, `".fasta"`, `".fq"`. |
|*runs_suffix*| A list (python format) with the two reads suffixes corresponding to the read1 and read2 for the paired-end sequencing. Example: `["_R1", "_R2"]`. |
|*blacklist_file*| The full path to a BED-formatted file containing the regions corresponding to blacklisted regions (regions masked for normalization and peak calling). Blacklisted regions for the most common genome assemblies can be downloaded from the [Boyle lab's git-page](https://github.com/Boyle-Lab/Blacklist/). |
|*genome_fasta*| The full path to a fasta-formatted file containing the sequence of there reference genome into which the reads will be aligned. If the index (.fai) file is not present in the same folder, it will be build by the pipeline during its execution. The reference genomes can be downloaded, for instance, from the [UCSC golden path](https://hgdownload.soe.ucsc.edu/goldenPath/). |
|*perform_HMCan_correction*| Default `"False"`. If set to `"True"` the signal will be corrected for the presence of CNVs (Copy Number Variations) and the peak calling will be performed by [HMCan](https://academic.oup.com/bioinformatics/article/29/23/2979/246862?login=false). This correction could be useful in the case of heterogeneous cancer samples. Instead, if set to `"False"`, the CNVs correction is skipped and the peak calling performed by MACS3. **The CNV correction implementation is still in beta-testing phase.**
|*effective_genomeSize*| Effective genome size used to normalize the ATAC-seq signal; e.g. for Hg38: `2900338458`. A table for the most common genome assemblies is available in this repository at [snakeATAC/resources/chromosome_sizes_for_normalization.yaml](https://github.com/sebastian-gregoricchio/snakeATAC/blob/main/resources/chromosomes_sizes_for_normalization.yaml). |
|*ignore_for_normalization*| A string of space separated chromosome names that should be excluded for ATAC-seq signal normalization, e.g. `"X Y MT M KI270728.1"`.  A table for the most common genome assemblies is available in this repository at [snakeATAC/resources/chromosome_sizes_for_normalization.yaml](https://github.com/sebastian-gregoricchio/snakeATAC/blob/main/resources/chromosomes_sizes_for_normalization.yaml). |
|*call_variants*| If `true`, variant calling by [GATK4](https://gatk.broadinstitute.org/hc/en-us) will be performed. **Variant calling is still in beta-test phase.** |
|*call_SNPs*| If `true`, Single Nucleotide Variation (SNP) calling by [GATK4](https://gatk.broadinstitute.org/hc/en-us) will be performed. **Variant calling is still in beta-test phase.** |
|*call_indels*| If `true`, Insertion/Deletion (indel) calling by [GATK4](https://gatk.broadinstitute.org/hc/en-us) will be performed. **Variant calling is still in beta-test phase.** |
|*dbsnp_file*| SNP database file (.dbsnp) for base recalibration required by [GATK4](https://gatk.broadinstitute.org/hc/en-us). It could happen that your .bam files contain the 'chr' prefix in the chromosome names while your dbSNP file does not (or viceversa). This can be fixed in the .dbsnp file with the [`bcftools annotate --rename-chrs`](http://samtools.github.io/bcftools/bcftools.html#annotate) command. For Hg38, for istance, the dbSNP file can be downloaded from the [broad institute cloud storage](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/). Do not forget to download the INDEX as well!|

<br></br>

#### Common section
| Parameter   |   Description   |
|------------:|:----------------|
|*fastQC_threads*| Default: `2`. Number of CPUs to use for [fastQC](https://github.com/s-andrews/FastQC) (fastq quality control) |
|*bwa_threads*| Default: `8`. Number of CPUs to use for the mapping performed by [bwa-mem](http://bio-bwa.sourceforge.net/bwa.shtml). |
|*mapQ_cutoff*| Default: `20`. All reads with a mapping quality (MAPQ) score lower than this value will be filtered out from the bam files. |
|*SAMtools_threads*| Default: `8`. Number of CPUs used by [samtools](http://www.htslib.org/doc/samtools.html) for bam indexing and filtering. |
|*remove_duplicates*| Default: `"true"`. Logical value to indicate whether the optical duplicates should be removed from the bams by [PICARD](https://broadinstitute.github.io/picard/). If set as `"true"` the suffix of output bam files will be '*_dedup*' and the duplicates will be removed, otherwise the suffix will be '*_mdup*' and the duplicates only marked when set to `"false"`.
|*PICARD_max_records_in_ram*| Default: `250000` This will specify the number of records stored in RAM before spilling to disk. The higher the number the higher the amount of RAM needed. |
|*PICARD_max_file_handles_for_read_ends_map*| Default: `4000` Maximum number of file handles to keep open when spilling read ends to disk. This number should be slightly lower than the per-process maximum number of file that may be open in your system (`ulimit -n`). |
|*minFragmentLength*| Default: `0`. Minimum fragment length needed for pair inclusion during Tn5 tagmentation bias correction by [deeptools alignmentSieve](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html) tool. |
|*maxFragmentLength*| Default: `0`. Maximum fragment length needed for pair inclusion during Tn5 tagmentation bias correction by [deeptools alignmentSieve](https://deeptools.readthedocs.io/en/develop/content/tools/alignmentSieve.html) tool. A value of 0 indicates no limit. |
|*plot_format*| Default: `"pdf"`. File format of the images of the fragment size distribution generated by [deeptools bamPEFragmentSize](https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html). |
|*bamPEFragmentSize_threads*| Default: `0`. Number of CPUs [deeptools bamPEFragmentSize](https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html) should use to compute the fragment size distribution. |
|*max_fragment_length*| Default: `2000`. Maximum fragment size length to be plotted by [deeptools bamPEFragmentSize](https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html). A value of 0 indicates to use twice the mean fragment length. |
|*window_length*| Default: `1000`. Size, in bp, of the sliding window used to compute the fragment size distribution by [deeptools bamPEFragmentSize](https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html). |
|*plotFingerprint_threads*| Default: `8`. Number of CPUs to be used by [deeptools plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html) to compute the Lorenz curves. |
|*plotFingerprint_binSize*| Default: `500`. Size of the bins, in bp, by which the genome should be subdivided in order to compute the Lorenz curves by [deeptools plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html). |
|*plotFingerprint_sampledRegions*| Default: `500000`. Number of regions to be samples in order to compute the Lorenz curves by [deeptools plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html). |
|*plotFingerprint_extra_parameters*| Default: `""` (empty). A string with any additional parameter to add for the Lorenz curve computation performed by [deeptools plotFingerprint](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html). |
|*binning_window_size*| Default: `10000`. Size of the bins, in bp, by which the genome should be subdivided in order to compute the average score matrix by [deeptools multiBigwigSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html). This matrix is used to compute sample correlation and Principal Component Analyses (PCA). |
|*multiBigwigSummary_threads*| Default: `4`. Number of CPUs to be used by [deeptools multiBigwigSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html) in order to compute an average signal over all the genome for each sample. The resulting matrix is used to compute sample correlation and Principal Component Analyses (PCA). |
|*heatmap_color*| Default: `"Blues"`. A string indicating the color gradient pattern to use for the correlation heatmaps. This value is passed to matplotlib/seaborn. Therefore, available options (see [matplotlib page](https://matplotlib.org/stable/tutorials/colors/colormaps.html) for examples) are the following: 'Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'Dark2', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired', 'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cividis', 'cool', 'coolwarm', 'copper', 'cubehelix', 'flag', 'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'icefire', 'inferno', 'jet', 'magma', 'mako', 'nipy_spectral', 'ocean', 'pink', 'plasma', 'prism', 'rainbow', 'rocket', 'seismic', 'spring', 'summer', 'tab10', 'tab20', 'tab20b', 'tab20c', 'terrain', 'twilight', 'twilight_shifted', 'viridis', 'vlag', 'winter'. |
|*zScore_heatmap_color*| Default: `"seismic"`. A string indicating the color gradient pattern to use for the peak score heatmaps. This value is passed to matplotlib/seaborn. Therefore, available options (see [matplotlib page](https://matplotlib.org/stable/tutorials/colors/colormaps.html) for examples) are the following: 'Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'Dark2', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired', 'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cividis', 'cool', 'coolwarm', 'copper', 'cubehelix', 'flag', 'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'icefire', 'inferno', 'jet', 'magma', 'mako', 'nipy_spectral', 'ocean', 'pink', 'plasma', 'prism', 'rainbow', 'rocket', 'seismic', 'spring', 'summer', 'tab10', 'tab20', 'tab20b', 'tab20c', 'terrain', 'twilight', 'twilight_shifted', 'viridis', 'vlag', 'winter'. |

<br></br>

*Copy Number Variation signal correction* (for details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/))
| Parameter   |   Description   |
|------------:|:----------------|
|*HMCan_path*| Full path to the folder containing the HMCan scripts. Instruction for the download and build of HMCan can be found at the the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*HMCan_threads*| Default: `12`. Number of CPUs to be used by [HMCan](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*chromosome_sizes_file*| Full path to the chromosome sizes file used to convert the CNV corrected signal into a bigWig. These files can be downloaded from the [UCSC golden path](https://hgdownload.soe.ucsc.edu/goldenPath/). |
|*reference_sample*| Default: `"NA"`. A string to define the reference sample (sample ID) to which all the samples should be normalized. `"NA"` indicates that the first samples in alphabetic order will be used as reference. |
|*format*| Default: `"BAM"` (used to build an HMCan configFile).
|*GCIndex*| (used to build an HMCan configFile) Full path to the GC index file provided at the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). Example: "/home/user/HMCan-master/data/GC_profile_25KbWindow_Mapp100bp_hg38.cnp". |
|*smallBinLength*| Default: `5` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*largeBinLength*| Default: `25000` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*pvalueThreshold| Default: `0.001` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*mergeDistance*| Default: `0` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*iterationThreshold*| Default: `5` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*finalThreshold*| Default: `0` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*maxIter*| Default: `20` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*PrintWig*| Default: `"False"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*PrintPosterior*| Default: `"True"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*PrintBedGraph*| Default: `"True"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*CallPeaks*| Default: `"True"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*pairedEnds*| Default: `"True"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*Type*| Default: `"ATAC-seq"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*GCmergeDistance*| Default: `1000` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*RemoveDuplicates*| Default: `"False"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*CNAnormalization*| Default: `"True"` (used to build an HMCan configFile). For details see the [HMCan page](https://bitbucket.org/pyminer/hmcan/src/master/). |
|*multiBamSummary_threads*| Default: `6`. Number of CPUs to be used by [deeptools multiBamSummary](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html) in order to compute an average signal over all the genome for each sample for the calculation of the scaling factors. This factors will be used to normalize the CNV corrected signal by sequencing depth.|

<br></br>

*Standard normalization and peak calling (without CNV correction)*
| Parameter   |   Description   |
|------------:|:----------------|
|*bigWig_binSize*| Default: `5`. Size, in bp, of the bins used to compute the normalized bigWig files. |
|*normalization_method*| Default: `"RPGC"`, reads per genomic content (1x normalization). Type of normalization to be used to generated the normalized bigWig files by [deeptools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html). |
|*bamCoverage_threads*| Default: `8`. Number of CPUs to be used by [deeptools bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) in order to compute the signal normalization and the bigWig generation. |
|*genome_size_MACS*| A string or a number indicating the genome size to use by [MACS3](https://github.com/macs3-project/MACS) for the peak calling. Example: `"hs"`. Some values: hs = 2.7e9, mm = 1.87e9, ce = 9e7, dm = 1.2e8. |
|*FDR_cutoff*| Deafult: `0.01`. False Discovery Ratio (FDR) cutoff used by [MACS3](https://github.com/macs3-project/MACS) to filter the significant peaks. |
|*call_summits*| Default: `"True"`. Logical value to define whether [MACS3](https://github.com/macs3-project/MACS) should also call the peak summits (position with the highest value). |
|*FRiP_threshold*| Default `20`. This value will be used to label the FRiP (Fraction of Reads in a Peak) score as "good" or "bad in the summary table for each single sample. A FRiP above the 20% (FRiP = 0.02) is considered a good score for ATAC-seq peaks by the [ENCODE guidelines](https://www.encodeproject.org/atac-seq/). |



<br></br>


## Results
Results interpretation and structure output

summmary file





<br></br>



# Changing names
```
for i in $(dir *f*.gz)
do
  R1=$(echo $i | sed 's/_1/_R1/')
  R2=$(echo $R1 | sed 's/_2/_R2/')
  mv $i $R2
done
```


<br></br>

-----------------
## Contact
For any suggestion, bug fixing, commentary please report it in the [issues](https://github.com/sebastian-gregoricchio/snakeATAC/issues)/[request](https://github.com/sebastian-gregoricchio/snakeATAC/pulls) tab of this repository.

## License
This repository is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE).

<br />

#### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
