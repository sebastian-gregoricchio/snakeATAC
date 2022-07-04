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
The snakeATAC configuration file is divided in two sections: a) 'experiment-specific', with al the parameters that most likely are changing from experiment to experiment; b) 'constant', with parameters that are quite stable independently of the experiments design. The latter should be changed only for very specific needs.
Hereafter, the meaning of the different parameters is described.

**Experiment-specific section**
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
|*dbsnp_file*| SNP database file (.dbsnp) for base recalibration required by [GATK4](https://gatk.broadinstitute.org/hc/en-us). It could happen that your .bam files contain the 'chr' prefix in the chromosome names while your dbsnp file does not (or viceversa). This can be fixed in the .dbsnp file with the [`bcftools annotate --rename-chrs`](http://samtools.github.io/bcftools/bcftools.html#annotate) command.


<br></br>


## Results
The snakemake pipeline r





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
For any suggestion, bug fixing, commentary please report it in the [issues](https://github.com/sebastian-gregoricchio/snakeATAC/issues)/(request)[https://github.com/sebastian-gregoricchio/snakeATAC/pulls] tab of this repository.

## License
This repository is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE).

<br />

#### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
