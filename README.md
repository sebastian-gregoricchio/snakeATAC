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

If no errors occur, the pipeline can be run with the same command line without the terminal `-n`. Notice that the absence of errors does not mean that the pipeline will run without any issues; the "dry-run" is only checking whether all the resources are available.


### Configuration file
The configuration file is a yaml-formatted file containing all the parameters that are passed to different steps of the pipelines such as the directory with the input files, reference genome, threads of the tools, etc.
The snakeATAC configuration file is divided in two sections: a) 'experiment-specific', with al the parameters that most likely are changing from experiment to experiment; b) 'constant', with parameters that are quite stable independently of the experiments design. The latter should be changed only for very specific needs.
Hereafter, the meaning of the different parameters is described.

**Experiment-specific section**
| Parameter   |   Description   |
|------------:|:----------------|
|runs_directory| "/home/user/ATAC_test/00_runs/"|
|output_directory| "/home/user/ATAC_test/" |
|fastq_extension| ".fastq.gz" |
|runs_suffix| ["_R1", "_R2"]|
|blacklist_file| "/home/user/annotations/blacklist/hg38-blacklist.v2.bed" # look at https://github.com/Boyle-Lab/Blacklist/ |
|genome_fasta|"/home/user/annotations/genomes/Hg38/Homo_sapiens.GRCh38_v102.dna.primary_assembly.fa" |
|perform_HMCan_correction| "False"
      # If 'True' the signal will be corrected for the presence of CNVs
      # and the peak calling will be performed by HMCan instead of MACS3 --> set the correct genomeSize paramater
      #
      # If 'False' the peak calling will be performed by MACS2 --> set the correct genomeSize paramater|
|effective_genomeSize| 2900338458 |
|ignore_for_normalization| "X Y MT M KI270728.1 KI270727.1 KI270442.1 KI270729.1 GL000225.1 KI270743.1 GL000008.2 GL000009.2 KI270747.1 KI270722.1 GL000194.1 KI270742.1 GL000205.2 GL000195.1 KI270736.1 KI270733.1 GL000224.1 GL000219.1 KI270719.1 GL000216.2 KI270712.1 KI270706.1 KI270725.1 KI270744.1 KI270734.1 GL000213.1 GL000220.1 KI270715.1 GL000218.1 KI270749.1 KI270741.1 GL000221.1 KI270716.1 KI270731.1 KI270751.1 KI270750.1 KI270519.1 GL000214.1 KI270708.1 KI270730.1 KI270438.1 KI270737.1 KI270721.1 KI270738.1 KI270748.1 KI270435.1 GL000208.1 KI270538.1 KI270756.1 KI270739.1 KI270757.1 KI270709.1 KI270746.1 KI270753.1 KI270589.1 KI270726.1 KI270735.1 KI270711.1 KI270745.1 KI270714.1 KI270732.1 KI270713.1 KI270754.1 KI270710.1 KI270717.1 KI270724.1 KI270720.1 KI270723.1 KI270718.1 KI270317.1 KI270740.1 KI270755.1 KI270707.1 KI270579.1 KI270752.1 KI270512.1 KI270322.1 GL000226.1 KI270311.1 KI270366.1 KI270511.1 KI270448.1 KI270521.1 KI270581.1 KI270582.1 KI270515.1 KI270588.1 KI270591.1 KI270522.1 KI270507.1 KI270590.1 KI270584.1 KI270320.1 KI270382.1 KI270468.1 KI270467.1 KI270362.1 KI270517.1 KI270593.1 KI270528.1 KI270587.1 KI270364.1 KI270371.1 KI270333.1 KI270374.1 KI270411.1 KI270414.1 KI270510.1 KI270390.1 KI270375.1 KI270420.1 KI270509.1 KI270315.1 KI270302.1 KI270518.1 KI270530.1 KI270304.1 KI270418.1 KI270424.1 KI270417.1 KI270508.1 KI270303.1 KI270381.1 KI270529.1 KI270425.1 KI270396.1 KI270363.1 KI270386.1 KI270465.1 KI270383.1 KI270384.1 KI270330.1 KI270372.1 KI270548.1 KI270580.1 KI270387.1 KI270391.1 KI270305.1 KI270373.1 KI270422.1 KI270316.1 KI270338.1 KI270340.1 KI270583.1 KI270334.1 KI270429.1 KI270393.1 KI270516.1 KI270389.1 KI270466.1 KI270388.1 KI270544.1 KI270310.1 KI270412.1 KI270395.1 KI270376.1 KI270337.1 KI270335.1 KI270378.1 KI270379.1 KI270329.1 KI270419.1 KI270336.1 KI270312.1 KI270539.1 KI270385.1 KI270423.1 KI270392.1 KI270394.1" |








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
