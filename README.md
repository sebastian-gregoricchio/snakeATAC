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


## Installation an dependencies
To install the pipeline it is required to download this repository and the installation of a conda environment is strongly suggested.
Follow the following steps for the installation:
* Place yourself in the directory where the repository should be downloaded with `cd </target/folder>`
* download the GitHub repository `git clone https://github.com/sebastian-gregoricchio/snakeATAC`
* install the conda environment from the env .yaml file contained in the repository:
```conda env create -f </target/folder>/snakeATAC/workflow/envs/snakeATAC_conda_env_stable.yamll```


# Changing names
```
for i in $(dir *f*.gz)
do
  R1=$(echo $i | sed 's/_1/_R1/')
  R2=$(echo $R1 | sed 's/_2/_R2/')
  mv $i $R2
done
```



-----------------
## Contact
For any suggestion, bug fixing, commentary please contact Sebastian Gregoricchio at [sebastian.gregoricchio@gmail.com](mailto:sebastian.gregoricchio@gmail.com).

## License
This repository is under a [GNU General Public License (version 3)](https://sebastian-gregoricchio.github.io/Rseb/LICENSE.md/LICENSE).

<br />

#### Contributors
[![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/Rseb?size=50&padding=5&bots=true)](https://sebastian-gregoricchio.github.io/)
