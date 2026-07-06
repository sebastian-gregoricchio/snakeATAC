[![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://sebastian-gregoricchio.github.io/snakeATAC/LICENSE.md/LICENSE)

# snakeATAC.utility
This R-package *snakeATAC.utility* is a companion R utilities for the `snakeATAC` pipeline,
a Snakemake workflow for the analysis of ATAC-seq data. Provides functions to import,
re-analyse and re-plot pipeline outputs — including sample-level quality control
(principal component analysis and correlation heatmaps from genome-wide coverage),
transcription-factor footprints and differential-binding volcanoes from TOBIAS results,
filtering of HOCOMOCO motif collections, and summaries of pipeline runtime and resource usage —
so that results can be inspected and figures customized directly in R.


## Installation
### Developmental versions
```r
## Install remotes from CRAN (if not already installed)
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install the snakeATAC package
remotes::install_github("sebastian-gregoricchio/snakeATAC",
                        subdir = "utility/snakeATAC.utility",
                        build_manual = TRUE,
                        build_vignettes = TRUE)
```


<br />

## Documentation
With the package a [web-manual](https://sebastian-gregoricchio.github.io/snakeATAC/utility/snakeATAC.utility/manual/index.html) and a [vignette](https://sebastian-gregoricchio.github.io/snakeATAC/utility/snakeATAC.utility/doc/snakeATAC.utility_overview.html) are available.
The vignette can be inspected on R as well by typing `browseVignettes("snakeATAC.utility")`.


<br />

## Package history and releases
A list of all releases and respective description of changes applied could be found [here](https://sebastian-gregoricchio.github.io/snakeATAC/utility/snakeATAC.utility/NEWS).

<br />

-----------------
## Contact
For any suggestion, bug fixing, commentary please report it in the [issues](https://github.com/sebastian-gregoricchio/snakeATAC/issues)/[request](https://github.com/sebastian-gregoricchio/snakeATAC/pulls) tab of this repository.

## License
This package is under a GNU General Public License (version 3).

<br />

#### Contributors
![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/snakeATAC?size=50&padding=5&bots=true)
