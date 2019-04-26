[![Travis-CI Build Status](https://travis-ci.org/David-J-R/BEclear.svg?branch=master)](https://travis-ci.org/David-J-R/BEclear)
[![Coverage Status](https://img.shields.io/codecov/c/github/David-J-R/BEclear/master.svg)](https://codecov.io/github/David-J-R/BEclear?branch=master)


# BEclear

## Description

This package provides functions to detect and correct for batch effects in
DNA methylation data. The core function for the data imputation is based on 
latent factor models and can also be used to predict missing values in any other 
matrix containing real numbers.

## Installation

```r
# Installation from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BEclear")
```

```r
# Installation of the development version from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("David-J-R/BEclear", build_opts = c())
```

## Usage example

For example code you can either read the [vignette](https://bioconductor.org/packages/devel/bioc/vignettes/BEclear/inst/doc/BEclear.html) on Bioconductor. Or you can open it local, after you installed the BEclear package with following command:

```r
browseVignettes("BEclear")
```

## Planned features

- [Treating Data-Sets Without Batches](https://github.com/David-J-R/BEclear/issues/22)
- [Dixon Test for Outlier Detection](https://github.com/David-J-R/BEclear/issues/20)
- [Make Saving Temporary Data Optional](https://github.com/David-J-R/BEclear/issues/21)
- [Bias Modelling](https://github.com/David-J-R/BEclear/issues/18)
- [ALS](https://github.com/David-J-R/BEclear/issues/19)
- [Test For Convergence Durin Epochs](https://github.com/David-J-R/BEclear/issues/23)
- [Testing BEclear on other data-sets](https://github.com/David-J-R/BEclear/issues/24)
- [After Merging Blocks Continue GD](https://github.com/David-J-R/BEclear/issues/25)

If you have a suggestion
you could either open an [issue](https://github.com/David-J-R/BEclear/issues) or 
write an email under the following address: [David.J.Rasp at gmail.com](mailto:David.J.Rasp@gmail.com)


## Citation

Akulenko, R., Merl, M., & Helms, V. (2016). BEclear: Batch effect detection and 
adjustment in DNA methylation data. PLoS ONE, 11(8), 1–17.
https://doi.org/10.1371/journal.pone.0159921
