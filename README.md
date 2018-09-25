# BEclear

This package provides functions to detect and correct for batch effects in
DNA methylation data. The core function "BEclear" is based on latent factor
models and can also be used to predict missing values in any other matrix
containing real numbers.

## Installation

```r
# Install the development version from GitHub:

# install.packages("devtools")
library(devtools)

devtools::install_github("David-J-R/BEclear", build_vignettes=TRUE)
```