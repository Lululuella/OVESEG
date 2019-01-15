# OVESEG

OVESEG-test (One Versus Everyone Subtype Exclusively-expressed Genes test) is a statistically-principled multiple-group comparison method that can detect tissue/cell-specific marker genes (MGs) among many subtypes. To assess the statistical significance of MGs, OVESEG-test uses a specifically designed test statistics that mathematically matches the definition of MGs, and employs a specifically designed novel permutation scheme to estimate the corresponding distribution under null hypothesis where the expression patterns of non-MGs can be highly complex.

This R package provides functions to compute OVESEG-test statistics and estimate p-values from weightedly aggregated permutations.

## Installation

You can install OVESEG, depending on R (>=3.5), from github by

``` r
devtools::install_github("Lululuella/OVESEG@pre-Bioconductor")
```

## Example

This is a basic example which shows you how to estimate p-values of MGs by OVESEG-test:

``` r
## specify the expression matrix, labels, number of permutations, and parallel option 
rtest <- OVESEGtest(y, group, NumPerm=999, BPPARAM=BiocParallel::SnowParam())
```

