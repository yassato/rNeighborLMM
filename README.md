# rNeighborLMM
rNeighborLMM: Fitting mixed models for oligogenic and polygenic indirect genetic effects


## Description
This R package serves as a wrapper of RAINBOWR (Hamazaki & Iwata 2020) and rNeighborGWAS (Sato et al. 2021), which solves a multi-kernel mixed model incorporating direct and indirect genetic effects (DGEs and IGEs) with their covariance. Genome-wide association study (GWAS) of oligogenic indirect genetic effects can be performed following the same way as rNeighborGWAS operates. See Sato & Hamazaki (2026) for theoretical details.

## Installation
Please install the rNeighborLMM package through GitHub. Note that this package is not available at CRAN.
```
devtools::install_github("yassato/rNeighborLMM", repo="master", build_vignettes = TRUE)
```

## Usage
See the package vignette after installation.
```
vignette("rNeighborLMM")
```
