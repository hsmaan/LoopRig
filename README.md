
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LoopRig <img src="man/figures/looprig_logo.png" height="180px" align="right"/>

[![Build
Status](https://travis-ci.com/hsmaan/LoopRig.svg?token=jBqxwnZzU1qwLZyzpxME&branch=master)](https://travis-ci.com/hsmaan/LoopRig)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/hsmaan/LoopRig?branch=master&svg=true)](https://ci.appveyor.com/project/hsmaan/LoopRig)
[![codecov](https://codecov.io/gh/hsmaan/LoopRig/branch/master/graph/badge.svg)](https://codecov.io/gh/hsmaan/LoopRig)

## Overview

LoopRig is an R package that aims to standardize complex
coordinate-based workflows utilizing chromatin loop and genomic element
data.

## Installation

LoopRig can be installed directly from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("hsmaan/LoopRig", build_vignettes = TRUE)
```

## Usage

Element data from BED4..n files and chromatin loop data from BEDPE files
is used as input for the `LoopsToRanges()` and `ElementsToRanges()`
functions, which create S3 containers for S4 *GRangesList* and *GRanges*
objects respectively. These containers are of class *LoopRanges* and
*ElementRanges*, and can be analyzed using the chromatin loop
manipulation and element linkage functions available in LoopRig.

An in-depth tutorial is available in the package vignettes:

``` r
browseVignettes("LoopRig")
vignette("LoopRig-Tutorial")
```

Complete package documentation is available at \[\]

## License

[GNU General Public
License 3.0](https://github.com/hsmaan/LoopRig/blob/master/LICENSE)
