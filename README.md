<!-- README.md is generated from README.Rmd. Please edit that file -->



# WeMix

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version-ago/WeMix)](https://www.r-pkg.org/badges/version-ago/WeMix)
[![CRAN weekly](https://cranlogs.r-pkg.org/badges/WeMix)](https://cranlogs.r-pkg.org/badges/WeMix)
[![CRAN grand total](https://cranlogs.r-pkg.org/badges/grand-total/WeMix)](https://cranlogs.r-pkg.org/badges/grand-total/WeMix)
<!-- badges: end -->

## Overview

Run mixed-effects models that include weights at every level. The WeMix package fits a weighted mixed model, also known as a multilevel, mixed, or hierarchical linear model (HLM). The weights could be inverse selection probabilities, such as those developed for an education survey where schools are sampled probabilistically, and then students inside of those schools are sampled probabilistically. Although mixed-effects models are already available in R, WeMix is unique in implementing methods for mixed models using weights at multiple levels. Both linear and logit models are supported. Models may have up to three levels.

## Pre-release Installation


``` r
# You can install the development version from GitHub with:
install.packages("devtools")
devtools::install_github("American-Institutes-for-Research/WeMix")
```


## CRAN-release Installation

``` r
# You can install the released version of EdSurvey from CRAN with:
install.packages("WeMix")
```

## Contributions

Contributions from external collaborators will be considered for inclusion in the package.
