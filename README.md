
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VARS

<!-- badges: start -->

<!-- badges: end -->

The goal of VARS is to implement various VAR estimation and
identification approaches, as well as inference tools.

## Progress

- [x] basic VAR model and VAR tools (IRF, FEVD, HDC)
- [x] classic VAR, i.e. OLS estimation, IV and recursive identification,
  bootstrap
- [x] bayesian VAR with (narrative) sign restrictions
- [x] add test data and code to the R package
- [x] add function and data documentation

## Installation

You can install the development version of VARS from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("GoldenBaozi/VARS")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(VARS)
## basic example code
## see /tests/testthat/test-bVAR.R
```

## Replic

Replication of Antolín-Díaz and Rubio-Ramírez (2018), figure 2

<figure>
<img src="./IRF_oil.png" alt="IRFs of global oil market" />
<figcaption aria-hidden="true">IRFs of global oil market</figcaption>
</figure>
