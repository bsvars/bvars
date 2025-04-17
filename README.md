
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bvarNWish

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bvarNWish/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bvarNWish/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

An **R** package for Bayesian Vector Autoregressions with Normal-Wishart
Priors

## Installation

#### The first time you install the package

You must have a **cpp** compiler. Follow the instructions from [Section
1.3. by Eddelbuettel & François
(2023)](https://cran.r-project.org/package=Rcpp/vignettes/Rcpp-FAQ.pdf).
In short, for **Windows:** install
[RTools](https://CRAN.R-project.org/bin/windows/Rtools/), for **macOS:**
install [Xcode Command Line
Tools](https://www.freecodecamp.org/news/install-xcode-command-line-tools/),
and for **Linux:** install the standard development packages.

#### Once that’s done:

The developer’s version of the package with the newest features can be
installed by typing:

    devtools::install_github("bsvars/bvarNWish")
