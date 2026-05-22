
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bvars

An **R** package for Bayesian Forecasting with Vector Autoregressions

<!-- badges: start -->

[![R-CMD-check](https://github.com/bsvars/bvars/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bsvars/bvars/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Provides fast and efficient procedures for Bayesian estimation and
forecasting using state-of-the-art Vector Autoregressions. This package
includes the model proposed by [Chan
(2020)](https://doi.org/10.1080/07350015.2018.1451336), that is, a
Bayesian Vector Autoregression with Minnesota priors and a flexible
structure of the error term specification. The latter includes:
conditional multivariate normal or Student’s t distributions, as well as
homoskedastic or heteroskedastic specifications with a common volatility
modelled by centred or non-centred Stochastic Volatility. Additionally,
the package facilitates predictive analyses using density forecasting
and forecast-error variance decompositions. All this is complemented by
simple workflows, useful plots and summary functions, and comprehensive
documentation. The ‘bvars’ package aligns with R packages ‘bsvars’ by
[Woźniak (2024)](https://doi.org/10.32614/CRAN.package.bsvars),
‘bsvarSIGNs’ by [Wang & Woźniak
(2025)](https://doi.org/10.32614/CRAN.package.bsvarSIGNs), and ‘bpvars’
by [Woźniak (2025)](https://doi.org/10.32614/CRAN.package.bpvars)
regarding objects, workflows, and code structure, and they constitute an
integrated toolset.

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

    devtools::install_github("bsvars/bvars")
