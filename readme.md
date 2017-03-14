# bvar

## Overview

bvar is a collection of R routines for estimating Linear and Nonlinear Bayesian Vector Autoregressive models in R.

Functionalities include:

* Linear VAR-models:
    * Minnesota and independent Normal-Wishart Prior
    * Impulse-Response Functions with Cholesky Decomposition
* Nonlinear VAR-models:
    * Threshold VAR with independent Normal-Wishart Prior and Generalized Impulse-Response Functions and Cholesky Decomposition
    
## Installation

To install the package you need the devtools package. If you don't have the devtools package, you can isntall it with

    install.packages("devtools")

Once you have installed the devtools package you can install the bvar package with

    library('devtools')
    devtools::install_github('joergrieger/bvar')
