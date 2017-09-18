# bvar

## Overview

bvar is a collection of R routines for estimating Linear and Nonlinear Bayesian Vector Autoregressive models in R. The R code is based on the Matlab Code by [Blake and Mumtaz (2012)](http://www.bankofengland.co.uk/education/Pages/ccbs/technical_handbooks/techbook4.aspx) and Koop and Koribilis (2009)

Models include:

* Linear VAR-models:
    * Minnesota and independent Normal-Wishart Prior
    * Impulse-Response Functions with Cholesky Decomposition or Sign Restrictions
* Nonlinear VAR-models:
    * Threshold VAR with independent Normal-Wishart Prior and Generalized Impulse-Response Functions and Cholesky Decomposition
* Other Models:
	* Factor-Augmented VARs with independent Normal-Wishart Prior 
	* Regime Switching Models
    
## Installation

To install the package you need the devtools package. If you don't have the devtools package, you can install it with

    install.packages("devtools")

Once you have installed the devtools package you can install the bvar package with

    library('devtools')
    devtools::install_github('joergrieger/bvar')
