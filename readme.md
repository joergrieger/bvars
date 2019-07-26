# bvar

## Overview

bvar is a collection of R routines for estimating Linear and Nonlinear Bayesian Vector Autoregressive models in R. The R code is based on the Matlab Code by [Blake and Mumtaz (2012)](http://www.bankofengland.co.uk/education/Pages/ccbs/technical_handbooks/techbook4.aspx) and Koop and Koribilis (2009)

Models and functionalities include:

* Linear VAR-models:
    * Minnesota and independent Normal-Wishart Prior
* Nonlinear VAR-models:
    * Threshold VAR with generalized impulse response-functions
* Other Models:
	* Factor-Augmented VARs
	* Regime Switching Models
* Prior for VARs
    * Independent Normal-Wishart
	* Minnesota Prior
	* Natural Conjugate prior
	* Uninformative prior 
* Identification of structural shocks
    * Cholesky decomposition
	* Sign restrictions
    
## Installation

To install the package you need the devtools package. If you don't have the devtools package, you can install it with

    install.packages("devtools")

Once you have installed the devtools package you can install the bvar package with

    library('devtools')
    devtools::install_github('joergrieger/bvars')
	
## To-do list

* improve numerical stability of Threshold-models
* speed up generalized impulse-response functions
* Documentation
* add functions for plotting impulse-respone functions, summary of inference, diagnostics, forecasting
* regime switching models with time-varying transition probabilities
* dummy observation prior
* ssvs - prior 


## Known issues and bugs

* generalized impulse functions work only for lags>1


