[![Build Status](https://travis-ci.com/joergrieger/bvars.svg?branch=master)](https://travis-ci.com/joergrieger/bvars) [![codecov](https://codecov.io/gh/joergrieger/bvars/branch/master/graph/badge.svg)](https://codecov.io/gh/joergrieger/bvars) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/bvars)]() [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![stability-stable](https://img.shields.io/badge/stability-stable-green.svg)](https://github.com/joethorley/stability-badges#stable)

# bvars

## Overview

bvar is a collection of R routines for estimating Linear and Nonlinear Bayesian Vector Autoregressive models in R. The R code is based on the Matlab Code by [Blake and Mumtaz (2012)](http://www.bankofengland.co.uk/education/Pages/ccbs/technical_handbooks/techbook4.aspx) and Koop and Koribilis (2009)

Models and functionalities include:

* VAR Models
  * Linear VARs
  * Regime Switching VARs
  * Threshold VARs
  * Factor-Augmented Models
* Identification of Structural Models
  * Cholesky decomposition
  * Sign Restrictions
* Functionalities to further analyze VARs
  * Impulse-Response Functions
  * Forecast error variance decomposition <not yet implemented>
  * Forecasting <not yet implemented>
  * historical decomposition <not yet implemented>
  
* [Project Homepage](https://joergrieger.github.io/bvars/)
    
## Installation

To install the package you need the devtools package. If you don't have the devtools package, you can install it with

    install.packages("devtools")

Once you have installed the devtools package you can install the bvar package with

    library('devtools')
    devtools::install_github('joergrieger/bvars')


## Known bugs and issues

* generalized impulse-response functions for threshold VARs work only for lags > 1
