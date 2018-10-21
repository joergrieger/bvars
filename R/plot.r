#
#
# Plot Impulse - Response functions or posteriors
#
#

plot.bvar <- function(bvarObj,type="irf",lag=1){
  if(type == "irf"){
    pltBvarIrf(bvarObj)
  }
  else if(type == "posterior"){
    print("posterior")

    if(lag > bvarObj$NoLags){
      stop("lag greater than laglength of model.")

    }
    pltBvarPosterior(bvarObj,lag=lag)
  }
  else{
    stop("type of plot must be irf or posterior")
  }
}

plot.tvar <- function(tvarObj, type = "irf",lag = 1){
  if(type == "irf"){
    pltTvarIrf(tvarObj)
  }
  else if(type == "posterior"){
    pltTvarPosterior(tvarObj)
  }
}

plot.msvar <- function(msvarObj, type = "irf"){
  if(type == "irf"){
    pltmsvarirf(msvarObj)
  }
}


#
# plot forecasts
#

plot.fcbvar <- function(fcobj){

  print("plotting forecasts from Bayesian VAR")
  pltfcbvar(fcobj)

}

plot.fctvar <- function(fcobj){

  print("Plotting forecast from Threshold VAR")
  pltfcbvar(fcobj)

}


