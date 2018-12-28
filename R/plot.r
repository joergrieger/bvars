#
#
# Plot Impulse - Response functions or posteriors
#
#

plot.bvirf <- function(irfObj){

  pltBvarIrf(irfObj)

}

plot.tvirf <- function(irfObj){

  pltTvarIrf(irfObj)

}

plot.bvar <- function(bvarObj,lag=0){

    if(lag > bvarObj$NoLags){

      stop("lag greater than laglength of model.")

    }
  if(lag == 0){
    for(ii in 1:bvarObj$NoLags){

      pltBvarPosterior(bvarObj,lag=ii)
      print("Plotting posterior density for lag")
      readline("Press [Enter] to continue")
    }

  }
  else{

    pltBvarPosterior(bvarObj,lag=lag)

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


