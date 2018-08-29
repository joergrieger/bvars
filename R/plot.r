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
