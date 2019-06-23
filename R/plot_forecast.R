pltfcbvar <- function(fcobj){

  cnames <- colnames(fcobj$forecast)
  nDim   <- ncol(fcobj$forecast)
  nLength <- floor(sqrt(nDim))

  if(is.ts(fcobj$forecast)){

    tsStart     = start(fcobj$forecast)
    tsFrequency = frequency(fcobj$forecast)

  }

  pltList <- list()

  for(ii in 1:nDim){

    tempDf <- data.frame(Upper = fcobj$Upper[,ii],forecast = fcobj$forecast[,ii],Lower = fcobj$Lower[,ii],Original=fcobj$Original[,ii])
    if(is.ts(fcobj$forecast)){
      tempts <- ts(tempDf,start=tsStart,frequency=tsFrequency)
    }
    else{
      tempts <- ts(tempDf,start=c(1))
    }

    p1 <- autoplot(tempts, facets=FALSE) +
      theme(legend.position = "none") +
      labs(title=cnames[ii])
    #p1
    #readline("Press [Enter] to continue")

    pltList[[ii]] <- p1

  }
  suppressWarnings(do.call("grid.arrange",c(pltList,ncol=nLength)))

}
