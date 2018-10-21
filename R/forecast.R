forecast.bvar <- function(bvarObj,forecastHorizon = 16,interval = c(0.95,0.05)){

  # Preliminary Calculations
  nVariables       <- dim(bvarObj$mydata)[2] # Number of variables
  nLength          <- dim(bvarObj$mydata)[1] # Length of time series
  nLags            <- bvarObj$NoLags  # Number of lags
  nForecasts       <- dim(bvarObj$betadraws)[3] # Number of forecasts, depends on the number sampled posteriors
  mForecastStorage <- array(NA,dim=c(forecastHorizon+nLags,nVariables,nForecasts)) # Matrix to storage forecasts

  if(is.ts(bvarObj$mydata)){
    tsStart          <- start(bvarObj$mydata)
    tsFrequency      <- frequency(bvarObj$mydata)
  }

  tsColnames <- colnames(bvarObj$mydata)

  # Check if mydata is a time series object
  # If it is a time series object get frequency, etc.
  if(is.ts(bvarObj$mydata)){
    nFreq <- frequency(bvarObj$mydata)
    nFirstDate <- min(time(bvarObj$mydata))
    nFirstYear <- floor(nFirstDate)
    nFirstMonthQuarter <- (nFirstDate-nFirstYear) * nFreq
    nLastDate  <- max(time(bvarObj$mydata))
  }

  # Do the forecasts
  for(ii in 1:nForecasts){
    mForecastStorage[1:nLags,,ii] <- (bvarObj$mydata[nLength:(nLength-nLags+1),])

    for(jj in 1:forecastHorizon){

      nStart <- jj
      nEnd   <- jj+nLags-1
      y <- mForecastStorage[nStart:nEnd,,ii]

      randDraw <- rnorm(nVariables) %*% t(chol(bvarObj$sigmadraws[,,ii]))

      if(bvarObj$intercept == TRUE){

        y <- c(1,t(y))

        tempForecast <- y%*%bvarObj$betadraws[,,ii] + randDraw

      }
      else{

        y <- t(y)
        tempForecast <- y%*%bvarObj$betadraws[,,ii] + randDraw

      }
      # Storing the forecast
      mForecastStorage[jj+nLags,,ii] <- tempForecast
    }

  }

  # Remove initial values
  mForecast <- mForecastStorage[-c(1:nLags),,]
  forecastMean  <- array(0,dim=c(forecastHorizon,nVariables))
  forecastUpper <- array(0,dim=c(forecastHorizon,nVariables))
  forecastLower <- array(0,dim=c(forecastHorizon,nVariables))


  for(ii in 1:nVariables){
    for(jj in 1:forecastHorizon){

      forecastMean[jj,ii]  <- mean(mForecast[jj,ii,])
      forecastUpper[jj,ii] <- quantile(mForecast[jj,ii,],probs=min(interval))
      forecastLower[jj,ii] <- quantile(mForecast[jj,ii,],probs=max(interval))

    }
  }

  #forecastFinalMean  <- rbind(as.matrix(bvarObj$mydata),forecastMean)
  #forecastFinalUpper <- rbind(as.matrix(bvarObj$mydata),forecastUpper)
  #forecastFinalLower <- rbind(as.matrix(bvarObj$mydata),forecastLower)

  OriginalPath       <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalMean  <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalUpper <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalLower <- array(NA,dim=c(nLength + forecastHorizon, nVariables))

  OriginalPath[1:nLength,] <- bvarObj$mydata
  forecastFinalMean[(nLength + 1):(nLength + forecastHorizon),] <- forecastMean
  forecastFinalUpper[(nLength + 1):(nLength + forecastHorizon),] <- forecastUpper
  forecastFinalLower[(nLength + 1):(nLength + forecastHorizon),] <- forecastLower


  if(is.ts(bvarObj$mydata)){

    forecastFinalMean  <- ts(forecastFinalMean,start=tsStart,frequency=tsFrequency)
    forecastFinalUpper <- ts(forecastFinalUpper,start=tsStart,frequency=tsFrequency)
    forecastFinalLower <- ts(forecastFinalLower,start=tsStart,frequency=tsFrequency)
    OriginalPath       <- ts(OriginalPath, start = tsStart, frequency = tsFrequency)

  }

  colnames(forecastFinalMean)  <- colnames(bvarObj$mydata)
  colnames(forecastFinalUpper) <- colnames(bvarObj$mydata)
  colnames(forecastFinalLower) <- colnames(bvarObj$mydata)
  colnames(OriginalPath)       <- colnames(bvarObj$mydata)

  retList <- structure(list(forecast = forecastFinalMean, Upper = forecastFinalUpper, Lower = forecastFinalLower,
                            Original = OriginalPath),class="fcbvar")
  return(retList)

}

forecast.tvar <- function(tvarObj, forecastHorizon = 4, interval =c(0.05,0.95)){

  # Preliminary Calculations
  nVariables       <- dim(tvarObj$mydata)[2] # Number of variables
  nLength          <- dim(tvarObj$mydata)[1] # Length of time series
  nLags            <- tvarObj$NoLags  # Number of lags
  nForecasts       <- dim(tvarObj$Alphadraws)[4] # Number of forecasts, depends on the number sampled posteriors
  mForecastStorage <- array(NA,dim=c(forecastHorizon+nLags,nVariables,nForecasts)) # Matrix to storage forecasts

  if(is.ts(tvarObj$mydata)){

    tsStart          <- start(tvarObj$mydata)
    tsFrequency      <- frequency(tvarObj$mydata)

  }

  nl1 <- max(tvarObj$thMax,tvarObj$NoLags)
  for(ii in 1:nForecasts){

    mForecastStorage[1:nl1,,ii] <- (tvarObj$mydata[nLength:(nLength-nl1+1),])

    for(jj in 1:forecastHorizon){

      nStart <- jj
      nEnd   <- jj+nLags-1
      y <- mForecastStorage[nStart:nEnd,,ii]
      thCheck <- mForecastStorage[jj + nl1 - tvarObj$deldraws[ii] - 1,tvarObj$thVar, ii]

      if(thCheck < tvarObj$tardraws[ii]){
        # First Regime
        Alpha <- tvarObj$Alphadraws[,,1,ii]
        Sigma <- tvarObj$Sigmadraws[,,1,ii]
      }
      else{
        # Second Regime
        Alpha <- tvarObj$Alphadraws[,,2,ii]
        Sigma <- tvarObj$Sigmadraws[,,2,ii]
      }

      randDraw <- rnorm(nVariables) %*% t(chol(Sigma))

      if(tv$Intercept){

        y <- c(1,t(y))
        tempForecast <- y%*%Alpha + randDraw

      }
      else{
        y <- t(y)
        tempForecast <- y%*%Alpha + randDraw
      }

      # Storing the forecast
      mForecastStorage[jj+nl1,,ii] <- tempForecast

    }

  }
  # Remove initial values
  mForecast <- mForecastStorage[-c(1:nl1),,]
  forecastMean  <- array(0,dim=c(forecastHorizon,nVariables))
  forecastUpper <- array(0,dim=c(forecastHorizon,nVariables))
  forecastLower <- array(0,dim=c(forecastHorizon,nVariables))
  for(ii in 1:nVariables){
    for(jj in 1:forecastHorizon){

      forecastMean[jj,ii]  <- mean(mForecast[jj,ii,])
      forecastUpper[jj,ii] <- quantile(mForecast[jj,ii,],probs=max(interval))
      forecastLower[jj,ii] <- quantile(mForecast[jj,ii,],probs=min(interval))

    }
  }

  #forecastFinalMean  <- rbind(as.matrix(tvarObj$mydata),forecastMean)
  #forecastFinalUpper <- rbind(as.matrix(tvarObj$mydata),forecastUpper)
  #forecastFinalLower <- rbind(as.matrix(tvarObj$mydata),forecastLower)
  OriginalPath       <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalMean  <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalUpper <- array(NA,dim=c(nLength + forecastHorizon, nVariables))
  forecastFinalLower <- array(NA,dim=c(nLength + forecastHorizon, nVariables))

  OriginalPath[1:nLength,] <-tvarObj$mydata
  forecastFinalMean[(nLength + 1):(nLength + forecastHorizon),] <- forecastMean
  forecastFinalUpper[(nLength + 1):(nLength + forecastHorizon),] <- forecastUpper
  forecastFinalLower[(nLength + 1):(nLength + forecastHorizon),] <- forecastLower

  if(is.ts(tvarObj$mydata)){

    forecastFinalMean  <- ts(forecastFinalMean,start=tsStart,frequency=tsFrequency)
    forecastFinalUpper <- ts(forecastFinalUpper,start=tsStart,frequency=tsFrequency)
    forecastFinalLower <- ts(forecastFinalLower,start=tsStart,frequency=tsFrequency)
    OriginalPath       <- ts(OriginalPath,start=tsStart,frequency=tsFrequency)

  }

  colnames(forecastFinalMean)  <- colnames(tvarObj$mydata)
  colnames(forecastFinalUpper) <- colnames(tvarObj$mydata)
  colnames(forecastFinalLower) <- colnames(tvarObj$mydata)
  colnames(OriginalPath)       <- colnames(tvarObj$mydata)

  retList <- structure(list(forecast = forecastFinalMean, Upper = forecastFinalUpper, Lower = forecastFinalLower,
                            Original = OriginalPath),class="fctvar")
  return(retList)

}
