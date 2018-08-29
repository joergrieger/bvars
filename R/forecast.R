forecast.bvar <- function(bvarObj,forecastHorizon = 16){

  # Preliminary Calculations
  nVariables       <- dim(bvarObj$mydata)[2] # Number of variables
  nLength          <- dim(bvarObj$mydata)[1] # Length of time series
  nLags            <- bvarObj$NoLags  # Number of lags
  nForecasts       <- dim(bvarObj$betadraws)[3] # Number of forecasts, depends on the number sampled posteriors
  mForecastStorage <- array(NA,dim=c(forecastHorizon+nLags,nVariables,nForecasts)) # Matrix to storage forecasts

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
    print(ii)
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
      forecastUpper[jj,ii] <- quantile(mForecast[jj,ii,],probs=0.95)
      forecastLower[jj,ii] <- quantile(mForecast[jj,ii,],probs=0.05)

    }
  }

  forecastFinalMean  <- rbind(as.matrix(bvarObj$mydata),forecastMean)
  forecastFinalUpper <- rbind(as.matrix(bvarObj$mydata),forecastUpper)
  forecastFinalLower <- rbind(as.matrix(bvarObj$mydata),forecastLower)

  retList <- list(forecast = forecastFinalMean, Upper = forecastFinalUpper, Lower = forecastFinalLower)

}
