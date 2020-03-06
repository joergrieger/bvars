#' @export
#' @title Geweke Convergence diagnostics
#' @param obj an estimated bvar model
#' @param frac1 first fraction
#' @param frac2 second fraction
#' @param ... not used
#' @rdname diag_geweke

diag_geweke <- function(obj,frac1,frac2,...) UseMethod("diag_geweke")


#' @export
#' @rdname diag_geweke
diag_geweke.bvar <- function(obj,frac1,frac2,...){

  # input check
  # run diagnostics
  nvar <- dim(obj$mcmc_draws$Alpha)[2]
  npar <- dim(obj$mcmc_draws$Alpha)[1]
  test_stats <- array(NA,dim=c(npar,nvar))
  for(ii in 1:nvar){
    for(jj in 1:npar){

      test_stats[jj,ii] <- univariate_geweke(obj$mcmc_draws$Alpha[jj,ii,],frac1,frac2)

    }
  }

  return(test_stats)
}

univariate_geweke <- function(data,frac1,frac2){

  nlength <- length(data)

  nrun1 <- floor(nlength * frac1)
  nrun2 <- ceiling(nlength * (1 - frac2))
  runlength2 <- nlength - nrun2 + 1

  mean1 <- mean(data[1:nrun1])
  mean2 <- mean(data[nrun2:nlength])

  # To Do: replace with spectral variance to account for autocorrelation
  bandwidth <- min(500,nlength/2)
  sd1 <- sqrt(spectral_variance(data[1:nrun1],bandwidth))
  sd2 <- sqrt(spectral_variance(data[nrun2:nlength],bandwidth))

  geweke.stat <- (mean1-mean2)/(sqrt(sd1/nrun1+sd2/runlength2))
  return(2 * (1 - pnorm(geweke.stat)))
}

spectral_variance <- function(data,bandwidth){

  nlength = length(data)
  bi = (1:bandwidth)/bandwidth
  vk1 <- 1 - 6 * bi^2 * (1 - bi)
  vk2 <- 2 * (1 - bi)^3

  vx <- data-mean(data) # demean data

  svar <- 0
  for(ii in 1:bandwidth){

    savar <- (t(vx[1:(nlength - ii)]) %*% vx[1:(nlength - ii)]) / nlength
    if(bi[ii] < 0.5){

      svar <- svar + vk1[ii] * savar

    }
    else{

      svar <- svar + vk2[ii] * savar

    }

  }

  res <- sum(vx^2)/nlength + 2 * (nlength)/(nlength - 1) * svar


}
