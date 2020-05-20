#' @export
#' @title Geweke Convergence diagnostics
#' @param obj an estimated bvar model
#' @param frac1 first fraction
#' @param frac2 second fraction
#' @param ... not used
#' @return returns a matrix with test results for each coefficient.
#' @details
#' Estimate convergence diagnostics of the mcmc sampler as in Geweke (1992). The diagnostics split the sample into two subsample, frac1 and frac2, and compares the mean and variance of these two subsample
#' \deqn{z=\frac{\hat{\theta}_1-\hat{\theta}_2}{\sqrt{Var(\hat\theta_1)+Var(\hat\theta_2)}}}
#' with \eqn{\theta_1} and \eqn{\theta_2} being the means of the coefficients of the two  subsamples. The variance of the coefficients of the subsample are estimated using a spectral approach to adjust for the dependence between the two subsamples.
#' This function uses the univariate version of the test in the sense that it performs a univariate test for each coefficient and calculates the probability that a coefficient has converged.
#' @reference Geweke, John, Evaluating the Accuracy of Sampling-Based approaches to the Calculation of Posterior Moments, In: J.M. Bernardo, J.O. Berger, A. P. Dawid and A.F.M. Smith, Eds., Bayesian Statistics, Vol. 4, Clarendon Press, Oxford, 1992, pp. 169-183
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
