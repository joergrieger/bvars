#' @export
#' Function to draw one single path for the impulse-response functions.
#' @param Alpha the estimate of the VAR-coefficients
#' @param Sigma the estimate of the Variance-covariance matrix
#' @param id_obj S3-object for identifying structural shocks
#' @param nolags Number of lags in the model
#' @param intercept Whether the model has an intercept or not
#' @param nhor the horizon of the impulse-response functions

compirf <- function(Alpha, Sigma,id_obj, nolags, intercept = TRUE, nhor){

  # Preliminaries

  K <- nrow(Sigma)
  bigj <- matrix(0, K, K * nolags)
  bigj[1:K,1:K] <- diag(K)

  if(intercept == TRUE){

    B <- Alpha[-c(1),]

  }
  else{

    B <- Alpha

  }

  PHI_Mat <- companionmatrix(B,nolags)
  biga <- PHI_Mat
  bigai <- biga

  # identify the model

  shock <- structural(id_obj, Alpha, Sigma)

  impresp <- matrix(0,K,K * nhor)
  impresp[1:K,1:K] <- shock

  for(ii in 1:(nhor-1)){

    impresp[,(ii * K + 1):(( ii + 1) * K)] <- (bigj %*% bigai %*% t(bigj) %*% shock)
    bigai <- bigai %*% biga

  }

  imp <- array(0,dim=c(K,K,nhor))
  jj <- 0

  for(ii in 1:K){
    for(ij in 1:nhor){

      jj <- ii + K * (ij - 1)
      imp[ii,,ij] <- impresp[,jj]

    }
  }

  return(imp)

}
