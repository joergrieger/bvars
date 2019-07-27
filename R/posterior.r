#' @title function to draw the posterior
#' @param priorObj object of class cnw
#' @param yLagged lhs of the VAR model
#' @param xLagged rhs of the VAR model
#' @param previous the previous draw
#' @param stabletest tests for stability of the draw of VAR-coefficients
#' @param ... currently not used
draw_posterior.cnw <- function(priorObj, yLagged, xLagged, previous, stabletest = TRUE,...){

  K <- ncol(yLagged)

  # preliminary calculations

  betaols <- solve( t( xLagged ) %*% xLagged ) %*% t( xLagged ) %*% yLagged
  residuals <- yLagged - xLagged %*% betaols
  SSE       <- t( residuals ) %*% residuals

  # Calculate posterior for coefficients
  Vpost <- solve( solve( priorObj$coefpriorvar ) + ( t( xLagged ) %*% xLagged ) )
  Apost <- Vpost %*% ( solve( priorObj$coefpriorvar ) %*% priorObj$coefprior + t( xLagged ) %*% xLagged %*% betaols )

  # Calculate posterior for variance - covariance matrix


  Spost     <- SSE + priorObj$varprior + t( betaols ) %*% t(xLagged) %*% xLagged %*% betaols +
    t(priorObj$coefprior) %*% solve( priorObj$coefpriorvar ) %*% priorObj$coefprior -
    t( Apost ) %*% ( solve( priorObj$coefpriorvar ) + t( xLagged) %*% xLagged ) %*% Apost

  vpost <- nrow(yLagged) + priorObj$varpriordf
  cova  <- previous$Sigma %x% Vpost

  # Draw posterior for coefficients and variance-covariance matrix


  # draw posterior for coefficients and test if it is stable

  stable <- 2
  while(stable > 1){

    alpha <- MASS::mvrnorm(mu = as.vector( Apost ), Sigma = cova)
    Alpha <- matrix(alpha,ncol = K)

    if( stabletest ){

      if(priorObj$intercept){

        Alphatest <- Alpha[-c(1),]

      }
      else{

        Alphatest <- Alpha

      }

      stable <- stability( betadraw = Alphatest, nolags = priorObj$nolags )

    }
    else{

      stable <- 0

    }

  }

  # Draw Variance-Covariance Matrix

  Sigma <- solve( stats::rWishart(1, vpost, solve( Spost ))[,,1] )

  # Return values
  return(list(Alpha = Alpha, Sigma = Sigma, addinfo = NULL ))


}

#' @title function to draw the posterior for an uninformative prior
#' @param priorObj object of class cnw
#' @param yLagged lhs of the VAR model
#' @param xLagged rhs of the VAR model
#' @param previous the previous draw
#' @param stabletest tests for stability of the draw of VAR-coefficients
#' @param ... currently not used

draw_posterior.unf <- function(priorObj,yLagged,xLagged,previous,stabletest = TRUE,...){

  # uninformative prior
  K   <- ncol(yLagged)
  obs <- nrow(yLagged)

  #variables
  Alpha     <- previous$Alpha
  Sigma     <- previous$Sigma
  intercept <- priorObj$intercept
  nolags    <- priorObj$nolags

  # Calculate posterior
  Vpost   <- Sigma %x% solve(t(xLagged) %*% xLagged)
  betaols <- solve(t(xLagged) %*% xLagged) %*% t(xLagged) %*% yLagged
  sse     <- t(yLagged - xLagged %*% betaols)%*%(yLagged - xLagged %*% betaols)

  # Draw posterior
  stable <-2
  while(stable>1){

    alpha <- MASS::mvrnorm(mu=as.vector(betaols),Sigma=Vpost)
    Alpha <- matrix(alpha,ncol=K)

    # Check if coefficients are good
    if(stabletest==TRUE){

      if(intercept==TRUE){

        Alphatest <- Alpha[-c(1),]

      }
      else{

        Alphatest <- Alpha

      }

      stable <- stability(betadraw = Alphatest, nolags = nolags)

    }
    else{

      stable <- 0

    }
  }

  Sigma <- solve(stats::rWishart(1,obs,solve(sse))[,,1])

  return(list(Alpha = Alpha,
              Sigma = Sigma,
              addInfo = NULL))

}


#' @title function to draw the posterior for a Minnesota Prior
#' @param priorObj object of class cnw
#' @param yLagged lhs of the VAR model
#' @param xLagged rhs of the VAR model
#' @param previous the previous draw
#' @param stabletest tests for stability of the draw of VAR-coefficients
#' @param ... currently not used

draw_posterior.minnesota <- function(priorObj,yLagged,xLagged,previous,stabletest = TRUE,...){

  # Declare variables
  K   <- ncol(yLagged)
  obs <- nrow(yLagged)

  # Previous draws

  Sigma <- previous$Sigma
  Alpha <- previous$Alpha

  # prior variables
  Vprior  <- priorObj$Vprior
  aprior  <- as.vector(priorObj$Aprior)
  betaols <- solve(t(xLagged) %*% xLagged) %*% t(xLagged) %*% yLagged
  betaols <- as.vector(betaols)


  Vpost <- solve(solve(Vprior) + solve(Sigma) %x% (t(xLagged) %*% xLagged))
  apost <- Vpost %*% (solve(Vprior) %*% aprior + (solve(Sigma) %x% (t(xLagged) %*% xLagged)) %*% betaols)


  # Draw coefficients

  stable <- 2
  while(stable>1){

    alpha <- MASS::mvrnorm(mu=apost,Sigma=Vpost)
    Alpha <- matrix(alpha,ncol=K)

    if(stabletest == TRUE){

      if(priorObj$intercept == TRUE){

        Alphatest <- Alpha[-c(1),]

      }
      else{

        Alphatest <- Alpha

      }

      stable <- stability(betadraw = Alphatest,nolags = priorObj$nolags)

    }

    else{

      stable <- 0

    }

  }

  # Draw Sigma

  residuals  <- yLagged - xLagged %*% Alpha

  wishart_scale        <- t(residuals) %*% residuals / obs
  Sigma <- solve(stats::rWishart(1,obs,wishart_scale)[,,1])

  return(list(Alpha = Alpha,
              Sigma = Sigma,
              addInfo = NULL))

}

#' @title draw posterior for VAR Model with SSVS-Prior
#' @param priorObj object of class cnw
#' @param yLagged lhs of the VAR model
#' @param xLagged rhs of the VAR model
#' @param previous the previous draw
#' @param stabletest tests for stability of the draw of VAR-coefficients
#' @param ... currently not used

draw_posterior.ssvs <- function(priorObj,yLagged,xLagged,previous,stabletest = TRUE,...){

  # Preliminaries
  kappa0 <- priorObj$kappa0
  kappa1 <- priorObj$kappa1
  tau0   <- priorObj$tau0
  tau1   <- priorObj$tau1

  intercept <- priorObj$intercept
  nolags    <- priorObj$nolags

  aprior   <- priorObj$aprior

  bi=0.01
  ai=0.01
  qij=0.5
  p_i=0.5

  K   <- ncol(yLagged)
  obs <- nrow(yLagged)
  norest <- ncol(yLagged) * ncol(xLagged)

  # OLS-draw
  betaols <- solve( t(xLagged) %*% xLagged ) %*% t(xLagged) %*% yLagged
  aols    <- as.vector(betaols)

  # get info from previous ddraw

  Alpha    <- previous$Alpha
  Sigma    <- previous$Sigma

  SSEGibbs <- previous$addInfo$SSEGibbs
  omega    <- previous$addInfo$omega
  gammas   <- previous$addInfo$gammas

  #
  # First: Draw Sigma
  #

  S <- array(list(),dim=c(K))

  for(kk2 in 1:K){

    S[[kk2]] <- SSEGibbs[1:kk2,1:kk2]

  }

  s <- array(list(),dim=(K-1))
  for(kk3 in 2:K){

    s[[kk3-1]] <- SSEGibbs[1:(kk3-1),kk3]

  }

  hh <- array(list(),dim=c(K-1))

  for(kk4 in 1:(K-1)){

    omeg <- omega[[kk4]]
    het <- hh[[kk4]]

    for(kkk in 1:dim(omeg)[1]){

      if(omeg[kkk]==0){

        het[kkk] <- kappa0

      }
      else{

        het[kkk] <- kappa1

      }
    }

    hh[[kk4]] <- het
  }

  Dj <- array(list(),dim=c(K-1))

  for(kk5 in 1:(K-1)){

    if(kk5==1){

      Dj[[kk5]] <- diag(hh[[kk5]],1)

    }
    else{

      Dj[[kk5]] <- diag(hh[[kk5]])

    }
  }

  DDj <- array(list(),dim=c(K-1))

  for(kk6 in 1:(K-1)){

    DD <- Dj[[kk6]]
    DDj[[kk6]] <- DD %*% DD

  }

  # Create B[i] matrix
  B <- array(list(),dim=c(K))
  for(rr in 1:K){
    if(rr==1){

      B[[rr]] <- bi+0.5*SSEGibbs[rr,rr]

    }
    else if(rr>1){

      si <- s[[rr-1]]
      Si <- S[[rr-1]]
      DiDi <- DDj[[rr-1]]
      B[[rr]] <- bi+0.5*(SSEGibbs[rr,rr]-t(si)%*%solve(Si+solve(DiDi))%*%si)

    }
  }

  psi_ii_sq = array(0,dim=c(K))

  for(kk7 in 1:K){

    psi_ii_sq[kk7] <- stats::rgamma(1,(ai+0.5*T),(B[[kk7]]))

  }

  #
  # Draw eta
  #

  eta <- array(list(),dim=c(K-1))
  for(kk8 in 1:(K-1)){

    si <- s[[kk8]]
    Si <- S[[kk8]]
    DiDi <- DDj[[kk8]]
    miuj <- -sqrt(psi_ii_sq[kk8+1])*(solve(Si+solve(DiDi))%*%si)
    Deltaj <- solve(Si+solve(DiDi))
    eta[[kk8]] <- miuj+t(chol(Deltaj))%*%stats::rnorm(kk8)

  }

  #
  # Draw Omega
  #

  for(kk9 in 1:(K-1)){

    omegg <- omega[[kk9]]
    etag  <- eta[[kk9]]
    omegavec <- matrix(nrow=0,ncol=1)

    for(nn in 1:dim(omegg)[1]){

      uij1 <- 1/(kappa0)*exp(-0.5*((etag[nn])^2)/(kappa0^2))*qij
      uij2 <- 1/(kappa1)*exp(-0.5*((etag[nn])^2)/(kappa1^2))*(1-qij)
      ost  <- max(0,uij1/(uij1+uij2))
      ost  <- min(1,ost)
      omegg <- stats::rbinom(1,1,ost)
      omegavec <- rbind(omegavec,omegg)

    }

    omega[[kk9]] <- omegavec

  }

  # Create PSI Matrix

  PSI_ALL <- array(0,dim=c(K,K))

  for(nn1 in 1:K){

    PSI_ALL[nn1,nn1] <- sqrt(psi_ii_sq[nn1])

  }
  for(nn2 in 1:(K-1)){

    etagg <- eta[[nn2]]

    for(nnn in 1:dim(etagg)[1]){

      PSI_ALL[nnn,nn2+1] <- etagg[nnn]

    }
  }

  Sigma <- solve(PSI_ALL %*% t(PSI_ALL))

  #
  # End Drawing Sigma
  # Draw Alphas now
  #

  hi <- array(0,dim=c(norest))
  for(nn3 in 1:norest){

    if(gammas[nn3]==0){

      hi[nn3] <- tau0#tau0[nn3]

    }
    else if(gammas[nn3]==1){

      hi[nn3] <- tau1#tau1[nn3]

    }
  }

  D <- diag(as.vector(t(hi) %*% diag(1,norest)))
  #D <- diag(hi)

  DD <- D %*% D

  isig <- solve(Sigma)
  psixx <- isig%x%(t(xLagged)%*%xLagged)
  Vpost <- solve(psixx+solve(DD))
  apost <- Vpost%*%(psixx%*%aols+solve(DD)%*%aprior)
  stable <- 2

  while(stable>1){

    alpha <- MASS::mvrnorm(mu=apost,Sigma=Vpost)
    Alpha <- matrix(alpha,ncol=K)

    # Check if coefficients are good
    if(stabletest==TRUE){

      if(priorObj$intercept == TRUE){

        Alphatest <- Alpha[-c(1),]

      }
      else{

        Alphatest <- Alpha

      }
      stable <- stability(betadraw = Alphatest,nolags = nolags)

    }
    else{
      stable <- 0
    }
  }
  for(nn6 in 1:norest){

    ui1 <- 1/tau0*exp(-0.5*(alpha[nn6]/tau0)^2)*p_i
    ui2 <- 1/tau1*exp(-0.5*(alpha[nn6]/tau1)^2)*(1-p_i)

    gst <- min(1,ui1/(ui1+ui2))
    gst <- max(0,gst)
    if(ui1==Inf){gst=1}
    if(ui1+ui2<1e-7){gst=0.5}
    gammas[nn6] <- stats::rbinom(1,1,gst)

  }

  SSEGibbs <- t(yLagged - xLagged %*% Alpha) %*% (yLagged - xLagged %*% Alpha)
  addInfo <- list(gammas = gammas,omega = omega, SSEGibbs = SSEGibbs)

  retlist <- list(Sigma=Sigma,Alpha=Alpha,addInfo = addInfo)
  return(retlist)


}

