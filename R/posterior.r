
#
# Function to draw from a posterior distribution using an independent Normal-Wishart Prior
#

postni <- function(y,x,aprior,Vprior,vprior,Sprior,Sigma){
  # declare variables
  T <- nrow(y)
  K <- ncol(y)
  z <- diag(K)%x%x
  
  # draw coefficients
  variance <- solve(Sigma)%x%diag(T)
  Vpost <- solve(solve(Vprior)+t(z)%*%variance%*%z)
  apost <- Vpost%*%(solve(Vprior)%*%aprior+t(z)%*%variance%*%as.vector(y))
  alpha <- mvrnorm(mu=apost,Sigma=Vpost)
  Alpha <- matrix(alpha,ncol=K)
  
  # draw variance
  
  vpost <- vprior+T
  res <- y-x%*%Alpha
  Spost <- solve(Sprior+t(res)%*%res)
  Sigma <- solve(rWishart(1,vpost,Spost)[,,1])

  return(list(Alpha=Alpha,Sigma=Sigma))
}

postmb <- function(y,x,aprior,Vprior,Sigma,betaols){
  
  # Declare variables
  K <- ncol(y)
  obs <- nrow(y)

  # Draw coefficienets
  Vpost <- solve(solve(Vprior)+solve(Sigma)%x%(t(x)%*%x))
  apost <- Vpost%*%(solve(Vprior)%*%aprior+(solve(Sigma)%x%(t(x)%*%x))%*%betaols)
  alpha <- mvrnorm(mu=apost,Sigma=Vpost)
  Alpha <- matrix(alpha,ncol=K)
  
  # Draw Sigmas
  e <- y-x%*%Alpha
  scale <- t(e)%*%e/obs
  Sigma <- solve(rWishart(1,obs,scale)[,,1])
  
  # return values in a list
  return(list(Alpha=Alpha,Sigma=Sigma))

                                 
}