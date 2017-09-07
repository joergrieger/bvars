extract <- function(x,K){
  n <- ncol(x)
  x <- as.matrix(x)
  x.x <- t(x)%*%x
  evectors <- eigen(x.x)$vectors
  ret.evectors <- sqrt(n)*evectors[,1:K]
  fac <- x%*%ret.evectors/n
  return(list(lam=ret.evectors,fac=fac))
}