facrot <- function(F,Ffast,Fslow){
  
  Ffast <- as.matrix(Ffast)
  k1    <- ncol(Ffast)
  fm    <- cbind(matrix(1,nrow(Ffast),1), Ffast, Fslow)
  b     <- olssvd(F,fm)
  Fr <- F-Ffast%*%b[2:(k1+1),]
  return(Fr)
}