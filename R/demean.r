demean <- function(x){
  x.nc <- ncol(x)
  for(i in 1:x.nc){
    x.mean <- mean(x[,i])
    x[,i]<-x[,i]-x.mean
  }
  return(x)
}