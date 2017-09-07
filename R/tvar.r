tvar <- function(mydata,lags=1,thMax,thresh=1,tarscale=0.5,tarstandard=NULL,intercept=TRUE,RandomWalk=TRUE,coefprior=NULL,coefpriorvar=1,varprior=10,varpriordof=10,irfhorizon=16,irfquantiles=c(0.05,0.95),reps=500,burnin=100,Restrictions=NULL){
  y <- as.matrix(mydata)
  T <- nrow(y)
  K <- ncol(y)
  obs <- T-max(thMax,lags)
  .tierror(mydata,lags,thDelay,thresh,tarscale,tarstandard,intercept,coefprior,coefpriorvar,varprior,varpriordof,irfhorizon,irfquantiles,reps,burnin,stabletest)
  prior <- .tiprior(y=y,lags=lags,thMax=thMax,intercept=intercept,coefprior=coefprior,coefpriorvar=coefpriorvar,varprior=varprior,RandomWalk=RandomWalk)
  results <- .tigibbs(y,lags,thMax,thresh,tarscale,tarstandard,intercept,prior$Aprior1,prior$Aprior2,prior$Vprior1,prior$Vprior2,prior$Sprior1,prior$Sprior2,varpriordof,irfhorizon=irfhorizon,irfquantiles=irfquantiles,reps=reps,burnin=burnin,Restrictions=Restrictions)
  retresults <-  structure(list(prior=prior,results=results),class="tvar")
  return(retresults)
}
.tierror <-function(mydata,lags=1,thDelay=1,thresh=1,tarscale=0.5,tarstandard=NULL,intercept=TRUE,coefprior=NULL,coefpriorvar=1,varprior=1,varpriordof=1,irfhorizon=16,irfquantiles=c(0.05,0.95),reps=300,burnin=100,Restrictions=Restricions){

  y <- as.matrix(mydata)
  T <- nrow(y)
  K <- ncol(y)

  constant=0
  if(intercept==TRUE) constant=1
  # check prior for coefficients
  if(!is.null(coefprior)){

    ndim1 <- K*lags+constant
    ndim2 <- K
    dim1 <- dim(coefprior)
    if(dim1[1]!=ndim1){
      print("Prior for coefficients provided by the user has an incorrect number of rows")
	  stop
    }
    if(dim1[2]!=ndim2){
      print("Prior for coefficients provided by the user has an incorrect number of columns")
	  stop
    }
  }
  if(!.isscalar(coefprior) && !is.null(coefprior)){
	 n <- K*lags+constant
     dim1 <- dim(coefpriorvar)
	 if(dim1[1]!=n){
	   print("Prior has to be an nxn matrix")
	   stop
	 }
	 if(dim1[2]!=n){
	   print("Prior has to be an nxn matrix")
	   stop
	 }
  }


}
