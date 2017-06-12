.SplitVariables <- function(y,lags,thDelay,thresh,tart,intercept){
    startest <- max(thDelay,lags)
	T <- nrow(y)
	K <- ncol(y)
    ytest <- y[(startest+1-thDelay):(T-thDelay),thresh]
	ystar <- y[(startest+1):T,]
	xstar <- embed(y,dimension=lags+1)[,-(1:K)]
	if(thDelay>lags){
	    diff1=(thDelay-lags)+1
		xnr <- nrow(xstar)
		xstar <- xstar[diff1:xnr,]
	}
	
	e1 <- ytest < tart
	e2 <- ytest >=tart
	y1 <- ystar[e1,]
	y2 <- ystar[e2,]
	x1 <- xstar[e1,]
	x2 <- xstar[e2,]
	if(intercept==TRUE){
	   x1 <- cbind(1,x1)
	   x2 <- cbind(1,x2)
	}
	return(list(y1=y1,y2=y2,x1=x1,x2=x2,xstar=xstar,ytest=ytest,ystar=ystar,e1=e1))
	   
}