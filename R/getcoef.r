.getcoef <- function(K,lags,apost01,Vpost01,intercept,maxeval=10000){
     stable=2
	 runs=1
    while(runs<=maxeval && stable>1){
	    
        alpha01 <- mvrnorm(mu=apost01,Sigma=Vpost01)
		Alpha01 <- matrix(alpha01,ncol=K)
		if(intercept==TRUE){
		    Alphatest <- Alpha01[2:(K*lags+1),]
		 }
		 else{
		     Alphatest <- Alpha01
		 }
		 stable <- .stability(Alphatest,lags,K)
		 runs <- runs+1
	 }
	 Problem=TRUE
	 if(runs >= maxeval ){
	     Problem=FALSE
	 }
	 return(list(alpha01=alpha01,Alpha01=Alpha01,Problem=Problem))
}