.tigibbs <- function(y,lags,thMax,thresh,tarscale=1,tarstandard=NULL,intercept=TRUE,Aprior1,Aprior2,Vprior1,Vprior2,Sprior1,Sprior2,varpriordof,irfhorizon=16,irfquantiles=c(0.1,0.9),reps=210,burnin=10,stabletest=FALSE){

    # Vectorize Prior
	aprior1 <- matrix(Aprior1,ncol=1)
	aprior2 <- matrix(Aprior2,ncol=1)
	#some preliminary calculations
	thDelay=thMax
	tard <- seq(1:thMax)
	startest <- max(thMax,lags)
	
	T <- nrow(y)
	K <- ncol(y)
	NoRegimes <- T-(startest+1)
	regimes  <- array(0,dim=c(reps-burnin,NoRegimes))
	constant=0
	if(intercept==TRUE) constant=1
	ytest <- y[(startest+1-thDelay):(T-thDelay),thresh]
	tart <- mean(ytest)
	tarmean <- tart
	ystar <- y[(startest+1):T,]
	xstar <- embed(y,dimension=lags+1)[,-(1:K)]
	# Initialize the Gibbs sampler
	xsplit <- .SplitVariables(y,lags,thDelay,thresh,tart,intercept)
	beta1 <- .estimBeta(xsplit$y1,xsplit$y2,xsplit$x1,xsplit$x2)
	beta01 <- beta1$beta01
	beta02 <- beta1$beta02
	SIGMA1 <- beta1$SIGMA1
	SIGMA2 <- beta1$SIGMA2

	# Define Variables
	ALPHAdraws1 <- array(0,dim=c(K*lags+constant,K,reps-burnin))
	ALPHAdraws2 <- array(0,dim=c(K*lags+constant,K,reps-burnin))
	SIGMAdraws1 <- array(0,dim=c(K,K,reps-burnin))
	SIGMAdraws2 <- array(0,dim=c(K,K,reps-burnin))
	tardraws    <- array(0,dim=c(reps-burnin))
	regdraws    <- array(0,dim=c(reps-burnin))
	irf1draws   <- array(0,dim=c(K,irfhorizon,K,reps-burnin))
	irf2draws   <- array(0,dim=c(K,irfhorizon,K,reps-burnin))
	if(is.null(tarstandard)){
	    tarstandard <- tarscale*sqrt(var(ytest))
	}
	# Start the Gibbs sampling
	ii=1
	while(ii<=reps){
	  cat(ii)
		cat("\n")
	    # Step 1: Split the data
		xsplit <- .SplitVariables(y,lags,thDelay,thresh,tart,intercept)
		obs1 <- nrow(xsplit$y1)
		obs2 <- nrow(xsplit$y2)
		x1 <- xsplit$x1
		x2 <- xsplit$x2
		y1 <- xsplit$y1
		y2 <- xsplit$y2
		e1 <- xsplit$e1
		xstar <- xsplit$xstar
		ystar <- xsplit$ystar
		ytest <- xsplit$ytest
		# Step 2: Get coefficients
		# 2a: Coefficients for the first regime
		z1 <- .id(K)%x%x1
		variance1 <- solve(SIGMA1)%x%.id(obs1)
		Vpost01   <- solve(solve(Vprior1)+t(z1)%*%variance1%*%z1)
		apost01   <- Vpost01%*%(solve(Vprior1)%*%aprior1+t(z1)%*%variance1%*%matrix(y1,ncol=1))

		alpha1  <- .getcoef(K,lags,apost01,Vpost01,intercept)
		if(alpha1$Problem){
		  alpha01 <- alpha1$Alpha01
		  Alpha01 <- alpha1$Alpha01
		}
		
		stable1 <- alpha1$Problem
		vpost01    <- obs1+varpriordof
		Spost01    <- Sprior1+t(y1-x1%*%Alpha01)%*%(y1-x1%*%Alpha01)
		SIGMA1     <- solve(rWishart(1,vpost01,Spost01)[,,1])
		ALPHA1 <- Alpha01

		# 2b: Coefficients for the second regime
		z2 <- .id(K)%x%x2
		variance2 <- solve(SIGMA2)%x%.id(obs2)
		Vpost02   <- solve(solve(Vprior2)+t(z2)%*%variance2%*%z2)
		apost02   <- Vpost02%*%(solve(Vprior2)%*%aprior2+t(z2)%*%variance2%*%matrix(y2,ncol=1))

		alpha2  <- .getcoef(K,lags,apost01,Vpost01,intercept)
		if(alpha2$Problem){
		  alpha02 <- alpha2$Alpha01
		  Alpha02 <- alpha2$Alpha01
		}

		stable2 <- alpha2$Problem
		vpost02 <- obs2+varpriordof
		Spost02 <- Sprior2+t(y2-x2%*%Alpha02)%*%(y2-x2%*%Alpha02)
		SIGMA2  <- solve(rWishart(1,vpost02,Spost02)[,,1])
		ALPHA2  <- Alpha02

		# Step 3: Sample new threshold using a Random Walk Metropolis-Hastings Algorithm
		tarnew <- tart+rnorm(1,sd=tarstandard)
		l1post <- .tarpost(xstar,Ystar=ystar,Ytest=ytest,ALPHA1,ALPHA2,SIGMA1,SIGMA2,tarnew,lags,intercept=intercept,tarmean,tarstandard,ncrit=0.15)
		l2post <- .tarpost(xstar,Ystar=ystar,Ytest=ytest,ALPHA1,ALPHA2,SIGMA1,SIGMA2,tart,lags,intercept=intercept,tarmean,tarstandard,ncrit=0.15)
		acc=min(1,exp(l1post$post-l2post$post))
		u <- runif(1)
		if(u<acc){
		  tart=tarnew
		}
		tarmean=tart
		# Step 4: Sample delay parameter
		prob <- matrix(0,nrow=thMax)
		for(jj in 1:thMax){
		  split1 <- .SplitVariables(y,lags,jj,thresh,tart,intercept)
			x <- .exptarpost(split1$xstar,split1$ystar,split1$ytest,ALPHA1,ALPHA2,SIGMA1,SIGMA2,tart,lags,intercept,tarmean,tarstandard,ncrit=0.15)
			prob[jj,1] <- x$post
		}
		prob <- prob/sum(prob)
		thDelay<-sample(thMax,1,replace=FALSE,prob)
		# Only use sample data if there were no problems
		problem <- stable1 && stable2
		if(ii>burnin && problem){
		  # Saving Values
		  ALPHAdraws1[,,ii-burnin] <- ALPHA1
			ALPHAdraws2[,,ii-burnin] <- ALPHA2
		  SIGMAdraws1[,,ii-burnin] <- SIGMA1
			SIGMAdraws2[,,ii-burnin] <- SIGMA2
			tardraws[ii-burnin] <- tart
			regdraws[ii-burnin] <- thDelay
			start <- thMax-thDelay
			nT <- length(e1)
			a<-nT-NoRegimes
			regimes[ii-burnin,]<-e1[(1+a):nT]
			#print(regimes[ii-burnin,])
			#readline(prompt="Press [enter] to continue")
			# Calculating Impulse-Response functions
			if(intercept==TRUE){
			  beta1 <- ALPHA1[2:(K*lags+1),]
			  beta2 <- ALPHA2[2:(K*lags+1),]
			}
			else{
			  beta1 <- ALPHA1
			  beta2 <- ALPHA2
			}
			for(ll in 1:K){
			  shockvar=ll
			  tirf1 <- .tirf(K,irfhorizon,lags,ytest,xstar,beta1,beta2,SIGMA1,SIGMA2,tart,thresh,thDelay,shockvar)
			  irf1draws[,,ll,ii-burnin]<-tirf1$irf1
			  irf2draws[,,ll,ii-burnin]<-tirf1$irf2
			}


		}
		if(problem){
		  ii <- ii+1
		}

	}
	upperquantile <- max(irfquantiles)
	lowerquantile <- min(irfquantiles)
	irffinal1 <- array(0,dim=c(K,irfhorizon,K,3))
	irffinal2 <- array(0,dim=c(K,irfhorizon,K,3))
	for(ii in 1:K){
	  for(jj in 1:K){
	    for(kk in 1:irfhorizon){
	      irffinal1[jj,kk,ii,1]<-quantile(irf1draws[jj,kk,ii,],probs=0.5)
	      irffinal1[jj,kk,ii,2]<-quantile(irf1draws[jj,kk,ii,],probs=lowerquantile)
	      irffinal1[jj,kk,ii,3]<-quantile(irf1draws[jj,kk,ii,],probs=upperquantile)
	      irffinal2[jj,kk,ii,1]<-quantile(irf2draws[jj,kk,ii,],probs=0.5)
	      irffinal2[jj,kk,ii,2]<-quantile(irf2draws[jj,kk,ii,],probs=lowerquantile)
	      irffinal2[jj,kk,ii,3]<-quantile(irf2draws[jj,kk,ii,],probs=upperquantile)
	    }
	  }
	}
	probRegimes <- array(0,dim=c(NoRegimes))
	for(ii in 1:NoRegimes){
	  probRegimes[ii] <- mean(regimes[,ii])
	}
	return(list(tardraws=tardraws,regdraws=regdraws,irf1=irffinal1,irf2=irffinal2,sigma1=SIGMAdraws1,sigma2=SIGMAdraws2,alpha1=ALPHAdraws1,alpha2=ALPHAdraws2,probRegimes=probRegimes))
}
