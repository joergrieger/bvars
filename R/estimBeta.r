.estimBeta <- function(y1,y2,x1,x2){
    beta01 <- solve(t(x1)%*%x1)%*%t(x1)%*%y1
	beta02 <- solve(t(x2)%*%x2)%*%t(x2)%*%y2
	SSE1 <- t(y1-x1%*%beta01)%*%(y1-x1%*%beta01)
	SSE2 <- t(y2-x2%*%beta02)%*%(y2-x2%*%beta02)
	obs1 <- nrow(y1)
	obs2 <- nrow(y2)
	SIGMA1 <- SSE1/obs1
	SIGMA2 <- SSE2/obs2
	return(list(beta01=beta01,beta02=beta02,SSE1=SSE1,SSE2=SSE2,SIGMA1=SIGMA1,SIGMA2=SIGMA2))
}