#' @export
#' @title plot residuals
#' @param obj a fitted VAR model
#' @param ... currently not used
#'

plot_residuals.bvar <- function(obj,...){

  nreps <- dim(obj$mcmc_draws$Alpha)[3]
  no_lags <- obj$general_info$nolags
  intercept <- obj$general_info$intercept

  lg <- lagdata(obj$data_info$data,intercept = obj$general_info$intercept,nolags=obj$general_info$nolags)
  nlength <- dim(lg$y)[1]

  resid <- array(0,dim=c(nlength,obj$data_info$no_variables,nreps))

  # Get residuals
  for(ii in 1:nreps){
    print(ii)

    resid[,,ii] <- lg$y - lg$x %*% obj$mcmc_draws$Alpha[,,ii]

  }
  residuals <- array(NA,dim=c(nlength,obj$data_info$no_variables,3))
  for(ii in 1:nlength){
    for(jj in 1:obj$data_info$no_variables){
      residuals[ii,jj,1] <- mean(resid[ii,jj,])
      residuals[ii,jj,2] <- quantile(resid[ii,jj,],prob=c(0.05))
      residuals[ii,jj,3] <- quantile(resid[ii,jj,],prob=c(0.95))
    }
  }

  if(stats::is.ts(obj$data_info$data)){

    timestamps <- time(obj$data_info$data)
    timestamps <- timestamps[-c(1:no_lags)]
  }
  else if(sum(class(obj$data_info$data) == "xts") > 0){

    timestamps <- time(obj$data_info$data)
    timestamps <- timestamps[-c(1:no_lags)]

  }
  else{

    timestamps <- seq(1:nlength)

  }
  pltList <- list()
  for(ii in 1:obj$data_info$no_variables){

    plot_prep_df <- data.frame(value = residuals[,ii,1],Time = timestamps,Lower = residuals[,ii,2], Upper = residuals[,ii,3])
    p1 <- ggplot2::ggplot(data = plot_prep_df) +
      ggplot2::geom_line(mapping = ggplot2::aes_(x=~Time,y=~value)) +
      ggplot2::ylab(obj$data_info$var_names[ii]) +
      ggplot2::theme(legend.position = "none")

    if(ii < obj$data_info$no_variables){
      p1 <- p1 + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    }

    pltList[[ii]]<- p1

  }
  # Create the plot
  do.call("grid.arrange",c(pltList,ncol=1,top="Residuals"))

}
