#' @export
#' @title plot residuals
#' @param obj a fitted VAR model (bvar,msvar or tvar)
#' @param ... currently not used
#' @rdname plotresids
plot_residuals <- function(obj,...) UseMethod("plot_residuals")

#' @export
#' @rdname plotresids
plot_residuals.bvar <- function(obj,...){

  nreps <- dim(obj$mcmc_draws$Alpha)[3]
  no_lags <- obj$general_info$nolags
  intercept <- obj$general_info$intercept

  lg <- lagdata(obj$data_info$data,intercept = obj$general_info$intercept,nolags=obj$general_info$nolags)
  nlength <- dim(lg$y)[1]

  resid <- array(0,dim=c(nlength,obj$data_info$no_variables,nreps))

  # Get residuals
  for(ii in 1:nreps){
    if(ii %% 10 == 0){

      print(ii)

    }


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
      ggplot2::geom_ribbon(ggplot2::aes_(ymin=~Lower,ymax=~Upper,x=~Time,alpha=0.01),fill="steelblue1")+
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

#' @export
#' @rdname plotresids
plot_residuals.msvar <- function(obj,...){

  # Preliminary stuff
  regimes <- obj$mcmc_draws$regimes
  nreps <- dim(obj$mcmc_draws$Alpha)[4]
  nT <- dim(regimes)[1]

  # Variable declarations
  residuals_tmp  <- array(NA,dim=c(nT,obj$data_info$no_variables,nreps))
  residuals <- array(NA,dim=c(nT,obj$data_info$no_variables,nreps))

  # lag data
  lg <- lagdata(obj$data_info$data,intercept = obj$general_info$intercept,nolags=obj$general_info$nolags)

  # Loop over all draws
  for(ii in 1:nreps){

    if(ii %% 10 == 0){
      print(ii)
    }

    # get residuals
    resid    <- array(0,dim=c(nT,obj$data_info$no_variables))
    for(jj in 1:nT){

      tmp_Alpha <- obj$mcmc_draws$Alpha[,,regimes[jj],ii]
      resid[jj,] <- lg$y[jj,] - lg$x[jj,] %*% tmp_Alpha
    }
    residuals_tmp[,,ii] <- resid
  }

  # Get mean and quantiles for residuals
  for(ii in 1:nT){
    for(jj in 1:obj$data_info$no_variables){
      residuals[ii,jj,1] <- mean(residuals_tmp[ii,jj,])
      residuals[ii,jj,2] <- quantile(residuals_tmp[ii,jj,],prob=c(0.05))
      residuals[ii,jj,3] <- quantile(residuals_tmp[ii,jj,],prob=c(0.95))
    }
  }

  # Create timestamps for the plot
  if(stats::is.ts(obj$data_info$data)){

    timestamps <- time(obj$data_info$data)
    timestamps <- timestamps[-c(1:obj$general_info$nolags)]
  }
  else if(sum(class(obj$data_info$data) == "xts") > 0){

    timestamps <- time(obj$data_info$data)
    timestamps <- timestamps[-c(1:obj$general_info$nolags)]

  }
  else{

    timestamps <- seq(1:nT)

  }
  # Plot residuals (maybe make an extra function for it?)
  pltList <- list()
  for(ii in 1:obj$data_info$no_variables){

    plot_prep_df <- data.frame(value = residuals[,ii,1],Time = timestamps,Lower = residuals[,ii,2], Upper = residuals[,ii,3])
    p1 <- ggplot2::ggplot(data = plot_prep_df) +
      ggplot2::geom_line(mapping = ggplot2::aes_(x=~Time,y=~value)) +
      ggplot2::geom_ribbon(ggplot2::aes_(ymin=~Lower,ymax=~Upper,x=~Time,alpha=0.01),fill="steelblue1")+
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

#' @export
#' @rdname plotresids
plot_residuals.tvar <- function(obj,...){
  nreps <- dim(obj$mcmc_draws$Alpha)[4]
  reg_length <- dim(obj$mcmc_draws$regimes)[1]
  regimes <- obj$mcmc_draws$regimes
  lg <- lagdata(obj$data_info$data, intercept = obj$general_info$intercept, nolags = obj$general_info$nolags)

  # cut off some of the data from to make the regime draws and lagged data the same length
  cut_off <- dim(lg$x)[1] - reg_length
  xdata <- lg$x[-c(1:cut_off),]
  ydata <- lg$y[-c(1:cut_off),]

  residuals <- array(0,dim=c(reg_length,obj$data_info$no_variables,3))
  residuals_tmp <- array(0,dim=c(reg_length,obj$data_info$no_variables,nreps))

  # Calculate residuals
  for(ii in 1:nreps){

    if(ii %% 10 == 0){
      print(ii)
    }
    resid    <- array(0,dim=c(reg_length,obj$data_info$no_variables))

    for(jj in 1:reg_length){
      tmp_Alpha <- obj$mcmc_draws$Alpha[,,regimes[jj]+1,ii]
      resid[jj,] <- ydata[jj,] - xdata[jj,] %*% tmp_Alpha
    }
    residuals_tmp[,,ii] <- resid
  }
  # Get mean and quantiles for residuals
  for(ii in 1:reg_length){
    for(jj in 1:obj$data_info$no_variables){
      residuals[ii,jj,1] <- mean(residuals_tmp[ii,jj,])
      residuals[ii,jj,2] <- quantile(residuals_tmp[ii,jj,],prob=c(0.05))
      residuals[ii,jj,3] <- quantile(residuals_tmp[ii,jj,],prob=c(0.95))
    }
  }
  # Create timestamps for the plot
  if(stats::is.ts(obj$data_info$data)){

    timestamps <- time(obj$data_info$data)
    nT <- dim(obj$data_info$data)[1]
    cut_off <- nT-reg_length
    timestamps <- timestamps[-c(1:cut_off)]
  }
  else if(sum(class(obj$data_info$data) == "xts") > 0){

    timestamps <- time(obj$data_info$data)
    nT <- dim(obj$data_info$data)[1]
    cut_off <- nT-reg_length
    timestamps <- timestamps[-c(1:cut_off)]

  }
  else{

    timestamps <- seq(1:reg_length)

  }
  # Plot residuals (maybe make an extra function for it?)
  pltList <- list()
  for(ii in 1:obj$data_info$no_variables){

    plot_prep_df <- data.frame(value = residuals[,ii,1],Time = timestamps,Lower = residuals[,ii,2], Upper = residuals[,ii,3])
    p1 <- ggplot2::ggplot(data = plot_prep_df) +
      ggplot2::geom_line(mapping = ggplot2::aes_(x=~Time,y=~value)) +
      ggplot2::geom_ribbon(ggplot2::aes_(ymin=~Lower,ymax=~Upper,x=~Time,alpha=0.01),fill="steelblue1")+
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
