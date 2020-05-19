#' @export
#' @title posterior density plots
#' @param obj an estimated BVAR-model
#' @param lag lag for parameters (0=intercept)
#' @param hpd level for credible interval
#' @param ... not used
#' @rdname plot_density
plot_density  <- function(obj,lag=1,hpd=NULL,...) UseMethod("plot_density")

#' @export
#' @rdname plot_density

plot_density.bvar <- function(obj,lag=1,hpd=NULL,...){
  # input check

  # want to plot density of intercept but no intercept -> stop
  intercept = obj$general_info$intercept
  if((lag == 0) && (intercept == FALSE)){

    stop("Model has no intercept")

  }
  # want to plot density of lag greater than used in the model -> stop
  if(lag > obj$general_info$nolags){

    stop("Number of lags is lower!")

  }
  # Sanity check: is lag-number of positive
  if(lag<0){

    stop("Number of lags must be positive")

  }

  coefs <- obj$mcmc_draws$Alpha
  ndim  <- dim(coefs)
  nruns <- ndim[3] # number of draws
  nvar  <- ndim[2] # number of variables

  pltlist <- list()

  if(lag == 0){
    trace_coefs <- coefs[1,,]
    for(ii in 1:nvar){

      tmp.df <- data.frame(draws = trace_coefs[ii,])
      p1 <- ggplot2::ggplot(data = tmp.df) + ggplot2::geom_density(mapping = ggplot2::aes_(x=~draws),fill="lightblue",alpha=0.75,size=1)
      p1 <- p1 + ggplot2::ylab(obj$data_info$var_names[ii])

      # Compute highest posterior density
      if(!is.null(hpd)){

        upper <- 1-hpd/2
        lower <- hpd/2
        lower_quantile <- stats::quantile(tmp.df$draws,probs=lower)
        upper_quantile <- stats::quantile(tmp.df$draws,probs=upper)

        p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=lower_quantile),color="dodgerblue1",size=1.1,linetype="dashed")
        p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=upper_quantile),color="dodgerblue1",size=1.1,linetype="dashed")

      }

      if(ii == nvar){

        p1 <- p1 + ggplot2::xlab("Value")

      }
      else{

        p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

      }
      pltlist[[ii]] <- p1

    }
  }
  if(lag > 0){
    # get coefficients
    const <- 0
    if(intercept == TRUE) constant <- 1
    nlow <- (lag - 1) * nvar + 1 + const
    nhigh <- lag * nvar + const
    trace_coefs <- coefs[nlow:nhigh,,]
    for(ii in 1:nvar){
      for(jj in 1:nvar){

        tmp.df <- data.frame(draws=trace_coefs[ii,jj,])
        p1 <- ggplot2::ggplot(data = tmp.df) + ggplot2::geom_density(mapping = ggplot2::aes_(x=~draws),fill="lightblue",alpha=0.75,size=1)

        if(!is.null(hpd)){

          upper <- 1-hpd/2
          lower <- hpd/2
          lower_quantile <- stats::quantile(tmp.df$draws,probs=lower)
          upper_quantile <- stats::quantile(tmp.df$draws,probs=upper)

          p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=lower_quantile),color="dodgerblue1",size=1.1,linetype="dashed")
          p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=upper_quantile),color="dodgerblue1",size=1.1,linetype="dashed")

        }

        if(ii == 1){

          p1 <- p1 + ggplot2::ylab(obj$data_info$var_names[jj])

        }
        else{

          p1 <- p1 + ggplot2::theme(axis.title.y = ggplot2::element_blank())

        }

        if(jj == 1){

          p1 <- p1 + ggplot2::ggtitle(obj$data_info$var_names[ii]) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

        }
        if(jj == nvar){

          p1 <- p1 + ggplot2::xlab("Value")

        }
        else{

          p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

        }

        pltlist[[(jj - 1) * nvar + ii]] <- p1

      }
    }
  }
  # Draw plots
  do.call("grid.arrange",c(pltlist,nrow=nvar))
}

#' @export
#' @param regime Regime for the density plot
#' @rdname plot_density

plot_density.msvar <- function(obj,lag=1,hpd=NULL,regime=1,...){

  # want to plot density of intercept but no intercept -> stop
  intercept = obj$general_info$intercept
  if((lag == 0) && (intercept == FALSE)){

    stop("Model has no intercept")

  }
  # want to plot density of lag greater than used in the model -> stop
  if(lag > obj$general_info$nolags){

    stop("Number of lags is lower!")

  }
  # Sanity check: is lag-number of positive
  if(lag<0){

    stop("Number of lags must be positive")

  }

  # Possible regime?
  if(regime > obj$general_info$noregimes){
    stop("Regime selected is higher than the number of regimes in the model")
  }
  if(regime < 1){

    stop("Selected Regime must be strictly positive")

  }

  coefs <- obj$mcmc_draws$Alpha[,,regime,]
  ndim  <- dim(coefs)
  nruns <- ndim[3] # number of draws
  nvar  <- ndim[2] # number of variables

  pltlist <- list()

  if(lag == 0){
    trace_coefs <- coefs[1,,]
    for(ii in 1:nvar){

      tmp.df <- data.frame(draws = trace_coefs[ii,])
      p1 <- ggplot2::ggplot(data = tmp.df) + ggplot2::geom_density(mapping = ggplot2::aes_(x=~draws),fill="lightblue",alpha=0.75,size=1)
      p1 <- p1 + ggplot2::ylab(obj$data_info$var_names[ii])

      # Compute highest posterior density
      if(!is.null(hpd)){

        upper <- 1-hpd/2
        lower <- hpd/2
        lower_quantile <- stats::quantile(tmp.df$draws,probs=lower)
        upper_quantile <- stats::quantile(tmp.df$draws,probs=upper)

        p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=lower_quantile),color="dodgerblue1",size=1.1,linetype="dashed")
        p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=upper_quantile),color="dodgerblue1",size=1.1,linetype="dashed")

      }

      if(ii == nvar){

        p1 <- p1 + ggplot2::xlab("Value")

      }
      else{

        p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

      }
      pltlist[[ii]] <- p1

    }
  }
  if(lag > 0){
    # get coefficients
    const <- 0
    if(intercept == TRUE) constant <- 1
    nlow <- (lag - 1) * nvar + 1 + const
    nhigh <- lag * nvar + const
    trace_coefs <- coefs[nlow:nhigh,,]
    for(ii in 1:nvar){
      for(jj in 1:nvar){

        tmp.df <- data.frame(draws=trace_coefs[ii,jj,])
        p1 <- ggplot2::ggplot(data = tmp.df) + ggplot2::geom_density(mapping = ggplot2::aes_(x=~draws),fill="lightblue",alpha=0.75,size=1)

        if(!is.null(hpd)){

          upper <- 1-hpd/2
          lower <- hpd/2
          lower_quantile <- stats::quantile(tmp.df$draws,probs=lower)
          upper_quantile <- stats::quantile(tmp.df$draws,probs=upper)

          p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=lower_quantile),color="dodgerblue1",size=1.1,linetype="dashed")
          p1 <- p1 + ggplot2::geom_vline(mapping = ggplot2::aes_(xintercept=upper_quantile),color="dodgerblue1",size=1.1,linetype="dashed")

        }

        if(ii == 1){

          p1 <- p1 + ggplot2::ylab(obj$data_info$var_names[jj])

        }
        else{

          p1 <- p1 + ggplot2::theme(axis.title.y = ggplot2::element_blank())

        }

        if(jj == 1){

          p1 <- p1 + ggplot2::ggtitle(obj$data_info$var_names[ii]) + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))

        }
        if(jj == nvar){

          p1 <- p1 + ggplot2::xlab("Value")

        }
        else{

          p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

        }

        pltlist[[(jj - 1) * nvar + ii]] <- p1

      }
    }
  }
  # Draw plots
  do.call("grid.arrange",c(pltlist,nrow=nvar))
}

