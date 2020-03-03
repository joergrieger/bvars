#' @export
#' @title autocorrelation of posterior draws
#' @param obj an estimated BVAR-model
#' @param lag lag for parameters (0=intercept)
#' @param maxlag maximum lag for autocorrelation
#' @param ... not used
#' @rdname plot_autocorr
plot_autocorr  <- function(obj,lag=1,maxlag=20,...) UseMethod("plot_autocorr")

#' @export
#' @rdname plot_autocorr
#'
plot_autocorr.bvar  <- function(obj,lag=1,maxlag=20,...){


  # input check

  # want to plot traces of intercept but no intercept -> stop
  intercept = obj$general_info$intercept
  if((lag == 0) && (intercept == FALSE)){

    stop("Model has no intercept")

  }
  # want to plot traces of lag greater than used in the model -> stop
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
  draw_no <- seq(0:maxlag) - 1

  pltlist <- list()

  if(lag == 0){

    trace_coefs <- coefs[1,,]
    for(ii in 1:nvar){
      # estimate autocorrelation
      coefs <- trace_coefs[ii,]
      tmp <- stats::acf(coefs,lag.max=maxlag)
      tmp.df <- data.frame(lag = draw_no, acf = tmp$acf)

      # plot autocorrelation
      p1 <- ggplot2::ggplot(data=tmp.df) +ggplot2::geom_point(mapping = ggplot2::aes_(x=~lag,y=~acf))

      p1 <- p1 + ggplot2::ylab(obj$data_info$var_names[ii])

      if(ii == nvar){

        p1 <- p1 + ggplot2::xlab("Lag")

      }
      else{

        p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

      }
      # Store plots in a list
      pltlist[[ii]] <- p1

    }

  }
  if(lag > 0 ){

    # get coefficients
    const <- 0
    if(intercept == TRUE) constant <- 1
    nlow <- (lag - 1) * nvar + 1 + const
    nhigh <- lag * nvar + const
    trace_coefs <- coefs[nlow:nhigh,,]

    for(ii in 1:nvar){
      for(jj in 1:nvar){

        coefs <- trace_coefs[ii,jj,]
        tmp <- stats::acf(coefs,lag.max=maxlag)
        tmp.df <- data.frame(lag = draw_no, acf = tmp$acf)
        p1 <- ggplot2::ggplot(data=tmp.df) + ggplot2::geom_point(mapping = ggplot2::aes_(x=~lag,y=~acf))

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

          p1 <- p1 + ggplot2::xlab("Lag")

        }
        else{

          p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

        }

        pltlist[[(jj - 1) * nvar + ii]] <- p1

      }
    }
  }
  do.call("grid.arrange",c(pltlist,nrow=nvar))

}
