#' @export
#' @title posterior trace plots
#' @param obj an estimated BVAR-model
#' @param lag lag for parameters (0=intercept)
#' @param ... not used
#'
#' @rdname plot_trace
plot_trace  <- function(obj,lag=1,...) UseMethod("plot_trace")

#' @export
#' @rdname plot_trace

plot_trace.bvar <- function(obj,lag=1,...){

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
  draw_no <- seq(1:nruns)

  pltlist <- list()

  # plot traces for intercept
  if(lag == 0){

    trace_coefs <- coefs[1,,]
    for(ii in 1:nvar){
      tmp.df <- data.frame(draw_no = draw_no, draw = trace_coefs[ii,])
      p1 <- ggplot2::ggplot(data=tmp.df) + ggplot2::geom_line(mapping = ggplot2::aes_(x=~draw_no,y=~draw))
      p1 <- p1 + ggplot2::ylab(obj$data_info$var_names[ii])

      if(ii == nvar){

        p1 <- p1 + ggplot2::xlab("Draw No.")

      }
      else{

        p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

      }
      # Store plots in a list
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

    # plot
    for(ii in 1:nvar){
      for(jj in 1:nvar){

        tmp.df <- data.frame(draw_no=draw_no,draw=trace_coefs[ii,jj,])
        p1 <- ggplot2::ggplot(data = tmp.df) +
          ggplot2::geom_line(mapping = ggplot2::aes_(x=~draw_no,y=~draw))

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

          p1 <- p1 + ggplot2::xlab("Draw No.")

        }
        else{

          p1 <- p1 + ggplot2::theme(axis.title.x=ggplot2::element_blank())

        }

        pltlist[[(jj - 1) * nvar + ii]] <- p1

      }
    }

  }
  # Create the final plot
  do.call("grid.arrange",c(pltlist,nrow=nvar))
}
