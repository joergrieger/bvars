draw_posterior  <- function(obj,...) UseMethod("draw_posterior")
structural <- function(id_obj,...) UseMethod("structural")
initialize_mcmc <- function(obj,...) UseMethod("initialize_mcmc")

#' @export
#' @title plot residuals
#' @param obj a fitted VAR model
#' @param ... currently not used
plot_residuals <- function(obj,...) UseMethod("plot_residuals")
