draw_posterior  <- function(obj,...) UseMethod("draw_posterior")
structural <- function(id_obj,...) UseMethod("structural")

initialize_mcmc <- function(obj,...) UseMethod("initialize_mcmc")
#' @export
#' @title Generic method for calculating Impulse-Response-Functions
#' @param obj an S3-object with an estimated VAR model
#' @param ... currently not used
#' @return returns an S3 object with an estimate of the impulse-response function
irf             <- function(obj,...) UseMethod("irf")

#' @export
#' @title compute forecasts
#' @param obj an S3-object with an estimated VAR model
#' @param ... currently not used
#' @return returns an S3 object with the forecasts
forecast <- function (obj,...) UseMethod("forecast")

#' @export
#' @title Historical decompositions
#' @param obj a fitted VAR model
#' @param ... currently not used
#' @returns returns an S3 object with the historical decompositions
#'
hd <- function(obj,...) UseMethod("hd")


#' @export
#' @title Forecast error variance decomposition
#' @param obj a fitted VAR model
#' @param ... currently not used
#' @returns returns an S3 object with the historical decompositions
#'
fevd <- function(obj,...) UseMethod("fevd")
