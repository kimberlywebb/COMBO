#' Test data for the COMBO_EM function
#'
#' A dataset for testing the \code{COMBO_EM} function, generated from the
#' \code{COMBO_data} function.
#'
#'
#' @format A list containing 6 variables for 1000 observations:
#' \describe{
#'   \item{Y}{The true outcome variable}
#'   \item{Ystar}{The observed outcome variable}
#'   \item{x_matrix}{A matrix of predictor values in the true outcome mechanism}
#'   \item{z_matrix}{A matrix of predictor values in the observed outcome mechanism}
#'   \item{true_beta}{Beta parameter values used for data generation in the true outcome mechanism}
#'   \item{true_gamma}{Gamma parameter values used for data generation in the observed outcome mechanism}
#' }
#'
#'
#' @examples
#' data("COMBO_EM_data")
#' head(COMBO_EM_data)
#'
"COMBO_EM_data"
