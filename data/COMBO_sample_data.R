#' Sample data set with misclassified binary outcomes
#'
#' A list containing true outcomes, observed (and potentially misclassified) outcomes,
#' predictors in the true outcome mechanism, and predictors in the observation mechanism.
#' Data is generated from the \code{COMBO_data} function.
#'
#' @format A list with 7 elements, containing the following characteristics for 500 observations:
#' \describe{
#'   \item{obs_Y}{A vector of observed outcomes.}
#'   \item{true_Y}{A vector of true outcomes.}
#'   \item{obs_Y_matrix}{A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row contains
#'   exactly one 0 entry and exactly one 1 entry.}
#'   \item{x}{A vector of generated predictor values in the true outcome
#'   mechanism, from the Normal distribution.}
#'   \item{z}{A vector of generated predictor values in the observation
#'   mechanism from the Gamma distribution.}
#'   \item{x_design_matrix}{The design matrix for the \code{x} predictor.}
#'   \item{z_design_matrix}{The design matrix for the \code{z} predictor.}
#' }
