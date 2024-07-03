#' Compute Conditional Probability of Each Second-Stage Observed Outcome Given Each True Outcome and First-Stage Observed Outcome, for Every Subject
#'
#'  Compute the conditional probability of observing second-stage outcome \eqn{Y^{*(2)} \in \{1, 2 \}} given
#'  the latent true outcome \eqn{Y \in \{1, 2 \}} and the first-stage outcome \eqn{Y^{*(1)} \in \{1, 2\}} as
#'  \eqn{\frac{\text{exp}\{\gamma^{(2)}_{\ell kj0} + \gamma^{(2)}_{\ell kjZ^{(2)}} Z^{(2)}\}}{1 + \text{exp}\{\gamma^{(2)}_{\ell kj0} + \gamma^{(2)}_{\ell kjZ^{(2)}} Z^{(2)}_i\}}}
#'  for each of the \eqn{i = 1, \dots,} \eqn{n} subjects.
#'
#' @param gamma2_array A numeric array of estimated regression parameters for the
#'   observation mechanism, \eqn{Y^{*(2)}| Y^{*(1)}, Y} (second-stage observed outcome,
#'   given the first-stage observed outcome and the true outcome)
#'   ~ \eqn{Z^{(2)}} (second-stage misclassification predictor matrix). Rows of the array
#'   correspond to parameters for the \eqn{Y^{*(2)} = 1} observed outcome, with the
#'   dimensions of \code{z2_matrix}. Columns of the array correspond to the first-stage
#'   outcome categories \eqn{k = 1, \dots,} \code{n_cat}. The third stage of the array
#'   corresponds to the true outcome categories \eqn{j = 1, \dots,} \code{n_cat}.
#'   The array should be obtained by \code{COMBO_EM} or \code{COMBO_MCMC}.
#' @param z2_matrix A numeric matrix of covariates in the second-stage observation mechanism.
#'   \code{z2_matrix} should not contain an intercept.
#'
#' @return \code{misclassification_prob2} returns a dataframe containing five columns.
#'   The first column, \code{Subject}, represents the subject ID, from \eqn{1} to \code{n},
#'   where \code{n} is the sample size, or equivalently, the number of rows in \code{z2_matrix}.
#'   The second column, \code{Y}, represents a true, latent outcome category \eqn{Y \in \{1, 2 \}}.
#'   The third column, \code{Ystar1}, represents a first-stage observed outcome category \eqn{Y^{*(1)} \in \{1, 2 \}}.
#'   The fourth column, \code{Ystar2}, represents a second-stage observed outcome category \eqn{Y^{*(2)}  \in \{1, 2 \}}.
#'   The last column, \code{Probability}, is the value of the equation
#'   \eqn{\frac{\text{exp}\{\gamma^{(2)}_{\ell kj0} + \gamma^{(2)}_{\ell kjZ^{(2)}} Z^{(2)}\}}{1 + \text{exp}\{\gamma^{(2)}_{\ell kj0} + \gamma^{(2)}_{\ell kjZ^{(2)}} Z^{(2)}_i\}}}
#'   computed for each subject, first-stage observed outcome category, second-stage
#'   observed outcome category, and true, latent outcome category.
#'
#' @include pitilde_compute.R
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' sample_size <- 1000
#' cov1 <- rnorm(sample_size)
#' cov2 <- rnorm(sample_size, 1, 2)
#' z2_matrix <- matrix(c(cov1, cov2), nrow = sample_size, byrow = FALSE)
#' estimated_gamma2 <- array(c(1, -1, .5, .2, -.6, 1.5,
#'                             -1, .5, -1, -.5, -1, -.5), dim = c(3,2,2))
#' P_Ystar2_Ystar1_Y <- misclassification_prob2(estimated_gamma2, z2_matrix)
#' head(P_Ystar2_Ystar1_Y)
misclassification_prob2 <- function(gamma2_array,
                                    z2_matrix){

  delta_array <- gamma2_array
  v_matrix <- z2_matrix

  n_cat = 2
  sample_size = nrow(v_matrix)

  if (is.data.frame(v_matrix))
    v_matrix <- as.matrix(v_matrix)
  if (!is.numeric(v_matrix))
    stop("'z2_matrix' should be a numeric matrix.")

  if (is.vector(v_matrix))
    v_matrix <- as.matrix(v_matrix)
  if (!is.matrix(v_matrix))
    stop("'z2_matrix' should be a matrix or data.frame.")

  V = matrix(c(rep(1, sample_size), c(v_matrix)),
             byrow = FALSE, nrow = sample_size)

  subject = rep(1:sample_size, n_cat * n_cat)
  Y_categories = rep(1:n_cat, each = sample_size * n_cat * n_cat)
  Ystar_categories = rep(c(1:n_cat, 1:n_cat), each = sample_size * n_cat)
  Ytilde_categories = rep(rep(1:n_cat, each = sample_size), n_cat * n_cat)
  pitilde_array = pitilde_compute(delta_array, V, sample_size, n_cat)
  pitilde_df = data.frame(Subject = subject,
                          Y = Y_categories,
                          Ystar1 = Ystar_categories,
                          Ystar2 = Ytilde_categories,
                          Probability = c(pitilde_array))

  return(pitilde_df)
}
