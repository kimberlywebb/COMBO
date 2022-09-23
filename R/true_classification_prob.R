#' Compute Probability of Each True Outcome, for Every Subject
#'
#'  Compute the probability of the latent true outcome \eqn{Y \in \{1, 2 \}} as
#'  \eqn{P(Y_i = j | X_i) = \frac{\exp(X_i \beta)}{1 + \exp(X_i \beta)}}
#'  for each of the \eqn{i = 1, \dots,} \code{n} subjects.
#'
#' @param beta_matrix A numeric column matrix of estimated regression parameters for the
#'   true outcome mechanism, \code{Y} (true outcome) ~ \code{X} (predictor matrix of interest),
#'   obtained from \code{COMBO_EM} or \code{COMBO_MCMC}.
#' @param x_matrix A numeric matrix of covariates in the true outcome mechanism.
#'   \code{x_matrix} should not contain an intercept.
#'
#' @return \code{true_classification_prob} returns a dataframe containing three columns.
#'   The first column, \code{Subject}, represents the subject ID, from \eqn{1} to \code{n},
#'   where \code{n} is the sample size, or equivalently, the number of rows in \code{x_matrix}.
#'   The second column, \code{Y}, represents a true, latent outcome category \eqn{Y \in \{1, 2 \}}.
#'   The last column, \code{Probability}, is the value of the equation
#'   \eqn{P(Y_i = j | X_i) = \frac{\exp(X_i \beta)}{1 + \exp(X_i \beta)}} computed
#'   for each subject and true, latent outcome category.
#'
#' @include pi_compute.R
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
#' x_matrix <- matrix(c(cov1, cov2), nrow = sample_size, byrow = FALSE)
#' estimated_betas <- matrix(c(1, -1, .5), ncol = 1)
#' P_Y <- true_classification_prob(estimated_betas, x_matrix)
#' head(P_Y)
true_classification_prob <- function(beta_matrix,
                                     x_matrix){

  n_cat = 2
  sample_size = nrow(x_matrix)

  if (!is.null(x_matrix)) {
    if (is.data.frame(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.numeric(x_matrix))
      stop("'x_matrix' must be numeric.")
    if (is.vector(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.matrix(x_matrix))
      stop("'x_matrix' must be a data.frame or matrix.")
  }

  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)

  subject = rep(1:sample_size, n_cat)
  Y_categories = rep(1:n_cat, each = sample_size)
  pi_matrix = pi_compute(beta_matrix, X, sample_size, n_cat)
  pi_df = data.frame(Subject = subject,
                     Y = Y_categories,
                     Probability = c(pi_matrix))

  return(pi_df)
}
