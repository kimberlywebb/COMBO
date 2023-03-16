#' EM-Algorithm Estimation of the Binary Outcome Misclassification Model while Assuming Perfect Sensitivity
#'
#' @param Ystar A numeric vector of indicator variables (1, 0) for the observed
#'   outcome \code{Y*}. The reference category is 0.
#' @param Z A numeric matrix of covariates in the true outcome mechanism.
#'   \code{Z} should not contain an intercept.
#' @param X A numeric matrix of covariates in the observation mechanism.
#'   \code{X} should not contain an intercept.
#' @param start Numeric vector of starting values for parameters in the true
#'   outcome mechanism (\eqn{\theta}) and the observation mechanism (\eqn{\beta}), respectively.
#' @param beta0_fixed Optional numeric vector of values of the observation mechanism
#'   intercept to profile over. If a single value is entered, this corresponds to
#'   fixing the intercept at the specified value. The default is \code{NULL}.
#' @param weights Optional vector of row-specific weights used for selection bias
#'   adjustment. The default is \code{NULL}.
#' @param expected A logical value indicating whether or not to calculate the
#'   covariance matrix via the expected Fisher information matrix. The default is \code{TRUE}.
#' @param tolerance A numeric value specifying when to stop estimation, based on
#'   the difference of subsequent log-likelihood estimates. The default is \code{1e-7}.
#' @param max_em_iterations An integer specifying the maximum number of
#'   iterations of the EM algorithm. The default is \code{1500}.
#'
#' @return \code{perfect_sensitivity_EM} returns a list containing nine elements.
#'   The elements are detailed in \code{?SAMBA::obsloglikEM} documentation. Code
#'   is adapted from the \code{SAMBA::obsloglikEM} function.
#'
#' @references Beesley, L. and Mukherjee, B. (2020).
#'   Statistical inference for association studies using electronic health records:
#'   Handling both selection bias and outcome misclassification.
#'   Biometrics, 78, 214-226.
#'
#' @include expit.R
#'
#' @importFrom stats rnorm rgamma rmultinom glm predict coef
#'
#' @examples \dontrun{
#' set.seed(123)
#' n <- 1000
#' x_mu <- 0
#' x_sigma <- 1
#' z_shape <- 1
#'
#' true_beta <- matrix(c(1, -.5), ncol = 1)
#' true_gamma <- matrix(c(1, 10, .5, -3), nrow = 2, byrow = FALSE)
#'
#' my_data <- COMBO_data(sample_size = n,
#'                       x_mu = x_mu, x_sigma = x_sigma,
#'                       z_shape = z_shape,
#'                       beta = true_beta, gamma = true_gamma)
#'
#' Ystar = ifelse(my_data[["obs_Y"]] == 1, 1, 0)
#' x_matrix = matrix(my_data[["x"]], ncol = 1)
#' z_matrix = matrix(my_data[["z"]], ncol = 1)
#'
#' starting_values <- rep(1,6)
#' beta_start <- matrix(starting_values[1:2], ncol = 1)
#' gamma_start <- matrix(starting_values[3:4], ncol = 1, nrow = 2, byrow = FALSE)
#'
#' perfect_sensitivity_results <- perfect_sensitivity_EM(Ystar,
#'                                                       Z = x_matrix, X = z_matrix,
#'                                                       start = c(beta_start, gamma_start),
#'                                                       beta0_fixed = NULL)
#' perfect_sensitivity_results$param
#' }
perfect_sensitivity_EM <- function(Ystar, Z, X, start, beta0_fixed = NULL,
                                   weights = NULL, expected = TRUE,
                                   tolerance = 1e-7, max_em_iterations = 1500)
{
  if (is.data.frame(Z))
    Z <- as.matrix(Z)
  if (!is.numeric(Z))
    stop("'Z' should be a numeric matrix.")

  if (is.vector(Z))
    Z <- as.matrix(Z)
  if (!is.matrix(Z))
    stop("'Z' should be a matrix or data.frame.")

  if (!is.null(X)) {
    if (is.data.frame(X))
      X <- as.matrix(X)
    if (!is.numeric(X))
      stop("'X' must be numeric.")
    if (is.vector(X))
      X <- as.matrix(X)
    if (!is.matrix(X))
      stop("'X' must be a data.frame or matrix.")
  }

  if (!is.numeric(Ystar) || !is.vector(Ystar))
    stop("'Ystar' must be a numeric vector.")
  if (length(setdiff(0:1, unique(Ystar))) != 0)
    stop("'Ystar' must be coded 0/1.")

  n <- length(Ystar)
  if (nrow(Z) != n)
    stop("The number of rows of 'Z' must match the length of 'Ystar'.")
  if (!is.null(X) && nrow(X) != n)
    stop("The number of rows of 'X' must match the length of 'Ystar'.")

  if (!is.logical(expected) || length(expected) > 1)
    stop("'expected' must be a length one logical.")

  # initialise p for EM
  theta <- start[1:(1 + ncol(Z))]
  beta  <- start[-(1:(1 + ncol(Z)))]
  pred1 <- expit(cbind(1, Z) %*% theta)
  pred2 <- expit(cbind(1, X) %*% beta)

  calculate.p <- function(pred1, pred2)
  {
   (Ystar) * (pred1 / (pred1 + (pred2 * (1 - pred1))))
  }
  p <- calculate.p(pred1, pred2)

  it <- 1
  converged <- F

  w <- 1

  param.seq  <- matrix(c(theta, beta), 1)
  loglik.seq <- -10 ^ 9
  while (!converged && it < max_em_iterations) {
    if (is.null(beta0_fixed)) {
      suppressWarnings({
        fit.beta <- stats::glm(Ystar ~ X, weights = (1 - p) * w,
                               family = stats::binomial())
      })
    } else {
      suppressWarnings({
        fit.beta <- stats::glm(Ystar ~ 0 + X, weights = (1 - p) * w,
                               offset = rep(beta0_fixed, length(p)),
                               family = stats::binomial())
      })
    }
    suppressWarnings({
      fit.theta <- stats::glm(p ~ Z, family = stats::binomial(),
                              weights = weights)
    })
    pred1 <- stats::predict(fit.theta, type = 'response')
    pred2 <- stats::predict(fit.beta, type = 'response')
    p     <- calculate.p(pred1, pred2)

    loglik <- sum(w * Ystar * log(pred1 + (pred2 * (1 - pred1))) +
                    w * (1 - Ystar) * log((1 - pred1) * (1 - pred2)))
    loglik.seq <- c(loglik.seq, loglik)

    it <- it + 1
    if (abs(loglik.seq[it] - loglik.seq[it - 1]) < tolerance)
      converged <- TRUE

    par <- c(stats::coef(fit.theta), beta0_fixed, stats::coef(fit.beta))
    param.seq <- rbind(param.seq, par)
  }

  param <- c(stats::coef(fit.theta), beta0_fixed, stats::coef(fit.beta))
  theta <- param[1:(ncol(Z) + 1)]
  beta  <- param[-(1:(ncol(Z) + 1))]

  if (is.null(weights)) {
    var <- SAMBA_obsloglik_var(Ystar, Z, X, theta, beta, beta0_fixed, expected)
  } else {
    var <- SAMBA_obsloglik_var_weighted(Ystar, Z, X, theta, beta, beta0_fixed,
                                  weights, expected)
  }

  results <- list(param = param, variance = var, param.seq = param.seq,
                  loglik.seq = loglik.seq, Ystar = Ystar, X = X, Z = Z,
                  weights = weights, beta0_fixed = beta0_fixed)

  return(results)
}
