#' Compute Conditional Probability of Each Second-Stage Observed Outcome Given Each True Outcome and First-Stage Observed Outcome, for Every Subject
#'
#' @param delta A numeric array of regression parameters for the second-stage observed
#'   outcome mechanism, \eqn{\tilde{Y} | Y^*, Y}
#'   (second-stage observed outcome, given the first-stage observed outcome and the true outcome) ~ \code{V} (misclassification
#'   predictor matrix). Rows of the matrix correspond to parameters for the \eqn{\tilde{Y} = 1}
#'   observed outcome, with the dimensions of \code{V}.
#'   Columns of the matrix correspond to the first-stage observed outcome categories
#'   \eqn{k = 1, \dots,} \code{n_cat}. The third dimension of the array
#'   corresponds to the true outcome categories \eqn{j = 1, \dots,} \code{n_cat}
#' @param V A numeric design matrix.
#' @param n An integer value specifying the number of observations in the sample.
#'   This value should be equal to the number of rows of the design matrix, \code{V}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcomes can take.
#'
#' @return \code{pitilde_compute} returns an array of conditional probabilities,
#'   \eqn{P(\tilde{Y}_i = \ell | Y^*_i = k, Y_i = j, V_i) = \frac{\text{exp}\{\delta_{\ell kj0} + \delta_{\ell kjV} V_i\}}{1 + \text{exp}\{\delta_{\ell kj0} + \delta_{\ell kjV} V_i\}}}
#'   for each of the \eqn{i = 1, \dots,} \code{n} subjects. Rows of the matrix
#'   correspond to each subject and second-stage observed outcome. Specifically, the probability
#'   for subject \eqn{i} and observed category $1$ occurs at row \eqn{i}. The probability
#'   for subject \eqn{i} and observed category $2$ occurs at row \eqn{i +} \code{n}.
#'   Columns of the matrix correspond to the first-stage outcome categories, \eqn{k = 1, \dots,} \code{n_cat}.
#'   The third dimension of the array corresponds to the true outcome categories,
#'   \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include sum_every_n.R
#' @include sum_every_n1.R
#'
#' @importFrom stats rnorm
#'
pitilde_compute <- function(delta, V, n, n_cat){

  exp_vd1 = exp(V %*% delta[,,1])
  exp_vd2 = exp(V %*% delta[,,2])

  pi_denominator1 = apply(exp_vd1, FUN = sum_every_n1, n, MARGIN = 2)
  pi_result1 = exp_vd1 / rbind(pi_denominator1)

  pi_denominator2 = apply(exp_vd2, FUN = sum_every_n1, n, MARGIN = 2)
  pi_result2 = exp_vd2 / rbind(pi_denominator2)

  pitilde_matrix1 = rbind(pi_result1,
                          1 - apply(pi_result1,
                                    FUN = sum_every_n, n = n,
                                    MARGIN = 2))

  pitilde_matrix2 = rbind(pi_result2,
                          1 - apply(pi_result2,
                                    FUN = sum_every_n, n = n,
                                    MARGIN = 2))


  pitilde_array = array(c(pitilde_matrix1, pitilde_matrix2),
                        dim = c(dim(pitilde_matrix1), 2))

  return(pitilde_array)
}
