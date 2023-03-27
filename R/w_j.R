#' Compute E-step for Binary Outcome Misclassification Model Estimated With the EM-Algorithm
#'
#' @param ystar_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param pistar_matrix A numeric matrix of conditional probabilities obtained from
#'   the internal function \code{pistar_compute}. Rows of the matrix correspond to
#'   each subject and to each observed outcome category. Columns of the matrix
#'   correspond to each true, latent outcome category.
#' @param pi_matrix A numeric matrix of probabilities obtained from the internal
#'   function \code{pi_compute}. Rows of the matrix correspond to each subject.
#'   Columns of the matrix correspond to each true, latent outcome category.
#' @param sample_size An integer value specifying the number of observations in
#'   the sample. This value should be equal to the number of rows of the observed
#'   outcome matrix, \code{ystar_matrix}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcome, \code{Y*}, can take.
#'
#' @return \code{w_j} returns a matrix of E-step weights for the EM-algorithm,
#'   computed as follows:
#'   \eqn{\sum_{k = 1}^2 \frac{y^*_{ik} \pi^*_{ikj} \pi_{ij}}{\sum_{\ell = 1}^2 \pi^*_{i k \ell} \pi_{i \ell}}}.
#'   Rows of the matrix correspond to each subject. Columns of the matrix correspond
#'   to the true outcome categories \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom
#'
w_j <- function(ystar_matrix, pistar_matrix, pi_matrix, sample_size, n_cat){

  pi_ij_1_allj_repped = do.call(rbind,
                                list(pi_matrix, pi_matrix))
  pistar_pi_1 = pistar_matrix * pi_ij_1_allj_repped
  suml_pistar_pi_1 = rowSums(pistar_pi_1)

  suml_pistar_pi_denominator_1 <- matrix(rep(suml_pistar_pi_1, n_cat),
                                         nrow = n_cat * sample_size,
                                         byrow = FALSE)
  obs_Y_matrix_repped_1 <- matrix(rep(c(ystar_matrix), each = n_cat),
                                  nrow = n_cat * sample_size, byrow = TRUE)
  weight_not_summed_1 <- obs_Y_matrix_repped_1 * (pistar_pi_1 / ifelse(suml_pistar_pi_denominator_1 == 0, .00000001, suml_pistar_pi_denominator_1))

  weight_1 <- matrix(NA, nrow = sample_size, ncol = n_cat)
  for(i in 1:sample_size){
    for(j in 1:n_cat){
      k_set = c(i, sample_size + i)
      sum_terms = weight_not_summed_1[c(k_set), j]
      weight_1[i, j] = sum(sum_terms)
    }
  }

  return(weight_1)
}
