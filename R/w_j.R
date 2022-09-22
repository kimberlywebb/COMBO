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
#'   \eqn{\sum_{k = 1}^2 \frac{y^*_{ik} \pi^*_{ikj} \pi_{ij}}{\sum_{\ell = 1}^2 \pi^*_{i \ell j} \pi_{ij}}}.
#'   Rows of the matrix correspond to each subject. Columns of the matrix correspond
#'   to the true outcome categories \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom
#'
#' @examples \dontrun{
#' set.seed(123)
#' n <- 1000
#' n_cat <- 2
#'
#' ones <- rep(1, n)
#' x <- rnorm(n)
#' X <- matrix(c(ones, x), nrow = n, byrow = FALSE)
#' z <- rgamma(n, shape = 1)
#' Z <- matrix(c(ones, z), nrow = n, byrow = FALSE)
#'
#' beta <- matrix(c(1, -2), ncol = 1)
#' gamma <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)
#'
#' probabilities <- pi_compute(beta, X, n, n_cat = 2)
#' conditional_probabilities <- pistar_compute(gamma, Z, n, n_cat = 2)
#'
#' true_Y <- rep(NA, n)
#' for(i in 1:n){
#'   true_Y[i] = which(rmultinom(1, 1, probabilities[i,]) == 1)
#' }
#' pistar_matrix_labels <- rep(1:n_cat, each = n)
#' obs_Y <- rep(NA, n)
#' for(i in 1:n){
#'  true_j = true_Y[i]
#'  obs_Y[i] = which(rmultinom(1, 1,
#'                             conditional_probabilities[c(i, n + i), true_j]) == 1)
#'  }
#'
#' obs_Y_reps <- matrix(rep(obs_Y, n_cat), nrow = n, byrow = FALSE)
#' category_matrix <- matrix(rep(1:n_cat, each = n), nrow = n, byrow = FALSE)
#' obs_Y_matrix <- 1 * (obs_Y_reps == category_matrix)
#'
#' e_step_weights <- w_j(ystar_matrix = obs_Y_matrix,
#'                       pistar_matrix = conditional_probabilities,
#'                       pi_matrix = probabilities,
#'                       sample_size = n, n_cat = n_cat)
#' head(e_step_weights)

#' }
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
  weight_not_summed_1 <- obs_Y_matrix_repped_1 * (pistar_pi_1 / suml_pistar_pi_denominator_1)

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
