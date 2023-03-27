#' Compute E-step for Two-Stage Binary Outcome Misclassification Model Estimated With the EM-Algorithm
#'
#' @param ystar_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \code{Y*}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param ytilde_matrix A numeric matrix of indicator variables (0, 1) for the observed
#'   outcome \eqn{\tilde{Y}}. Rows of the matrix correspond to each subject. Columns of
#'   the matrix correspond to each observed outcome category. Each row should contain
#'   exactly one 0 entry and exactly one 1 entry.
#' @param pitilde_array A numeric array of conditional probabilities obtained from
#'   the internal function \code{pitilde_compute}. Rows of the matrices correspond
#'   to each subject and to each second-stage observed outcome category. Columns of the matrix correspond
#'   to each first-stage observed outcome category. The third dimension of the array
#'   corresponds to each true, latent outcome category.
#' @param pistar_matrix A numeric matrix of conditional probabilities obtained from
#'   the internal function \code{pistar_compute}. Rows of the matrix correspond to
#'   each subject and to each first-stage observed outcome category. Columns of the matrix
#'   correspond to each true, latent outcome category.
#' @param pi_matrix A numeric matrix of probabilities obtained from the internal
#'   function \code{pi_compute}. Rows of the matrix correspond to each subject.
#'   Columns of the matrix correspond to each true, latent outcome category.
#' @param sample_size An integer value specifying the number of observations in
#'   the sample. This value should be equal to the number of rows of the observed
#'   outcome matrices, \code{ystar_matrix} and \code{ytilde_matrix}.
#' @param n_cat The number of categorical values that the true outcome, \code{Y},
#'   and the observed outcomes can take.
#'
#' @return \code{w_j} returns a matrix of E-step weights for the EM-algorithm,
#'   computed as follows:
#'   \eqn{\sum_{k = 1}^2 \sum_{\ell = 1}^2 \frac{y^*_{ik} \tilde{y}_{i \ell} \tilde{\pi}_{i \ell kj} \pi^*_{ikj} \pi_{ij}}{\sum_{h = 1}^2 \tilde{\pi}_{i \ell kh} \pi^*_{ikh} \pi_{ih}}}.
#'   Rows of the matrix correspond to each subject. Columns of the matrix correspond
#'   to the true outcome categories \eqn{j = 1, \dots,} \code{n_cat}.
#'
#' @include pi_compute.R
#' @include pistar_compute.R
#' @include pitilde_compute.R
#'
#' @importFrom stats rnorm rgamma rmultinom
#'
w_j_2stage <- function(ystar_matrix, ytilde_matrix,
                       pitilde_array, pistar_matrix, pi_matrix,
                       sample_size, n_cat){

  big_Ystar_matrix = array(c(rep(ystar_matrix[,1], 2),
                             rep(ystar_matrix[,2], 2),
                             rep(ystar_matrix[,1], 2),
                             rep(ystar_matrix[,2], 2)),
                           dim = c(2 * sample_size, 2, 2))

  big_Ytilde_matrix = array(rep(c(ytilde_matrix[,1],
                                  ytilde_matrix[,2]), 4),
                            dim = c(2 * sample_size, 2, 2))

  big_pi_array = array(c(rep(pi_matrix[,1], 4),
                         rep(pi_matrix[,2], 4)),
                       dim = c(2 * sample_size, 2, 2))

  big_pistar_array = array(c(rep(pistar_matrix[1:sample_size, 1], 2),
                             rep(pistar_matrix[(sample_size + 1):(2 * sample_size), 1], 2),
                             rep(pistar_matrix[1:sample_size, 2], 2),
                             rep(pistar_matrix[(sample_size + 1):(2 * sample_size), 2], 2)),
                           dim = c(2 * sample_size, 2, 2))

  pi_multiply <- pitilde_array * big_pistar_array * big_pi_array
  weight_denominator <- pi_multiply[,,1] + pi_multiply[,,2]

  pi_y_multiply <- big_Ystar_matrix * big_Ytilde_matrix * pi_multiply

  pi_y_before_sum <- pi_y_multiply / array(c(weight_denominator, weight_denominator),
                                           dim = c(dim(weight_denominator), 2))

  pi_y_suml <- suml_function(pi_y_before_sum)
  pi_y_sumk <- apply(pi_y_suml, 3, rowSums)

  weight_matrix <- pi_y_sumk

  return(weight_matrix)
}

suml_function <- function(my_array){
  my_n <- dim(my_array[,,1])[1] / 2

  sum1 <- matrix(c(my_array[1:my_n, 1, 1] + my_array[(my_n + 1):(2 * my_n), 1, 1],
                   my_array[1:my_n, 2, 1] + my_array[(my_n + 1):(2 * my_n), 2, 1]),
                 nrow = my_n)

  sum2 <- matrix(c(my_array[1:my_n, 1, 2] + my_array[(my_n + 1):(2 * my_n), 1, 2],
                   my_array[1:my_n, 2, 2] + my_array[(my_n + 1):(2 * my_n), 2, 2]),
                 nrow = my_n)

  return_array <- array(c(sum1, sum2), dim = c(my_n, 2, 2))
  return(return_array)
}
