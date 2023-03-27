#' Compute the Mean Conditional Probability of Correct Classification, by True Outcome Across all Subjects
#'
#' @param pistar_matrix A numeric matrix of conditional probabilities obtained from
#'   the internal function \code{pistar_compute_for_chains}. Rows of the matrix correspond to
#'   each subject and to each observed outcome category. Columns of the matrix
#'   correspond to each true, latent outcome category.
#' @param j An integer value representing the true outcome category to compute
#'   the average conditional probability of correct classification for.
#'   \code{j} can take on values \code{1} and \code{2}.
#' @param sample_size An integer value specifying the number of observations in the sample.
#'
#' @return \code{mean_pistarjj_compute} returns a numeric value equal to the average
#'   conditional probability \eqn{P(Y^* = j | Y = j, Z)} across all subjects.
#'
#' @importFrom stats rnorm
#'
mean_pistarjj_compute <- function(pistar_matrix, j, sample_size){
  k_index = ifelse(j == 1, 1:sample_size,
                   ifelse(j == 2, (sample_size + 1):(sample_size*2), NA))
  mean_pistar = mean(pistar_matrix[k_index, j])
  return(mean_pistar)
}
