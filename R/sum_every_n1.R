#' Sum Every "n"th Element, then add 1
#'
#' @param x A numeric vector to sum over
#' @param n A numeric value specifying the distance between the reference index and the next index to be summed
#'
#' @return \code{sum_every_n1} returns a vector of sums of every \code{n}th element of the vector \code{x}, plus 1.
#'
sum_every_n1 <- function(x, n){
  vector_groups = split(x,
                        ceiling(seq_along(x) / n))
  sum_x = Reduce(`+`, vector_groups) + 1

  return(sum_x)
}
