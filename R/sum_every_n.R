#' Sum Every "n"th Element
#'
#' @param x A numeric vector to sum over
#' @param n A numeric value specifying the distance between the reference index and the next index to be summed
#'
#' @return \code{sum_every_n} returns a vector of sums of every \code{n}th element of the vector \code{x}.
#'
sum_every_n <- function(x, n){
  vector_groups = split(x,
                        ceiling(seq_along(x) / n))
  sum_x = Reduce(`+`, vector_groups)

  return(sum_x)
}
