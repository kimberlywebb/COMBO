
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


w_j <- function(ystar_matrix, ytilde_matrix, pitilde_array, pistar_matrix, pi_matrix, sample_size, n_cat){
  
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
