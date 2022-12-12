
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
