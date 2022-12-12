q_delta_f <- function(delta_v, V, obs_Ystar_matrix, obs_Ytilde_matrix, w_mat,
                      sample_size, n_cat){
  
  delta_mat = array(delta_v, dim = c(ncol(V), 2, 2))
  
  pitilde_terms_array = pitilde_compute(delta_mat, V, sample_size, n_cat)
  
  big_Ystar_matrix = array(c(rep(obs_Ystar_matrix[,1], 2),
                             rep(obs_Ystar_matrix[,2], 2),
                             rep(obs_Ystar_matrix[,1], 2),
                             rep(obs_Ystar_matrix[,2], 2)),
                           dim = c(2 * sample_size, 2, 2))
  
  big_Ytilde_matrix = array(rep(c(obs_Ytilde_matrix[,1],
                                  obs_Ytilde_matrix[,2]), 4),
                            dim = c(2 * sample_size, 2, 2))
  
  big_w_matrix = array(c(rep(w_mat[,1], 4),
                         rep(w_mat[,2], 4)),
                       dim = c(2 * sample_size, 2, 2))
  
  summand = big_Ystar_matrix * big_Ytilde_matrix * big_w_matrix * log(pitilde_terms_array)
  result = -sum(summand)
  return(result)
}
