
COMBO_data <- function(sample_size,
                       x_mu, x_sigma,
                       z_shape, v_shape,
                       beta, gamma, delta){

  n_cat <- 2
  x <- rnorm(sample_size, x_mu, x_sigma)
  x_matrix <- matrix(c(rep(1, sample_size),
                       x),
                     nrow = sample_size, byrow = FALSE)

  z <- rgamma(sample_size, z_shape)
  z_matrix <- matrix(c(rep(1, sample_size),
                       z),
                     nrow = sample_size, byrow = FALSE)
  
  v <- rgamma(sample_size, v_shape)
  v_matrix <- matrix(c(rep(1, sample_size),
                       v),
                     nrow = sample_size, byrow = FALSE)

  pi_matrix <- pi_compute(beta, x_matrix, sample_size, n_cat)

  true_Y <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_Y[i] = which(rmultinom(1, 1, pi_matrix[i,]) == 1)
  }

  pistar_matrix <- pistar_compute(gamma, z_matrix, sample_size, n_cat)

  obs_Ystar <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_j = true_Y[i]
    obs_Ystar[i] = which(rmultinom(1, 1,
                                   pistar_matrix[c(i,sample_size + i),
                                                 true_j]) == 1)
  }

  obs_Ystar_reps <- matrix(rep(obs_Ystar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Ystar_matrix <- 1 * (obs_Ystar_reps == category_matrix)
  
  pitilde_array <- pitilde_compute(delta, v_matrix, sample_size, n_cat)
  
  obs_Ytilde <- rep(NA, sample_size)
  for(i in 1:sample_size){
    true_j = true_Y[i]
    obs_k = obs_Ystar[i]
    obs_Ytilde[i] = which(rmultinom(1, 1,
                                   pitilde_array[c(i,sample_size + i),
                                                 obs_k, true_j]) == 1)
  }
  
  obs_Ytilde_reps <- matrix(rep(obs_Ytilde, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Ytilde_matrix <- 1 * (obs_Ytilde_reps == category_matrix)

  data_output <- list(obs_Ystar = obs_Ystar,
                      obs_Ytilde = obs_Ytilde,
                      true_Y = true_Y,
                      obs_Ystar_matrix = obs_Ystar_matrix,
                      obs_Ytilde_matrix = obs_Ytilde_matrix,
                      x = x,
                      z = z,
                      v = v,
                      x_design_matrix = x_matrix,
                      z_design_matrix = z_matrix,
                      v_design_matrix = v_matrix)

  return(data_output)

}
