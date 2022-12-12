
loglik <- function(param_current,
                   obs_Ystar_matrix, obs_Ytilde_matrix,
                   X, Z, V,
                   sample_size, n_cat){

  beta_current = matrix(param_current[1:ncol(X)], ncol = 1)
  gamma_current = matrix(c(param_current[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))]),
                         ncol = n_cat, byrow = FALSE)
  delta_current = array(c(param_current[(ncol(X) + (n_cat * ncol(Z)) + 1):length(param_current)]),
                         dim = c(ncol(Z), 2, 2))

  pi_terms_v = pi_compute(beta_current, X, sample_size, n_cat)
  pistar_terms_v = pistar_compute(gamma_current, Z, sample_size, n_cat)
  pitilde_terms_v = pitilde_compute(delta_current, V, sample_size, n_cat)

  weights = w_j(obs_Ystar_matrix, pistar_terms_v, pi_terms_v, sample_size, n_cat)

  loglikelihood = sum(
    (q_beta_f(beta_current, X = X, w_mat = weights,
              sample_size = sample_size, n_cat = n_cat)) +
      (q_gamma_f(c(gamma_current), Z = Z,
                 obs_Y_matrix = obs_Y_matrix,
                 w_mat = weights,
                 sample_size = sample_size, n_cat = n_cat)) + 
      (q_delta_f(c(delta_current), V = V,
                 obs_Ystar_matrix, obs_Ytilde_matrix,
                 w_mat = weights,
                 sample_size = sample_size, n_cat = n_cat)))

  return(loglikelihood)
}

