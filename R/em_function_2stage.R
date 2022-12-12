
em_function <- function(param_current,
                        obs_Ystar_matrix, obs_Ytilde_matrix,
                        X, Z, V,
                        sample_size, n_cat){

  beta_current = matrix(param_current[1:ncol(X)], ncol = 1)
  gamma_current = matrix(c(param_current[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))]),
                         ncol = n_cat, byrow = FALSE)
  delta_current = array(c(param_current[(ncol(X) + (n_cat * ncol(Z)) + 1):length(param_current)]),
                        dim = c(ncol(Z), 2, 2))

  probabilities = pi_compute(beta_current, X, sample_size, n_cat)
  conditional_probabilities = pistar_compute(gamma_current, Z, sample_size, n_cat)
  conditional_probabilities2 = pitilde_compute(delta_current, V, sample_size, n_cat)

  weights = w_j(ystar_matrix = obs_Ystar_matrix,
                ytilde_matrix = obs_Ytilde_matrix,
                pitilde_array = conditional_probabilities2,
                pistar_matrix = conditional_probabilities,
                pi_matrix = probabilities,
                sample_size = sample_size, n_cat = n_cat)

  Ystar01 = obs_Ystar_matrix[,1]
  fit.gamma1 <- suppressWarnings( stats::glm(Ystar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,1],
                           family = "binomial"(link = "logit")) )
  gamma1_new <- unname(coefficients(fit.gamma1))

  fit.gamma2 <- suppressWarnings( stats::glm(Ystar01 ~ . + 0, as.data.frame(Z),
                           weights = weights[,2],
                           family = "binomial"(link = "logit")) )
  gamma2_new <- unname(coefficients(fit.gamma2))

  fit.beta <- suppressWarnings( stats::glm(weights[,1] ~ . + 0, as.data.frame(X),
                         family = stats::binomial()) )
  beta_new <- unname(coefficients(fit.beta))
  
  Ytilde01 <- obs_Ytilde_matrix[,1]
  
  Y_star1_tilde <- ifelse(Ystar01 == 1 & Ytilde01 == 1, 1, 0)
  fit.delta11 <- suppressWarnings( stats::glm(Y_star1_tilde ~ . + 0, as.data.frame(V),
                                             weights = weights[,1],
                                             family = "binomial"(link = "logit")) )
  delta11_new <- unname(coefficients(fit.delta11))
  
  Y_star2_tilde <- ifelse(Ystar01 == 0 & Ytilde01 == 1, 1, 0)
  fit.delta12 <- suppressWarnings( stats::glm(Y_star2_tilde ~ . + 0, as.data.frame(V),
                                              weights = weights[,1],
                                              family = "binomial"(link = "logit")) )
  delta12_new <- unname(coefficients(fit.delta12))
  
  fit.delta21 <- suppressWarnings( stats::glm(Y_star1_tilde ~ . + 0, as.data.frame(V),
                                              weights = weights[,2],
                                              family = "binomial"(link = "logit")) )
  delta21_new <- unname(coefficients(fit.delta21))
  
  fit.delta22 <- suppressWarnings( stats::glm(Y_star2_tilde ~ . + 0, as.data.frame(V),
                                              weights = weights[,2],
                                              family = "binomial"(link = "logit")) )
  delta22_new <- unname(coefficients(fit.delta22))
  
  delta_estimates <- optim(par = rep(1, 8),
                           fn = q_delta_f,
                           V = V,
                           obs_Ystar_matrix = obs_Ystar_matrix,
                           obs_Ytilde_matrix = obs_Ytilde_matrix,
                           w_mat = weights,
                           sample_size = sample_size, n = n,
                           method = "BFGS")
  
  delta_new <- delta_estimates$par

  param_new = c(beta_new, gamma1_new, gamma2_new, delta_new)
  param_current = param_new
  return(param_new)

}
