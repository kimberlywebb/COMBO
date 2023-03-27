#' Select a Naive Two-Stage Regression Model for a Given Prior
#'
#' @param prior A character string specifying the prior distribution for the naive
#'   \eqn{\beta} parameters. Options are \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"} (double Exponential, or Weibull).
#'
#' @return \code{naive_model_picker_2stage} returns a character string specifying the
#'   logistic regression model to be turned into a .BUG file and used
#'   for MCMC estimation with \code{rjags}.
#'
naive_model_picker_2stage <- function(prior) {

  unif_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi2_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

    for(k in 1:n_cat){
      for(j in 1:n_cat){

      log(phitilde[i, k, j]) <- delta[k, j, 1:dim_v] %*% v[i, 1:dim_v]
      pi2_obs2[i, k, j] <- phitilde[i, k, j] / (sum(phitilde[i, 1:n_cat, j]))

      }

    pi2_obs[i, k] <- sum(pi2_obs2[i, k, 1:n_cat] * pi_obs[i, 1:n_cat])
    }

  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ dunif(unif_l_beta[1, l], unif_u_beta[1, l])
    beta[2, l] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_v){

      delta[1, m, n] ~ dunif(unif_l_delta[1, m, n], unif_u_delta[1, m, n])
      delta[2, m, n] <- 0
    }

  }

  }"

  t_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi2_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

    for(k in 1:n_cat){
      for(j in 1:n_cat){

      log(phitilde[i, k, j]) <- delta[k, j, 1:dim_v] %*% v[i, 1:dim_v]
      pi2_obs2[i, k, j] <- phitilde[i, k, j] / (sum(phitilde[i, 1:n_cat, j]))

      }

    pi2_obs[i, k] <- sum(pi2_obs2[i, k, 1:n_cat] * pi_obs[i, 1:n_cat])
    }

  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ dt(t_mu_beta[1, l], t_tau_beta[1, l], t_df_beta[1, l])
    beta[2, l] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_v){

      delta[1, m, n] ~ dt(t_mu_delta[1, m, n], t_tau_delta[1, m, n], t_df_delta[1, m, n])
      delta[2, m, n] <- 0
    }

  }

  }"

  normal_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi2_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

    for(k in 1:n_cat){
      for(j in 1:n_cat){

      log(phitilde[i, k, j]) <- delta[k, j, 1:dim_v] %*% v[i, 1:dim_v]
      pi2_obs2[i, k, j] <- phitilde[i, k, j] / (sum(phitilde[i, 1:n_cat, j]))

      }

    pi2_obs[i, k] <- sum(pi2_obs2[i, k, 1:n_cat] * pi_obs[i, 1:n_cat])
    }

  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ dnorm(normal_mu_beta[1, l], normal_sigma_beta[1, l])
    beta[2, l] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_v){

      delta[1, m, n] ~ dunif(normal_mu_delta[1, m, n], normal_sigma_delta[1, m, n])
      delta[2, m, n] <- 0
    }

  }

  }"

  dexp_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi2_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

    for(k in 1:n_cat){
      for(j in 1:n_cat){

      log(phitilde[i, k, j]) <- delta[k, j, 1:dim_v] %*% v[i, 1:dim_v]
      pi2_obs2[i, k, j] <- phitilde[i, k, j] / (sum(phitilde[i, 1:n_cat, j]))

      }

    pi2_obs[i, k] <- sum(pi2_obs2[i, k, 1:n_cat] * pi_obs[i, 1:n_cat])
    }

  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ ddexp(dexp_mu_beta[1, l], dexp_b_beta[1, l])
    beta[2, l] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_v){

      delta[1, m, n] ~ dunif(dexp_mu_delta[1, m, n], dexp_b_delta[1, m, n])
      delta[2, m, n] <- 0
    }

  }

  }"

  selected_model = ifelse(prior == "t", t_modelstring,
                          ifelse(prior == "uniform", unif_modelstring,
                                 ifelse(prior == "normal", normal_modelstring,
                                        ifelse(prior == "dexp", dexp_modelstring,
                                               print("Please select a prior distribution.")))))

  return(selected_model)
}
