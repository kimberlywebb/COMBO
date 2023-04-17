#' Select a Logisitic Regression Model for a Given Prior
#'
#' @param prior A character string specifying the prior distribution for the naive
#'   \eqn{\beta} parameters. Options are \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"} (double Exponential, or Weibull).
#'
#' @return \code{naive_model_picker} returns a character string specifying the
#'   logistic regression model to be turned into a .BUG file and used
#'   for MCMC estimation with \code{rjags}.
#'
naive_model_picker <- function(prior) {

  unif_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    obs_Y[i] ~ dcat(pi_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }
  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ dunif(unif_l_beta[1, l], unif_u_beta[1, l])
    beta[2, l] <- 0
  }

  }"

  t_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    obs_Y[i] ~ dcat(pi_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }
  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ dt(t_mu_beta[1, l], t_tau_beta[1, l], t_df_beta[1, l])
    beta[2, l] <- 0
  }

  }"

  normal_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    obs_Y[i] ~ dcat(pi_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }
  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ dnorm(normal_mu_beta[1, l], normal_sigma_beta[1, l])
    beta[2, l] <- 0
  }

  }"

  dexp_modelstring = "model{

# likelihood
  for(i in 1:sample_size){

    obs_Y[i] ~ dcat(pi_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi_obs[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }
  }

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ ddexp(dexp_mu_beta[1, l], dexp_b_beta[1, l])
    beta[2, l] <- 0
  }

  }"

  selected_model = ifelse(prior == "t", t_modelstring,
                          ifelse(prior == "uniform", unif_modelstring,
                                 ifelse(prior == "normal", normal_modelstring,
                                        ifelse(prior == "dexp", dexp_modelstring,
                                               stop("Please select a prior distribution.")))))

  return(selected_model)
}
