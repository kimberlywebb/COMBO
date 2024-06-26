#' Select a Two-Stage Binary Outcome Misclassification Model for a Given Prior
#'
#' @param prior A character string specifying the prior distribution for the
#'   \eqn{\beta}, \eqn{\gamma}, and \eqn{\delta} parameters. Options are \code{"t"},
#'   \code{"uniform"}, \code{"normal"}, or \code{"dexp"} (double Exponential, or Weibull).
#'
#' @return \code{model_picker} returns a character string specifying the two-stage binary
#'   outcome misclassification model to be turned into a .BUG file and used
#'   for MCMC estimation with \code{rjags}.
#'
model_picker_2stage <- function(prior){

  unif_modelstring = "
  model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_star_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi_tilde_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

     for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phistar[i, k, j]) <- gamma1[k, j, 1:dim_z] %*% z[i,1:dim_z]
        pistar[i, k, j] <- phistar[i, k, j] / (sum(phistar[i, 1:n_cat, j]))

      }

      pi_star_obs[i, k] <- sum(pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

    }

  for(l in 1:n_cat){
     for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phitilde[i, l, k, j]) <- gamma2[l, k, j, 1:dim_v] %*% v[i, 1:dim_v]
        pitilde[i, l, k, j] <- phitilde[i, l, k, j] / (sum(phitilde[i, 1:n_cat, k, j]))

      }

      pi_tilde_obs1[i, l, k] <- sum(pitilde[i, l, k, 1:n_cat] * pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

     }

      pi_tilde_obs[i, l] <- sum(pi_tilde_obs1[i, l, 1:n_cat])
    }

  }

# reference categories

  #beta[n_cat, 1:dim_x] <- 0
  #gamma[n_cat, 1:n_cat, 1:dim_z] <- 0

# priors
  for(p in 1:dim_x){
    beta[1, p] ~ dunif(unif_l_beta[1, p], unif_u_beta[1, p])
    beta[2, p] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_z){

      gamma1[1, m, n] ~ dunif(unif_l_gamma[1, m, n], unif_u_gamma[1, m, n])
      gamma1[2, m, n] <- 0
    }

  }

  for(q in 1:n_cat){
  for(r in 1:n_cat){
    for(s in 1:dim_v){

      gamma2[1, q, r, s] ~ dunif(unif_l_delta[1, q, r, s], unif_u_delta[1, q, r, s])
      gamma2[2, q, r, s] <- 0
    }

  }

  }
}
"

t_modelstring = "
  model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_star_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi_tilde_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

    for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phistar[i, k, j]) <- gamma1[k, j, 1:dim_z] %*% z[i,1:dim_z]
        pistar[i, k, j] <- phistar[i, k, j] / (sum(phistar[i, 1:n_cat, j]))

      }

      pi_star_obs[i, k] <- sum(pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

    }

    for(l in 1:n_cat){
     for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phitilde[i, l, k, j]) <- gamma2[l, k, j, 1:dim_v] %*% v[i, 1:dim_v]
        pitilde[i, l, k, j] <- phitilde[i, l, k, j] / (sum(phitilde[i, 1:n_cat, k, j]))

      }

      pi_tilde_obs1[i, l, k] <- sum(pitilde[i, l, k, 1:n_cat] * pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

     }

      pi_tilde_obs[i, l] <- sum(pi_tilde_obs1[i, l, 1:n_cat])
    }

  }

# reference categories

  #beta[n_cat, 1:dim_x] <- 0
  #gamma[n_cat, 1:n_cat, 1:dim_z] <- 0

# priors
  for(p in 1:dim_x){
    beta[1, p] ~ dt(t_mu_beta[1,p], t_tau_beta[1,p], t_df_beta[1,p])
    beta[2, p] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_z){
      gamma1[1, m, n] ~ dt(t_mu_gamma[1,m,n], t_tau_gamma[1,m,n], t_df_gamma[1,m,n])
      gamma1[2, m, n] <- 0
    }

  }

  for(q in 1:n_cat){
    for(r in 1:n_cat){
      for(s in 1:dim_v){

      gamma2[1, q, r, s] ~ dt(t_mu_delta[1, q, r, s], t_tau_delta[1, q, r, s], t_df_delta[1, q, r, s])
      gamma2[2, q, r, s] <- 0
    }

  }

  }

}
"

normal_modelstring = "
  model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_star_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi_tilde_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

    for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phistar[i, k, j]) <- gamma1[k, j, 1:dim_z] %*% z[i,1:dim_z]
        pistar[i, k, j] <- phistar[i, k, j] / (sum(phistar[i, 1:n_cat, j]))

      }

      pi_star_obs[i, k] <- sum(pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

    }

    for(l in 1:n_cat){
     for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phitilde[i, l, k, j]) <- gamma2[l, k, j, 1:dim_v] %*% v[i, 1:dim_v]
        pitilde[i, l, k, j] <- phitilde[i, l, k, j] / (sum(phitilde[i, 1:n_cat, k, j]))

      }

      pi_tilde_obs1[i, l, k] <- sum(pitilde[i, l, k, 1:n_cat] * pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

     }

      pi_tilde_obs[i, l] <- sum(pi_tilde_obs1[i, l, 1:n_cat])
    }

  }

# reference categories

  #beta[n_cat, 1:dim_x] <- 0
  #gamma[n_cat, 1:n_cat, 1:dim_z] <- 0

# priors
  for(p in 1:dim_x){
    beta[1, p] ~ dnorm(normal_mu_beta[1, p], normal_sigma_beta[1, p])
    beta[2, p] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_z){

      gamma1[1, m, n] ~ dnorm(normal_mu_gamma[1, m, n], normal_sigma_gamma[1, m, n])
      gamma1[2, m, n] <- 0
    }

  }

for(q in 1:n_cat){
  for(r in 1:n_cat){
    for(s in 1:dim_v){

      gamma2[1, q, r, s] ~ dnorm(normal_mu_delta[1, q, r, s], normal_sigma_delta[1, q, r, s])
      gamma2[2, q, r, s] <- 0
    }

  }

  }

}
"

dexp_modelstring = "
  model{

# likelihood
  for(i in 1:sample_size){

    Y_star[i] ~ dcat(pi_star_obs[i, 1:n_cat])
    Y_tilde[i] ~ dcat(pi_tilde_obs[i, 1:n_cat])

  # regression
    for(j in 1:n_cat){

      log(phi[i, j]) <- beta[j,1:dim_x] %*% x[i,1:dim_x]
      pi[i, j] <- phi[i, j] / sum(phi[i, 1:n_cat])

    }

    for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phistar[i, k, j]) <- gamma1[k, j, 1:dim_z] %*% z[i,1:dim_z]
        pistar[i, k, j] <- phistar[i, k, j] / (sum(phistar[i, 1:n_cat, j]))

      }

      pi_star_obs[i, k] <- sum(pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

    }

    for(l in 1:n_cat){
     for(k in 1:n_cat){
      for(j in 1:n_cat){

        log(phitilde[i, l, k, j]) <- gamma2[l, k, j, 1:dim_v] %*% v[i, 1:dim_v]
        pitilde[i, l, k, j] <- phitilde[i, l, k, j] / (sum(phitilde[i, 1:n_cat, k, j]))

      }

      pi_tilde_obs1[i, l, k] <- sum(pitilde[i, l, k, 1:n_cat] * pistar[i, k, 1:n_cat] * pi[i, 1:n_cat])

     }

      pi_tilde_obs[i, l] <- sum(pi_tilde_obs1[i, l, 1:n_cat])
    }

  }

# reference categories

  #beta[n_cat, 1:dim_x] <- 0
  #gamma[n_cat, 1:n_cat, 1:dim_z] <- 0

# priors
  for(l in 1:dim_x){
    beta[1, l] ~ ddexp(dexp_mu_beta[1, l], dexp_b_beta[1, l])
    beta[2, l] <- 0
  }

  for(m in 1:n_cat){
    for(n in 1:dim_z){

      gamma1[1, m, n] ~ ddexp(dexp_mu_gamma[1, m, n], dexp_b_gamma[1, m, n])
      gamma1[2, m, n] <- 0
    }

  }

  for(q in 1:n_cat){
  for(r in 1:n_cat){
    for(s in 1:dim_v){

      gamma2[1, q, r, s] ~ ddexp(dexp_mu_delta[1, q, r, s], dexp_b_delta[1, q, r, s])
      gamma2[2, q, r, s] <- 0
    }

  }

  }

}
"

selected_model = ifelse(prior == "t", t_modelstring,
                        ifelse(prior == "uniform", unif_modelstring,
                               ifelse(prior == "normal", normal_modelstring,
                                      ifelse(prior == "dexp", dexp_modelstring,
                                             stop("Please select a prior distribution.")))))

return(selected_model)

}
