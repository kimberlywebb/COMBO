
COMBO_EM <- function(Ystar, Ytilde,
                     x_matrix, z_matrix, v_matrix,
                     beta_start, gamma_start, delta_start,
                     tolerance = 1e-7, max_em_iterations = 1500,
                     em_method = "squarem"){

  if (is.data.frame(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.numeric(z_matrix))
    stop("'z_matrix' should be a numeric matrix.")

  if (is.vector(z_matrix))
    z_matrix <- as.matrix(z_matrix)
  if (!is.matrix(z_matrix))
    stop("'z_matrix' should be a matrix or data.frame.")

  if (!is.null(x_matrix)) {
    if (is.data.frame(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.numeric(x_matrix))
      stop("'x_matrix' must be numeric.")
    if (is.vector(x_matrix))
      x_matrix <- as.matrix(x_matrix)
    if (!is.matrix(x_matrix))
      stop("'x_matrix' must be a data.frame or matrix.")
  }

  if (!is.numeric(Ystar) || !is.vector(Ystar))
    stop("'Ystar' must be a numeric vector.")
  if (length(setdiff(1:2, unique(Ystar))) != 0)
    stop("'Ystar' must be coded 1/2, where the reference category is 2.")

  n_cat = 2
  sample_size = length(Ystar)

  if (nrow(z_matrix) != sample_size)
    stop("The number of rows of 'z_matrix' must match the length of 'Ystar'.")
  if (!is.null(x_matrix) && nrow(x_matrix) != sample_size)
    stop("The number of rows of 'x_matrix' must match the length of 'Ystar'.")

  X = matrix(c(rep(1, sample_size), c(x_matrix)),
             byrow = FALSE, nrow = sample_size)
  Z = matrix(c(rep(1, sample_size), c(z_matrix)),
             byrow = FALSE, nrow = sample_size)
  V = matrix(c(rep(1, sample_size), c(v_matrix)),
             byrow = FALSE, nrow = sample_size)

  obs_Ystar_reps = matrix(rep(Ystar, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix = matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                           byrow = FALSE)
  obs_Ystar_matrix = 1 * (obs_Y_reps == category_matrix)
  
  obs_Ytilde_reps <- matrix(rep(Ytilde, n_cat), nrow = sample_size, byrow = FALSE)
  category_matrix <- matrix(rep(1:n_cat, each = sample_size), nrow = sample_size,
                            byrow = FALSE)
  obs_Ytilde_matrix <- 1 * (obs_Ytilde_reps == category_matrix)

  control_settings = list(convtype = "parameter", tol = tolerance,
                          stoptype = "maxiter", maxiter = max_em_iterations)

  results = turboEM::turboem(par = c(c(beta_start), c(gamma_start), c(delta_start)),
                             fixptfn = em_function, #objfn = loglik,
                             method = c(em_method),
                             obs_Ystar_matrix = obs_Ystar_matrix,
                             obs_Ytilde_matrix = obs_Ytilde_matrix,
                             X = X, Z = Z, V = V,
                             sample_size = sample_size, n_cat = n_cat,
                             control.run = control_settings)
  
  ##############################################################################
  # NOTHING AFTER THIS IS EDITED FOR THE TWO-STAGE MODEL
  ##############################################################################

  Ystar01 = ifelse(Ystar == 1, 1, ifelse(Ystar == 2, 0, NA))
  log_reg = stats::glm(Ystar01 ~ . + 0, as.data.frame(X),
                       family = "binomial"(link = "logit"))

  SAMBA_start <- c(beta_start, c(gamma_start)[1:(1 + ncol(z_matrix))])
  SAMBA_i <- SAMBA::obsloglikEM(Ystar01, Z = x_matrix,
                                X = z_matrix, start = SAMBA_start,
                                tol = tolerance,
                                maxit = max_em_iterations)

  perfect_sens_start <- c(beta_start, c(gamma_start)[(2 + ncol(z_matrix)):length(c(gamma_start))])
  perfect_sens_i <- perfect_sensitivity_EM(Ystar01, Z = x_matrix,
                                           X = z_matrix, start = perfect_sens_start,
                                           tolerance = tolerance,
                                           max_em_iterations = max_em_iterations)

  # Do label switching correction within the EM algorithm simulation
  results_i_gamma <- matrix(turboEM::pars(results)[(ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))],
                            ncol = n_cat, byrow = FALSE)
  results_i_pistar_v <- pistar_compute(results_i_gamma, Z, sample_size, n_cat)
  pistar_11 <- mean(results_i_pistar_v[1:sample_size, 1])
  pistar_22 <- mean(results_i_pistar_v[(sample_size + 1):(2*sample_size), 2])
  estimates_i <- if ((pistar_11 > .50 | pistar_22 > .50) |
                     (is.na(pistar_11) & is.na(pistar_22))) {
    # If turboem cannot estimate the parameters they will be NA.
    turboEM::pars(results)
  } else {
    gamma_index = (ncol(X) + 1):(ncol(X) + (n_cat * ncol(Z)))
    n_gamma_param = length(gamma_index) / n_cat
    gamma_flip_index = ncol(X) + c((n_gamma_param + 1):length(gamma_index), 1:n_gamma_param)
    c(-1*turboEM::pars(results)[1:ncol(X)], turboEM::pars(results)[gamma_flip_index])
  }

  #sigma_EM = tryCatch(solve(turboEM::hessian(results)[[1]]), silent = TRUE,
  #                    error = function(e) NA)
  #SE_EM = tryCatch(sqrt(diag(matrix(Matrix::nearPD(sigma_EM)$mat,
  #                                  row = length(c(c(beta_start), c(gamma_start))),
  #                                  byrow = FALSE))),
  #                 silent = TRUE,
  #                 error = function(e) rep(NA, ncol(X) + (n_cat * ncol(Z))))

  sigma_EM = solve(turboEM::hessian(results)[[1]])
  SE_EM = sqrt(diag(matrix(Matrix::nearPD(sigma_EM)$mat,
                           nrow = length(c(c(beta_start), c(gamma_start))),
                           byrow = FALSE)))

  beta_param_names <- paste0(rep("beta", ncol(X)), 1:ncol(X))
  gamma_param_names <- paste0(rep("gamma", (n_cat * ncol(Z))),
                              rep(1:ncol(Z), n_cat),
                              rep(1:n_cat, each = ncol(Z)))
  SAMBA_beta_param_names <- paste0(rep("SAMBA_beta", ncol(X)), 1:ncol(X))
  SAMBA_gamma_param_names <- paste0(rep("SAMBA_gamma", (ncol(Z))),
                              rep(1:ncol(Z), 1),
                              rep(1, ncol(Z)))
  PSens_beta_param_names <- paste0(rep("PSens_beta", ncol(X)), 1:ncol(X))
  PSens_gamma_param_names <- paste0(rep("PSens_gamma", (ncol(Z))),
                                    rep(1:ncol(Z), 1),
                                    rep(2, ncol(Z)))
  naive_beta_param_names <- paste0("naive_", beta_param_names)
  estimates <- data.frame(Parameter = c(beta_param_names,
                                        gamma_param_names,
                                        SAMBA_beta_param_names,
                                        SAMBA_gamma_param_names,
                                        PSens_beta_param_names,
                                        PSens_gamma_param_names,
                                        naive_beta_param_names),
                          Estimates = c(estimates_i,
                                        unname(SAMBA_i$param),
                                        unname(perfect_sens_i$param),
                                        unname(log_reg$coefficients)),
                          SE = c(SE_EM, sqrt(diag(SAMBA_i$variance)),
                                 sqrt(diag(perfect_sens_i$variance)),
                                 unname(summary(log_reg)$coefficients[,2])),
                          Convergence = c(rep(results$convergence,
                                              length(c(beta_param_names,
                                                       gamma_param_names))),
                                          rep(NA, length(c(SAMBA_beta_param_names,
                                                           SAMBA_gamma_param_names))),
                                          rep(NA, length(c(PSens_beta_param_names,
                                                           PSens_gamma_param_names))),
                                          rep(log_reg$converged, ncol(X))))

  return(estimates)
}
