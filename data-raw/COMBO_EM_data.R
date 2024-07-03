## code to prepare `COMBO_EM_data` dataset goes here

set.seed(123)
n <- 1000
x_mu <- 0
x_sigma <- 1
z_shape <- 1

true_beta <- matrix(c(1, -2), ncol = 1)
true_gamma <- matrix(c(.5, 1, -.5, -1), nrow = 2, byrow = FALSE)

x_matrix = matrix(rnorm(n, x_mu, x_sigma), ncol = 1)
X = matrix(c(rep(1, n), x_matrix[,1]), ncol = 2, byrow = FALSE)
z_matrix = matrix(rgamma(n, z_shape), ncol = 1)
Z = matrix(c(rep(1, n), z_matrix[,1]), ncol = 2, byrow = FALSE)

exp_xb = exp(X %*% true_beta)
pi_result = exp_xb[,1] / (exp_xb[,1] + 1)
pi_matrix = matrix(c(pi_result, 1 - pi_result), ncol = 2, byrow = FALSE)

true_Y <- rep(NA, n)
for(i in 1:n){
  true_Y[i] = which(stats::rmultinom(1, 1, pi_matrix[i,]) == 1)
}

exp_zg = exp(Z %*% true_gamma)
pistar_denominator = matrix(c(1 + exp_zg[,1], 1 + exp_zg[,2]), ncol = 2, byrow = FALSE)
pistar_result = exp_zg / pistar_denominator

pistar_matrix = matrix(c(pistar_result[,1], 1 - pistar_result[,1],
                         pistar_result[,2], 1 - pistar_result[,2]),
                       ncol = 2, byrow = FALSE)

obs_Y <- rep(NA, n)
for(i in 1:n){
  true_j = true_Y[i]
  obs_Y[i] = which(rmultinom(1, 1,
                             pistar_matrix[c(i, n + i),
                                           true_j]) == 1)
}

Ystar <- obs_Y

COMBO_EM_data = list(
  Y = true_Y,
  Ystar = Ystar,
  x_matrix = x_matrix,
  z_matrix = z_matrix,
  true_beta = true_beta,
  true_gamma = true_gamma
)
usethis::use_data(COMBO_EM_data, overwrite = TRUE)
