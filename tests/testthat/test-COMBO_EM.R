testthat::test_that("COMBO_EM converges and gives fixed answer", {
  starting_values <- rep(1,6)
  beta_start <- matrix(starting_values[1:2], ncol = 1)
  gamma_start <- matrix(starting_values[3:6], ncol = 2, nrow = 2, byrow = FALSE)

  d = COMBO_EM_data
  EM_results <- COMBO_EM(d$Ystar,
                         x_matrix = d$x_matrix,
                         z_matrix = d$z_matrix,
                         beta_start = beta_start,
                         gamma_start = gamma_start)

  set.seed(20231129)
  testthat::expect_named(EM_results,
                         c("Parameter", "Estimates", "SE", "Convergence"))

  testthat::expect_equal(
    EM_results$Parameter,
    c("beta1", "beta2", "gamma11", "gamma21", "gamma12", "gamma22",
      "SAMBA_beta1", "SAMBA_beta2", "SAMBA_gamma11", "SAMBA_gamma21",
      "PSens_beta1", "PSens_beta2", "PSens_gamma12", "PSens_gamma22",
      "naive_beta1", "naive_beta2")
  )

  testthat::expect_equal(
    EM_results$Estimates,
    c(0.948901200733174, -2.35565018127966, 0.537732671019593, 0.948995256168159,
      -0.0314394214697733, -1.42142132233738, 1.02579268978015, -1.13762209307613,
      1.22380032952442, 0.578072505296911, 0.369582202204925, -0.722315049144411,
      -0.755327827106668, -43.3997446544131, 0.383157793606271, -0.715413926060254
    ),
    tolerance = 1e-7
  )

  testthat::expect_identical(
    EM_results$SE,
    c(0.275893850397179, 0.165185172609195, 0.172067103928561, 0.210221989482556,
      0.128817744504468, 0.258025405640404, 0.29776853156754, 0.201198942967142,
      0.303545982852963, 0.35267048809751, 7.02366053029974, 2.81551719137552,
      3.28479253270602, 22.1436043722328, 0.0681163450300179, 0.0753607738803868
    ),
    tolerance = 1e-7
  )

  testthat::expect_true(
    all(EM_results$Convergence, na.rm = TRUE)
  )

})
