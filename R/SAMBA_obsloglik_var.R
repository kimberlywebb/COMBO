# Note that this code is from the SAMBA R package, created by Lauren Beesley and
## Bhramar Mukherjee

# Beesley, L. and Mukherjee, B. (2020).
#'   Statistical inference for association studies using electronic health records:
#'   Handling both selection bias and outcome misclassification.
#'   Biometrics, 78, 214-226.

SAMBA_obsloglik_var <- function(Dstar, Z, X, theta, beta, beta0_fixed,
                          expected = TRUE)
{
  fixed <- !is.null(beta0_fixed) & length(beta0_fixed) < 2

  n <- length(Dstar)

  X1 <- cbind(1, X)
  Z1 <- cbind(1, Z)

  XBeta <- X1 %*% beta
  ZTheta <- Z1 %*% theta

  #fixed logit(1-specificity). Currently, we only support specificity = 1
  YAlpha <- matrix(-Inf, n)

  expit.zt <- expit(ZTheta)
  expit.xb <- expit(XBeta)
  expit.ya <- expit(YAlpha)

  ### Calculate Derivatives ###
  K1 <- as.vector(expit.xb * expit.zt + expit.ya * (1 - expit.zt))

  exp.xb <- exp(XBeta)
  exp.zt <- exp(ZTheta)

  dK1.dB  <- expit.xb / (1 + exp.xb) * expit.zt
  dK1.dT <- (exp.zt / ((1 + exp.zt) ^ 2)) * (expit.xb -  expit.ya)

  dK1.dBdB <- (1 - exp.xb) * exp.xb * ((1 / (1 + exp.xb)) ^ 3) * expit.zt
  dK1.dBdT <- (exp.xb / ((1 + exp.xb) ^ 2)) * (exp.zt / ((1 + exp.zt) ^ 2))
  dK1.dTdT <- (1 - exp.zt) * exp.zt * ((1 / (1 + exp.zt)) ^ 3)
  dK1.dTdT <- dK1.dTdT * (expit.xb - expit.ya)

  if (expected) {
    tmp <- as.vector(1 / (K1 * (1 - K1)))
    meat.bb <- -dK1.dB * dK1.dB * tmp
    meat.bt <- -dK1.dB * dK1.dT * tmp
    meat.tt <- -dK1.dT * dK1.dT * tmp

  } else {
    tmp <- Dstar / K1 ^ 2
    meat.bb <- tmp * (K1 * dK1.dBdB - dK1.dB * dK1.dB)
    meat.bt <- tmp * (K1 * dK1.dBdT - dK1.dB * dK1.dT)
    meat.tt <- tmp * (K1 * dK1.dTdT - dK1.dT * dK1.dT)

    tmp <- (1 - Dstar) / (1 - K1) ^ 2
    meat.bb <- meat.bb - tmp * ((1 - K1) * dK1.dBdB + dK1.dB * dK1.dB)
    meat.bt <- meat.bt - tmp * ((1 - K1) * dK1.dBdT + dK1.dB * dK1.dT)
    meat.tt <- meat.tt - tmp * ((1 - K1) * dK1.dTdT + dK1.dT * dK1.dT)
  }

  I.betabeta   <- t(apply(X1, 2, function(x1) x1 * as.vector(meat.bb))) %*% X1
  I.betatheta  <- t(apply(X1, 2, function(x1) x1 * as.vector(meat.bt))) %*% Z1
  I.thetatheta <- t(apply(Z1, 2, function(x1) x1 * as.vector(meat.tt))) %*% Z1

  info <- rbind(cbind(I.thetatheta, t(I.betatheta)),
                cbind(I.betatheta, I.betabeta))

  if (fixed) {
    k <- ncol(Z) + 2
    var <- -solve(info[-k , -k])
    var <- cbind(var[, 1:(k - 1)], NA, var[, k:ncol(var)])
    var <- rbind(var[1:(k - 1),], NA, var[k:nrow(var),])
  } else {
    var <- -solve(info)
  }

  var
}
