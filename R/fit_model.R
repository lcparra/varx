library(MASS) # For generalized inverse functions
library(Matrix) # For dealing with sparse matrices

fit_model <- function(Rxx, Rxy, ryy, gamma, base) {
  # Apply basis functions if available
  if (!is.null(base[[1]])) {
    B <- do.call(bdiag, base)  # Construct a block diagonal matrix from the list of basis matrices
    Rxx <- t(B) %*% Rxx %*% B
    Rxy <- t(B) %*% Rxy
  }
  
  # Regularization
  Gamma <- gamma * diag(diag(Rxx)) # Tikhonov regularization scaled by diagonal elements of Rxx
  
  # Least squares estimate with regularization
  h <- solve(Rxx + Gamma, Rxy)
  
  # Mean error square
  Rxyest <- Rxx %*% h
  s2 <- rowSums(h * Rxyest) - 2 * rowSums(h * Rxy) + ryy
  
  # Bias term for ridge regression bias
  Bias <- if (gamma > 0) {
    residual <- Rxy - Rxyest
    rowSums((residual * solve(Rxx, residual)) / s2) / 2
  } else {
    rep(0, length(ryy))
  }
  
  return(list(h = h, s2 = s2, Bias = Bias))
}