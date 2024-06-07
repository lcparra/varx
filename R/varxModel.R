varx <- function(Y, na, X = NULL, nb = NULL, gamma = 0) {
  # Check if X and nb are provided and handle defaults
  if (is.null(X) || all(nb == 0)) {
    X <- NULL
    nb <- 0
  }
  
  # Ensure Y and X are lists (to handle multiple datasets)
  if (!is.list(Y)) {
    Y <- list(Y)
  }
  if (!is.null(X) && !is.list(X)) {
    X <- list(X)
  }
  
  # Dimensions
  ydim <- ncol(Y[[1]])
  xdim <- if (!is.null(X)) ncol(X[[1]]) else 0
  
  # Initialize basis functions and parameters
  base <- vector("list", ydim + (if (!is.null(X)) xdim else 0))
  if (is.matrix(nb)) {
    m <- list(base = nb)
    base <- rep(list(diag(na)), ydim)
    base <- c(base, list(nb))
    nb <- nrow(nb)
    bparams <- ncol(nb)
  } else {
    m <- list(base = NULL)
    bparams <- nb
  }
  
  # Setup lags and parameters
  lags <- matrix(rep(na, ydim + xdim), ncol = 1)
  params <- lags
  
  # Calculate correlations
  Rxx <- 0
  Rxy <- 0
  ryy <- 0
  T <- 0
  for (i in seq_along(Y)) {
    x <- cbind(Y[[i]][-nrow(Y[[i]]), ], if (!is.null(X)) X[[i]][-nrow(X[[i]]), ])
    y <- Y[[i]][-1, ]
    
    # Compute correlations using myxcorr
    res <- myxcorr(x, y, lags)
    Rxx <- Rxx + res$Rxx
    Rxy <- Rxy + res$Rxy
    ryy <- ryy + res$ryy
    T <- T + res$T
  }
  
  # Regularization adjustment
  if (!is.null(gamma)) {
    gamma <- gamma / sqrt(T - sum(lags))
  }
  
  # Fit model
  fit <- fit_model(Rxx, Rxy, ryy, gamma, base)
  
  # Extract coefficients
  A <- array(dim = c(na, ydim, ydim))
  B <- array(dim = c(nb, xdim, ydim))
  
  AB <- fit$h
  s2 <- fit$s2
  Bias <- fit$Bias
  
  # Setup model outputs
  A <- aperm(AB[1:(ydim * na), , drop = FALSE], c(1, 3, 2))
  if (xdim > 0) {
    B <- array(aperm(AB[(ydim * na + 1):length(AB), , drop = FALSE], c(1, 3, 2)), dim = c(nb, ydim, xdim))
  }
  
  m$A <- A
  m$B <- B
  m$B_coeff <- B
  if (!is.null(base[[1]])) {
    m$B <- tensor::ttm(base, B, c(2, 1))
  }
  
  return(m)
}

library(pracma)  # For toeplitz and other numerical functions
library(stats)   # For filter

myxcorr <- function(x, y, lags) {
  if (is.null(lags)) {
    lags <- nrow(x) - 2
  } else {
    lags <- matrix(lags, nrow = ncol(x), ncol = 1)
  }
  
  Q <- max(lags)
  
  # Handling NaN values
  z1 <- rep(NaN, Q)
  filtered <- filter(rep(1, Q), 1, colSums(x, na.rm = TRUE), method = "recursive", init = z1)
  valid <- !is.na(filtered + colSums(y, na.rm = TRUE))
  
  # Number of valid samples
  T <- sum(valid)
  
  # Remove mean
  x[valid, ] <- sweep(x[valid, ], 2, colMeans(x[valid, ], na.rm = TRUE))
  y[valid, ] <- sweep(y[valid, ], 2, colMeans(y[valid, ], na.rm = TRUE))
  
  # Compute correlations using block-Toeplitz matrices
  X <- matrix(0, nrow = nrow(x), ncol = sum(lags))
  for (i in seq_len(ncol(x))) {
    startidx <- sum(lags[1:(i-1)])
    endidx <- startidx + lags[i, 1]
    toeplitz_matrix <- toeplitz(x[, i], c(x[1, i], rep(0, lags[i, 1] - 1)))
    X[, (startidx+1):endidx] <- toeplitz_matrix
  }
  
  Rxx <- t(X[valid, ]) %*% X[valid, ]
  Rxy <- t(X[valid, ]) %*% y[valid, ]
  
  # Power of y
  ryy <- colSums(abs(y[valid, ])^2)
  
  return(list(Rxx = Rxx, Rxy = Rxy, ryy = ryy, T = T))
}

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

library(pracma)  # For the hann window and matrix operations

basis <- function(T, n, type) {
  r <- T / n
  b <- matrix(0, nrow = T, ncol = n)
  
  if (type == 'hanning') {
    for (i in 1:n) {
      b[1:round(r*4), i] <- hann(round(r*4))
    }
    b <- rev(b[(floor(r) + 1):T, ])
    
  } else if (type == 'normal') {
    t <- 1:T
    for (i in 1:n) {
      b[, i] <- exp(-(t - (i * r))^2 / r^2)
    }
    
  } else {
    b <- diag(1, T, n)
  }
  
  return(b)
}