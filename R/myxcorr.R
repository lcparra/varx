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
