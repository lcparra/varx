
tensor_product <- function(A, B, dimA, dimB) {
  # Check the dimensions
  if (length(dim(A)) < dimA || length(dim(B)) < dimB) {
    stop("Invalid dimensions.")
  }
  
  # Get the dimensions of the input matrices
  dim_A <- dim(A)
  dim_B <- dim(B)
  
  # Initialize the output matrix
  C <- matrix(0, nrow = dim_A[1], ncol = dim_B[2])
  
  # Perform the tensor product
  for (i in 1:dim_A[1]) {
    for (j in 1:dim_B[2]) {
      C[i, j] <- sum(A[i, ] * B[, j])
    }
  }
  
  return(C)
}
library(tensor)
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
    bparams <- ncol(nb)
    nb <- nrow(nb)
  } else {
    m <- list(base = NULL)
    bparams <- nb
  }
  
  # Setup lags and parameters
  lags <- matrix(rep(na, ydim), ncol = 1)
  params <- lags
  
  # ... unless using basis function and need only including when modeling MA of external input
  if (nb != 0) {
    # Append 'xdim' number of rows with all elements being 'nb' to 'lags'
    lags <- rbind(lags, matrix(rep(nb, xdim), ncol = 1))
    
    # Append 'xdim' number of rows with all elements being 'bparams' to 'params'
    params <- rbind(params, matrix(rep(bparams, xdim), ncol = 1))
  }

  # Calculate correlations
  Rxx <- 0
  Rxy <- 0
  ryy <- 0
  T <- 0
  for (i in seq_along(Y)) {
    x <- Y[[i]][-nrow(Y[[i]]), ]
    y <- Y[[i]][-1, ]

    if (!is.null(X)) {
      x <- cbind(x, X[[i]][-1, ])
    }
    
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
  A <- aperm(array(AB[1:(ydim*na),], dim = c(na, ydim, ydim)), c(1, 3, 2))
  if (xdim > 0) {
    B <- aperm(array(AB[((1 + ydim*na):nrow(AB)),], dim = c(params[length(params)], xdim, ydim)), c(1, 3, 2))
  }
  
  m$A <- A
  m$B <- B
  m$B_coeff <- B
  if (!is.null(base[[1]])) {
    dims <- dim(B); dims <- dims[dims != 1]; B <- array(B, dim = dims)

    m$B <- tensor::tensor(base[[ydim+1]], B, alongA = 2, alongB = 1)
  }
  
  

  # Granger Causal test for all inputs (external and recurrent)
  #initialize deviance and pval matrices size is 2,xdim
  xdim <- ncol(x)
  Deviance <- matrix(0, nrow = 2, ncol = xdim)
  pval <- matrix(0, nrow = 2, ncol = xdim)
  iterarray = 1:xdim
  for (i in xdim:1) { # same as above but with reduced model removing i-th input
    ii <- 1:sum(lags)
    idxaddition <- if (i == 1) 0 else sum(lags[1:(i-1)])
    ii <- ii[-((1:lags[i]) + idxaddition)] # use inputs excluding i-th
    fit <- fit_model(Rxx[ii, ii], Rxy[ii, ], ryy, gamma, base[setdiff(iterarray,i)])
    df <- T - sum(params) # degrees of freedom of the full model
    Deviance[, i] <- df * log(fit$s2 / s2) - T * fit$Bias + T * Bias # not the exact formula, but I calibrated and seems to work well for small T
    pval[, i] <- 1 - pchisq(Deviance[, i], params[i])
  }

  # store additional outputs in model structure
  m$A_pval <- pval[, 1:ydim]
  m$B_pval <- pval[, (ydim+1):ncol(pval)]
  m$A_Deviance <- Deviance[, 1:ydim]
  m$B_Deviance <- Deviance[, (ydim+1):ncol(Deviance)]
  m$A_Rvalue <- sqrt(1 - exp(-m$A_Deviance / T))
  m$B_Rvalue <- sqrt(1 - exp(-m$B_Deviance / T))
  m$T <- T
  m$s2 <- s2 / T


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
  z <- rep(NA, Q-1) 
  filtered <- filter(rowSums(x, na.rm = TRUE), rep(1, Q), sides = 1, init = z)
  sum_y <- rowSums(y, na.rm = TRUE)
  valid <- !is.na(filtered + sum_y)
  
  # Number of valid samples
  T <- sum(valid)
  
  # Remove mean
  x<- sweep(x, 2, colMeans(x[valid, ], na.rm = TRUE))
  y <- sweep(y, 2, colMeans(y[valid, ], na.rm = TRUE))
  
  # Compute correlations using block-Toeplitz matrices
  X <- matrix(0, nrow = nrow(x), ncol = sum(lags))
  for (i in ncol(x):1) {
    idxaddition <- if (i == 1) 0 else sum(lags[1:(i-1)])
    X[,(1:lags[i]) + idxaddition] <- Toeplitz(x[, i], c(x[1, i], rep(0, lags[i] - 1)))
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
  s2 <- t(colSums(h * Rxyest) - 2 * colSums(h * Rxy) + ryy)


  # Bias term for ridge regression bias
  Bias <- if (gamma > 0) {
    residual <- Rxy - Rxyest
    t(colSums((residual * solve(Rxx, residual))) / s2 / 2)
  } else {
    0
  }
  
  return(list(h = h, s2 = t(s2), Bias = Bias))
}

library(pracma)  # For the matrix operations
library(gsignal, include.only = c("hann"))
# hann <- function(N) {
#   0.5 * (1 - cos(2 * pi * (0:(N-1)) / (N-1)))
# }

basis <- function(T, n, type) {
  r <- T / n
  b <- matrix(0, nrow = T*2, ncol = n)
  
  if (type == 'hanning') {
    for (i in 1:n) {
      b[(1:ceil(r*4)) + ceil((i-1)*r), i] <- hann(round(r*4))
    }
    b <- b[(floor(r) + (1:T)), ]
    b <- b[rev(seq_len(nrow(b))),]
    
  } else if (type == 'normal') {
    t <- 1:T
    for (i in 1:n) {
      b[, i] <- exp(-(t-(i-1)*r)^2/r^2)
    }
    b <- b[1:T,]
    
  } else {
    b <- diag(1, T, n)
  }
  
  return(b)
}