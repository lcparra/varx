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
