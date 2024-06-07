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