# Function to check if a package is installed
is_package_installed <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Install and load required packages
required_packages <- c("pracma", "MASS", "Matrix", "R.matlab","gsignal","tensor")
for (pkg in required_packages) {
  if (!is_package_installed(pkg)) {
    install.packages(pkg)
  }
}

# Load the model
source("R/varxModel.R")  # Assuming all your functions are saved here

# # Load necessary libraries
library(R.matlab)

# Load data from .mat files
x_data <- readMat('R/testdata/x.mat')
y_data <- readMat('R/testdata/y.mat')

# Extract matrices from list
x <- x_data$x
y <- y_data$y

# Parameters for the VARX model
L <- 10
na <- 10
nb <- 20
gamma <- 0

model <- varx(y, na, x, nb, gamma)

print(model)



## Example with basis functions

# Load data from .mat files
x_data <- readMat('R/testdata/x_basis.mat')
y_data <- readMat('R/testdata/y_basis.mat')

# Extract matrices from list
x <- x_data$x
y <- y_data$y

# Parameters for the VARX model
L <- 10
na <- 3; nb <- 4; gamma <- 0.1

# Basis functions
base = basis(30, nb, 'hanning')
# base = basis(30,nb,'hanning')
model = varx(y,na,x,base,gamma)

print(model)