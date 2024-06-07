# Load necessary libraries
library(R.matlab)

# Function to source if saved separately
source("R/varxModel.R")  # Assuming all your functions are saved here

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

trace(varx, browser, at = 2) 
# Call the varx function from your sourced R scripts
model <- varx(y, na, x, nb, gamma)

# Optionally print the model to check outputs
print(model)
