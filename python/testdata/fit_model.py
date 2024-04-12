from scipy.linalg import block_diag
from scipy.linalg import solve
import numpy as np
# fit_model
def fit_model(Rxx, Ryy, ryy, gamma, base):
    # apply basis functions, if available 
    if base[0] is not None:
        B = block_diag(*base)
        Rxx = B.T @ Rxx @ B
        Rxy = B.T @ Rxy
        
    # Regularizer
    Gamma = gamma * np.diag(np.diag(Rxx)) # Tikhonov, scaled for all variables to be regularized equally, regardless of magnitude

    # Least squares estimate with regularization
    h = solve(Rxx + Gamma, Rxy)

    # mean error square
    Rxyest = Rxx @ h
    s2 = (np.sum(h * Rxyest) - 2 * np.sum(h * Rxy) + ryy).T

    # Bias term for ridge regression bias -- see Babadi derivation 
    Bias = np.sum((Rxy - Rxyest) * solve(Rxx, Rxy - Rxyest)).T / s2 / 2 if gamma > 0 else 0

    return h, s2, Bias