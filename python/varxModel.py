import numpy as np
from scipy.linalg import toeplitz
from scipy.signal import lfilter

def varx(Y, na, X=None, nb=0, lambda_=0):
    # If X is not provided or nb is 0, set X to None and nb to 0
    if X is None or nb == 0:
        X = None
        nb = 0

    # Make Y and X into lists if they are not already, to handle multiple data records
    if not isinstance(Y, list):
        Y = [Y]
    if not isinstance(X, list):
        X = [X]

    # Get dimensions
    ydim = Y[0].shape[1]
    xdim = X[0].shape[1] if X[0] is not None else 0

    # Initialize basis functions
    if isinstance(nb, list) or isinstance(nb, np.ndarray):
        base = [np.eye(na) for _ in range(ydim)]
        for _ in range(xdim):
            base.append(nb)
        nb, bparams = nb.shape[0], nb.shape[1]
    else:
        base = [None] * (ydim + xdim)
        bparams = nb

    lags = np.ones((ydim, 1)) * na
    params = lags.copy()

    if nb:
        lags = np.vstack((lags, np.ones((xdim, 1)) * nb))
        params = np.vstack((params, np.ones((xdim, 1)) * bparams))

    # Initialize correlation matrices
    Rxx = np.zeros((na * (ydim + xdim), na * (ydim + xdim)))
    Rxy = np.zeros((na * (ydim + xdim), ydim))
    ryy = np.zeros(ydim)
    T = 0

    for i in range(len(Y)):
        x = Y[i][:-1]
        y = Y[i][1:]

        if nb:
            x = np.hstack((x, X[i][1:])) if X[i] is not None else x

        # Compute auto and cross correlations
        Rxx_, Rxy_, ryy_, T_ = myxcorr(x, y, lags)
        Rxx += Rxx_
        Rxy += Rxy_
        ryy += ryy_
        T += T_

    # Regularization
    gamma = lambda_ / np.sqrt(T - np.sum(lags))
    Gamma = gamma * np.diag(np.diag(Rxx))

    # Fit model
    AB, s2, Bias = fit_model(Rxx, Rxy, ryy, Gamma, base)

    # Reshape and store filter values
    A = np.transpose(np.reshape(AB[:ydim * na], (na, ydim, ydim)), (0, 2, 1))
    B = np.transpose(np.reshape(AB[ydim * na:], (params[-1], xdim, ydim)), (0, 2, 1))

    # Apply basis functions if available
    B_coeff = B.copy()
    if base[0] is not None:
        B = tensorprod(base, B_coeff, 2, 1)

    # Granger Causal test
    Deviance = np.zeros((na + nb, xdim))
    pval = np.zeros((na + nb, xdim))
    for i in range(xdim, 0, -1):
        ii = np.delete(np.arange(np.sum(lags)), np.arange((i - 1) * nb, i * nb))
        AB_r, s2_r, Bias_r = fit_model(Rxx[ii][:, ii], Rxy[ii], ryy, Gamma, base[:i - 1] + base[i:])
        df = T - np.sum(params)
        Deviance[:, i - 1] = df * np.log(s2_r / s2) - T * Bias_r + T * Bias
        pval[:, i - 1] = 1 - chi2.cdf(Deviance[:, i - 1], params.flatten())

    # Store additional outputs
    model = {
        "A": A,
        "B": B,
        "B_coeff": B_coeff,
        "A_pval": pval[:, :ydim],
        "B_pval": pval[:, ydim:],
        "A_Deviance": Deviance[:, :ydim],
        "B_Deviance": Deviance[:, ydim:],
        "T": T
    }

    return model


def myxcorr(x,y,lags):
    '''
    Computes auto- and corr-correlation matrices Rxx, Rxy and ryy defined as
    (after the mean is subtracted from x and y):

        Rxx(l)= sum_n x'(n)*x(n+l)
        Rxy(l)= sum_n x'(n)*y(n+l)
        ryy   = sum_n |y(n)|^2


    This is computed in the frequency domain, returning only zero and
    positive delays l = 0 ... lags-1. Note that these definitions differ from
    the standard definition of the cross correlation and from matlab's
    conventional xcorr() implementation. Here y is delayed, not x, which is
    equivalent to keeping the negative delays in the standard definition.

    Rxx and Rxy are arranged as required to implement conventional MIMO
    system identification (see example below). Rxx is a square block-Toeplitz
    auto-correlation matrix of size sum(lags) and Rxy is a cross-correlation
    block-matrix of size [sum(lags), ydim]. lags is a vector of lags used for
    each of xdim columns in x. 

    x and y are required input, all others are optional. lags defaults to
    size(x,1)-1 and is converted into a vector of size xdim if it is given as
    a scalar.  x and y should be arranged as time by dimensions. y should
    not be longer in time than x; if necessary, pad x with zeros before
    calling.

    Correlations are computed using Toeplitz matrices. x or y may contain NaN
    which are removed from the sum over n which have NaN in any row of y, or
    do not have a max(lags) history in the rows of x. T is returned to report
    how many samples were used. The max(lags) samples at the start are also
    removed (this is known as the covariance method, see Ljung, System
    Identification, Ch. 10.1). Mean subtraction is done with the mean of
    valid samples only.

    This function can be used for MIMO FIR identification as follows: 

    [Rxx,Rxy] = myxcorr(x,y,lags);
    h = inv(Rxx)*Ry;
    yest = filterMIMO(h,x,lags);
    error = y-yest;

    (c) (matlab) April 21, 2021, Lucas C Parra
                 12/12/2023, last version, Lucas Parra
    (c) (python) April 12, 2024, Aimar Silvan, based on matlab version 12/12/2023
    '''
    
    # Deal with inputs
    if lags is None:
        lags = x.shape[0] - 1 -1
    else:
        lags = lags*np.ones((1,x.shape[1]))
    
    # Compute correlations up to largest possible lag
    Q = int(np.max(lags))
    
    # Find valid samples without NaN in all rows of Y and Q history of all X rows
    z = np.empty((Q-1))
    z[:] = np.nan
    valid = ~np.isnan(lfilter(np.ones((Q)),1,np.sum(x,1)) + np.sum(y,1))
    
    # number of valid samples
    T = sum(valid)
    
    # remove offset
    x = x - np.mean(x[valid,:],0)
    y = y - np.mean(y[valid,:],0)   

    # #compute correlations with block-Toeplitz matrices
    X = np.zeros((x.shape[0], int(np.sum(lags))))
    for i in range(x.shape[1]-1, -1, -1):
        startidx = i*int(lags[0,i])
        endidx = int(startidx + (lags[0,i]))
        toeplitz_matrix = toeplitz(x[:, i], np.concatenate(([x[0, i]], np.zeros(int(lags[0, i])-1))))
        X[:, startidx:endidx] = toeplitz_matrix

    Rxx = X[valid, :].T @ X[valid, :]
    Rxy = X[valid, :].T @ y[valid, :]
    
    # Power of y
    ryy = np.sum(np.abs(y[valid,:])**2,0)
    
    return Rxx, Rxy, ryy, T 
    
    #testcode
    #data is same from 'testcode fir MIMO  MA identification' in 'matlab/myxcorr.m'
    import scipy.io as sio
    x = sio.loadmat('testdata/x.mat')['x']
    y = sio.loadmat('testdata/y.mat')['y']
    L = 10
    [Rxx,Rxy,ryy,T] = myxcorr(x,y,L)
    




def fit_model(Rxx, Rxy, ryy, gamma, base):
    if base[0] is not None:
        B = np.block_diag(*base)
        Rxx = np.matmul(np.matmul(B.T, Rxx), B)
        Rxy = np.matmul(B.T, Rxy)

    Gamma = gamma * np.diag(np.diag(Rxx))
    h = np.linalg.solve(Rxx + Gamma, Rxy)
    Rxy_est = np.matmul(Rxx, h)
    s2 = np.mean(np.multiply(h, Rxy_est) - 2 * np.multiply(h, Rxy) + ryy, axis=0)

    if np.any(s2 < 0):
        print("Square error of linear regression can not be negative. Please debug ...")
        # Handle error appropriately

    if gamma > 0:
        Bias = np.sum(np.multiply(Rxy - Rxy_est, np.linalg.solve(Rxx, Rxy - Rxy_est)), axis=0) / s2 / 2
    else:
        Bias = 0

    return h, s2, Bias


def tensorprod(base, B_coeff, axis1, axis2):
    B = np.zeros((B_coeff.shape[0], B_coeff.shape[2], base[0].shape[0]))
    for i in range(B_coeff.shape[0]):
        B[i] = np.matmul(base[i].T, B_coeff[i])

    return B
