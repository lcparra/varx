import numpy as np
from scipy.linalg import toeplitz
from scipy.signal import lfilter, lfilter_zi
#myxcorr
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
        lags = lags*np.ones((x.shape[1],1))
    
    # Compute correlations up to largest possible lag
    Q = int(np.max(lags))
    
    # Find valid samples without NaN in all rows of Y and Q history of all X rows
    # z = np.empty((Q-1))
    # z[:] = np.nan
    z1 = lfilter_zi(np.ones((Q)),1)
    z1 = z1*np.nan
    filtered, _ = lfilter(np.ones((Q)),1,np.sum(x,1),zi = z1)
    valid = ~np.isnan(filtered + np.sum(y,1))
    
    # number of valid samples
    T = sum(valid)
    
    # remove offset
    x = x - np.mean(x[valid,:],0)
    y = y - np.mean(y[valid,:],0)   

    # #compute correlations with block-Toeplitz matrices
    X = np.zeros((x.shape[0], int(np.sum(lags))))
    for i in range(x.shape[1]-1, -1, -1):
        startidx = int(np.sum(lags[0:i]))
        endidx = startidx + int(lags[i,0])
        toeplitz_matrix = toeplitz(x[:, i], np.concatenate(([x[0, i]], np.zeros(int(lags[i, 0])-1))))
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