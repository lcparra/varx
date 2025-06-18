import numpy as np
from scipy.linalg import toeplitz
from scipy.signal import lfilter, lfilter_zi
# MYXCORR
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


# BASIS
import numpy as np
from scipy.signal import windows
def basis(T,n,type):
    '''
    b=basis(T,n,type) makes basis functions of type 'hanning' or 'normal'. T
    is the length of the basis. n is how many to use. Typically T>>n, if the
    goal is to represent a filter with fewer parameters. Omit output argument
    to see how the basis functions look.
    '''
    
    r = T / n
    b = np.zeros((T, n))

    if type == 'hanning':
        for i in range(n):
            b[:int(round(r*4)), i] = windows.hann(int(round(r*4)))
        b = np.flipud(b[int(np.floor(r)):T, :])

    elif type == 'normal':
        t = np.arange(1, T+1)
        for i in range(n):
            b[:, i] = np.exp(-np.power(t - (i * r), 2) / np.power(r, 2))

    else:
        b = np.eye(T)
        
    return b

from scipy.linalg import block_diag
from scipy.linalg import solve
# FIT_MODEL
def fit_model(Rxx, Rxy, ryy, gamma, base):
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
    s2 = (np.sum(h * Rxyest, 0) - 2 * np.sum(h * Rxy,0) + ryy).T

    # Bias term for ridge regression bias -- see Babadi derivation 
    Bias = np.sum((Rxy - Rxyest) * solve(Rxx, Rxy - Rxyest),0).T / s2 / 2 if gamma > 0 else 0

    return h, s2, Bias

# VARX
import numpy as np
from scipy.stats import chi2
def varx(Y, na, X=None, nb=0, gamma=0):
    '''
    model = varx(Y,na,X,nb,gamma) fits an vectorial ARX model to the MIMO
    system output Y with input X by minimizing the equation error e(t), i.e.
    equation error model:

    Y(t) = A*Y(t-1) + B*X(t) + e(t)

    where * represents a convolution.  The model contains the following
    variables, stored as stucture elements:

    model = A, B, A_pval, B_pval, A_Deviance,B_Deviance, A_Rvalue,B_Rvalue,s2,T

    A and B are filter model parameters found with conventional least squares
    with ridge regression. They are stored as tensor of size [na,ydim,ydim]
    and [nb,ydim,xdim] respectively. na and nb are the legth of the filters.
    gamma is the regularization for the ridge (shrinkage) regularization and
    defaults to 0 and should not be selected larger than 1. Note that x(t)
    represents the history including the current sample in the input. Thus,
    we are allowing for instant effects. This is the norm in the signal
    processing literature but no in the Granger Causality VAR models,
    although there is no theoretical reason not to include instant effect in
    the external input. To avoid instant effects, the user can simply delay
    the input by one sample.

    A_pval,B_pval are  P-values for each channel (for all delays together)
    using the Deviance formalism.

    A_Deviance, B_Deviance ,T are Deviance and number of sample used in the
    estimation of p-values, and Deviance/T can serve as a measure of effect
    size, and can be used to compute generalized R-square: R2 = 1 -
    exp(-Devinace/T). exp(-Devinace/T). These are returned as A_Rvalue, B_Rvalue.
    s2 is the mean squared error.

    varx(Y,na,X,base,gamma) If base is not a scalar, it is assumed that it
    represent basis functions for filters B of size [filter length, number of
    basis functions]. B will have size [size(base,2),ydim,xdim], i.e. as many
    parameters for each path as basis functions. The actual filters can be
    obtained as tensorprod(base,B,2,1);

    varx(Y,na,X,nb,gamma) If Y is a cell array, then the model is fit on all
    data records in X and Y. All elements in the cell arrays X and Y have to
    have the same xdim and ydim, but may have different numer of rows (time
    samples).

    varx(Y,na) Only fitst the AR portion. To provide gamma, set set x=[] and
    nb=0.

    If the intention is to only fit a MA model, then the Granger formalism
    requires at least an AR portion for each output channel, without the
    interaction between ouput channels. If that is the intention, then one
    should call this function for each output channel separatelly, e.g.
    varx(y(:,i),na,x,nb)

    model can be used by varx_display(model) for display. 

    (c) (matlab) July 10, 2023 Lucas C Parra
        04/11/2024, last version, Lucas Parra
    (c) (python) April 12, 2024, Aimar Silvan, based on matlab version 04/11/2024

    '''
    # If not simulating eXternal MA channel then xdim=0
    if X is None or np.all(nb == 0):
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
    if isinstance(nb, np.ndarray):
        m = {'base': nb}  # save for output
        # internally, base variable is a list with one base for each dimension
        base = [np.eye(na) for _ in range(ydim)]
        base.extend([nb for _ in range(xdim)])
        nb, bparams = nb.shape  # lags according and number of parameters according to basis
    else:
        m = {'base': None}  # save for output
        base = [None for _ in range(ydim + xdim)]  # empty bases
        bparams = nb  # number of parameters same as number of lags
    
    # number of lags and number of parameters equal ...
    lags = np.ones((ydim, 1)) * na
    params = lags.copy()
    # ... unless using basis function and need only including when modeling MA of external input
    if nb:
        lags = np.concatenate((lags, np.ones((xdim, 1)) * nb))
        params = np.concatenate((params, np.ones((xdim, 1)) * bparams))

    # calculate correlations
    Rxx = 0
    Rxy = 0
    ryy = 0
    T = 0
    for i in range(len(Y)):
        # Set preceding output and input both as input to the LS problem
        x = Y[i][:-1, :]
        y = Y[i][1:, :]
        
        # if modeling also the MA of eXternal input
        if nb:
            x = np.concatenate((x, X[i][1:, :]), axis=1)

        # Compute auto and cross correlations
        Rxx_, Rxy_, ryy_, T_ = myxcorr(x, y, lags)

        # accumulate over all data records
        Rxx += Rxx_
        Rxy += Rxy_
        ryy += ryy_
        T += T_
        
    if gamma is None:
        gamma = 0
    else:
        gamma = gamma/np.sqrt(T-np.sum(lags)) # regularization decreasing with degrees of freedom
        
    AB, s2, Bias = fit_model(Rxx, Rxy, ryy, gamma, base)
    
    
    A = np.transpose(AB[0:ydim*na, :].reshape(na, ydim, ydim, order='F'), (0, 2, 1)) # F so it follows Fortran-style order (consistent with matlab)
    B = np.squeeze(np.transpose(AB[ydim*na:].reshape(int(params[-1,0]), xdim, ydim, order='F'), (0, 2 ,1)))

    m['A'] = A
    m['B'] = B
    
    # if we used a base, return filters B with base applied
    m['B_coeff'] = m['B']
    if base[0] is not None:
        m['B'] = np.tensordot(m['base'], m['B_coeff'], axes=([1], [0]))
    
    # Granger Causal test for all inputs (external and recurrent)
    xdim = x.shape[1]
    Deviance = np.zeros((ydim, xdim))
    pval = np.zeros((ydim, xdim))
    for i in range(xdim-1, -1, -1):
        ii = np.arange(0, np.sum(lags)).astype(int)
        startidx = int(np.sum(lags[0:i]))
        endidx = startidx + int(lags[i,0])
        ii = np.delete(ii, np.arange(startidx, endidx).astype(int))
        
        _, s2r, Biasr = fit_model(Rxx[ii, :][:, ii], Rxy[ii, :], ryy, gamma, [base[j] for j in range(xdim) if j != i])
        
        df = T - np.sum(params)  # degrees of freedom of the full model
        
        Deviance[:, i] = np.maximum(df * np.log(s2r / s2) - T * Biasr + T * Bias,0)  # not the exact formula, but I calibrated and seems to work well for small T. For large parameter count this approximation may be negative, hence the maximum()
        
        pval[:, i] = 1 - chi2.cdf(Deviance[:, i], params[i])
    
    # store additional outputs in model dictionary
    m['A_pval'] = pval[:, :ydim]
    m['B_pval'] = pval[:, ydim:]
    m['A_Deviance'] = Deviance[:, :ydim]
    m['B_Deviance'] = Deviance[:, ydim:]
    m['T'] = T
    
    # Compute effect values and add to output
    m['A_Rvalue'] = np.sqrt(1 - np.exp(-m['A_Deviance'] / m['T']))
    m['B_Rvalue'] = np.sqrt(1 - np.exp(-m['B_Deviance'] / m['T']))
    
    return m