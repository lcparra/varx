# varx
import numpy as np
def varx(Y, na, X, nb, gamma):
    '''
    model = varx(Y,na,X,nb,gamma) fits an vectorial ARX model to the MIMO
    system output Y with input X by minimizing the equation error e(t), i.e.
    equation error model:

    Y(t) = A*Y(t-1) + B*X(t) + e(t)

    where * represents a convolution.  The model contains the following
    variables, stored as stucture elements:

    model = A, B, A_pval, B_pval, A_Deviance,B_Deviance, T

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
    exp(-Devinace/T).

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

    (c) July 10, 2023 Lucas C Parra
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
        
    Ab, s2, Bias = fit_model(Rxx, Rxy, ryy, gamma, base)