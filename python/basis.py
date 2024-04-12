# basis
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