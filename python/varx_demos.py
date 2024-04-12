# # test varx (without basis)
import varxModel
import scipy.io as sio
x = sio.loadmat('testdata/x.mat')['x']
y = sio.loadmat('testdata/y.mat')['y']
L = 10
na = 10; nb = 20; gamma = 0

model = varxModel.varx(y,na,x,nb,gamma)
varxModel.varx()


    