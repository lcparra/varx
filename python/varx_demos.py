# # test varx (without basis)
import varxModel
import scipy.io as sio
x = sio.loadmat('testdata/x.mat')['x']
y = sio.loadmat('testdata/y.mat')['y']
L = 10
na = 10; nb = 20; gamma = 0

model = varxModel.varx(y,na,x,nb,gamma)

# # test varx (with basis)
x = sio.loadmat('x_basis.mat')['x']
y = sio.loadmat('y_basis.mat')['y']
L = 10
na = 3; nb = 4; gamma = 0.1
base = basis(30, nb, 'normal')
model = varxModel.varx(y,na,x,base,gamma)
    
