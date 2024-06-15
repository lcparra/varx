from neurolib.utils.loadData import Dataset
from neurolib.models.aln import ALNModel
import matplotlib.pyplot as plt
import neurolib.utils.functions as func
import numpy as np
import scipy.io as sio
from scipy import signal


ds = Dataset("gw")
fs = 1000 # in Hz
model = ALNModel(Cmat = ds.Cmat, Dmat = ds.Dmat)
model.params['duration'] = 5*60*fs # in ms
model.params['dt'] = 0.1 # in ms


# model.run(bold=True)
model.run()

# R = func.fc(model.BOLD.BOLD[:, 5:])
R = func.fc(model.output, )

plt.figure(1) 
plt.clf()

# donwsample x by a factor of 100 with antialiasing using scipy functions
q = 100
output = signal.decimate(model.output,q)
fs = fs/q

# plt.subplot(2,2,2), plt.imshow(model.BOLD.BOLD[:, 5:],aspect='auto')
plt.subplot(2,1,1), plt.imshow(output,aspect='auto')
plt.subplot(2,2,3), plt.imshow(ds.Cmat), plt.colorbar()
plt.subplot(2,2,4), plt.imshow(R-np.eye(R.shape[0])), plt.colorbar()
plt.show()

# save model.output and ds.Cmat a single matlab file
sio.savemat('../../data/neurolib_5min_rest_model_output.mat', {'output': output, 'Cmat': ds.Cmat, 'fs': fs})

# show the size of model.output
## show the size of model.output
print(model.output.shape)



## 

