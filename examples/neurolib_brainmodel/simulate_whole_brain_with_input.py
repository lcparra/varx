from neurolib.utils.loadData import Dataset
from neurolib.models.aln import ALNModel
import matplotlib.pyplot as plt
import neurolib.utils.functions as func
import numpy as np


ds = Dataset("gw")

model = ALNModel(Cmat = ds.Cmat, Dmat = ds.Dmat)
model.params['duration'] = 5*60*10 # in ms

# model.run(bold=True)
model.run()


## why is this not working
# R = func.fc(model.BOLD.BOLD[:, 5:])
R = func.fc(model.output)

plt.figure(1) 
plt.clf()

# plt.subplot(2,2,2), plt.imshow(model.BOLD.BOLD[:, 5:],aspect='auto')
plt.subplot(2,1,1), plt.imshow(model.output,aspect='auto')
plt.subplot(2,2,3), plt.imshow(ds.Cmat), plt.colorbar()
plt.subplot(2,2,4), plt.imshow(R-np.eye(R.shape[0])), plt.colorbar()
plt.show()

##
# save model.output and ds.Cmat a single matlab file
import scipy.io as sio
sio.savemat('model_output.mat', {'output': model.output, 'Cmat': ds.Cmat})

# show the size of model.output
## show the size of model.output
print(model.output.shape)



## 

