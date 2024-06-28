from neurolib.utils.loadData import Dataset
from neurolib.models.aln import ALNModel
import matplotlib.pyplot as plt
import neurolib.utils.functions as func
import neurolib.utils.stimulus as stim
import numpy as np
import scipy.io as sio
from scipy import signal

# adding varx to the system path
import sys
sys.path.insert(0, '../../python')
from varxModel import varx

add_stimulus = False

ds = Dataset("gw")
model = ALNModel(Cmat = ds.Cmat, Dmat = ds.Dmat)

# we chose a parameterization in which the brain network oscillates slowly
# between up- and down-states

model.params['dt'] = 0.1 # in ms
model.params['duration'] = 5*60*1000 # in ms

if add_stimulus:
    # now we create multi-d input of 10Hz
    ac_stimulus = stim.SinusoidalInput(amplitude=1, frequency=15.0).to_model(model)

    # We set the input to a bunch of nodes to zero.
    ac_stimulus[20:, :] = 0

    # assign to input current
    model.params["ext_exc_current"] = ac_stimulus

# model.run(chunkwise=True,bold=True)
model.run()

# donwsample x by a factor of 100 with antialiasing using scipy functions
q = 100
fs = 1/model.params.dt*1000/q # in Hz
output = signal.decimate(model.output,q)

# output = output[:,round(fs*5):] # cut out the first 5 seconds of the data

varx_model = varx(output.T,2)
R = varx_model['A_Rvalue'] - np.diag(np.diag(varx_model['A_Rvalue']))

plt.figure(1), plt.clf()
plt.subplot(2,1,2),
frs, powers = func.getMeanPowerSpectrum(output, dt=1/fs*1000, spectrum_windowsize=0.5)
plt.plot(frs, powers, c="k")
plt.subplot(2,1,1), plt.imshow(output,aspect='auto')
plt.show()

plt.figure(2), plt.clf()
plt.subplot(2,2,1), plt.imshow(ds.Cmat)
plt.subplot(2,2,2), plt.imshow(R)
plt.subplot(2,2,3), plt.scatter(np.sqrt(ds.Cmat.flatten()), R.flatten(), s=1)
plt.xlabel("True Connectivity")
plt.ylabel("VARX estimate R")
plt.show()

# compute Spearman's correlation between A_pval and ds.Cmat
from scipy.stats import spearmanr
corr, pval = spearmanr(R.flatten(), ds.Cmat.flatten())
print(f"Spearman's correlation between A_Rvalue and ds.Cmat (no input): {corr:.4f}")

if add_stimulus:
    input = signal.decimate(ac_stimulus[0:1, :], q)  # weird trick to keep it a matrix
    varx_model = varx(output.T,2,input.T,10,0)
    R_input = varx_model['A_Rvalue'] - np.diag(np.diag(varx_model['A_Rvalue']))
    plt.figure(2), plt.clf()
    plt.subplot(2, 2,1), plt.imshow(ds.Cmat)
    plt.subplot(2,2,2), plt.imshow(R_input)
    plt.subplot(2,2,3), plt.scatter(np.sqrt(ds.Cmat.flatten()), R_input.flatten(), s=1)
    plt.xlabel("True Connectivity")
    plt.ylabel("VARX estimate R")
    plt.show()
    corr, pval = spearmanr(R_input.flatten(), ds.Cmat.flatten())
    print(f"Spearman's correlation between A_Rvalue and ds.Cmat (input): {corr:.4f}")


# save model.output and ds.Cmat a single matlab file
if add_stimulus:
    sio.savemat('../../data/neurolib_5min_stimulus_model_output.mat', {'output': output, 'Cmat': ds.Cmat, 'fs': fs})
else:
    sio.savemat('../../data/neurolib_5min_rest_model_output.mat', {'output': output, 'Cmat': ds.Cmat, 'fs': fs})



#
