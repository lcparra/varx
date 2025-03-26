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

add_stimulus = False # add an external stimlus
asymmetry = 1.0 # make the connectivity matrix asymmetric by this factor

ds = Dataset("gw")

# multiply the upper triangle of matrix ds.Cmat with asymmetry
# ds.Cmat[np.triu_indices(ds.Cmat.shape[0], k=2)] *= asymmetry
# multiply a few select rows with asymmetry
ds.Cmat[0:10,:] *= asymmetry
ds.Cmat[64:66,:] *= asymmetry
ds.Cmat[32:36,:] *= asymmetry
ds.Cmat[54:56,:] *= asymmetry

plt.figure(1), plt.clf()
plt.imshow(ds.Cmat)
plt.show()


model = ALNModel(Cmat = ds.Cmat, Dmat = ds.Dmat)

# we chose a parameterization in which the brain network oscillates slowly
# between up- and down-states

model.params['dt'] = 0.1 # in ms
model.params['duration'] = 5*60*1000 # in ms

if add_stimulus:
    # creat a 1/f stimulus in all channels
    L = int(model.params['duration']/model.params['dt'])
    freq = np.fft.rfftfreq(L); freq[0]=1
    noise = np.fft.rfft(np.random.randn(model.params['N'],L))
    noise = np.divide(noise,freq)
    noise = np.fft.irfft(noise)
    noise = (noise - np.mean(noise, axis=1, keepdims=True)) / np.std(noise, axis=1, keepdims=True)
    plt.clf(); plt.plot(noise.T); plt.show()
    noise[20:, :] = 0
    model.params["ext_exc_current"] = noise
    del noise

#    ac_stimulus = stim.SinusoidalInput(amplitude=1, frequency=15.0).to_model(model)
#    ac_stimulus[20:, :] = 0
#    model.params["ext_exc_current"] = ac_stimulus

# model.run(chunkwise=True,bold=False)
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
#    input = signal.decimate(ac_stimulus[0:1, :], q)  # weird trick to keep it a matrix
    input = signal.decimate(model.params["ext_exc_current"], q)
    varx_model = varx(output.T,2,input.T,10,0)
    R_input = varx_model['A_Rvalue'] - np.diag(np.diag(varx_model['A_Rvalue']))
    plt.figure(2), plt.clf()
    plt.subplot(2, 2,1), plt.imshow(ds.Cmat)
    plt.subplot(2,2,2), plt.imshow(R_input)
    plt.subplot(2,2,3), plt.scatter(np.sqrt(ds.Cmat.flatten()), R_input.flatten(), s=1)
    plt.xlabel("True Connectivity")
    plt.ylabel("VARX estimate R")
    corr, pval = spearmanr(R_input.flatten(), ds.Cmat.flatten())
    print(f"Spearman's correlation between A_Rvalue and ds.Cmat (input): {corr:.4f}")

if asymmetry != 1.0:
    plt.clf()
    diff_Cmat = ds.Cmat - ds.Cmat.T
    plt.subplot(2, 2, 1), plt.imshow(diff_Cmat)
    diff_R = R - R.T
    plt.subplot(2, 2, 2), plt.imshow(diff_R)
    plt.subplot(2, 2, 3),
    plt.scatter(diff_Cmat.flatten(), diff_R.flatten(), s=1)
    plt.xlabel("True Connectivity C-C^T")
    plt.ylabel("VARX estimate R-R^T")
    plt.subplot(2, 2, 4), plt.imshow(ds.Cmat)
    plt.show()

# save model.output and ds.Cmat a single matlab file
filename = '../../data/neurolib_model_output.mat'
if asymmetry != 1.0:
    filename = filename[:-4] + f'_asymmetry_{asymmetry:.2f}.mat'  # add asymmetry as string to filename
if add_stimulus:
    sio.savemat(filename[:-4] +  '_stimulus.mat', {'output': output,  'input':input, 'Cmat': ds.Cmat, 'fs': fs})
else:
    sio.savemat(filename[:-4] +  '_rest.mat', {'output': output, 'Cmat': ds.Cmat, 'fs': fs})

