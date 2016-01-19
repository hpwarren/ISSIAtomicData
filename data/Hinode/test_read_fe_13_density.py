
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# read file
file = 'eis_l1_20130708_002042.fe_density.h5'
data = h5py.File(file,'r')

# extract some variables
intensities = np.array(data['intensities'])
intensities_error = np.array(data['intensities_error'])
index = np.array(data['index'])
logn = np.array(data['logn'])
emissivity = np.array(data['emissivity'])

# ratio from CHIANTI
ratio_chianti = emissivity[4, :]/emissivity[2, :]
logn_obs_fn = interp1d(ratio_chianti, logn, kind='cubic')

# compute logn for some observed ratios
ratio_obs = []
logn_obs = []
for idx in range(0, 99):
  x = intensities[4, idx]/intensities[2, idx]
  y = logn_obs_fn(x)  
  ratio_obs.append(x)
  logn_obs.append(y)

plt.plot(logn, ratio_chianti)
plt.title('203.826/202.044')
plt.xlabel('Log density')
plt.ylabel('Ratio (energy units)')
plt.plot(logn_obs, ratio_obs, marker='o', linestyle='')
plt.show()







