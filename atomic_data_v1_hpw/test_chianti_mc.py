'''
reads hdf5 file output from 'chianti_ratio_mc.pro'
plots line ratios as a function of density
'''

import h5py
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# read file
file = 'fe_13.monte_carlo.h5'
mc = h5py.File(file,'r')

# assign datasets to variables
ioneq = mc['ioneq_file']
chianti_perturb_aval = mc['chianti_perturb_aval']
chianti_perturb_spl1 = mc['chianti_perturb_spl1']
chianti_perturb_spl2 = mc['chianti_perturb_spl2']
time_stamp = np.array(mc['time_stamp'])
print time_stamp

# get data arrays
logn = np.array(mc['logn'])
wvl = np.array(mc['wavelength'])
trans = np.array(mc['transition'])
em = np.array(mc['emissivity'])
text = np.array(mc['text'])

print text

for t in trans:
  print t

# specify line ratios
l1 = [203.82600403,200.02099609,201.12100220,203.15199280]
l2 = [202.04400635,202.04400635,202.04400635,202.04400635]

# get index for line ratios
l1ndx = [np.where(wvl==i)[0][0] for i in l1]
l2ndx = [np.where(wvl==i)[0][0] for i in l2]

# get ratios
ratio = em[l1ndx,:,:]/em[l2ndx,:,:]

# make plots
plt.figure(figsize=(10, 8), dpi=100, facecolor='w', edgecolor='k')
plt.suptitle('Chianti Monte Carlo',fontsize=12)
plt.subplot(221)
plt.plot(logn,ratio[0,:])
plt.title('203.826/202.044')
plt.ylabel('Line Ratio',fontsize=10)
plt.subplot(222)
plt.plot(logn,ratio[1,:])
plt.title('200.021/202.044')
plt.subplot(223)
plt.plot(logn,ratio[2,:])
plt.title('201.121/202.044')
plt.xlabel('log(n$_{e}$)',fontsize=10)
plt.ylabel('Line Ratio',fontsize=10)
plt.subplot(224)
plt.plot(logn,ratio[3,:])
plt.title('203.152/202.044')
plt.xlabel('log(n$_{e}$)',fontsize=10)
plt.show()
