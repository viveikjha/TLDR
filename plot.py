import pickle
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import UnivariateSpline as UVS
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import minimize

'''IMPORT DATA FILES'''
inputfile = 'vdm.csv'
tdf = np.transpose(np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0))
print('TDF: ',tdf.shape)
nlines = len(tdf) # length first axis

inputfile = '../../data/rvm_wavelengths.csv'
rvm_spectra_lam = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
print('lam: ',len(rvm_spectra_lam))
nlines2 = len(rvm_spectra_lam)



#must correct for redshift to arp 151: z = 0.021091
lamHB = 5500.0
vmin = -4000.0
c = 299792.458
lammin = (vmin/c*lamHB)+lamHB

vmax = 4000.0
lammax = (vmax/c*lamHB)+lamHB
indexHB = np.argmin(abs(rvm_spectra_lam - lamHB))
print 'indexHB: ', indexHB
lamarr = rvm_spectra_lam[np.where(rvm_spectra_lam >= lammin)]
lamarr = rvm_spectra_lam[np.where(lamarr <= lammax)]

arr = np.empty([len(lamarr),50])
for i in range(0,nlines):
	for j in range(0,len(lamarr)):
		if lamarr[j] == rvm_spectra_lam[i]:
			#print '!!'
			for k in range(0,50):
				arr[j][k]=tdf[i+24][k]
#print arr		
print lammin
print lammax

print tdf[324]


fig, ax = plt.subplots(figsize=(6,8))
ax.imshow(np.rot90((np.sqrt(arr)),1), cmap='Reds',extent=[-4000,4000,0, 20])
Scale=1.5
ax.set_aspect(Scale*8000.0/(20.0))
ax.set_xlabel('V(km/s)')
ax.set_ylabel(r'$\tau$ (days)')
plt.title(r'Velocity Delay Map of $H\beta$')

fig, ax = plt.subplots(figsize=(12,8))
ax.imshow(np.rot90((np.sqrt(tdf)),1), cmap='Reds',extent=[min(rvm_spectra_lam),max(rvm_spectra_lam),0, 50],aspect="auto")
plt.title('TDF Reconstruction')
ax.set_xlabel(r'$\lambda$ (nm)')
ax.set_ylabel(r'$\tau$ (days)')
print len(tdf[0])
ax.set_aspect(10)




plt.show()


