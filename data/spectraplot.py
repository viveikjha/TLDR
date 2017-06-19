import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty
from numpy import loadtxt


flux=np.loadtxt("rvm_fluxes.csv",delimiter=",")
errflux=np.loadtxt("rvm_errfluxes.csv",delimiter=",")
waves=np.loadtxt("rvm_wavelengths.csv",delimiter=",")
fs=flux[:,0]
efs=errflux[:,0]
PlotPretty.pp('white')

halpha=656.3*10.0
hbeta=486.1*10.0
hgamma=434.0*10.0
Heii=4686.0
Hei=5876.0
ymin=-0.5
ymax=3.5
plt.figure()
plt.minorticks_on()
plt.ylim(ymin,ymax)
plt.plot(waves,fs)

plt.plot([halpha,halpha],[ymin,ymax],'r--')
plt.annotate(r'$H\alpha$', xy=(halpha+3.0, 3.0))
plt.plot([hbeta,hbeta],[ymin,ymax],'r--')
plt.annotate(r'$H\beta$', xy=(hbeta+3.0, 3.0))
plt.plot([hgamma,hgamma],[ymin,ymax],'r--')
plt.annotate(r'$H\gamma$', xy=(hgamma+3.0, 3.0))

plt.plot([Hei,Hei],[ymin,ymax],'g--')
plt.annotate(r'$HeI$', xy=(Hei-250.0, 2.75))
plt.plot([Heii,Heii],[ymin,ymax],'g--')
plt.annotate(r'$HeII$', xy=(Heii-300.0, 2.75))


plt.xlabel(r"$Wavelength$ $(\AA)$")
plt.ylabel(r"$Flux$")
plt.title(r"$Broad$ $Line$ $Spectra$")
#plt.savefig("BLRspectralabelled.jpg")

fs=flux[:,10]
plt.plot(waves,fs, 'y')





plt.show()
