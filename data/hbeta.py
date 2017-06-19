import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty
from numpy import loadtxt


flux=np.loadtxt("rvm_fluxes_trimmed.csv",delimiter=",")
errflux=np.loadtxt("rvm_errfluxes_trimmed.csv",delimiter=",")
waves=np.loadtxt("rvm_wavelengths_trimmed.csv",delimiter=",")
dates=np.loadtxt("rvm_dates.csv",delimiter=",")
fs=flux[:,0]
efs=errflux[:,0]
PlotPretty.pp('white')

halpha=656.3*10.0
hbeta=486.1*10.0
hgamma=434.0*10.0
Heii=4686.0
Hei=5876.0
ymin=-0.5
ymax=2.0
plt.figure()
plt.minorticks_on()
plt.ylim(ymin,ymax)
plt.plot(waves,fs,label=dates[0])
lab=r"$HJD-$"+str(int(np.around(245000+dates[0])))
plt.annotate(lab,xy=(5030,ymax-0.05-(1.0*0.1)),color='b')
#plt.plot([hbeta,hbeta],[ymin,ymax],'r--')
#plt.annotate(r'$H\beta$', xy=(hbeta+3.0, 3.0))



plt.xlabel(r"$Wavelength$ $(\AA)$")
plt.ylabel(r"$Flux$")
plt.title(r"$Broad$ $Line$ $Spectra$")

fs=flux[:,50]
lab=r"$HJD-$"+str(int(np.around(245000+dates[50])))
plt.plot(waves,fs, 'y',label=np.around(dates[50]))
plt.annotate(lab,xy=(5030,ymax-0.05-(2.0*0.1)),color='y')


fs=flux[:,60]
plt.plot(waves,fs, 'g',label=np.around(dates[60]))
lab=r"$HJD-$"+str(int(np.around(245000+dates[60])))
plt.annotate(lab,xy=(5030,ymax-0.05-(3.0*0.1)),color='g')

fs=flux[:,80]
plt.plot(waves,fs, 'r',label=dates[80])
lab=r"$HJD-$"+str(int(np.around(245000+dates[80])))
plt.annotate(lab,xy=(5030,ymax-0.05-(4.0*0.1)),color='r')

plt.savefig("hbetaD.jpg")
plt.show()
