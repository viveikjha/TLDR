import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
#import PlotPretty
from numpy import loadtxt

Flux = loadtxt("rvm_fluxes.csv",delimiter=",")
ErrFlux = loadtxt("rvm_errfluxes.csv",delimiter=",")
SNR = np.ravel(Flux / ErrFlux)


print "Max Flux: ", np.max(Flux)
print "Max Flux Error: ", np.max(ErrFlux)
print "Max SNR: ", np.max(Flux/ErrFlux)
print "------------------------"
print "Median Flux: ", np.median(Flux)
print "Median Flux Error: ", np.median(ErrFlux)
print "Median SNR: ", np.median(Flux/ErrFlux)
print "------------------------"
print "Mean Flux: ", np.mean(Flux)
print "Mean Flux Error: ", np.mean(ErrFlux)
print "Mean SNR: ", np.mean(Flux/ErrFlux)
#SNR Distribution
#S=np.around(SNR,decimals=0)
#histogram = np.histogram(S,bins=(np.max(S)-np.min(S)))
#h = histogram[0]
#b = histogram[1]
#b2=b[0:(np.size(b)-1)]
#print np.shape(h)
#print np.shape(b2)

#plt.figure(figsize=(10,5))
#plt.bar(b2,np.log10(h),width=1.0)
#plt.ylabel(r"$log_{10}$(Number of Pixels)")
#plt.xlabel(r"Signal to Noise")
#plt.xlim([-9,34])
#plt.show()
