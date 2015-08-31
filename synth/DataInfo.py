import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import UnivariateSpline as UVS
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import minimize

'''IMPORT DATA FILES'''
inputfile = '../data/rvm_fluxes.csv'
data_flux = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
print 'rvm_flux data: ',data_flux.shape

#print('Flux Array: ',data_flux.shape)
#nlines = len(data_flux) # length first axis

inputfile = '../data/rvm_errfluxes.csv'
data_errflux = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
print 'rvm_errflux data: ', data_errflux.shape
#print('Error Flux Array: ',data_errflux.shape)
#nlines = len(data_errflux) # length first axis

inputfile = '../data/rvm_wavelengths.csv'
rvm_spectra_lam = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
print 'rvm_spectra_lam data: ', rvm_spectra_lam.shape

#print('lam: ',len(rvm_spectra_lam))
#nlines2 = len(rvm_spectra_lam)

#if (nlines != nlines2):

#   print("The number of lines is inconsistent")

inputfile = '../data/rvm_dates.csv'
data_flux_dates = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
ndataperline = len(data_flux_dates)
print 'data_flux_dates data: ', data_flux_dates.shape

#print('date: ',ndataperline)


inputfile = '../data/arp151.b.dat'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
print 'data data: ', data.shape

continuum_date = data[:,0]
continuum_scale = 1.0e3
continuum_flux = data[:,1]*continuum_scale
continuum_errflux = data[:,2]
#print('cont: ', data.shape)
'''DONE IMPORTING DATA FILES'''



