'''THIS FILE CREATES A DECENT HBETA IMAGE FOR ARP151 DATA'''
'''KEEP THIS FILE AS IS FOR REFERENCE FOR ALGORITHIM CHANGES'''
'''-------------------DO NOT ALTER!-------------------------'''

import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as UVS
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import minimize

'''IMPORT DATA FILES'''
inputfile = 'rvm_flxS.txt'
data_flux = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
print('Flux Array: ',data_flux.shape)
nlines = len(data_flux) # length first axis

inputfile = 'rvm_flx_errS.txt'
data_errflux = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
print('Error Flux Array: ',data_errflux.shape)
#nlines = len(data_errflux) # length first axis

inputfile = '../data/rvm_wavelengths.csv'
rvm_spectra_lam = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
print('lam: ',len(rvm_spectra_lam))
nlines2 = len(rvm_spectra_lam)
if (nlines != nlines2):
   print("The number of lines is inconsistent")

inputfile = '../data/rvm_dates.csv'
data_flux_dates = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
ndataperline = len(data_flux_dates)
print('date: ',ndataperline)


inputfile = 'continuumS.txt'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
continuum_date = data[:,0]
continuum_scale = 1.0e3
continuum_flux = data[:,1]*continuum_scale
continuum_errflux = data[:,2]
print('cont: ', data.shape)
'''DONE IMPORTING DATA FILES'''


'''SETTING UP TIME DELAY FUNCTION'''
ntdftimes = 50         #Number of points in the time delay function
tdf_values = np.zeros(ntdftimes)
tdf_times  = np.linspace(0 , 20, ntdftimes)        #Time delay times

'''DONE SETTING UP TIME DELAY FUNCTION'''

'''COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS'''
interpolation_points = np.empty([ndataperline,ntdftimes])
for idate in range(0,ndataperline):
        for idelay in range(0,ntdftimes):
                       interpolation_points[idate,idelay] = data_flux_dates[idate]-tdf_times[idelay]

#f= interp1d(continuum_date,continuum_flux)                                                    #Creates interpolation function
f= UVS(continuum_date,continuum_flux,k=2)                                                      #Creates interpolation function
int_continuum_flux = f(np.ravel(interpolation_points))						#Interpolates continuum over points
int_continuum_flux  = np.reshape(int_continuum_flux,(ndataperline,ntdftimes))/1000000.
int_continuum_dates = interpolation_points                                                      #only to make a matching set

g= UVS(continuum_date,continuum_errflux,k=2)                                                      #Creates interpolation function
int_continuum_errflux = g(np.ravel(interpolation_points))						#Interpolates continuum over points
int_continuum_errflux  = np.reshape(int_continuum_errflux,(ndataperline,ntdftimes))/1000000.


#print('Continuum points\n', int_continuum_flux)
#print(int_continuum_flux.shape)

'''DONE COMPUTING CONTINUUM FUNCTION FOR REQUIRED POINTS'''

'''COMPUTE THE HESSIAN MATRIX'''
'''try:
        inputfile = 'Hessian.npy'
        Hessian = np.load(inputfile)
        print('Hessian Found')
except:
        Hessian = np.empty([nlines,ntdftimes,ntdftimes])
        for nline in range(0,nlines):
                for jouter in range(0,ntdftimes):
	                for inner in range(0,ntdftimes):
                                Hessian[nline,jouter,inner]=np.sum(2.0/data_errflux[nline,:]*int_continuum_flux[:,jouter]*int_continuum_flux[:,inner])
        np.save('Hessian',Hessian)
        print('Hessian saved for future use.')
'''
'''DONE COMPUTING HESSIAN'''


store = np.zeros([nlines,ntdftimes])
chi2 = np.zeros(nlines)
alpha = 500.0
def fct(tdf_values,data_flux, data_errflux):
        model_flux = np.empty(ndataperline)
        for idate in range(0,ndataperline):
            model_flux[idate] = np.sum(tdf_values*int_continuum_flux[idate,:])
        chi2 = np.sum( ((model_flux-data_flux)/data_errflux)**2)
#        reg =  np.sum( np.abs( tdf_values[1:ntdftimes-1] -tdf_values[0:ntdftimes-2] ) )
        reg =  np.sum( (tdf_values[1:ntdftimes-1] -tdf_values[0:ntdftimes-2])**2 )
#        print(chi2, alpha*reg)
        return (chi2 + alpha*reg)

def grad(tdf_values, data_flux, data_errflux):
        model_flux = np.empty(ndataperline)
        for idate in range(0,ndataperline):
            model_flux[idate] = np.sum(tdf_values*int_continuum_flux[idate,:])

        chi2grad = np.zeros(ntdftimes)
        reggrad = np.zeros(ntdftimes)
        for idelay in range(0,ntdftimes):
            chi2grad[idelay] = 2.0*np.sum( (model_flux - data_flux)/data_errflux**2 * int_continuum_flux[:,idelay])
            if((idelay > 0) and (idelay < (ntdftimes-1))):
                 reggrad[idelay] = 4.*tdf_values[idelay]-2.*tdf_values[idelay-1]-2.*tdf_values[idelay+1]
#                reggrad[idelay] = np.sign(tdf_values[idelay+1]-tdf_values[idelay])+np.sign(tdf_values[idelay]-tdf_values[idelay-1])
        return (chi2grad+alpha*reggrad)
        #Done Gradient

for iline in range(0, nlines):
        flx = data_flux[iline,:]
        errflx = data_errflux[iline,:]
        bnd = ntdftimes*[(0.0,None)]
        #H = Hessian[iline]

        #value = np.average(flx)
        #tdf_values.fill(value/50)
        #res=minimize(fct, tdf_values, args=[flx,errflx], method='TNC', jac=grad, bounds=bnd,options={'maxiter': 200})
        #res=minimize(fct, tdf_values, args=[flx,errflx], method='Newton-CG', jac=grad,hess=H, bounds=bnd)
        #res=minimize(fct, tdf_values, args=[flx,errflx], method='trust-ncg', jac=grad,hess=H)
        
        res=minimize(fct, tdf_values, args=[flx,errflx], method='TNC', jac=grad, bounds=bnd)
        
        store[iline,:] = res.x
        chi2[iline] = fct(res.x,flx,errflx)/ndataperline
        print('iline = ', iline, '/', nlines, 'chi2r = ',chi2[iline],' status = ',res.success)

np.savetxt('tdf.txt',store,delimiter=' ',newline='\n')
np.savetxt('chi2.txt',chi2,delimiter=' ',newline='\n')

#fig, ax = plt.subplots(figsize=(12,8))
#ax.imshow(store, cmap='Reds', extent=[0,20,nlines-1, 0])
#ax.imshow(np.rot90(np.sqrt(store),1), cmap='Reds',extent=[0,20,nlines-1, 0])
#plt.title('fit')
#ax.set_aspect(ntdftimes/(4.0*nlines))

#figb,axb = plt.subplots(figsize=(6,4))
#axb.plot(chi2)

#plt.title('chi2')



plt.show()
