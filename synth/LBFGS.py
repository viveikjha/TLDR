import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import UnivariateSpline as UVS
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import minimize

'''IMPORT DATA FILES'''
inputfile = 'rvm_flxS.txt'
data_flux = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
spectral_scale = 1.0e0
data_flux=spectral_scale*data_flux
print('Flux Array: ',data_flux.shape)
print data_flux
nlines = len(data_flux) # length first axis

inputfile = 'rvm_flx_errS.txt'
data_errflux = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
data_errflux=spectral_scale*data_errflux
#data_errflux = 5.0/100*data_flux
print('Error Flux Array: ',data_errflux.shape)
print data_errflux
nlines = len(data_errflux) # length first axis

inputfile = '../data/rvm_wavelengths.csv'
rvm_spectra_lam = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
rvm_spectra_lam = np.array([1])
print('lam: ',len(rvm_spectra_lam))           
nlines2 = len(rvm_spectra_lam)
#if (nlines != nlines2):
#   print("The number of lines is inconsistent")

inputfile = '../data/rvm_dates.csv'
data_flux_dates = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
ndataperline = len(data_flux_dates)
print('date: ',ndataperline)


inputfile = 'continuumS.txt'
#continuum_scale = 1.0
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
continuum_date = data[:,0]


continuum_scale = 1.0

alpha = 5.0e2

continuum_flux = data[:,1]*continuum_scale
print 'continuum mean ', np.mean(continuum_flux)

continuum_errflux = data[:,2]*continuum_scale
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

#f= interp1d(continuum_date,continuum_flux)                                                     #Creates interpolation function
f= UVS(continuum_date,continuum_flux,k=2)                                                       #Creates interpolation function
int_continuum_flux = f(np.ravel(interpolation_points))                               			#Interpolates continuum over points
int_continuum_flux  = np.reshape(int_continuum_flux,(ndataperline,ntdftimes))
int_continuum_dates = interpolation_points                                                      #only to make a matching set
#print('Continuum points\n', int_continuum_flux)



store = np.zeros([nlines,ntdftimes])
chi2 = np.zeros(nlines)


def model(tdf_values,data_flux,data_errflux):
		model_flux = np.empty(ndataperline)
		for idate in range(0,ndataperline):
			model_flux[idate] = np.sum(tdf_values*int_continuum_flux[idate,:])
		return model_flux

def fct(tdf_values,data_flux, data_errflux):
        model_flux = np.empty(ndataperline)
        for idate in range(0,ndataperline):
            model_flux[idate] = np.sum(tdf_values*int_continuum_flux[idate,:])
        chi2 = np.sum( ((model_flux-data_flux)/data_errflux)**2)
        #print len(model_flux),len(data_errflux)
        #reg =  np.sum( np.abs( tdf_values[1:ntdftimes-1] -tdf_values[0:ntdftimes-2] ) )#TV
        reg =  np.sum( (tdf_values[1:ntdftimes-1] -tdf_values[0:ntdftimes-2])**2 ) #L2^2
#        print(chi2, alpha*reg)
        return (chi2 + alpha*reg)
        #return(chi2)

def grad(tdf_values, data_flux, data_errflux):
        model_flux = np.empty(ndataperline)
        for idate in range(0,ndataperline):
            model_flux[idate] = np.sum(tdf_values*int_continuum_flux[idate,:])

        chi2grad = np.zeros(ntdftimes)
        reggrad = np.zeros(ntdftimes)
        for idelay in range(0,ntdftimes):
            chi2grad[idelay] = 2.0*np.sum( (model_flux - data_flux)/data_errflux**2 * int_continuum_flux[:,idelay])
            if((idelay > 0) and (idelay < (ntdftimes-1))):
                #reggrad[idelay] = 4.*tdf_values[idelay]-2.*tdf_values[idelay-1]-2.*tdf_values[idelay+1]
                reggrad[idelay] = np.sign(tdf_values[idelay+1]-tdf_values[idelay])+np.sign(tdf_values[idelay]-tdf_values[idelay-1])
        return (chi2grad+alpha*reggrad)
        #return (chi2grad)
        
        
#def hes(tdf_values,data_flux,data_errflux):
#		return(Hessian)        
        
        #Done Gradient
for iline in range(0, 1):
        flx = data_flux
        errflx = data_errflux
        bnd = ntdftimes*[(0.0,None)]

        #res=minimize(fct, tdf_values, args=[flx,errflx], method='TNC', jac=grad, bounds=bnd,options={'maxiter': 200})
        res=minimize(fct, tdf_values, args=[flx,errflx], method='TNC', jac=grad, bounds=bnd)
        #res=minimize(fct, tdf_values, args=[flx,errflx], method='L-BFGS-B', jac=grad, bounds=bnd)
        
        store[iline] = res.x
        chi2[iline] = fct(res.x,flx,errflx)/ndataperline
        print'iline = ', iline, '/', nlines, 'chi2r = ',chi2[iline],' status = ',res.success 

#np.savetxt('tdf.txt',store,delimiter=' ',newline='\n')
#np.savetxt('chi2.txt',chi2,delimiter=' ',newline='\n')

inputfile = 'PsiS.txt'
Psi_act = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)

#print Psi_act[1]
#print res.x

real = model(Psi_act,data_flux,data_errflux)
chi2 = np.sum( ((real-data_flux)/data_errflux)**2)
#print 'real: ', chi2
plt.plot(Psi_act[1],'r')
plt.plot(store[0],'k*')



plt.show()
