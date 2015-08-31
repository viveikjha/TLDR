import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def Model(Psi,C):
	L = np.empty([len(C[:,0])])
	for date in range(0,len(C[0])):
		L[date] = np.sum(Psi*C[date,:])
	return L			


def convolve(tdf_values,data_flux):
		val = len(data_flux)
		model_flux = np.empty(val)
		for idate in range(0,val):
			model_flux[idate] = np.sum(tdf_values*data_flux[idate])
		return model_flux

#TDF function from Blandford and McKee
def RAN(t,r,s):
	temp = np.empty([len(t)])
	for i in range(0,len(t)):
		if t[i] > 0 and t[i] <= (2.0*r):
			temp[i] = t[i]/(2.0*r)+np.log(2.0)
		elif t[i] > (2.0*r) and t[i] <= (4.0*r):
			temp[i] = 2.0 - t[i]/(2.0*r)+np.log(4.0*r/t[i])
		else:
			temp[i]=0
	return (1.0/(2.0*r)*temp)*s
	
	'''IMPORT DATA FILES'''
inputfile = '../data/rvm_fluxes.csv'
data_flux = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print('Flux Array: ',data_flux.shape)
nlines = len(data_flux) # length first axis

rvm_flux = data_flux[0,:]
print "spectra flux: ",rvm_flux.shape

inputfile = '../data/rvm_errfluxes.csv'
data_errflux = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print('Error Flux Array: ',data_errflux.shape)

rvm_errflux = data_errflux[0,:]
print "spectra error: ",rvm_errflux.shape

inputfile = '../data/rvm_wavelengths.csv'
rvm_spectra_lam = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
lam = rvm_spectra_lam[0]

inputfile = '../data/rvm_dates.csv'
rvm_dates = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
ndataperline = len(rvm_dates)

inputfile = '../data/arp151.b.dat'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
continuum_date = data[:,0]
continuum_scale = 1.0
continuum_flux = data[:,1]*continuum_scale
continuum_errflux = data[:,2]

print "continuum flux: ",continuum_flux.shape


'''SETTING UP TIME DELAY FUNCTION'''
ntdftimes = 50         #Number of points in the time delay function
tdf_values = np.zeros(ntdftimes)
tdf_times  = np.linspace(0 , ntdftimes, ntdftimes)        #Time delay times

'''DONE SETTING UP TIME DELAY FUNCTION'''

'''COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS'''
interpolation_points = np.empty([ndataperline,ntdftimes])
for idate in range(0,ndataperline):
        for idelay in range(0,ntdftimes):
                       interpolation_points[idate,idelay] = rvm_dates[idate]-tdf_times[idelay]
                       #print interpolation_points[idate,idelay]
f= interp1d(continuum_date,continuum_flux,kind='linear',bounds_error=False,fill_value=continuum_flux[continuum_flux.size-1])                                                    #Creates interpolation function
#f= UVS(continuum_date,continuum_flux,k=2)                                                      #Creates interpolation function
int_continuum_flux = f(np.ravel(interpolation_points))						#Interpolates continuum over points
int_continuum_flux  = np.reshape(int_continuum_flux,(ndataperline,ntdftimes))
int_continuum_dates = interpolation_points                                                      #only to make a matching set

icf = np.empty([ndataperline,ntdftimes])
for date in range(0,ndataperline):
	icf[date]=f(interpolation_points[date,:])

g= interp1d(continuum_date,continuum_errflux,bounds_error=False,fill_value=continuum_errflux[continuum_errflux.size-1])                                                      #Creates interpolation function
int_continuum_errflux = g(np.ravel(interpolation_points))						#Interpolates continuum over points
int_continuum_errflux  = np.reshape(int_continuum_errflux,(ndataperline,ntdftimes))

print "interpolated continuum: ",int_continuum_flux.shape

PsiSyn = RAN(tdf_times,15.0,10.0)

fig,(ax1,ax2,ax3) = plt.subplots(3)
fig.subplots_adjust(hspace=0.5)
ax1.set_title("Continuum")
ax1.plot(continuum_date,continuum_flux)
ax1.plot(continuum_date,continuum_flux,'b*')
#ax1.plot(int_continuum_dates[20],int_continuum_flux[20],'r*')
#ax1.plot(int_continuum_dates[80],int_continuum_flux[80],'ko')

ax2.set_title("Time Delay Function")
ax2.plot(tdf_times,PsiSyn)


synth = Model(PsiSyn,int_continuum_flux)
ax3.set_title("Spectral Line")
ax3.plot(synth)

plt.show()
