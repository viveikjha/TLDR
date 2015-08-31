import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline as UVS



def Model(X,ICF):
	L = len(ICF)
	MF=np.zeros((L),dtype='float')
	for i in range(0,L):
		MF[i]=np.sum(X*ICF[i,:])
	print 'Model: ',MF	
	return MF


def convolve(tdf_values,data_flux):
		val = len(data_flux)
		model_flux = np.empty(val)
		for idate in range(0,val):
			model_flux[idate] = np.sum(tdf_values*data_flux[idate])
		return model_flux

#TDF function from Blandford and McKee
def RAN(t,r):
	temp = np.empty([len(t)])
	for i in range(0,len(t)):
		if t[i] > 0 and t[i] <= (2.0*r):
			temp[i] = t[i]/(2.0*r)+np.log(2.0)
		elif t[i] > (2.0*r) and t[i] <= (4.0*r):
			temp[i] = 2.0 - t[i]/(2.0*r)+np.log(4.0*r/t[i])
		else:
			temp[i]=0
		#print 't: ', t[i], 'psi: ', temp[i]
	return 1.0/(2.0*r)*temp

'''IMPORT DATA FILES'''
#inputfile = '../Ver1/tdfL.txt'
#tdfL = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
#tdfL=tdfL[330]

inputfile = '../data/rvm_fluxes.csv'
rvm_spectra_flux = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print 'Flux Array: ',rvm_spectra_flux.shape
spectra_scale=1e0
rvm_spectra_flux=rvm_spectra_flux*spectra_scale
#print rvm_spectra_flux.shape

inputfile = '../data/rvm_fluxes.csv'
data_flux = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print('Flux Array: ',data_errflux.shape)
#nlines = len(data_errflux) # length first axis

inputfile = '../data/rvm_errfluxes.csv'
data_errflux = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print('Error Flux Array: ',data_errflux.shape)
#nlines = len(data_errflux) # length first axis

inputfile = '../data/rvm_wavelengths.csv'
rvm_spectra_lam = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print 'lam: ',len(rvm_spectra_lam)

inputfile = '../data/rvm_dates.csv'
rvm_spectra_date = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
ndataperline = len(rvm_spectra_date)
#print 'date: ',len(rvm_spectra_date)

inputfile = '../data/arp151.b.dat'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
continuum_date = data[:,0]
continuum_scale = 1.0e3
continuum_flux = data[:,1]*continuum_scale


inputfile = '../data/arp151.b.dat'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
continuum_date = data[:,0]
continuum_scale = 1.0
continuum_flux = data[:,1]*continuum_scale
av = np.mean(continuum_flux)
print 'average: ', av
print continuum_flux
ax = np.max(continuum_flux)
continuum_errflux = data[:,2]

print('cont: ', data.shape)
'''DONE IMPORTING DATA FILES'''


'''----------------------------------------------------'''
'''===================================================='''
'''----------------------------------------------------'''
tau = 50

c = 1.0
r = 15.0
tau = np.linspace(0,50,49)
Psi = RAN(tau,r)	#TDF from Blanford & McKee
print tau
#plt.figure()
#plt.plot(tau,Psi)
#plt.show()



'''COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS'''
pts_required=len(rvm_spectra_date)*len(tau)
interpolation_points = np.empty([len(rvm_spectra_date),len(tau)])
q=0
for idate in range(0,len(rvm_spectra_date)):
	for jdelay in range(0,len(tau)):
		interpolation_points[idate][jdelay] = rvm_spectra_date[idate]-tau[jdelay]
#f= interp1d(continuum_date,continuum_flux)					#Creates interpolation function
f= UVS(continuum_date,continuum_flux,k=2)					#Creates interpolation function

ICF = f(np.ravel(interpolation_points))			#Interpolates continuum over points
#print len(ICF)
#print len(rvm_spectra_date)*len(tau)
ICF=np.reshape(ICF,(len(rvm_spectra_date),len(tau)))
int_continuum_dates = interpolation_points					#only to make a matching set

ICF=ICF/1.0e6
print ICF.shape






#INTERPOLATION!!!
#f= interp1d(continuum_date,continuum_flux)  
#f = interp1d(data[0,:],data[1,:],bounds_error=False,fill_value=data[65,2])
#ICF = f(data_flux_dates)
#ICF[impulse_position]=impulse_value
print 'ICF: ', ICF.shape

L = Model(Psi,ICF)
LR=L
percent_err = 5.0
for j in range(0,len(data_errflux)):
	for h in range(0,len(L)):
		data_flux[j][h]=np.random.normal(L[h],0.02*L[h])
	#data_errflux[j] = percent_err/100*L

print data.shape

TDF = np.array([tau,Psi])
#print TDF
np.savetxt('rvm_flxS.csv',data_flux,delimiter=',')
print "Wrote: rvm_flxS.csv"
np.savetxt('rvm_flx_errS.csv',data_errflux,delimiter=',')	
print "Wrote: rvm_flx_errS.csv"
np.savetxt('continuumS.csv',data,delimiter=',')
print "Wrote: continuumS.csv"
np.savetxt('PsiS.txt',TDF,delimiter=',')

plt.figure(1)
#plt.plot(data_flux_dates,L)
plt.plot(data_flux[0])
plt.plot(L,'r--')
#plt.figure()
#plt.plot(data[:,1])
plt.show()


#print data[:,1]
