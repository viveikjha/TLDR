import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.fftpack import fft,fftshift
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline as UVS
from scipy.integrate import simps
from scipy.optimize import fmin_l_bfgs_b as minimize
import scipy
import time

def ell2norm(X):
	X = np.square(X)
	norm = np.sqrt(np.sum(X))
	#print norm
	return norm
	
def  ell1norm(PsiB,U,rho,alpha):
	#Psib and U must be for t-1
	PsiA_s = PsiB-(U/rho)
	PsiA = np.zeros(len(PsiA_s))
	for f in range(0,len(PsiA_s)):
		if PsiA_s[f] > alpha:
			PsiA[f] = PsiA_s[f]
	return PsiA
	
		

'''IMPORT DATA FILES'''
inputfile = 'rvm_flxS.txt'
rvm_spectra_flux = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
#print 'Flux Array: ',rvm_spectra_flux.shape
rvm_spectra_flux=rvm_spectra_flux
#print rvm_spectra_flux.shape

inputfile = 'rvm_flxS.txt'
data_flux = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
#print('Flux Array: ',data_errflux.shape)
#nlines = len(data_errflux) # length first axis

inputfile = 'rvm_flx_errS.txt'
data_errflux = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
#print('Error Flux Array: ',data_errflux.shape)
#nlines = len(data_errflux) # length first axis

inputfile = '../data/rvm_wavelengths.csv'
rvm_spectra_lam = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print 'lam: ',len(rvm_spectra_lam)

inputfile = '../data/rvm_dates.csv'
rvm_spectra_date = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)
#print 'date: ',len(rvm_spectra_date)

inputfile = 'continuumS.txt'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
continuum_date = data[:,0]
#continuum_flux = data[:,1]

continuum_scale = 1.0e3
continuum_flux = data[:,1]*continuum_scale

continuum_errflux = data[:,2]
#print 'cont: ', data.shape
print 'conterr: ', continuum_errflux.shape

'''DONE IMPORTING DATA FILES'''

'''SETTING UP TIME DELAY FUNCTION'''
tdfpoints = 50 	#Number of points in the time delay function
tdfvalue = 0.0	#Initial Value of the timedelay function
tdf_values = np.zeros([tdfpoints])
tdf_values.fill(tdfvalue) 
tdf_times=np.linspace(0,20,tdfpoints)	#Time delay times
'''DONE SETTING UP TIME DELAY FUNCTION'''

'''COMPUTING THE CONTINUUM FUNCTION FOR REQUIRED POINTS'''
pts_required=len(rvm_spectra_date)*len(tdf_times)
interpolation_points = np.empty([len(rvm_spectra_date),len(tdf_times)])
q=0
for idate in range(0,len(rvm_spectra_date)):
	for jdelay in range(0,len(tdf_times)):
		interpolation_points[idate][jdelay] = rvm_spectra_date[idate]-tdf_times[jdelay]
#f= interp1d(continuum_date,continuum_flux)					#Creates interpolation function
f= UVS(continuum_date,continuum_flux,k=2)					#Creates interpolation function

int_continuum_flux = f(np.ravel(interpolation_points))			#Interpolates continuum over points
#print len(int_continuum_flux)
#print len(rvm_spectra_date)*len(tdf_times)
int_continuum_flux=np.reshape(int_continuum_flux,(len(rvm_spectra_date),len(tdf_times)))
int_continuum_dates = interpolation_points					#only to make a matching set

def chi2(tik,tdf_values,data_flux, data_errflux):
	ndataperline = len(data_flux)
	model_flux = np.empty(ndataperline)
	for idate in range(0,ndataperline):
		model_flux[idate] = np.sum(tdf_values*int_continuum_flux[idate,:])
		chi2 = np.sum( ((model_flux-data_flux)/data_errflux)**2)
	reg =  np.sum( (tdf_values[1:tdfpoints-1] -tdf_values[0:tdfpoints-2])**2 ) #L2^2
	#print len(model_flux),len(data_errflux)
	#return (chi2+tik*reg)
	return (chi2)


'''DONE COMPUTING CONTINUUM FUNCTION FOR REQUIRED POINTS'''
'''INITIALIZING ADMM PARAMETERS'''
nits = 100 #NUMBER OF ADMM ITERATIONS
initial = 0.00 #INITIAL VALUE FOR PSIs

#REGULARIZING WEIGHTS!!!!
rho=1.0e0
mu = 0.0 #PRIOR WEIGHT
#tik = rho/2
#alpha = mu/rho	#REGULARIZATION WEIGHT

alpha = 1.0e3
tik = 1.0e3
#CONVERGENCE PARAMETERS
gamma = 1.25
epsA = 0.0
epsR =  1.0e-3



TDFA = np.empty((len(rvm_spectra_lam),len(tdf_times)))
TDFB = TDFA
TTDF = np.empty((len(rvm_spectra_lam),tdfpoints))

starttime = time.time()

'''Precomputing Matrices for Tikhonov in Step 2.'''
continuum=int_continuum_flux
n_tdf_pts=tdfpoints
n_line_pts=len(rvm_spectra_date)

#Build mapping matrix
CMAT = np.zeros((n_line_pts,n_tdf_pts))
for lns in range(0,n_line_pts):
	for tdf in range(0,n_tdf_pts):
		if lns >= tdf:
			CMAT[lns][tdf]=continuum[lns,tdf]

#Building Covariance matrix

P = np.empty((len(rvm_spectra_lam),len(data_errflux[0]),len(data_errflux[0])))
print data_errflux.shape
for lam in range(0,len(rvm_spectra_lam)):
	T=np.identity(len(data_errflux[0]))
	for i in range(0,len(data_errflux[0])):
		for j in range(0,len(data_errflux[0])):
			if i == j:
				T[i][j] = data_errflux[lam,i]
	P[lam]=T
print 'shape of P: ',P.shape				
						
						
						
CMATT=np.transpose(CMAT)
M=np.dot(np.dot(CMATT,P),CMAT)

Gamma = np.identity(n_tdf_pts)
GammaT = np.transpose(Gamma)
Q = np.dot(GammaT,Gamma)

#Compute the possible np.dot(CMAT,L)
a='shape of C^T: ',CMATT.shape
b = 'shape of C: ',CMAT.shape
#print b
#print a

#CALCULATE C^T P C
B = np.empty((len(rvm_spectra_flux),tdfpoints,tdfpoints))
G = np.empty((len(rvm_spectra_flux),tdfpoints))
for lam in range(0,len(rvm_spectra_flux)):
	L = rvm_spectra_flux[lam]
	A = np.dot(CMATT,P[lam])
	B[lam] = np.dot(A,CMAT)
	G[lam] = np.dot(A,L)
	
#print 'shape of L: ',len(L)
#print 'shape of G: ',G.shape
print 'Done precomputing matrices'	
	

#Initialize matrices for ADMM:
value = 1.0
U = np.zeros((nits+1,len(rvm_spectra_flux),tdfpoints))
U[0].fill(value)
PsiA = np.zeros((nits+1,len(rvm_spectra_flux),tdfpoints))
PsiA_s = np.zeros((nits+1,len(rvm_spectra_flux),tdfpoints))
PsiB = np.zeros((nits+1,len(rvm_spectra_flux),tdfpoints))
PsiB_s = np.zeros((nits+1,len(rvm_spectra_flux),tdfpoints))

TDFARRA = np.zeros((nits+1,len(rvm_spectra_flux),tdfpoints))
TDFARRB = TDFARRA

'''ADMM LOOP'''

t=0
PsiA[0][:].fill(initial)
PsiA_s[0][:].fill(initial)
PsiB[0][:].fill(initial)
PsiB_s[0][:].fill(initial)
U[0][:].fill(value)

for iter in range(1,nits+1):
	#Computing some Tikhonov Matrices:
	C2 = tik*Q

	
	#for i in range(0,len(rvm_spectra_flux)):
	for i in range(0,1):
	
		#Step 1 Proximal Solution
		#L1 norm proximal Operator
		#PsiA[iter][i]=ell1norm(PsiB[iter-1][i],U[iter-1][i],rho,alpha)
		
		#L2 norm proximal Operator
		for f in range(0,len(PsiA[iter][i])):
			if PsiA[iter][i][f] > 0.0:
				PsiA[iter][i][f]=PsiA_s[iter][i][f]/(2.0*alpha+1.0)
			else: 
				PsiA[iter][i][f]=0.0
	
		#Step 2 Tikhonov Regularization
		PsiB_s[iter][i] = PsiA[iter][i] + U[iter-1][i]/rho
		
		
		#PsiB[iter][i] = np.dot(np.linalg.inv(B[i]+tik*Q),G[i]+np.dot(tik*Q,PsiB_s[iter][i]))
		
		C1 = B[i]
		#print 'C1:', C1.shape
		#C2 = tik*Q
		#print 'C2:', C2.shape
		C3 = G[i]
		#print 'C3:', C3.shape
		C4 = np.dot(C2,PsiB_s[iter][i])
		#print 'C4:', C4.shape
		X = np.dot(np.linalg.inv(C1+C2),(C3+C4))
		#print 'result shape: ',X.shape
		PsiB[iter][i]=X
		#Step 3 Update Multipliers
		U[iter][i] = U[iter-1][i]+rho*(PsiA[iter][i]-PsiB[iter][i])
		
		TDFARRA[iter][i]=PsiA[iter][i]
		TDFARRA[iter][i]=PsiB[iter][i]
	chilam=50
	tdf_values = TDFARRA[iter][chilam]
	print 'chi2: ', chi2(alpha,tdf_values,data_flux[chilam], data_errflux[chilam]), 'iteration: ', iter,' of ', nits
	
#print (time.time()-starttime)
#np.savetxt('TDFA.txt',TDFARRA[nits],delimiter=' ')	
#np.savetxt('TDFB.txt',TDFARRB[nits],delimiter=' ')

#fig, ax = plt.subplots(figsize=(12,8))
#ax.imshow(store, cmap='Reds', extent=[0,20,nlines-1, 0])
#ax.imshow(np.rot90(np.sqrt(TDFARRA[nits]),1), cmap='Reds',extent=[4000,7000,0, 3000])
#plt.title('Velocity Delay Map')
#ax.set_aspect(ntdftimes/(4.0*nlines))	
#plt.show()

inputfile = 'PsiS.txt'
Psi_act = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
#print res.x

plt.plot(Psi_act[1],'r')
plt.plot(PsiB[nits][0],'b')
plt.plot(PsiA[nits][0],'g')

plt.show()