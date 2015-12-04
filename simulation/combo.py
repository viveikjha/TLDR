import numpy as np
import matplotlib.pyplot as plt
lams = np.loadtxt('../data/rvm_wavelengths.csv',delimiter=',')

LOSvel = np.loadtxt('LOSvel.csv',delimiter=',')	
Fim = np.loadtxt('Fim.csv',delimiter=',')	
X = np.loadtxt('X.csv',delimiter=',')	
Y = np.loadtxt('Y.csv',delimiter=',')	
radii = np.loadtxt('radii.csv',delimiter=',')	
delay = np.loadtxt('delay.csv',delimiter=',')	
angle = np.loadtxt('angle.csv',delimiter=',')	

ndelays = 50
nlams = len(lams)
vdm = np.zeros([ndelays,nlams])

sol = 299792.458 #km/s

#Center Wavelengths
Ha =6563.0    
Hb =4861.0
Hg =4341.0

#Relative Strength (By eye)
Ra = 1.0
Rb = 1.0/3.0
RG = 1.0/6.0

#wavelengths = np.empty([])
delays = delay
#H-Alpha:
DL_a = Ha*LOSvel/sol
newl_a = Ha+DL_a 
wavelengths = newl_a

#H-Beta:
DL_b = Hb*LOSvel/sol
newl_b = Hb+DL_b 
wavelengths = np.append(wavelengths,newl_b)
delays = np.append(delays,delay)

#H-Gamma:
DL_g = Hg*LOSvel/sol
newl_g = Hg+DL_g 
wavelengths = np.append(wavelengths,newl_g)
delays = np.append(delays,delay)

'''Need to put velocities and delays into the vdm grid'''
#Wavelength Bins
lam_min = min(lams)
lam_max = max(lams)

del_min = 0.0
del_max = 20.0
delay_bins = np.linspace(del_min,del_max,ndelays)

init_sig = 1.0

#Binning & Scaling
for k in range(0,len(wavelengths)):
#for k in range(0,3):
	lam_binset=0
	del_binset=0
	delay_bin_index=0
	lam_bin_index = 0
	for l in range(1,nlams-1):
		bin_min = lams[l]-(lams[l]-lams[l-1])/2.0
		bin_max = lams[l]+(lams[l+1]-lams[l])/2.0
		if bin_min < wavelengths[k] <= bin_max:
			#print bin_min, " < ",wavelengths[k], " <= ", bin_max 
			lam_bin_index = l
			lam_binset=1
		elif wavelengths[k] > lams[nlams-1]-(lams[nlams-1]-lams[nlams-2])/2.0:
			lam_bin_index = nlams-1
		#else:
			#print "L: Problem!"
	for t in range(1,ndelays-1):
		#Get delay index
		bin_min = delay_bins[t]-(delay_bins[t]-delay_bins[t-1])/2.0
		bin_max = delay_bins[t]+(delay_bins[t+1]-delay_bins[t])/2.0
		if delays[k] > bin_min and delays[k] <= bin_max:
			delay_bin_index = t
		elif delays[k] <= delay_bins[1] - (delay_bins[1]-delay_bins[0])/2.0:
			delay_bin_index=0
		elif delays[k] >  delay_bins[ndelays-1]-(delay_bins[ndelays-1]-delay_bins[ndelays-2])/2.0:
			delay_bin_index=ndelays-1
	vdm[delay_bin_index,lam_bin_index] = vdm[delay_bin_index,lam_bin_index] + init_sig/(delays[k])**2

np.savetxt('simulated_vdm.csv',vdm,delimiter=',')


#lineax=np.empty(ndelays);lineax.fill(Ha)
#linebx=np.empty(ndelays);lineax.fill(Hb)
#linegx=np.empty(ndelays);lineax.fill(Hg)


liney=delay_bins

plt.figure()
plt.imshow(np.log(vdm),cmap='Reds',aspect='auto',origin='lower',extent=[min(lams),max(lams),0,20.0])
#plt.plot(lineax,liney,'k')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Delay (Days)')


plt.figure()
plt.imshow((vdm),cmap='Reds',aspect='auto',origin='lower',extent=[min(lams),max(lams),0,20.0])
#plt.plot(lineax,liney,'k')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Delay (Days)')

#plt.plot(linebx,liney)
#plt.plot(linegx,liney)
#plt.savefig('sim_vdm.png', format='png')
plt.show()


print "done"



