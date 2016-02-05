import numpy as np
import matplotlib.pyplot as plt
lams = np.loadtxt('../data/rvm_wavelengths.csv',delimiter='\n')

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

#H-Beta:
DL_b = Hb*LOSvel/sol	#Calculate wavelength shifts about hbeta
newl_b = Hb+DL_b 			#Puts wavelengths for the vdm about hbeta
Hb_wavelengths = newl_b
Hb_delays = delay

plt.figure()
plt.imshow(np.log(Fim),cmap='Reds',aspect='auto',origin='lower',extent=[min(Hb_wavelengths),max(Hb_wavelengths),0,50.0])
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Delay (Days)')

plt.show()

