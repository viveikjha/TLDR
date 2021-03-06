import numpy as np
import matplotlib.pyplot as plt
lams = np.loadtxt('data/rvm_wavelengths.csv',delimiter=',')

#vdm = np.loadtxt('simulation/simulated_vdm.csv',delimiter=',')
vdm = np.loadtxt('box_20x50/RevMapResult.csv',delimiter=',')


#lineax=np.empty(ndelays);lineax.fill(Ha)
#linebx=np.empty(ndelays);lineax.fill(Hb)
#linegx=np.empty(ndelays);lineax.fill(Hg)




#plt.figure()
#plt.imshow(np.log(vdm),cmap='Reds',aspect='auto',origin='lower',extent=[min(lams),max(lams),0,20.0])
#plt.xlabel(r'Wavelength ($\AA$)')
#plt.ylabel('Delay (Days)')


plt.figure()
plt.imshow((vdm),cmap='Greys',aspect='auto',origin='lower',extent=[min(lams),max(lams),0,20.0])
#plt.plot(lineax,liney,'k')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Delay (Days)')

plt.show()
