import numpy as np
import matplotlib.pyplot as plt
import PlotPretty
PlotPretty.pp('white')
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def spiral(phi,A,B,N):
	#print phi, " ", B*np.tan(phi/(2.0*N))
	r = A*1.0/(np.log(B*np.tan(phi/(2.0*N))))
	return r
weight = 1.0
nps = 500
phi = np.linspace(0.2,8.0*np.pi/4.0,nps)
r = np.array([])
phi2 = np.array([])

for i in range(0,len(phi)):
	try:
		r = np.append(r,spiral(phi[i],5.0,0.5,4.0))
		phi2 = np.append(phi2,phi[i])
	except:
		print phi[i], " failed."

apn = 10
r_n = np.array([])
phi_n = np.array([])
for i in range(0,len(r)):
	new_rs = np.random.normal(np.abs(r[i]),np.abs(r[i])*0.08,apn)
	new_phis = np.random.normal(phi[i],phi[i]*0.01,apn)
	r_n = np.append(r_n,new_rs)
	phi_n = np.append(phi_n,new_phis)

r_n = np.append(r_n,r_n)
phi_n = np.append(phi_n,phi_n+np.pi)+np.pi/4.0

#for i in range(0,len(r_n)):
#	print r_n[i], " ", phi_n[i]

G = 5.702*10.0**-11.0 #light days (solar masses)^-1 (speed of light)^2
M = 1.0e8 # solar masses
c = 1.0
X,Y = pol2cart(r_n,-phi_n)

delay =  r_n+(r_n*np.cos(phi_n))

velocity =np.sqrt(G*M/r_n)
LOSvel = velocity*np.cos((np.pi/2.0)-phi_n)

#fig = plt.figure(figsize=(4,12))

#ax = fig.add_subplot(121)
#ax = plt.subplot2grid((3,1),(0,0),rowspan=1)
#ax.set_title("BLR Geometry")
col=LOSvel
sizes = (delay/max(delay))*100.0
#ax.scatter(X,Y,c=col,s=sizes,edgecolor='k',cmap="seismic",vmin=min(LOSvel),vmax=max(LOSvel))


#ISODELAY CURVES
############################################
t = 0.2
theta = np.linspace(-0.935*np.pi,0.935*np.pi,200)
r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

#plt.plot(x,y,'k--',linewidth=weight)
st = r'$\tau = $'+ str(t)
#plt.annotate(st, xy=(x[199]-0.55, y[199]-0.05))
############################################
############################################
t = 1.0
theta = np.linspace(-0.85*np.pi,0.85*np.pi,200)
r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

#plt.plot(x,y,'w--',linewidth=weight)
st = r'$\tau = $'+ str(t)
#plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################
############################################
t = 4.0
theta = np.linspace(-0.67*np.pi,0.67*np.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

#plt.plot(x,y,'k--',linewidth=weight)
st = r'$\tau = $'+ str(t)
#plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################

############################################
t = 10.0
theta = np.linspace(-0.43*np.pi,0.43*np.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

#plt.plot(x,y,'k--',linewidth=weight)
st = r'$\tau = $'+ str(t)
#plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################
#LABELS
xmin = -15.0
xmax = 8.0
npts = 200
x = np.linspace(xmin,xmax,npts)
y = np.zeros([npts])
#plt.plot(x,y,'k',linewidth=weight)
dirstr = r"$\Longleftarrow Observer$"
#ax.annotate(dirstr, xy=(xmin, 0.02))

#plt.plot([0],[0],'ko',markersize=6)
#ax.set_xlim(xmin,xmax)





#ax=fig.add_subplot(122)
#ax = plt.subplot2grid((3,1),(0,1),rowspan=2)
#ax.set_title("Velocity Delay Map")
#ax.set_xlabel("L.O.S. Velocity (c)")
#ax.set_ylabel("Delay")
#ax.set_xlim(-0.5,0.5)
#ax.scatter(LOSvel,delay,c=-col,edgecolor='k',cmap="seismic",vmin=min(LOSvel),vmax=max(LOSvel))


###########################################
#RELATIVE SIGNAL STRENTGTH
IFE = 1.0
RF = IFE/1.0
Emissivity = delay/max(delay)
RF = Emissivity*IFE #ADD EMISSIVITY AS FUNCTION OF MAXIMUM DELAY.


pixlen =50
pixlen_waves=20
pixlen_times=50
Fim = np.zeros([pixlen_times,pixlen_waves],dtype='float')
pixbinT = np.linspace(min(delay),max(delay),pixlen_times+1)
pixbinV = np.linspace(min(LOSvel),max(LOSvel),pixlen_waves+1)

wt=0
wv=0
pts = 0.0
c = 0
for j in range(0,len(LOSvel)): #ind
	for m in range(0,pixlen_waves): #V
		for n in range(0,pixlen_times): #T
			if (LOSvel[j] > pixbinV[m]) and (LOSvel[j] <= pixbinV[m+1]) and (delay[j] > pixbinT[n]) and (delay[j] <=pixbinT[n+1]):
				if m == 0 and n == 0:
					c+=1
				Fim[n,m] = Fim[n,m]+RF[j]
				pts = pts+1
#print c
#print "Max Signal: ", np.max(Fim)



#ax = plt.subplot2grid((3,1),(2,0),colspan=1)
#ax.imshow((Fim),cmap='Reds',origin='lower',extent=[min(LOSvel),max(LOSvel),0,12],aspect='auto')
#plt.minorticks_on()
#ax.set_title("Binned Relative Signal Strength")
#ax.set_ylabel("Delay")
#ax.set_xlabel("Binned L.O.S. Velocity (km/s)" )



#plt.tight_layout()
#plt.savefig('disk.png', format='png')
#plt.savefig("Spiral/Spiral_Geometry.png")
#plt.show()

Ha =6563.0
sol = 299792.458 #km/s
#H-Alpha:
DL_a = Ha*LOSvel/sol
newl_a = Ha+DL_a
wavelengths = newl_a
print "shape of vdm: ", np.shape(Fim)
#print "min: ", min(newl_a)," max: ", max(newl_a)
plt.figure()
plt.imshow((Fim), origin="lower", aspect="auto",interpolation="None",cmap="Reds")
plt.show()


#np.savetxt('Spiral/LOSvel.csv',LOSvel,delimiter=',')
#np.savetxt('Spiral/Fim.csv',(Fim/np.max(Fim)),delimiter=',')
#np.savetxt('Spiral/X.csv',X,delimiter=',')
#np.savetxt('Spiral/Y.csv',Y,delimiter=',')
#np.savetxt('Spiral/radii.csv',r_n,delimiter=',')
#np.savetxt('Spiral/delay.csv',delay,delimiter=',')
#np.savetxt('Spiral/angle.csv',phi_n,delimiter=',')
#np.savetxt('Spiral/new_wavelengths.csv',newl_a,delimiter=',')
#np.savetxt('Spiral/simulated_vdm.csv',Fim,delimiter=',')

Fim=(Fim/np.max(Fim))*0.05
np.savetxt('../Paper/simulated_spiral_vdm.csv',Fim,delimiter=',')

#plt.show()
