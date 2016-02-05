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
		r = np.append(r,spiral(phi[i],1.0,0.5,4.0))
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

fig = plt.figure(figsize=(14,6))

#ax = fig.add_subplot(121)
ax = plt.subplot2grid((1,2),(0,0),rowspan=1)
ax.set_title("BLR Geometry")
col=LOSvel
sizes = (delay/max(delay))*100.0
ax.scatter(X,Y,c=col,s=sizes,edgecolor='k',cmap="seismic",vmin=min(LOSvel),vmax=max(LOSvel))


#ISODELAY CURVES
############################################
t = 0.2
theta = np.linspace(-0.935*np.pi,0.935*np.pi,200)
r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'k--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.55, y[199]-0.05))
############################################
############################################
t = 1.0
theta = np.linspace(-0.85*np.pi,0.85*np.pi,200)
r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'w--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################
############################################
t = 4.0
theta = np.linspace(-0.67*np.pi,0.67*np.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'k--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################

############################################
t = 10.0
theta = np.linspace(-0.43*np.pi,0.43*np.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'k--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################
#LABELS
xmin = -15.0
xmax = 8.0
npts = 200
x = np.linspace(xmin,xmax,npts)
y = np.zeros([npts])
plt.plot(x,y,'k',linewidth=weight)
dirstr = r"$\Longleftarrow Observer$"
ax.annotate(dirstr, xy=(xmin, 0.02))

plt.plot([0],[0],'ko',markersize=6)
ax.set_xlim(xmin,xmax)





#ax=fig.add_subplot(122)
ax = plt.subplot2grid((1,2),(0,1),rowspan=2)
ax.set_title("Velocity Delay Map")
ax.set_xlabel("L.O.S. Velocity (c)")
ax.set_ylabel("Delay")
#ax.set_xlim(-0.5,0.5)
ax.scatter(LOSvel,delay,c=-col,edgecolor='k',cmap="seismic",vmin=min(LOSvel),vmax=max(LOSvel))


###########################################
#RELATIVE SIGNAL STRENTGTH
IFE = 1.0

RF = IFE/1.0

pixlen =50
Fim = np.zeros([pixlen,pixlen],dtype='float')
pixbinT = np.linspace(min(delay),max(delay),pixlen+1)
pixbinV = np.linspace(min(LOSvel),max(LOSvel),pixlen+1)

wt=0
wv=0
pts = 0.0
c = 0
for j in range(0,len(LOSvel)): #ind
	for m in range(0,pixlen): #V
		for n in range(0,pixlen): #T
			if (LOSvel[j] > pixbinV[n]) and (LOSvel[j] <= pixbinV[n+1]) and (delay[j] > pixbinT[m]) and (delay[j] <=pixbinT[m+1]):
				if m == 0 and n == 0:
					c+=1		
				Fim[m,n] = Fim[m,n]+RF
				pts = pts+1
print c
print "Max Signal: ", np.max(Fim)


		
ax = plt.subplot2grid((3,1),(2,0),colspan=1)
ax.imshow((Fim),cmap='Reds',origin='lower',extent=[min(LOSvel),max(LOSvel),0,12],aspect='auto')
plt.minorticks_on()
ax.set_title("Binned Relative Signal Strength")
ax.set_ylabel("Delay")
ax.set_xlabel("Binned L.O.S. Velocity (km/s)" )



plt.tight_layout()
#plt.savefig('disk.png', format='png')
plt.show()

Ha =6563.0 
sol = 299792.458 #km/s
#H-Alpha:
DL_a = Ha*LOSvel/sol
newl_a = Ha+DL_a 
wavelengths = newl_a

print "min: ", min(newl_a)," max: ", max(newl_a)
plt.figure()
plt.imshow(np.log(Fim), origin="lower", extent=[min(newl_a),max(newl_a),0,12],aspect="auto",interpolation="None")
plt.show()


np.savetxt('LOSvel.csv',LOSvel,delimiter=',')	
np.savetxt('Fim.csv',(Fim/np.max(Fim)),delimiter=',')	
np.savetxt('X.csv',X,delimiter=',')	
np.savetxt('Y.csv',Y,delimiter=',')	
np.savetxt('radii.csv',radii,delimiter=',')	
np.savetxt('delay.csv',delay,delimiter=',')	
np.savetxt('angle.csv',angle,delimiter=',')	
np.savetxt('new_wavelengths.csv',newl_a,delimiter=',')
np.savetxt('simulated_vdm.csv',Fim,delimiter=',')



plt.show()
