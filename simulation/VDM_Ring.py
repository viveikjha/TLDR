import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty

PlotPretty.pp('white')


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)
weight = 1.0 #LINE WEIGHT!
r_o = 5.0
r_i = 1.0
r = 0.5*(r_o-r_i)+r_i
stdev=0.50
center = [0.0,0.0]

nps = 10000

#radii = np.random.normal(r,stdev,nps)
radii = np.random.uniform(r_i,r_o,nps)
angle = np.random.uniform(0.0,2.0*math.pi,nps)

X=np.zeros(nps)
Y=np.zeros(nps)
for i in range(0,nps):
	X[i],Y[i]=pol2cart(radii[i],angle[i])
delay =  radii+(radii*np.cos(angle))

G = 5.702*10.0**-11.0 #light days (solar masses)^-1 (speed of light)^2
M = 7.0e6 # solar masses


sol = 299792.458 #km/s
velocity =np.sqrt(G*M/radii)*sol
LOSvel = velocity*np.cos((math.pi/2.0)-angle)


fig = plt.figure(figsize=(4,12))

#ax = fig.add_subplot(121)
ax = plt.subplot2grid((3,1),(0,0),colspan=1)
ax.set_title("BLR Geometry")
col=LOSvel
sizes = (delay/max(delay))*100.0
ax.scatter(X,Y,c=-col,s=sizes,cmap="bwr",vmin=min(LOSvel),vmax=max(LOSvel))
plt.minorticks_on()
#ax.set_aspect('equal')


c = 1.0



#CIRCLE AT RADIUS
theta = np.linspace(0.0,2.0*math.pi,200)
x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'c',linewidth=2.0)

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

#ISODELAY CURVES
############################################
t = 0.2
theta = np.linspace(-0.935*math.pi,0.935*math.pi,200)
r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'c--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.55, y[199]-0.05))
############################################
############################################
t = 1.0
theta = np.linspace(-0.85*math.pi,0.85*math.pi,200)
r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'c--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################
############################################
t = 4.0
theta = np.linspace(-0.67*math.pi,0.67*math.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'c--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################

############################################
t = 10.0
theta = np.linspace(-0.43*math.pi,0.43*math.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'c--',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################

ax.set_xlim(xmin,xmax)
ax.set_ylim(-10.0,10.0)


#ax=fig.add_subplot(122)
ax = plt.subplot2grid((3,1),(1,0),colspan=1)
ax.set_title("Velocity Delay Map")
ax.set_xlabel("L.O.S. Velocity (km/s)")
ax.set_ylabel("Delay")
ax.scatter(LOSvel,delay,c=-col,cmap="bwr",vmin=min(LOSvel),vmax=max(LOSvel))
plt.minorticks_on()

###########################################
#RELATIVE SIGNAL STRENTGTH
IFE = 1.0

RF = IFE/1.0
Emissivity = delay/max(delay)
RF = Emissivity*IFE #ADD EMISSIVITY AS FUNCTION OF MAXIMUM DELAY.

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
				Fim[m,n] = Fim[m,n]+RF[j]
				pts = pts+1
print c
print "Max Signal: ", np.max(Fim)

vc = np.empty([pixlen])
for n in range(0,pixlen):
	vc[n]=(pixbinV[n+1]-pixbinV[n])/2.0+pixbinV[n] #Bin Centers


ax = plt.subplot2grid((3,1),(2,0),colspan=1)
ax.imshow((Fim),cmap='Reds',origin='lower',extent=[min(LOSvel),max(LOSvel),0,12],aspect='auto')
plt.minorticks_on()
ax.set_title("Binned Relative Signal Strength")
ax.set_ylabel("Delay")
ax.set_xlabel("Binned L.O.S. Velocity (km/s)" )



plt.tight_layout()
plt.savefig('disk.png', format='png')
#plt.show()

Ha =6563.0
sol = 299792.458 #km/s
#H-Alpha:
DL_a = Ha*LOSvel/sol
newl_a = Ha+(vc/sol*Ha)			#Wavelengths of the Bin Centers
print("!!: ",len(newl_a))
wavelengths = newl_a
print(newl_a)


print "min: ", min(newl_a)," max: ", max(newl_a)
plt.figure()
plt.imshow(np.log(Fim), origin="lower", extent=[min(newl_a),max(newl_a),0,12],aspect="auto",interpolation="None")
plt.show()


np.savetxt('Ring/LOSvel.csv',LOSvel,delimiter=',')
np.savetxt('Ring/Fim.csv',(Fim/np.max(Fim)),delimiter=',')
np.savetxt('Ring/X.csv',X,delimiter=',')
np.savetxt('Ring/Y.csv',Y,delimiter=',')
np.savetxt('Ring/radii.csv',radii,delimiter=',')
np.savetxt('Ring/delay.csv',delay,delimiter=',')
np.savetxt('Ring/angle.csv',angle,delimiter=',')
np.savetxt('Ring/new_wavelengths.csv',newl_a,delimiter=',')
np.savetxt('Ring/simulated_vdm.csv',Fim,delimiter=',')
