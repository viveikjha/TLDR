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

r_o = 50.0	#Outer Radius
r_i = 5.0	#Inner Radius
r = 0.5*(r_o-r_i)+r_i #Average Radius
stdev=0.50
center = [0.0,0.0]

nps = 500

#radii = np.random.normal(r,stdev,nps)
radii = np.random.uniform(r_i,r_o,nps)
angle = np.random.uniform(0.0,2.0*math.pi,nps)
i = 0.0
#radii = np.array([20.0,20.0,20.0,20.0])
#angle = np.array([0.0,math.pi/2,math.pi,3.0*math.pi/2])
X=np.zeros(nps)
Y=np.zeros(nps)
for i in range(0,nps):
	X[i],Y[i]=pol2cart(radii[i],angle[i])	
#delay =  radii+(radii*np.cos(angle))

G = 5.702*10.0**-11.0 #light days (solar masses)^-1 (speed of light)^2
M = 1.0e8 # solar masses


sol = 299792.458 #km/s
velocity =np.sqrt(G*M/radii)*sol

LOSvel = velocity*np.cos((math.pi/2.0)-angle)*np.sin(i)
delay = radii*(1.0+(np.cos(angle)*np.sin(i)))

fig = plt.figure(figsize=(6,6))

#ax = fig.add_subplot(121)
ax = plt.subplot2grid((1,1),(0,0),colspan=1)
ax.set_title("BLR Geometry")
col=LOSvel
sizes = (delay/max(delay))*100.0
ax.scatter(X,Y,c=-col,s=sizes,cmap="bwr",vmin=min(LOSvel),vmax=max(LOSvel))
#ax.scatter(X,Y,c=-col,cmap="bwr",vmin=min(LOSvel),vmax=max(LOSvel))
plt.minorticks_on()
#ax.set_aspect('equal')


#LABELS
xmin = -100.0
xmax = 60.0
npts = 200
x = np.linspace(xmin,xmax,npts)
y = np.zeros([npts])
plt.plot(x,y,'k',linewidth=weight)
dirstr = r"$\Longleftarrow Observer$"
ax.annotate(dirstr, xy=(xmin, 0.02))

plt.plot([0],[0],'ko',markersize=6)
ax.set_xlim(xmin,xmax)


ax.set_xlim(xmin,xmax)
ax.set_ylim(-60.0,60.0)


plt.show()




