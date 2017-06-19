import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty
from numpy import loadtxt

PlotPretty.pp('white')


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)
weight = 2.0 #LINE WEIGHT!
r_o = 5.0
r_i = 1.0
r = 0.5*(r_o-r_i)+r_i
stdev=0.50
center = [0.0,0.0]

nps = 10000


radii = loadtxt("radii.csv")
angle = loadtxt("angle.csv")

X=loadtxt("X.csv")
Y=loadtxt("Y.csv")
delay = loadtxt("delay.csv")

G = 5.702*10.0**-11.0 #light days (solar masses)^-1 (speed of light)^2
M = 7.0e6 # solar masses

LOSvel = loadtxt("LOSvel.csv")

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
fig = plt.figure(figsize=(5,10))
#fig=plt.figure()
ax = plt.subplot2grid((2,1),(0,0),rowspan=1)
#fig = plt.figure(1)
#ax = plt.subplot(111)
ax.set_title(r"$BLR$ $Geometry$")
col=LOSvel
sizes = (delay/max(delay))*100.0
#ax.scatter(X,Y,c=-col,s=sizes,cmap="bwr",vmin=min(LOSvel),vmax=max(LOSvel))
ax.scatter(X,Y,c=-col,cmap="bwr",vmin=min(LOSvel),vmax=max(LOSvel))
plt.minorticks_on()
ax.set_aspect('equal')
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

c = 1.0



#CIRCLE AT RADIUS
theta = np.linspace(0.0,2.0*math.pi,200)
x = r * np.cos(theta)
y = r * np.sin(theta)

#plt.plot(x,y,'c',linewidth=2.0)

#LABELS
xmin = -15.0
xmax = 8.0
npts = 200
x = np.linspace(xmin,xmax,npts)
y = np.zeros([npts])
plt.plot(x,y,'k',linewidth=1.0)
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

plt.plot(x,y,'y',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.55, y[199]-0.05))
############################################
############################################
t = 1.0
theta = np.linspace(-0.85*math.pi,0.85*math.pi,200)
r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'y',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################
############################################
t = 4.0
theta = np.linspace(-0.67*math.pi,0.67*math.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'y',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################

############################################
t = 10.0
theta = np.linspace(-0.43*math.pi,0.43*math.pi,200)

r = c * t / (1.0+np.cos(theta))

x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x,y,'y',linewidth=weight)
st = r'$\tau = $'+ str(t)
plt.annotate(st, xy=(x[199]-0.5, y[199]-0.05))
############################################

ax.set_xlim(xmin,xmax)
ax.set_ylim(-10.0,10.0)

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
#ax=fig.add_subplot(122)
ax = plt.subplot2grid((2,1),(1,0),rowspan=1)
#fig=plt.figure(2)
#ax = plt.subplot(111)
ax.set_title(r"$Velocity$ $Delay$ $Map$")
ax.set_xlabel(r"$L.O.S.$ $Velocity$ $(km/s)$")
ax.set_ylabel(r"$Delay$ $(\tau)$")
ax.scatter(LOSvel,delay,c=col,cmap="bwr",vmin=min(LOSvel),vmax=max(LOSvel))
ax.plot([-8000,8000],[4.0,4.0],'y',linewidth=weight)
ax.plot([-8000,8000],[1.0,1.0],'y',linewidth=weight)
ax.plot([-8000,8000],[10.0,10.0],'y',linewidth=weight)
ax.plot([-8000,8000],[0.2,0.2],'y',linewidth=weight)
ax.set_xlim(-8000,8000)

plt.minorticks_on()
############################################################################
############################################################################
############################################################################

plt.tight_layout()
plt.savefig("Geometry_VDM2.jpg")
plt.show()
