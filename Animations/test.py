import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.animation as animation
import math
import PlotPretty
PlotPretty.pp('white')

stepsize=1.0/5.0


def rndstep(value):
    steps = np.ceil((np.max(value)/stepsize))+1
    setvals = np.arange(0.0,steps*stepsize,stepsize)
    for i in range(0,np.size(value)):
        for j in range(0,np.size(setvals)-1):
            if value[i] > setvals[j] and value[i] < setvals[j+1]:
                if value[i] < setvals[j+1]-(0.5*stepsize):
                    value[i]=setvals[j]
                else:
                    value[i]=setvals[j+1]
    return np.around(value,5)

# UNDERLAY
X=np.loadtxt("X.csv",delimiter=',')
Y=np.loadtxt("Y.csv",delimiter=',')
LOSV=np.loadtxt("LOSvel.csv",delimiter=',')
D=np.loadtxt("delay.csv",delimiter=',')
print np.shape(LOSV)
R=rndstep(np.loadtxt("radii.csv",delimiter=','))
car=np.ones(np.size(X))
lr=len(R)

def ler(r,el):
    return el==r

def rcheck(r,shiftX,shiftY,LVSH,delays):
    for i in range(0,lr):
        if R[i] == r:
            shiftX = np.append(shiftX,X[i])
            shiftY = np.append(shiftY,Y[i])
            LVSH = np.append(LVSH,LOSV[i])
            delays = np.append(delays,D[i])
    return shiftX,shiftY,LVSH,delays

def dcheck(shiftX,LVSH,delays,stop):
    for i in range(0,np.size(shiftX)):
        if shiftX[i] <= stop:
            LVSH=np.append(LVSH,LOSV[i])
            delays = np.append(delays,D[i])
    return LVSH,delays

def pol2cart(r,phi):
    x=r*np.cos(phi)
    y=r*np.sin(phi)
    return x,y

def update_plot(frame, phi,r):
    global shiftX
    global shiftY
    global LVSH
    global delays
    radcur=frame*stepsize
    plt.clf()
    ax1 = plt.subplot2grid((1,2),(0,0),colspan=1)
    ax1.set_xlim(-30, 30)
    ax1.set_ylim(-30, 30)
    ax1.set_xlabel('x')
    ax2 = plt.subplot2grid((1,2),(0,1),colspan=1)

    stop=-18
    ax1.set_xlim(-20, 10)
    ax1.set_ylim(-10, 10)
    ax1.plot([-stop,-stop],[-10,10],color='#ffa500')
    titlestring='Radius = ' + str(radcur)
    ax1.set_title(titlestring)
    ax1.scatter(X,Y,c='grey')
    shiftX,shiftY,LVSH,delays=rcheck(radcur,shiftX,shiftY,LVSH,delays)
    if np.size(shiftX)>1:
        shiftX=shiftX-stepsize
        ax1.scatter(shiftX,shiftY,c=-LVSH,cmap="bwr")
        LVSH,delays=dcheck(shiftX,LVSH,delays,stop)
    x,y = pol2cart(frame/5.0*r,phi)
    ax1.plot(x,y)
    ax2.set_xlabel("LOS Velocity")
    ax2.set_ylabel("Time")

    if radcur >= abs(stop):
        ax2.scatter(LVSH,delays,c=-LVSH,cmap='bwr')

    return


ints=200
fig1 = plt.figure(figsize=(10,4))
npts=50
shiftX=np.empty([])
shiftY=np.empty([])
contflx=np.empty([])
LVSH=np.empty([])
delays=np.empty([])
phi = np.linspace(0.0,2.0*math.pi,npts)
data = np.random.rand(2, 25)
ax1 = plt.subplot2grid((1,2),(0,0),colspan=1)
ax1.set_xlim(-30, 30)
ax1.set_ylim(-30, 30)
ax1.set_xlabel('x')
ax2 = plt.subplot2grid((1,2),(0,1),colspan=1)

r=1.0

line_ani = animation.FuncAnimation(fig1, update_plot, ints, fargs=(phi, r), interval=10, blit=False)
#line_ani.save('lines.mp4')
plt.show()
