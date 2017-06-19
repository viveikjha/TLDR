import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.animation as animation
import math
import PlotPretty
PlotPretty.pp('white')

stepsize=1.0/5.0
stop=-18
npts = 50
# UNDERLAY
X=np.loadtxt("X.csv",delimiter=',')
Y=np.loadtxt("Y.csv",delimiter=',')
LOSV=np.loadtxt("LOSvel.csv",delimiter=',')
D=np.loadtxt("delay.csv",delimiter=',')
R=(np.loadtxt("radii.csv",delimiter=','))
phi = np.linspace(0.0,2.0*math.pi,npts)
shiftX=X
n=np.size(X)
rf = np.zeros(n) #radius flag
df = np.zeros(n) #delay flagdef pol2cart(r,phi):

def pol2cart(r,phi):
    x=r*np.cos(phi)
    y=r*np.sin(phi)
    return x,y

def update_plot(frame,phi,r):
    r=1.0
    global shiftX
    global shiftY
    global LVSH
    global delays
    global pd
    plt.clf()
    radcur=frame*stepsize #current radius
    ax1 = plt.subplot2grid((1,2),(0,0),colspan=1)
    ax1.set_xlim(-30, 30)
    ax1.set_ylim(-30, 30)
    ax1.set_xlabel('x')
    ax2 = plt.subplot2grid((1,2),(0,1),colspan=1)
    stop=-10
    ax1.set_xlim(-20, 10)
    ax1.set_ylim(-10, 10)
    ax1.plot([-stop,-stop],[-10,10],color='#ffa500')
    titlestring='Radius = ' + str(radcur)
    ax1.set_title(titlestring)
    ax1.scatter(X,Y,c='grey')
    for i in range(0,n):  #CHECK Eligibility
        if R[i] <= radcur:
            rf[i]=1
        if shiftX[i] <= stop:
            df[i]=1
        if rf[i] == 1: #Shift & Plot
            shiftX[i]=shiftX[i]-stepsize
            ax1.scatter(shiftX[i],Y[i],c=-LOSV[i],cmap='bwr')
        if df[i] == 1:  #Plot
            ax2.scatter(LVSH[i],delays[i])
    x,y = pol2cart(frame/5.0*r,phi)

    ax1.plot(x,y)
    ax2.set_xlabel("LOS Velocity")
    ax2.set_ylabel("Time")
    return


ints=200
r=1.0
fig1 = plt.figure(figsize=(10,4))
contflx=np.empty([])
LVSH=np.empty([])
delays=np.empty([])
data = np.random.rand(2, 25)
r=1.0

#line_ani = animation.FuncAnimation(fig1, update_plot, ints, fargs=(phi), interval=10, blit=False)
line_ani = animation.FuncAnimation(fig1, update_plot, ints, fargs=(phi, r), interval=10, blit=False)

line_ani.save('lines.mp4')
#plt.show()
