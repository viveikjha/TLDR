import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.animation as animation
import math
import PlotPretty
PlotPretty.pp('white')

stepsize=1.0/10.0
LOSvel=np.loadtxt('LOSvel.csv',delimiter=',')
Fim=np.loadtxt('Fim.csv',delimiter=',')
X=np.loadtxt('X.csv',delimiter=',')
Y=np.loadtxt('Y.csv',delimiter=',')
radii=np.loadtxt('radii.csv',delimiter=',')
delay=np.loadtxt('delay.csv',delimiter=',')
angle=np.loadtxt('angle.csv',delimiter=',')
lam=np.loadtxt('new_wavelengths.csv',delimiter=',')
X2=np.copy(X)
Y2=np.copy(Y)

ints = 200
cont = np.ones(ints)
time = np.around(np.linspace(0,ints*stepsize,ints),1)
stop=-15


def update_plot(frame, c,s):
    global X2
    global Y2
    global time
    global cont
    radcur=frame*stepsize
    time=time+stepsize
    for i in range(0,len(X)):
        if radii[i] < radcur:
            X2[i]-=stepsize
    x=radcur*c
    y=radcur*s
    ax1.cla()
    ax1.set_title("Geometry")
    ax1.set_xlim(-20,20)
    ax1.set_ylim(-20,20)
    ax1.plot(x,y)
    ax1.scatter(X,Y,c='grey')
    ax1.scatter(X2,Y2,c=-LOSvel,cmap="bwr")

    if radcur > 15 and radcur < 21:
        if radcur == 15:
            cont[0]=2.0
        elif radcur == 16:
            cont[0] = 1.0
            cont[1] = 2.0
        elif radcur == 17:
            cont[1] = 1.0
            cont[2] = 2.0
        elif radcur == 18:
            cont[2] = 1.0
            cont[3] = 2.0
        elif radcur == 19:
            cont[3] = 1.0
            cont[4] = 2.0
        elif radcur==20:
            cont[4] = 1.0
    ax2.cla()
    ax2.set_title("Continuum")
    ax2.plot(time,cont)
    plt.draw()


npts=50
phi = np.linspace(0.0,2.0*math.pi,npts)
c=np.cos(phi)
s=np.sin(phi)
fig1=plt.figure()
ax1 = plt.subplot2grid((2,2),(0,0),colspan=1)
ax1.set_title("Geometry")
ax2 = plt.subplot2grid((2,2),(0,1),colspan=1)
#ax2.set_title("Delay Map")
#ax3 = plt.subplot2grid((2,2),(1,0),colspan=1)
#ax3.set_title("Continuum")
#ax4 = plt.subplot2grid((2,2),(1,1),colspan=1)
#ax4.set_title("Spectra")

line_ani = animation.FuncAnimation(fig1, update_plot, ints, fargs=( c,s), interval=10, blit=False)
#line_ani.save('lines.mp4')
plt.show()
