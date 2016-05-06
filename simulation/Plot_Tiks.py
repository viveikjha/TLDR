import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty
from numpy import loadtxt

PlotPretty.pp('white')


arr = loadtxt("TikSol.csv",delimiter=",")

mu = arr[:,0]
residual = arr[:,1]
flux = arr[:,2]
chi2 = arr[:,3]

fig=plt.figure(figsize=(10,10))
ax = plt.subplot2grid((2,2),(0,0),colspan=1)
ax.set_title(r"L-Curve?")
plt.loglog(flux,residual)
plt.title("L-Curve?")
plt.xlabel(r"$\|\|Tikhonov$ $Solution\|\|_2$")
plt.ylabel(r"$\|\|Tikhonov$ $Solution-VDM_{True}\|\|_2$")

ax = plt.subplot2grid((2,2),(0,1),colspan=1)
ax.set_title(r"Residual vs $\mu$")
plt.loglog(mu,residual)
plt.xlabel(r"$\mu$")
plt.ylabel(r"$\|\|Tikhonov$ $Solution-VDM_{True}\|\|_2$")


ax = plt.subplot2grid((2,2),(1,0),colspan=1)
ax.set_title(r"Reduced $\chi^2$ vs $\mu$")
plt.loglog(mu,chi2)
plt.xlabel(r"$\mu$")
plt.ylabel(r"Reduced $\chi^2$")


ax = plt.subplot2grid((2,2),(1,1),colspan=1)
ax.set_title(r"Residual vs Reduced $\chi^2$")
plt.loglog(chi2,residual)
plt.ylabel(r"$\|\|Tikhonov$ $Solution-VDM_{True}\|\|_2$")
plt.xlabel(r"Reduced $\chi^2$")

plt.tight_layout()
plt.savefig("TikSols.png")

plt.show()
