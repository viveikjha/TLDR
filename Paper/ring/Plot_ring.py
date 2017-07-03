import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty

PlotPretty.pp('white')

actual=np.loadtxt("simulated_ring_vdm.csv",delimiter=",")
recovered=np.loadtxt("RevMapResult.csv",delimiter=",")

#fig = plt.figure(figsize=(6,4))


ax = plt.subplot2grid((1,2),(0,0))
ax.set_title(r"Actual TDF")
ax.imshow(actual,origin="lower",cmap="Greys",aspect="auto",interpolation='None')
plt.minorticks_on()

ax = plt.subplot2grid((1,2),(0,1))
ax.set_title(r"Recovered TDF")
ax.imshow(recovered,origin="lower",cmap="Greys",aspect="auto",interpolation='None')
plt.minorticks_on()

#plt.show()
plt.savefig("Ring_Recovery.png")
