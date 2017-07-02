import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty

PlotPretty.pp('white')

actual=np.loadtxt("UT_vdm.csv",delimiter=",")
recovered=np.loadtxt("test_box_result.csv",delimiter=",")

#fig = plt.figure(figsize=(6,4))


ax = plt.subplot2grid((1,2),(0,0))
ax.set_title(r"Actual TDF")
ax.imshow(actual,origin="lower",cmap="Greys",aspect="auto",interpolation="none")
plt.minorticks_on()

ax = plt.subplot2grid((1,2),(0,1))
ax.set_title(r"Recovered TDF")
ax.imshow(recovered,origin="lower",cmap="Greys",aspect="auto",interpolation="none")
plt.minorticks_on()

plt.show()
#plt.savefig("Box_Recovery.png")
