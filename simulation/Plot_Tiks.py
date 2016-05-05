import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
import PlotPretty
from numpy import loadtxt

PlotPretty.pp('white')

mu=loadtxt("TikMus.csv",delimiter=",")
residuals=loadtxt("TikRes.csv",delimiter=",")


plt.figure()
plt.loglog(mu,residuals)
plt.title("Tikhonov Solution")
plt.xlabel(r"$\mu$")
plt.ylabel(r"$Residual$")



plt.tight_layout()
plt.savefig("TikSols.png")

plt.show()
