import numpy as np
import matplotlib.pyplot as plt
#from matsemilogylib import gridspec
import math
import PlotPretty
PlotPretty.pp('white')

cf=np.loadtxt("Conflag.csv",delimiter=",");
chi2=np.loadtxt("Chi2.csv",delimiter=",");
J=np.loadtxt("J.csv",delimiter=",");
N=np.loadtxt("N_res.csv",delimiter=",");
P=np.loadtxt("P_res.csv",delimiter=",");
T=np.loadtxt("T_res.csv",delimiter=",");
V=np.loadtxt("V_res.csv",delimiter=",");
Z=np.loadtxt("Z_res.csv",delimiter=",");
Dchi2=[]
DJ=[]
DN=[]
DP=[]
DT=[]
DV=[]
DZ=[]
for i in range(0,len(cf)):
    if cf[i] == 0 or cf[i]==1:
        Dchi2=np.append(Dchi2,chi2[i])
        DJ=np.append(DJ,J[i])
        DN=np.append(DN,N[i])
        DT=np.append(DT,T[i])
        DV=np.append(DV,V[i])
        DP=np.append(DP,P[i])
        DZ=np.append(DZ,Z[i])
plt.figure()
plt.title(r"Reconstruction Residuals")
#plt.plot(Dchi2/np.max(Dchi2),label=r"$\chi^2_{reduced}$")
#plt.plot(DJ/np.max(DJ)label)
plt.plot(DZ/np.max(DZ),"-",label=r"$\|\|X-Z\|\|_2$")
plt.plot(DN/np.max(DN),"--",label=r"$\|\|X-N\|\|_2$")
plt.plot(DP/np.max(DP),"-.",label=r"$\|\|X-P\|\|_2$")
plt.plot(DT/np.max(DT),"-",label=r"$\|\|D_TX-T\|\|_2$")
plt.plot(DV/np.max(DV),"--",label=r"$\|\|XD_V-V\|\|_2$")

plt.xlabel(r"Iteration")
plt.ylabel(r"Normalized Residual")
plt.legend()
plt.minorticks_on()

plt.show()
#plt.savefig("Residuals.png")
