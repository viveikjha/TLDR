import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


inputfile = 'tdf.txt'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
rec = data[0,:]
print rec.shape


inputfile = 'PsiS.txt'
data = np.loadtxt(open(inputfile,'rb'),delimiter=' ',skiprows=0)
act = data[1,:]
print act.shape


plt.figure()
plt.plot(act,'r')
#plt.plot(rec,'b')

plt.show()