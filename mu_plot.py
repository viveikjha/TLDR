import numpy as np
import matplotlib.pyplot as plt

inputfile = 'mu_test_mu.csv'
mu_arr = np.loadtxt(open(inputfile,'rb'),delimiter=',',skiprows=0)

mu = mu_arr[:,0]
chi2 = mu_arr[:,1]
its = mu_arr[:,2]

plt.figure(1)
plt.plot(mu,chi2)
plt.xlabel('mu')
plt.ylabel('chi2')

#plt.figure(2)
#plt.plot(mu,its)
#plt.xlabel('mu')
#plt.ylabel('Number of iterations to Converge')
plt.show()
