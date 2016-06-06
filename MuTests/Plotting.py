import numpy as np
import matplotlib.pyplot as plt


I=np.loadtxt('Results.csv',delimiter=',')
mu =I[:,0]
chi2=I[:,1]
reg =I[:,2]
res = I[:,3]
con = I[:,4]

fig, ax1 = plt.subplots()
ax1.loglog(mu, chi2, 'b')
ax1.set_xlabel('res')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('chi2', color='b')
for i in range(1,np.size(mu)):
    if con[i] == 1.0:
        ax1.loglog(mu[i],chi2[i],"g.")
    else:
        ax1.loglog(mu[i],chi2[i],"r.")

for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()

ax2.loglog(mu, reg, 'r')
ax2.set_ylabel('reg', color='r')
for i in range(1,np.size(mu)):
    if con[i] == 1.0:
        ax2.loglog(mu[i],reg[i],"g.")
    else:
        ax2.loglog(mu[i],reg[i],"r.")

for tl in ax2.get_yticklabels():
    tl.set_color('r')
plt.show()
