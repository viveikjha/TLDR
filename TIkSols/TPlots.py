import matplotlib.pyplot as plt
import numpy as np
import PlotPretty
PlotPretty.pp('white')

S16=np.loadtxt("Res_L_0.1.csv",delimiter=',')
mini=np.argmin(S16[:,1])
plt.figure()
plt.semilogx(S16[:,0],S16[:,1],linestyle='--',color='k')

labstr=r'Best $\chi^2$ = '+ str(np.around(S16[mini,1],2)) + r' $\mu$ = ' + str(np.around(S16[mini,0],1))
plt.semilogx(S16[mini,0],S16[mini,1],'k.',markersize=15,label=labstr)
plt.axes().tick_params(axis='y', which='minor',left='on')
plt.ylim(0,200)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\chi^2$')
plt.title(r'Tikhonov Initialization With Various Smoothing Values')

plt.legend()
#plt.show()
plt.savefig("Tik_Init.png")
