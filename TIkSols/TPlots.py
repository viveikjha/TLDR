import matplotlib.pyplot as plt
import numpy as np
import PlotPretty
PlotPretty.pp('white')

S16=np.loadtxt("Res_L_0.1_SNR_1.6.csv",delimiter=',')
S18=np.loadtxt("Res_L_0.1_SNR_1.8.csv",delimiter=',')
S20=np.loadtxt("Res_L_0.1_SNR_2.0.csv",delimiter=',')
S25=np.loadtxt("Res_L_0.1_SNR_2.5.csv",delimiter=',')
S30=np.loadtxt("Res_L_0.1_SNR_3.0.csv",delimiter=',')
S35=np.loadtxt("Res_L_0.1_SNR_3.5.csv",delimiter=',')
S40=np.loadtxt("Res_L_0.1_SNR_4.0.csv",delimiter=',')

plt.figure()
plt.semilogx(S16[:,0],S16[:,1],linestyle='--',color='k',label=r'$\bar{n}=1.6$')

#plt.semilogx(S18[:,0],S18[:,1],linestyle=':')
plt.semilogx(S20[:,0],S20[:,1],linestyle='-.', color='k',label=r'$\bar{n}=2.0$')
plt.semilogx(S25[:,0],S25[:,1],linestyle='-',color='0.75',label=r'$\bar{n}=2.5$')
plt.semilogx(S30[:,0],S30[:,1],linestyle=':',color='k',label=r'$\bar{n}=3.0$')
#plt.semilogx(S35[:,0],S35[:,1],linestyle=':')
plt.semilogx(S40[:,0],S40[:,1],linestyle='-',color='k',label=r'$\bar{n}=4.0$')
plt.ylim(0,200)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\chi^2$')
plt.title(r'Tikhonov Initialization At Various Noise Levels')
plt.legend()
#plt.show()
plt.savefig("Tik_init_lvl_0.1.png")
