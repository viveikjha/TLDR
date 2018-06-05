import matplotlib.pyplot as plt
import numpy as np
import PlotPretty
PlotPretty.pp('white')

Result=np.loadtxt("Res_L_1.0.csv",delimiter=',')
mini=np.argmin(Result[:,1]) #Lowest Chi2
mxpsnr=np.argmax(Result[:,3]) #Highest PSNR
#plt.figure()
#plt.semilogx(Result[:,0],Result[:,1],linestyle='--',color='k')
#labstr=r'Best $\chi^2$ = '+ str(np.around(Result[mini,1],2)) + r' $\mu$ = ' + str(np.around(Result[mini,0],1))
#plt.semilogx(Result[mini,0],Result[mini,1],'k.',markersize=15,label=labstr)

#plt.axes().tick_params(axis='y', which='minor',left='on')
#plt.ylim(0,25)
#plt.xlabel(r'$\mu$')
#plt.ylabel(r'$\chi^2$')
#plt.title(r'Tikhonov Initialization With Various Smoothing Values')

#plt.legend()

#plt.show()



#plt.figure()
#plt.semilogx(Result[:,0],Result[:,3])
#plt.show()

fig=plt.figure()
plt.title(r'Tikhonov Initialization With Various Smoothing Values',y=1.02)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax=fig.add_subplot(111)

a,= ax.semilogx(Result[:,0],Result[:,1],linestyle='-',color='k',label=r'$\chi^2$')
labstrchi2=r'Best $\chi^2$ = '+ str(np.around(Result[mini,1],2)) + r' $\mu$ = ' + str(np.around(Result[mini,0],1))
b,= ax.semilogx(Result[mini,0],Result[mini,1],'k.',markersize=15,markerfacecolor='none',label=labstrchi2)
#ax.set_ylim(0,30000)
ax2=ax.twinx()
c,=ax2.semilogx(Result[:,0],Result[:,3],linestyle='--',color='k',label="PSNR")
labstrpsnr=r'Best $PSNR$ = '+ str(np.around(Result[mxpsnr,3],2)) + r' $\mu$ = ' + str(np.around(Result[mxpsnr,0],1))
d,=ax2.semilogx(Result[mxpsnr,0],Result[mxpsnr,3],'k.',markersize=15,label=labstrpsnr)
ax.set_xlabel(r'$\mu$')
ax.set_ylabel(r'$\chi^2$')
ax2.set_ylabel(r'$PSNR (dB)$')
lns=[a,c,b,d]
labs=[a.get_label(),c.get_label(),b.get_label(),d.get_label()]
ax.legend(lns,labs)
#plt.tight_layout()
#plt.show()

plt.savefig("Tik_Init.png")
