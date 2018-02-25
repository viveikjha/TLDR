import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from itertools import product
import numpy as np
import PlotPretty

PlotPretty.pp('white')
names=["checkerboard","circle","diagonal","diagonalinverted","halfbottom","halfleft","halfright","halftop","horizontalstripe","invertedbox","invertedhorizontalstripe","invertedverticalstripe","lowertri","reverseddiagonal","reverseddiagonalinverted","ring","uppertri","verticalstripe"]
names2=["checkerboard","circle","diagonal","diagonalinverted","halfbottom","halfleft","halfright","halftop","horizontalstripe","invertedbox","invertedhorizontalstripe","invertedverticalstripe","lowertri","reverseddiagonal","reverseddiagonalinverted","ring","uppertri","verticalstripe"]
#18

#i=1
#im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
#im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
fig=plt.figure(figsize=(14,8))
grid=[5, 17]
hratios=[8,1,8,1,8,1]
wratios=[4,4,1,4,4,1,4,4,1,4,4,1,4,4,1,4,4] #width_ratios=ratios
gs = gridspec.GridSpec(6, 17,hspace=0.05,wspace=0.1,width_ratios=wratios, height_ratios=hratios )
#Set1
row=0
place=row*17
i=0
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+1])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])


#Set2
i=1
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+3])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+4])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set3
i=2
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+6])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+7])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set4
i=3
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+9])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+10])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set5
i=4
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+12])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+13])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set6
i=5
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+15])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+16])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#labels
row=2
place=row*17
txt="#1"
ax=plt.subplot(gs[1,0:2])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#2"
ax=plt.subplot(gs[1,3:5])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#3"
ax=plt.subplot(gs[1,6:8])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#4"
ax=plt.subplot(gs[1,9:11])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#5"
ax=plt.subplot(gs[1,12:14])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#6"
ax=plt.subplot(gs[1,15:17])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)


#Set1
row=2
place=row*17
i=6
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+1])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])


#Set2
i=7
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+3])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+4])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set3
i=8
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+6])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+7])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set4
i=9
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+9])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+10])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set5
i=10
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+12])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+13])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set6
i=11
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+15])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+16])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

txt="#1"
ax=plt.subplot(gs[3,0:2])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#2"
ax=plt.subplot(gs[3,3:5])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#3"
ax=plt.subplot(gs[3,6:8])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#4"
ax=plt.subplot(gs[3,9:11])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#5"
ax=plt.subplot(gs[3,12:14])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#6"
ax=plt.subplot(gs[3,15:17])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)


#Set1
row=4
place=row*17
i=12
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+1])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])


#Set2
i=13
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+3])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+4])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set3
i=14
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+6])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+7])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set4
i=15
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+9])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+10])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set5
i=16
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+12])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+13])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

#Set6
i=17
im1=np.loadtxt((names[i]+"/"+names[i]+".csv"),delimiter=",")
im2=np.loadtxt((names[i]+"/"+names2[i]+".csv"),delimiter=",")
ax=plt.subplot(gs[place+15])
ax.imshow(im1,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])
ax=plt.subplot(gs[place+16])
ax.imshow(im2,interpolation="none",cmap="Greys")
ax.set_xticks([])
ax.set_yticks([])

txt="#1"
ax=plt.subplot(gs[5,0:2])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#2"
ax=plt.subplot(gs[5,3:5])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#3"
ax=plt.subplot(gs[5,6:8])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#4"
ax=plt.subplot(gs[5,9:11])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#5"
ax=plt.subplot(gs[5,12:14])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)

txt="#6"
ax=plt.subplot(gs[5,15:17])
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.5, 0.5, txt, ha='center', va='center', size=18)


#gs.update(left=0.125,right=0.9,top=0.9,bottom=0.1,wspace=0.2,hspace=0.2)
#plt.tight_layout()
plt.show()
