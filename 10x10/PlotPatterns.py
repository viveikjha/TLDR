import matplotlib.pyplot as plt
import numpy as np

checkerboard=np.loadtxt("checkerboard.csv",delimiter=",")
circle= np.loadtxt("circle.csv",delimiter=",")
diagonal= np.loadtxt("diagonal.csv",delimiter=",")
diagonal_inverted=np.loadtxt("diagonalinverted.csv",delimiter=",")
halfbottom=np.loadtxt("halfbottom.csv",delimiter=",")
halfleft=  np.loadtxt("halfleft.csv",delimiter=",")
halfright=np.loadtxt("halfright.csv",delimiter=",")
halftop=np.loadtxt("halftop.csv",delimiter=",")
horizontalstripe=np.loadtxt("horizontalstripe.csv",delimiter=",")
invertedbox=np.loadtxt("invertedbox.csv",delimiter=",")
invertedhorizontalstripe=np.loadtxt("invertedhorizontalstripe.csv",delimiter=",")
invertedverticalstripe=np.loadtxt("invertedverticalstripe.csv",delimiter=",")
lowertri=np.loadtxt("lowertri.csv",delimiter=",")
reverseddiagonal=np.loadtxt("reverseddiagonal.csv",delimiter=",")
reverseddiagonalinverted=np.loadtxt("reverseddiagonalinverted.csv",delimiter=",")
ring=np.loadtxt("ring.csv",delimiter=",")
uppertri=np.loadtxt("uppertri.csv",delimiter=",")
verticalstripe=np.loadtxt("verticalstripe.csv",delimiter=",")

#ROW 1
ax=plt.subplot2grid((3,6),(0,0));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(circle)

ax=plt.subplot2grid((3,6),(0,1));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(checkerboard)

ax=plt.subplot2grid((3,6),(0,2));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(diagonal)

ax=plt.subplot2grid((3,6),(0,3));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(diagonal_inverted)

ax=plt.subplot2grid((3,6),(0,4));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(halfbottom)

ax=plt.subplot2grid((3,6),(0,5));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(halfleft)

#ROW 2
ax=plt.subplot2grid((3,6),(1,0));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(halfright)

ax=plt.subplot2grid((3,6),(1,1));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(halftop)

ax=plt.subplot2grid((3,6),(1,2));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(horizontalstripe)

ax=plt.subplot2grid((3,6),(1,3));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(invertedbox)

ax=plt.subplot2grid((3,6),(1,4));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(invertedhorizontalstripe)

ax=plt.subplot2grid((3,6),(1,5));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(invertedverticalstripe)

#ROW 3 lowertri reverseddiagonal reverseddiagonalinverted ring uppertri verticalstripe
row=2
ax=plt.subplot2grid((3,6),(row,0));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(lowertri)

ax=plt.subplot2grid((3,6),(row,1));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(reverseddiagonal)

ax=plt.subplot2grid((3,6),(row,2));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(reverseddiagonalinverted)

ax=plt.subplot2grid((3,6),(row,3));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(ring)

ax=plt.subplot2grid((3,6),(row,4));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(uppertri)

ax=plt.subplot2grid((3,6),(row,5));ax.set_xticklabels([]);ax.set_yticklabels([]);ax.set_xticks([]);ax.set_yticks([]);
ax.imshow(verticalstripe)

plt.tight_layout()
plt.show()
