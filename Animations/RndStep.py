import numpy as np

def rndstep(value,step):
    steps = np.ceil((np.max(value)/step))+1
    setvals = np.arange(0.0,steps*step,step)
    for i in range(0,np.size(value)):
        for j in range(0,np.size(setvals)-1):
            if value[i] > setvals[j] and value[i] < setvals[j+1]:
                if value[i] < setvals[j+1]-(0.5*step):
                    value[i]=setvals[j]
                else:
                    value[i]=setvals[j+1]
    return np.around(value,5)


A=[1.0, 1.1, 1.6,1.9,2.0]
step=1.0/5.0
print rndstep(A,step)
