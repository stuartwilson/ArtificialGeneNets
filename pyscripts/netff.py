import h5py
import numpy as np
from conns import *

h5f = h5py.File('../configs/network.h5','w')

N = 12
n = np.array([N],dtype=int)
inputs = np.array([0,1,2,3],dtype=int)
outputs = np.array([4,5,6,7],dtype=int)
knockouts = np.array([8,9,10,11],dtype=int)
weightbounds = np.array([-2.0,+2.0],dtype=float)
pre = np.array([],dtype=int)
post = np.array([],dtype=int)

######## NETWORK SPEC

extra = np.arange(12,N)

pre, post = recur(pre,post,np.arange(N))
#Layer1 = inputs
#Layer2 = knockouts
#Layer3 = extra
#Layer3 = outputs

#pre,post = ffcon(pre, post, Layer1, Layer2)
#pre,post = ffcon(pre, post, Layer2, Layer3)
#pre,post = ffcon(pre, post, Layer3, Layer4)
#pre,post = recur(pre, post, Layer2)
#pre,post = recur(pre, post, Layer3)

######### NETWORK SPEC


h5f.create_dataset('N', data=n)
h5f.create_dataset('inputs', data=inputs)
h5f.create_dataset('outputs', data=outputs)
h5f.create_dataset('knockouts', data=knockouts)
h5f.create_dataset('pre', data=pre)
h5f.create_dataset('post', data=post)
h5f.create_dataset('weightbounds', data=weightbounds)

h5f.close()
