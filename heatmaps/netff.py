import h5py
import numpy as np
from conns import *

h5f = h5py.File('network.h5','w')

N = 5
n = np.array([N],dtype=int)
inputs = np.array([0,1],dtype=int)
outputs = np.array([2],dtype=int)
knockouts = np.array([],dtype=int)
weightbounds = np.array([-2.0,+2.0],dtype=float)
pre = np.array([],dtype=int)
post = np.array([],dtype=int)

######## NETWORK SPEC

#Layer1 = inputs
#Layer2 = np.arange(2,7)
#Layer3 = np.arange(7,N)
#Layer4 = [2]

pre, post = recur(pre,post,np.arange(N))
#Layer1 = inputs
#Layer2 = knockouts
#Layer3 = extra
#Layer3 = outputs

#pre,post = ffcon(pre, post, Layer1, Layer2)
#pre,post = ffcon(pre, post, Layer2, Layer3)
#pre,post = ffcon(pre, post, Layer3, Layer4)
#pre,post = ffcon(pre, post, Layer2, Layer3)
#pre,post = ffcon(pre, post, Layer3, Layer4)
#pre,post = recur(pre, post, Layer2)
#pre,post = recur(pre, post, Layer3)

######### NETWORK SPEC

print(n)

h5f.create_dataset('N', data=n)
h5f.create_dataset('inputs', data=inputs)
h5f.create_dataset('outputs', data=outputs)
h5f.create_dataset('knockouts', data=knockouts)
h5f.create_dataset('pre', data=pre)
h5f.create_dataset('post', data=post)
h5f.create_dataset('weightbounds', data=weightbounds)

h5f.close()
