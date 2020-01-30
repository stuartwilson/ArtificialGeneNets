import h5py
import numpy as np
from conns import *

h5f = h5py.File('network.h5','w')

N = 22
n = np.array([N],dtype=int)
inputs = np.array([0,1,2,3],dtype=int)
outputs = np.array([4,5,6,7],dtype=int)
knockouts = np.array([8,9,10,11],dtype=int)
pre = np.array([],dtype=int)
post = np.array([],dtype=int)

######## NETWORK SPEC
pre, post = recur(pre,post,np.arange(N))
######### NETWORK SPEC

h5f.create_dataset('N', data=n)
h5f.create_dataset('inputs', data=inputs)
h5f.create_dataset('outputs', data=outputs)
h5f.create_dataset('knockouts', data=knockouts)
h5f.create_dataset('pre', data=pre)
h5f.create_dataset('post', data=post)

h5f.close()