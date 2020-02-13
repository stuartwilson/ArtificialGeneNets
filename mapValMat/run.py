import h5py
import numpy as np
from conns import *
import subprocess
import time

h5f = h5py.File('network.h5','w')

N = 32
n = np.array([N],dtype=int)
inputs = np.array([0,1],dtype=int)
outputs = np.array([2],dtype=int)
context = np.array([3],dtype=int)
pre = np.array([],dtype=int)
post = np.array([],dtype=int)

######## NETWORK SPEC
pre, post = recur(pre,post,np.arange(N))
######### NETWORK SPEC

h5f.create_dataset('N', data=n)
h5f.create_dataset('inputs', data=inputs)
h5f.create_dataset('outputs', data=outputs)
h5f.create_dataset('context', data=context)
h5f.create_dataset('pre', data=pre)
h5f.create_dataset('post', data=post)

h5f.close()



######## RUN THE MODEL
p = []

t = 1000000
n = 10
nbatch = 3

k = 0
while(k<n):
    for i in range(nbatch):
        dst = 'data/expt'+str(k)
        subprocess.run('mkdir '+dst,shell=True)
        subprocess.run('cp inputs.h5 '+dst+'/inputs.h5',shell=True)
        subprocess.run('cp network.h5 '+dst+'/network.h5',shell=True)
        p = np.hstack([p,subprocess.Popen('./../build/sim/hmap config.json  '+dst+'  '+str(t)+' '+str(k),shell=True)])
        k = k+1

    ready=False
    while(not ready):
        ready=True
        for i in range(len(p)):
            if(p[i].poll()is None):
                ready=False



