import h5py
import numpy as np
from conns import *
import subprocess
import sys

# COMMON
dstRoot = 'data/expt'

if(len(sys.argv)<5):
    print("Specify: Nnodes Nbatch Nsims timesteps")

Nnodes = int(sys.argv[1])
Nbatch = int(sys.argv[2])
Nsims = int(sys.argv[3])
t = int(sys.argv[4])

Seeds = np.arange(Nsims)
Sizes = np.ones(Nsims)*Nnodes

cullMax = 0.95
toCull = np.zeros(Nsims,dtype=int)
for i in range(Nsims):
    toCull[i] = int(cullMax*(i/(Nsims-1))*Nnodes*Nnodes)

finErr = np.zeros(Nsims)
minErr = np.zeros(Nsims)

######## RUN THE MODEL
j = 0
k = 0
running = True
while(running):

    P = []
    for i in range(Nbatch):

        if not running: break

        dst = dstRoot+str(j)
        subprocess.run('mkdir '+dst,shell=True)
        subprocess.run('cp configs/configEnucleateOriginal.json '+dst+'/config.json',shell=True)
        subprocess.run('cp maps/sighted.h5 '+dst+'/sighted.h5',shell=True)
        subprocess.run('cp maps/enucleate.h5 '+dst+'/enucleate.h5',shell=True)

        ### NETWORK SPEC
        N = Sizes[j]
        Narr = np.array([N],dtype=int)
        pre = np.array([],dtype=int)
        post = np.array([],dtype=int)
        pre, post = recur(pre,post,np.arange(N))
        for c in range(toCull[j]):
            pre, post = cullRand(pre,post)

        ###
        h5f = h5py.File(dst+'/network.h5','w')
        h5f.create_dataset('pre', data=pre)
        h5f.create_dataset('post', data=post)
        h5f.close()

        ###
        #subprocess.run('cp inputs.h5 '+dst+'/inputs.h5',shell=True)
        P = np.hstack([P,subprocess.Popen('./../build/sim/hmap '+dst+'  '+str(t)+' '+str(Seeds[j]),shell=True)])

        ###
        running=j<(Nsims-1)
        j+=1

    waitUntilReady(P)


    ######## PERFORM ANALYSIS

    for i in range(Nbatch):
        if not running: break

        try:
            dst = dstRoot+str(k)
            h5f = h5py.File(dst + '/outputs.h5','r')
            err = h5f['error'][:]
            h5f.close()

            finErr[k] = err[-1]
            minErr[k] = np.min(err)
        except:
            print("no output"+str(k)+"\n")

        running=k<(Nsims-1)
        k+=1
        print(k)

    ### store results periodically
    h5f = h5py.File('summary.h5','w')
    h5f.create_dataset('finErr', data=finErr)
    h5f.create_dataset('minErr', data=minErr)
    h5f.close()
