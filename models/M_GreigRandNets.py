import h5py
import numpy as np
from conns import *
import subprocess
import sys
import pylab as pl

def getMap(fname):
    X = np.array([],dtype=int)
    Y = np.array([],dtype=int)
    M = np.array([],dtype=int)
    I = pl.imread(fname)
    Imax = 1./np.max(I.shape)
    for y in range(I.shape[0]):
        for x in range(I.shape[1]):
            if(I[y,x,3]==1):
                p = I[y,x,:3]
                X = np.hstack([X,x*Imax])
                Y = np.hstack([Y,(I.shape[0]-y-1)*Imax])
                if(np.array_equal(p,np.array([1.0,0.0,0.0],dtype='float32'))):
                    M = np.hstack([M,1]) # red = v1
                elif(np.array_equal(p,np.array([0.0,1.0,0.0],dtype='float32'))):
                    M = np.hstack([M,3]) # green = m1
                elif(np.array_equal(p,np.array([0.0,0.0,1.0],dtype='float32'))):
                    M = np.hstack([M,2]) # blue = s1
                elif(np.array_equal(p,np.array([1.0,0.0,1.0],dtype='float32'))):
                    M = np.hstack([M,4]) # pink = a1
                elif(np.array_equal(p,np.array([0.0,0.0,0.0],dtype='float32'))):
                    M = np.hstack([M,5]) # black is mixed V1/S1
                else:
                    M = np.hstack([M,0])
    return X, Y, M, len(M)

L = np.zeros([5])

x,y,w,L[0] = getMap('images/greig/wild.png')
x,y,ko1,L[1] = getMap('images/greig/coup.png')
x,y,ko2,L[2] = getMap('images/greig/pax6.png')
x,y,ko3,L[3] = getMap('images/greig/sp8.png')
x,y,ko4,L[4] = getMap('images/greig/emx2.png')

minL = int(np.min(L))
x = x[:minL]
y = y[:minL]
w = w[:minL]
ko1 = ko1[:minL]
ko2 = ko2[:minL]
ko3 = ko3[:minL]
ko4 = ko4[:minL]

inPatterns = np.array([],dtype=float)
inPatterns = np.hstack([inPatterns, x])
inPatterns = np.hstack([inPatterns, y])
inPatterns = np.hstack([inPatterns, max(x)-x])
inPatterns = np.hstack([inPatterns, max(y)-y])

outPatterns = np.array([],dtype=float)
outPatterns = np.hstack([outPatterns, [0.0,0.0,0.0,0.0]]) # no field
outPatterns = np.hstack([outPatterns, [1.0,0.0,0.0,0.0]]) # V1
outPatterns = np.hstack([outPatterns, [0.0,1.0,0.0,0.0]]) # S1
outPatterns = np.hstack([outPatterns, [0.0,0.0,1.0,0.0]]) # M1
outPatterns = np.hstack([outPatterns, [0.0,0.0,0.0,1.0]]) # A1
outPatterns = np.hstack([outPatterns, [1.0,1.0,0.0,0.0]]) # mixed V1/S1

maps = np.array([],dtype=int)
maps = np.hstack([maps, w])     # wild-type
maps = np.hstack([maps, ko1])   # coup knockout pattern
maps = np.hstack([maps, ko2])   # Pax6 knockout pattern
maps = np.hstack([maps, ko3])   # Sp8  knockout pattern
maps = np.hstack([maps, ko4])   # Emx2 knockout pattern


# COMMON
dstRoot = 'data/expt'

inputs = np.array([0,1,2,3],dtype=int)
outputs = np.array([4,5,6,7],dtype=int)
knockouts = np.array([8,9,10,11],dtype=int)

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
        h5f.create_dataset('inputs', data=inputs)
        h5f.create_dataset('outputs', data=outputs)
        h5f.create_dataset('knockouts', data=knockouts)
        h5f.create_dataset('pre', data=pre)
        h5f.create_dataset('post', data=post)
        h5f.create_dataset('x', data=x)
        h5f.create_dataset('y', data=y)
        h5f.create_dataset('inPatterns', data=inPatterns)
        h5f.create_dataset('outPatterns', data=outPatterns)
        h5f.create_dataset('maps', data=maps)
        h5f.close()

        ###
        #subprocess.run('cp inputs.h5 '+dst+'/inputs.h5',shell=True)
        P = np.hstack([P,subprocess.Popen('./../build/sim/greig '+dst+'  '+str(t)+' '+str(Seeds[j]),shell=True)])

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




