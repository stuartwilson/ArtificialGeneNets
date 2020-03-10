import h5py
import numpy as np
from conns import *
import subprocess
import sys

def getImg(fname):
    x = np.genfromtxt(fname, delimiter=',')
    x = x[1:,1:]
    y = x
    for i in range(x.shape[1]):
        n = len(np.where(np.isnan(x[:,i]))[0])
        line = np.hstack([np.ones(int(n/2))*np.nan,x[:-int(n/2),i]])
        if(line.shape[0]):
            y[:,i] = line
    return y

def getData(a):
    x=np.array([],dtype=float)
    y=np.array([],dtype=float)
    z=np.array([],dtype=float)

    norm1 = a.shape[0]/a.shape[1]
    norm2 = 1.0/np.fmax(a.shape[0],a.shape[1])
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if(~np.isnan(a[i,j])):

                y = np.hstack([y,(1.*i)*norm2])
                x = np.hstack([x,(1.*j)*norm1*norm2])
                z = np.hstack([z,a[i,j]])
    return x,y,z

'''
    SETUP INPUTS
'''

a = getImg('images/case_19_50_P36_pup_sighted_control.csv')
b = getImg('images/case_19_18_P36_pup_enucleated_at_P4.csv')

ax,ay,az = getData(a)
bx,by,bz = getData(b)
L = np.array([len(ax),len(bx)])
minL = int(np.min(L))
ax=ax[:minL]
ay=ay[:minL]
az=az[:minL]
bx=bx[:minL]
by=by[:minL]
bz=bz[:minL]

jointnorm = False

if jointnorm :

    minz = np.min([np.min(az),np.min(bz)])
    maxz = np.max([np.max(az),np.max(bz)])
    norm = 1./(maxz-minz)
    for i in range(len(az)):
        az[i] = (az[i]-minz)*norm
    for i in range(len(bz)):
        bz[i] = (bz[i]-minz)*norm

else:

    minz = np.min(az)
    norm = 1./(np.max(az)-minz)
    for i in range(len(az)):
        az[i] = (az[i]-minz)*norm

    minz = np.min(bz)
    norm = 1./(np.max(bz)-minz)
    for i in range(len(bz)):
        bz[i] = (bz[i]-minz)*norm

#az = 1./(1.+np.exp(-(az-0.3)*1000.))
#bz = 1./(1.+np.exp(-(bz-0.3)*1000.))

dims = np.array([2,minL],dtype=int)

inPatterns = np.array([],dtype=float)
inPatterns = np.hstack([inPatterns, ax])
inPatterns = np.hstack([inPatterns, ay])
inPatterns = np.hstack([inPatterns, bx])
inPatterns = np.hstack([inPatterns, by])

maps = np.array([],dtype=float)
maps = np.hstack([maps, az])     # control
maps = np.hstack([maps, bz])     # enucleate

#h5f.close()



# COMMON
dstRoot = 'data/expt'

inputs = np.array([0,1],dtype=int)
outputs = np.array([2],dtype=int)
context = np.array([3],dtype=int)
knockouts = np.array([],dtype=int)

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
        h5f.create_dataset('N', data=Narr)
        h5f.create_dataset('inputs', data=inputs)
        h5f.create_dataset('outputs', data=outputs)
        h5f.create_dataset('knockouts', data=knockouts)
        h5f.create_dataset('context', data=context)
        h5f.create_dataset('pre', data=pre)
        h5f.create_dataset('post', data=post)
        h5f.create_dataset('inPatterns', data=inPatterns)
        h5f.create_dataset('maps', data=maps)
        h5f.create_dataset('dims', data=dims)
        h5f.close()

        ###
        subprocess.run('cp inputs.h5 '+dst+'/inputs.h5',shell=True)
        P = np.hstack([P,subprocess.Popen('./../build/sim/hmap config.json  '+dst+'  '+str(t)+' '+str(Seeds[j]),shell=True)])

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
