import numpy as np
import pylab as pl
import h5py
from conns import *

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

h5f = h5py.File('inputs.h5','w')
h5f.create_dataset('inPatterns', data=inPatterns)
h5f.create_dataset('maps', data=maps)
h5f.create_dataset('dims', data=dims)

h5f.close()

'''
    SETUP NETWORK
'''


N = 12
n = np.array([N],dtype=int)
inputs = np.array([0,1],dtype=int)
outputs = np.array([2],dtype=int)
knockouts = np.array([],dtype=int)
weightbounds = np.array([-2.5,+2.5],dtype=float)
pre = np.array([],dtype=int)
post = np.array([],dtype=int)

######## NETWORK SPEC
pre, post = recur(pre,post,np.arange(N))

h5f = h5py.File('network.h5','w')
h5f.create_dataset('N', data=n)
h5f.create_dataset('inputs', data=inputs)
h5f.create_dataset('outputs', data=outputs)
h5f.create_dataset('knockouts', data=knockouts)
h5f.create_dataset('pre', data=pre)
h5f.create_dataset('post', data=post)
h5f.create_dataset('weightbounds', data=weightbounds)

h5f.close()
