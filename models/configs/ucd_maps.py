import h5py
import numpy as np
import pylab as pl
import sys

def interp(x):
    return (np.max(x)-np.min(x))*np.argsort(x)/(len(x)-1)+np.min(x)


#file = '18_159'
file = sys.argv[1]

A = h5py.File(file,'r')

n = A['nframes'][:][0]

X = np.array([],dtype=float)
Y = np.array([],dtype=float)
F = np.array([],dtype=float)

for i in range(n-1):
    a = A['Frame'+str(i+1).zfill(3)]
    f = np.array([a['means_autoscaled'][:]],dtype=float)[0]
    y = a['sbox_linear_distance'][:]
    x = np.ones(len(f))*(1.0)*i/(n-1)
    X = np.hstack([X,x])
    Y = np.hstack([Y,interp(y)])
    F = np.hstack([F,f])
A.close()

Y = Y - np.min(Y)
X = X*np.max(Y)

B = h5py.File(file[:-3]+'_c.h5','w')
B.create_dataset('X', data=X)
B.create_dataset('Y', data=Y)
B.create_dataset('Z', data=np.zeros(len(Y)))
B.create_dataset('F', data=F)
B.close()

'''
F = pl.figure(figsize=(10,10))
f = F.add_subplot(111)
f.plot(X,Y,'.')
pl.show()
'''
