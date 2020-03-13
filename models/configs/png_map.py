import numpy as np
import pylab as pl
import h5py
import sys

def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])

def getImgData(fname):
    I = pl.imread(fname)
    Q = rgb2gray(I)
    print(Q.shape)
    X = np.array([],dtype=int)
    Y = np.array([],dtype=int)
    Z = np.array([],dtype=int)
    Imax = 1./np.max(I.shape)
    for y in range(I.shape[0]):
        for x in range(I.shape[1]):
            if (I[y,x,3]):
                X = np.hstack([X,x*Imax])
                Y = np.hstack([Y,(I.shape[0]-y-1)*Imax])
                Z = np.hstack([Z,Q[y,x]])

    return X,Y,Z

X,Y,Z = getImgData(sys.argv[1])

h5f = h5py.File(sys.argv[2],'w')
h5f.create_dataset('X', data=X)
h5f.create_dataset('Y', data=Y)
h5f.create_dataset('Z', data=np.zeros([len(X)]))
h5f.create_dataset('F', data=Z)
h5f.close()

'''
if(len(sys.argv)==3):
    if (int(sys.argv[2])==1):
        F = pl.figure()
        f = F.add_subplot(111)
        for i in range(len(X)):
            f.plot(X[i],Y[i],'.',color=np.ones([3])*Z[i])
        f.set_aspect('equal')
        pl.show()
'''
