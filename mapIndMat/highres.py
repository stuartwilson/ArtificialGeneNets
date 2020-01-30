import h5py
import numpy as np
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

x,y,w,L[0] = getMap('images/wild.png')
x,y,ko1,L[1] = getMap('images/coup.png')
x,y,ko2,L[2] = getMap('images/pax6.png')
x,y,ko3,L[3] = getMap('images/sp8.png')
x,y,ko4,L[4] = getMap('images/emx2.png')

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

h5f = h5py.File('inputs.h5','w')

h5f.create_dataset('x', data=x)
h5f.create_dataset('y', data=y)
h5f.create_dataset('inPatterns', data=inPatterns)
h5f.create_dataset('outPatterns', data=outPatterns)
h5f.create_dataset('maps', data=maps)

h5f.close()
