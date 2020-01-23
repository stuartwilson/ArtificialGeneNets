import h5py
import numpy as np
import pylab as pl
import sys

logdir = sys.argv[1]
print(logdir)

h5f = h5py.File(logdir + '/outputs.h5','r')
err = h5f['error'][:]
h5f.close()


h5f = h5py.File(logdir + '/weightmat.h5','r')
W = h5f['weightmat'][:]
n = int(np.sqrt(W.shape[0]))
W = W.reshape([n,n])
h5f.close()


F = pl.figure(figsize=(10,8))
f = F.add_subplot(211)
f.plot(err)
f.set_xlabel('time')
f.set_ylabel('error')

f = F.add_subplot(212)
g = f.pcolor(W,cmap='jet')
F.colorbar(g)
f.set_xlabel('presynaptic')
f.set_ylabel('postsynaptic')
f.set_aspect(1.0)

pl.savefig('test.pdf')

pl.show()
