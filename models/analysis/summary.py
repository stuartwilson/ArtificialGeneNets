import h5py
import sys
import numpy as np
import pylab as pl

h5f = h5py.File(sys.argv[1],'r')
minErr = h5f['minErr'][:]
finErr = h5f['finErr'][:]

zeroErrID = np.where(minErr==0) # zero errors are spurious
minErrCp = minErr.copy()
minErrCp[zeroErrID] = 10000
best = np.argmin(minErrCp)
minErr[zeroErrID]=np.NaN
x = np.arange(len(minErr))/(len(minErr))

print("BEST: "+str(best))


### PLOTTING

fs = 15

F = pl.figure(figsize=(6,5))
f = F.add_subplot(111)
f.plot(x[best],minErr[best],'o',color='red',markersize=8)
f.plot(x,minErr,'.',color='k')
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_xlabel('x',fontsize=fs)
f.set_ylabel('error',fontsize=fs)
f.set_title("BEST: "+str(best))

F.savefig(sys.argv[1]+'.pdf')
pl.show()
