import h5py
import numpy as np
import pylab as pl

'''
scp -r spwilson@corebeast.shef.ac.uk:MODELS/gitprojects/ArtificialGeneNets/mapIndMat/summary.h5/ ~/MODELS/gitprojects/ArtificialGeneNets/sandbox/

scp -r spwilson@corebeast.shef.ac.uk:MODELS/gitprojects/ArtificialGeneNets/mapValMat/summary.h5/ ~/MODELS/gitprojects/ArtificialGeneNets/sandbox/

scp -r spwilson@corebeast.shef.ac.uk:MODELS/gitprojects/ArtificialGeneNets/model/summary.h5/ ~/MODELS/gitprojects/ArtificialGeneNets/sandbox/

'''


fs = 15
range = 0.95

h5f = h5py.File('summary.h5','r')
minErr = h5f['minErr'][:]
finErr = h5f['finErr'][:]

zeroErrID = np.where(minErr==0) # zero errors are spurious
minErrCp = minErr.copy()
minErrCp[zeroErrID] = 10000
best = np.argmin(minErrCp)

print("BEST: "+str(best))
print("scp -r spwilson@corebeast.shef.ac.uk:MODELS/gitprojects/ArtificialGeneNets/mapIndMat/data/expt"+str(best)+"/ ~/MODELS/gitprojects/ArtificialGeneNets/sandbox/tmp")


minErr[zeroErrID]=np.NaN
x = 1-range*np.arange(len(minErr))/(len(minErr))

F = pl.figure(figsize=(6,5))
f = F.add_subplot(111)
f.plot(x[best],minErr[best],'o',color='red',markersize=8)
f.plot(x,minErr,'.',color='k')
#f.axis([0,1.0,0,0.5])
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_xlabel('proportion of connections',fontsize=fs)
f.set_ylabel('error',fontsize=fs)
f.set_title('$N=12$')

'''
f = F.add_subplot(122)
f.plot(np.log(1-x),np.log(minErr),'.',color='k')
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_xlabel('ln(proportion of connections)')
f.set_ylabel('ln(error)')
'''
F.savefig('summary.pdf')
pl.show()
