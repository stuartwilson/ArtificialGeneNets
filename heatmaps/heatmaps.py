import numpy as np
import pylab as pl
import h5py

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
    minz = np.min(z)
    norm = 1./(np.max(z)-minz)
    for i in range(len(z)):
        z[i] = (z[i]-minz)*norm
    return x,y,z


a = getImg('case_19_50_P36_pup_sighted_control.csv')
b = getImg('case_19_18_P36_pup_enucleated_at_P4.csv')

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

az = 1./(1.+np.exp(-(az-0.3)*1000.))
bz = 1./(1.+np.exp(-(bz-0.3)*1000.))

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

h5f.close()


F = pl.figure()
f = F.add_subplot(121)
maxV = np.max(az)
minV = np.min(az)
C = (az-minV)/(maxV-minV)

print(np.max(C))
print(np.min(C))
for i in range(minL):
    c=C[i]
    f.plot(ax[i],ay[i],'o',color=(c,c,c))
f.set_aspect(1)
F.savefig('test.pdf')

pl.show()
