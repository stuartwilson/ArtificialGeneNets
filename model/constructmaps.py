import h5py
import numpy as np

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

a = getImg('images/hmaps/case_19_50_P36_pup_sighted_control.csv')
b = getImg('images/hmaps/case_19_18_P36_pup_enucleated_at_P4.csv')

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


h5f = h5py.File('maps/sighted.h5','w')
h5f.create_dataset('X', data=ax)
h5f.create_dataset('Y', data=ay)
h5f.create_dataset('Z', data=ay)
#h5f.create_dataset('Z', data=np.zeros(len(ax)))
h5f.create_dataset('F', data=az)
h5f.close()

h5f = h5py.File('maps/enucleate.h5','w')
h5f.create_dataset('X', data=bx)
h5f.create_dataset('Y', data=by)
#h5f.create_dataset('Z', data=np.zeros(len(bx)))
h5f.create_dataset('Z', data=by)
h5f.create_dataset('F', data=bz)
h5f.close()
