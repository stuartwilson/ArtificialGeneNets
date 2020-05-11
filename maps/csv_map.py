import h5py
import numpy as np
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
    # norm the map
    minz = np.min(z)
    norm = 1./(np.max(z)-minz)
    for i in range(len(z)):
        z[i] = (z[i]-minz)*norm

    return x,y,z

'''
    SETUP INPUTS
'''

#'images/hmaps/case_19_18_P36_pup_enucleated_at_P4.csv'
#'images/hmaps/case_19_50_P36_pup_sighted_control.csv'

x,y,f = getData(getImg(argv[1]))

h5f = h5py.File(argv[2],'w')
h5f.create_dataset('X', data=x)
h5f.create_dataset('Y', data=y)
h5f.create_dataset('Z', data=np.zeros(len(x)))
h5f.create_dataset('F', data=f)
h5f.close()
