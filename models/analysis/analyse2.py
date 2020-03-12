import h5py
import numpy as np
import pylab as pl
from scipy import signal
import sys

logdir = sys.argv[1]
print(logdir)

h5f = h5py.File(logdir + '/outputs.h5','r')
err = h5f['error'][:]
h5f.close()


h5f = h5py.File(logdir + '/weights.h5','r')
W = h5f['weightmat'][:]
n = int(np.sqrt(W.shape[0]))
W = W.reshape([n,n])
h5f.close()

noConn = np.where(W==0)




Wimg = W.copy()
Wimg[noConn] = np.NaN

F = pl.figure(figsize=(10,6))
f = F.add_subplot(121)
g = f.pcolor(Wimg,cmap='jet')
F.colorbar(g)
f.set_xlabel('pre')
f.set_ylabel('post')
f.set_aspect(1.0)


t = 2.*np.pi*np.arange(n)/(n-1)
offset = 0.0
head = 0.82
tail = 0.76
f = F.add_subplot(122)
for i in range(n):
    for j in range(n):
        if(~(W[i,j]==0)):
            if(W[i,j]>0):
                col = 'r'
            if(W[i,j]<0):
                col = 'b'

            xorg = np.cos(t[i])+offset*np.cos(t[i]+np.pi*0.5)
            yorg = np.sin(t[i])+offset*np.sin(t[i]+np.pi*0.5)
            xdst = np.cos(t[j])+offset*np.cos(t[i]+np.pi*0.5)
            ydst = np.sin(t[j])+offset*np.sin(t[i]+np.pi*0.5)
            phi = np.arctan2(ydst-yorg,xdst-xorg)
            rho = np.sqrt((ydst-yorg)**2+(xdst-xorg)**2)
            f.plot(
                [xorg,xorg+rho*np.cos(phi)],
                [yorg,yorg+rho*np.sin(phi)],'--',color='k',linewidth=0)
            f.plot(
                [xorg+head*rho*np.cos(phi),xorg+rho*np.cos(phi)],
                [yorg+head*rho*np.sin(phi),yorg+rho*np.sin(phi)],color=col,linewidth=3)
            f.plot(
                [xorg+tail*rho*np.cos(phi),xorg+rho*np.cos(phi)],
                [yorg+tail*rho*np.sin(phi),yorg+rho*np.sin(phi)],color=col,linewidth=1)

A = [r'$x$',r'$y$','output','context']
for i in range(n):
    col = 'k'
    if(W[-1,i]>0):
        col = 'r'
    if(W[-1,i]<0):
        col = 'b'
    f.plot(np.cos(t[i]),np.sin(t[i]),'o',color=col,markersize=20)
    f.plot(np.cos(t[i]),np.sin(t[i]),'o',color='w',markersize=15)
    if(i<len(A)):
        f.annotate(A[i],[np.cos(t[i])*1.15,np.sin(t[i])*1.15],fontsize=15)

f.axis(np.array([-1,1,-1,1])*1.3)
f.set_aspect(1.0)
f.set_xticks([])
f.set_yticks([])

pl.savefig(logdir+'/weights.png')

pl.show()
