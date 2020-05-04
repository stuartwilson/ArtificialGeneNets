import h5py
import numpy as np
import pylab as pl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy import signal
import sys

def getColor(x):
    if (x==0):
        return (0,1,0)
    if(x<0):
        return (0,0,1)
    if(x>0):
        return (1,0,0)


fs = 15

logdir = sys.argv[1]
print(logdir)

h5f = h5py.File(logdir + '/outputs.h5','r')
err = h5f['error'][:]
h5f.close()

h5f = h5py.File(logdir + '/weights.h5','r')
W = h5f['weightmat'][:]
n = int(np.sqrt(W.shape[0])) # Note that with bias weights this will be #nodes+1
W = W.reshape([n,n]).T
h5f.close()

noConn = np.where(W==0)

Wimg = W.copy()
Wimg[noConn] = np.NaN

#print(Wimg)

F = pl.figure(figsize=(15,5))
f = F.add_subplot(131)
f.plot(np.arange(len(err))/10.,err,color='k',linewidth=2)
#f.plot(err2,linewidth=3)
f.set_aspect(np.diff(f.get_xlim())/np.diff(f.get_ylim()))
f.set_xlabel('trials ($x10^4$)',fontsize=fs)
f.set_ylabel('error',fontsize=fs)

f = F.add_subplot(132)
g = f.pcolor(Wimg,cmap='jet')
c = F.colorbar(g,fraction=0.03, pad=0.04,ticks=[np.min(W),0,np.max(W)])
c.ax.set_yticklabels(['{0:.2f}'.format(np.min(W)),'0.00','{0:.2f}'.format(np.max(W))])

f.set_ylim(0,n-1)   # use only when bias included
f.set_yticks(np.arange(0,n-1)+0.5)
f.set_yticklabels(np.arange(1,n))
f.set_xticks(np.arange(0,n)+0.5)
f.set_xticklabels(np.arange(1,n+1))
labels = [item.get_text() for item in f.get_xticklabels()]
labels[-1]='bias'
f.set_xticklabels(labels)

f.set_xlabel('pre',fontsize=fs)
f.set_ylabel('post',fontsize=fs)
f.set_aspect(1.0)

#pl.savefig(logdir+'/figA.pdf')

P = []
t = 2.*np.pi*(np.arange(n)-2)/(n-1)
offset = 0.0
head = 0.82
tail = 0.76
#F = pl.figure(figsize=(6,6))
f = F.add_subplot(133)

for j in range(n-1):
    for i in range(n-1):
        if(~(W[j,i]==0)):
            xorg = np.cos(t[i])
            yorg = np.sin(t[i])
            xdst = np.cos(t[j])
            ydst = np.sin(t[j])
            f.plot([xorg,xdst],[yorg,ydst],'-',linewidth=1,color='k')


rho = 0.4
wid = 0.04

for i in range(n-1):
    for j in range(n):
        if(~(W[j,i]==0)):
            col = getColor(W[j,i])
            xorg = np.cos(t[i])
            yorg = np.sin(t[i])
            xdst = np.cos(t[j])
            ydst = np.sin(t[j])
            dy = ydst-yorg
            dx = xdst-xorg
            ang = np.arctan2(dy,dx)
            p0x = xdst
            p0y = ydst
            p1x = p0x-rho*np.cos(ang)+wid*np.cos(ang+np.pi*0.5)
            p1y = p0y-rho*np.sin(ang)+wid*np.sin(ang+np.pi*0.5)
            p2x = p0x-rho*np.cos(ang)-wid*np.cos(ang+np.pi*0.5)
            p2y = p0y-rho*np.sin(ang)-wid*np.sin(ang+np.pi*0.5)

            p = Polygon(np.array([[p0x,p0y],[p1x,p1y],[p2x,p2y]]),closed=False, color=col,alpha=1.0,fill=True, edgecolor=None)
            f.add_line(p)
            P.append(p)

A = [r'1. $x$',r'2. $y$','3. context','4. output','5.']
for i in range(n-1):
    col = getColor(W[i,-1])
    f.plot(np.cos(t[i]),np.sin(t[i]),'o',color=col,markersize=20)
    f.plot(np.cos(t[i]),np.sin(t[i]),'o',color='w',markersize=15)
    if(i<len(A)):
        f.annotate(A[i],[np.cos(t[i])*1.4,np.sin(t[i])*1.4],va='center',fontsize=fs)

f.axis(np.array([-1,1,-1,1])*1.4)
f.set_aspect(1.0)
f.set_xticks([])
f.set_yticks([])
f.axis('off')

pl.savefig(logdir+'/fig.pdf')

