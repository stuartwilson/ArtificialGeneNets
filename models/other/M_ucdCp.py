import sys
#from tools import *
from tools import recur, waitUntilReady
import h5py
import subprocess
import numpy as np


def run(Nnodes, Nbatch, Nsims, t, src, dst, loc='data'):

    Seeds = np.arange(Nsims)
    Sizes = np.linspace(Nnodes[0],Nnodes[1],Nsims)
    Sizes = np.ceil(Sizes)

    finErr = np.zeros(Nsims)
    minErr = np.zeros(Nsims)

    ######## RUN THE MODEL
    j = 0
    k = 0
    running = True
    while(running):
        P = []
        for i in range(Nbatch):
            if not running: break
            dir = loc+'/expt'+str(j)
            subprocess.run('mkdir '+dir,shell=True)
            for q in range(len(src)):
                subprocess.run('cp '+src[q]+' '+dir+'/'+dst[q],shell=True)
            ### NETWORK SPEC
            N = Sizes[j]
            Narr = np.array([N],dtype=int)
            pre = np.array([],dtype=int)
            post = np.array([],dtype=int)
            pre, post = recur(pre,post,np.arange(N))
            h5f = h5py.File(dir+'/network.h5','w')
            h5f.create_dataset('pre', data=pre)
            h5f.create_dataset('post', data=post)
            h5f.close()
            P = np.hstack([P,subprocess.Popen('./../build/recurrentnet '+dir+'  '+str(Seeds[j])+' '+str(t),shell=True)])
            running=j<(Nsims-1)
            j+=1
        waitUntilReady(P)

        ######## PERFORM ANALYSIS
        for i in range(Nbatch):
            if not running: break
            try:
                dir = loc+'/expt'+str(k)
                h5f = h5py.File(dir + '/outputs.h5','r')
                err = h5f['error'][:]
                h5f.close()
                finErr[k] = err[-1]
                minErr[k] = np.min(err)
            except:
                print("no output"+str(k)+"\n")

            running=k<(Nsims-1)
            k+=1
            print(k)

        ### store results periodically
        h5f = h5py.File(loc+'/summary.h5','w')
        h5f.create_dataset('finErr', data=finErr)
        h5f.create_dataset('minErr', data=minErr)
        h5f.close()

if(len(sys.argv)<6):
    print("Specify: MINnodes MAXnodes Nbatch Nsims timesteps")

minNodes = int(sys.argv[1])
maxNodes = int(sys.argv[2])
Nbatch = int(sys.argv[3])
Nsims = int(sys.argv[4])
t = int(sys.argv[5])

s = 'configs/maps/ucd/'

src = ['configs/config2out1ctxt.json',
        s+'18_147_id2_c.h5',    # CTL ID2
        s+'17_268_id2_c.h5',    # P04 ID2
        s+'18_147_rzrb_c.h5',   # CTL RZRb
        s+'17_268_rzrb_c.h5']   # P04 RZRb

'''
src = [ 'configs/config2out1ctxt.json',
        s+'17_280_id2_c.h5',
        s+'17_268_id2_c.h5',
        s+'17_280_rzrb_c.h5',
        s+'17_268_rzrb_c.h5']
'''

dst = [ 'config.json',
        'ctrl1.h5',
        'expt1.h5',
        'ctrl2.h5',
        'expt2.h5']

#runCull(Nnodes, Nbatch, Nsims, t, src, dst)
run([minNodes,maxNodes],Nbatch,Nsims,t,src,dst)

#18_160 and 19_47 were ok!

'''
??? 17_295.h5
CTL 17_291.h5
CTL 17_299.h5
CTL 18_160.h5
CTL 19_46.h5
CTL 19_50.h5
P04 18_159.h5
P04 19_18.h5
P04 19_47.h5
P12 17_290.h5
P12 17_292.h5
P12 17_294.h5

P04 17_268_id2_c.h5
P04 17_268_rzrb_c.h5
CTL 17_280_id2_c.h5
CTL 17_280_rzrb_c.h5
CTL 18_147_id2_c.h5
CTL 18_147_rzrb_c.h5
P04 18_148_id2_c.h5
P04 18_148_rzrb_c.h5
'''
