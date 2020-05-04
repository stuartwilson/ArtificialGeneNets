import h5py
import subprocess
import numpy as np

def ffcon(pre, post, input,output):
    for i in input:
        for j in output:
            pre = np.hstack([pre,i])
            post = np.hstack([post,j])
    return pre, post

def recur(pre, post, inds, selfconns=False):
    for i in inds:
        for j in inds:
            if (selfconns or ~(i==j)):
                pre = np.hstack([pre,i])
                post = np.hstack([post,j])
    return pre, post

def cullRand(pre, post):
    if(len(pre)):
        ind = np.floor(np.random.rand()*len(pre))
        pre = np.delete(pre, ind)
        post = np.delete(post, ind)
    return pre, post


def waitUntilReady(x):
    ready=False
    while(not ready):
        ready=True
        for i in range(len(x)):
            if(x[i].poll() is None):
                ready=False

def removeConnectionsTo(pre, post, i):
    k = np.where(post!=i)
    pre = pre[k]
    post = post[k]
    return pre, post

def removeConnectionsFrom(pre, post, i):
    k = np.where(pre!=i)
    pre = pre[k]
    post = post[k]
    return pre, post



def runCull(Nnodes, Nbatch, Nsims, t, src, dst, loc='data', cullMax=0.95):

    Seeds = np.arange(Nsims)
    Sizes = np.ones(Nsims)*Nnodes
    finErr = np.zeros(Nsims)
    minErr = np.zeros(Nsims)
    toCull = np.zeros(Nsims,dtype=int)
    for i in range(Nsims):
        toCull[i] = int(cullMax*(i/(Nsims-1))*Nnodes*Nnodes)

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

            N = Sizes[j]
            Narr = np.array([N],dtype=int)
            pre = np.array([],dtype=int)
            post = np.array([],dtype=int)
            pre, post = recur(pre,post,np.arange(N))

            for c in range(toCull[j]):
                pre, post = cullRand(pre,post)

            h5f = h5py.File(dir+'/network.h5','w')
            h5f.create_dataset('pre', data=pre)
            h5f.create_dataset('post', data=post)
            h5f.close()
            P = np.hstack([P,subprocess.Popen('./../build/src/model '+dir+'  '+str(t)+' '+str(Seeds[j]),shell=True)])
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
            P = np.hstack([P,subprocess.Popen('./../build/src/model '+dir+'  '+str(t)+' '+str(Seeds[j]),shell=True)])
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

def run2(Nnodes, Nbatch, Nsims, t, src, dst, loc='data'):

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
            pre, post = removeConnectionsTo(pre, post, 2)
            pre, post = removeConnectionsFrom(pre, post, 2)
            pre = np.hstack([pre,[2,2]])
            post = np.hstack([post,[3,4]])
            h5f = h5py.File(dir+'/network.h5','w')
            h5f.create_dataset('pre', data=pre)
            h5f.create_dataset('post', data=post)
            h5f.close()
            P = np.hstack([P,subprocess.Popen('./../build/src/model '+dir+'  '+str(t)+' '+str(Seeds[j]),shell=True)])
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
