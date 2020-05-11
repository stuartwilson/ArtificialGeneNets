import h5py
import json
import subprocess
import numpy as np

'''
    GENERAL TOOLS
'''

def waitUntilReady(x):
    ready=False
    while(not ready):
        ready=True
        for i in range(len(x)):
            if(x[i].poll() is None):
                ready=False


'''
    NETWORK TOPOLOGY TOOLS
'''

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


'''
    RUNNING METHODS - PARENT: ALGORITHMS
'''

class BatchRun:
    def __init__(self, Sizes, Nbatch, T, mapFiles, configfilecontents, params='', cmd='recurrentnet',dir='data'):
        self.cmd = cmd
        self.dir = dir
        self.Nbatch = Nbatch
        self.simIndex = 0
        self.Sizes = Sizes
        self.Nsims = len(Sizes)
        self.N = Sizes[0]
        self.T = T
        self.mapFiles = mapFiles
        self.configfilecontents = configfilecontents
        self.params = params
        self.j = 0
        self.k = 0

    def buildNet(self, pre, post):
        pre, post = recur(pre,post,np.arange(self.N))
        return pre, post

    def run(self):

        Seeds = np.arange(self.Nsims)
        finErr = np.zeros(self.Nsims)
        minErr = np.zeros(self.Nsims)

        ######## RUN THE MODEL
        running = True
        while(running):
            P = []
            for i in range(self.Nbatch):
                if not running: break
                self.N=self.Sizes[self.j]
                loc = self.dir+'/expt'+str(self.j)
                subprocess.run('mkdir '+loc,shell=True)
                with open(loc+'/config.json', 'w') as outfile: json.dump(self.configfilecontents, outfile)
                for q in range(len(self.mapFiles)):
                    subprocess.run('cp '+self.mapFiles[q]+' '+loc+'/map'+str(q)+'.h5',shell=True)
                pre = np.array([],dtype=int)
                post = np.array([],dtype=int)
                pre, post = self.buildNet(pre,post)
                h5f = h5py.File(loc+'/network.h5','w')
                h5f.create_dataset('pre', data=pre)
                h5f.create_dataset('post', data=post)
                h5f.close()
                P = np.hstack([P,subprocess.Popen('./../build/'+self.cmd+' '+loc+'  '+str(Seeds[self.j])+' '+str(self.T),shell=True)])
                running=self.j<(self.Nsims-1)
                self.j+=1
            waitUntilReady(P)

            ######## PERFORM ANALYSIS
            for i in range(self.Nbatch):
                if not running: break
                try:
                    loc = self.dir+'/expt'+str(self.k)
                    h5f = h5py.File(loc + '/outputs.h5','r')
                    err = h5f['error'][:]
                    h5f.close()
                    finErr[self.k] = err[-1]
                    minErr[self.k] = np.min(err)
                except:
                    print("no output"+str(self.k)+"\n")

                running=self.k<(self.Nsims-1)
                self.k+=1
                print(self.k)

            ### store results periodically
            h5f = h5py.File(self.dir+'/summary.h5','w')
            h5f.create_dataset('finErr', data=finErr)
            h5f.create_dataset('minErr', data=minErr)
            h5f.close()



'''
    RUNNING METHODS - CHILD: TOPOLOGY VARIATIONS
'''

class FullyRecurrent(BatchRun):
    def buildNet(self, pre, post):
        pre, post = recur(pre,post,np.arange(self.N))
        return pre, post


class RandomCull(BatchRun):
    def buildNet(self, pre, post):
        pre, post = recur(pre,post,np.arange(self.N))
        toCull = int(self.params[0]*(self.j/(self.Nsims-1))*self.N)
        for c in range(toCull):
            pre, post = cullRand(pre,post)
        return pre, post
