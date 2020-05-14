import sys
import networks
import numpy as np
import h5py

######## SPECIFY TRAINING MAPS

#mapDir = '../maps/ucd/maps/'
mapDir = 'tmp/'
mapFiles = [mapDir+'18_147_id2_c.h5',    # CTL ID2
            mapDir+'17_268_id2_c.h5',    # P04 ID2
            mapDir+'18_147_rzrb_c.h5',   # CTL RZRb
            mapDir+'17_268_rzrb_c.h5']   # P04 RZRb


######## CONFIGURE THE NETWORK

config = {}
config['maps'] = []
config['maps'].append({'filename': 'map0.h5', 'outputID': 3, 'contextID': 2, 'contextVal': 0.0})
config['maps'].append({'filename': 'map1.h5', 'outputID': 3, 'contextID': 2, 'contextVal': 1.0})
config['maps'].append({'filename': 'map2.h5', 'outputID': 4, 'contextID': 2, 'contextVal': 0.0})
config['maps'].append({'filename': 'map3.h5', 'outputID': 4, 'contextID': 2, 'contextVal': 1.0})



class EvolveTopology(networks.BatchRun):

    def initGenes(self):
        self.P = int(self.Nsims)
        self.Ngenes = int(self.N**2)
        self.X = np.random.rand(self.P,self.Ngenes)<0.5 # genome
        self.Error = np.zeros(self.P)
        self.Edges = np.zeros(self.P)

        preIndices,postIndices = np.meshgrid(np.arange(int(self.N)),np.arange(int(self.N)))
        self.preIndices = preIndices.flatten()
        self.postIndices = postIndices.flatten()
        self.IndicesMap = np.random.permutation(int(self.N)**2)


    def buildNet(self, pre, post):
        p = self.j
        for i in range(int(self.N)**2):
            if(self.X[p,i]==1):
                pre = np.hstack([pre,self.preIndices[self.IndicesMap[i]]])
                post = np.hstack([post,self.postIndices[self.IndicesMap[i]]])
        return pre, post

        '''
        c = 0
        for i in range(int(self.N)):
            for j in range(int(self.N)):
                if(self.X[p,c]==1):
                    pre = np.hstack([pre,i])
                    post = np.hstack([post,j])
                c=c+1
        '''

    def recombine(self, F):

        lower = np.zeros([self.P])
        upper = np.zeros([self.P])

        cumSum = 0.
        for p in range(self.P):
            lower[p] = cumSum
            cumSum += F[p]
            upper[p] = cumSum

        Mum = np.zeros([self.P,self.Ngenes])
        for p in range(self.P):
            r = np.random.rand()*cumSum
            for q in range(self.P):
                if ((lower[q]<=r) and (r<upper[q])):
                    Mum[p] = self.X[q]

        Dad = np.zeros([self.P,self.Ngenes])
        for p in range(self.P):
            r = np.random.rand()*cumSum
            for q in range(self.P):
                if ((lower[q]<=r) and (r<upper[q])):
                    Dad[p] = self.X[q]

        self.X = Mum
        for p in range(self.P):
            for n in range(self.Ngenes):
                if(n<self.Ngenes/2):
                    self.X[p,n] = Dad[p,n]


    def evolve(self, nGens):
        for i in range(nGens):

            Error,dummy = self.run(overwrite=True)

            Edges = np.sum(self.X*1.0,axis=1)

            #F = 1.0-(Err-np.min(Err))/(np.max(Err)-np.min(Err))
            F = 1.0-Error
            self.recombine(F)
            print('Generation: '+str(i))
            self.Error = np.vstack([self.Error,Error])
            self.Edges = np.vstack([self.Edges,Edges])

            h5f = h5py.File(self.dir+'/evo.h5','w')
            h5f.create_dataset('Error', data=self.Error[1:,:])
            h5f.create_dataset('Edges', data=self.Edges[1:,:])
            h5f.close()







######## CONFIGURE THE RUN

'''
if(len(sys.argv)<6): print("Specify: min-nodes max-nodes num-per-batch num-total-sims timesteps")
minNodes = int(sys.argv[1])
maxNodes = int(sys.argv[2])
Nbatch = int(sys.argv[3])
Nsims = int(sys.argv[4])
T = int(sys.argv[5])
'''

######## CHOOSE ALGORITHM

#Sizes = np.ceil(np.linspace(minNodes,maxNodes,Nsims))
#net = networks.FullyRecurrent(Sizes,Nbatch,T,mapFiles,config)
#net = networks.RandomCull(np.ones(Nsims)*maxNodes,Nbatch,T,mapFiles,config,params=[0.95])

popSize = 500
maxNumGenes = 10
nBatch = 10
nGens = 1000

net = EvolveTopology(np.ones(popSize)*maxNumGenes,nBatch,10000,mapFiles,config)
net.initGenes()
net.evolve(nGens)

######## DO THE RUN


