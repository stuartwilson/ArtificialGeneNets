import sys
import networks
import numpy as np

######## SPECIFY TRAINING MAPS

mapDir = '../maps/ucd/maps/'
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


######## CONFIGURE THE RUN

if(len(sys.argv)<6): print("Specify: min-nodes max-nodes num-per-batch num-total-sims timesteps")
minNodes = int(sys.argv[1])
maxNodes = int(sys.argv[2])
Nbatch = int(sys.argv[3])
Nsims = int(sys.argv[4])
T = int(sys.argv[5])


######## CHOOSE ALGORITHM

Sizes = np.ceil(np.linspace(minNodes,maxNodes,Nsims))
net = networks.FullyRecurrent(Sizes,Nbatch,T,mapFiles,config)
#net = networks.RandomCull(np.ones(Nsims)*maxNodes,Nbatch,T,mapFiles,config,params=[0.95])


######## DO THE RUN

net.run()
