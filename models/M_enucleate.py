import sys
from tools import *

if(len(sys.argv)<5):
    print("Specify: Nnodes Nbatch Nsims timesteps")

Nnodes = int(sys.argv[1])
Nbatch = int(sys.argv[2])
Nsims = int(sys.argv[3])
t = int(sys.argv[4])

ins = ['configs/config1out1ctxt.json','configs/maps/csv/sighted.h5','configs/maps/csv/enucleate.h5']
outs = ['config.json','ctrl.h5','expt.h5']

runCull(Nnodes, Nbatch, Nsims, t, ins, outs)


