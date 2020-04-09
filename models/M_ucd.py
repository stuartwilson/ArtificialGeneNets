import sys
from tools import *

if(len(sys.argv)<5):
    print("Specify: Nnodes Nbatch Nsims timesteps")

Nnodes = int(sys.argv[1])
Nbatch = int(sys.argv[2])
Nsims = int(sys.argv[3])
t = int(sys.argv[4])

s = 'configs/maps/ucd/'
src = [ 'configs/config1out1ctxt.json',
        s+'19_46_c.h5',
        s+'17_295_c.h5']

dst = [ 'config.json',
        'ctrl.h5',
        'expt.h5']

#runCull(Nnodes, Nbatch, Nsims, t, src, dst)
run([5,5],Nbatch,Nsims,t,src,dst)

#18_160 and 19_47 were ok!

'''
17_290_c.h5 #   P12enuc
17_291_c.h5 #   ????
17_295_c.h5 #   P12enuc
17_299_c.h5 #   control **      ##
18_159_c.h5 #   P4enuc
18_160_c.h5 #   control **
19_18_c.h5  #   P4enuc **       ##
19_47_c.h5  #   P4enuc **
'''
