import sys
from tools import *

if(len(sys.argv)<6):
    print("Specify: MINnodes MAXnodes Nbatch Nsims timesteps")

minNodes = int(sys.argv[1])
maxNodes = int(sys.argv[2])
Nbatch = int(sys.argv[3])
Nsims = int(sys.argv[4])
t = int(sys.argv[5])

s = 'configs/maps/ucd/'

src = [ 'configs/config2out1ctxt.json',
        s+'19_46_c.h5', #   ctrl
        s+'19_47_c.h5', #   p4
        s+'19_46_c.h5', #   ctrl
        s+'17_290_c.h5' #   p12
        ]

''' P4:
        s+'18_147_rzrb_c.h5',
        s+'17_268_rzrb_c.h5',
        s+'18_147_id2_c.h5',
        s+'17_268_id2_c.h5'
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

ID2

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
