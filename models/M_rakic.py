import sys
from tools import *

if(len(sys.argv)<5):
    print("Specify: Nnodes Nbatch Nsims timesteps")

Nnodes = int(sys.argv[1])
Nbatch = int(sys.argv[2])
Nsims = int(sys.argv[3])
t = int(sys.argv[4])

s = 'configs/maps/rakic/'
src = [ 'configs/configRakic.json',
        s+'15_Fgf8.h5',
        s+'8_Emx2.h5',
        s+'6_Couptf1.h5',
        s+'28_Sp8.h5',
        s+'23_Pax6.h5']

dst = [ 'config.json',
        'gene1.h5',
        'gene2.h5',
        'gene3.h5',
        'gene4.h5',
        'gene5.h5']

#runCull(Nnodes, Nbatch, Nsims, t, src, dst)
run([7,21],Nbatch,Nsims,t,src,dst)

'''
'1_BMPs.h5'
'2_Cadherin2.h5'
'3_Cadherin6.h5'
'4_Cadherin8.h5'
'5_Clim1a.h5'
'6_Couptf1.h5'
'7_EGFs.h5'
'8_Emx2.h5'
'9_EphA2.h5'
'10_EphA7.h5'
'11_EphrinA5.h5'
'12_Etv1.h5'
'13_Fgf15.h5'
'14_Fgf17.h5'
'15_Fgf8.h5'
'16_Foxg1.h5'
'17_Is2.h5'
'18_Lhx2.h5'
'19_Lmo3.h5'
'20_Lmo4.h5'
'20_Nt3.h5'
'21_Otx2.h5'
'22_P75.h5'
'23_Pax6.h5'
'24_RORbeta.h5'
'25_SCIP.h5'
'26_sfrps.h5'
'27_Sp5.h5'
'28_Sp8.h5'
'29_Tbr1.h5'
'30_Wnts.h5'
'''
