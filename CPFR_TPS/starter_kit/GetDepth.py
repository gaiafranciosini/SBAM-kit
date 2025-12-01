#!/usr/bin/env python3
########################################################################
########################################################################
# Nov 2025 A.Burattini
########################################################################
########################################################################

import argparse
import numpy as np
from mhd_io import *

parser = argparse.ArgumentParser(
    description='Evaluates isodoses and saves them on a mhd map'
)
parser.add_argument("-dose", help="Dose map", required=True)
#parser.add_argument("-xy", help="X and Y size", nargs=2, type=float, required=True)
parser.add_argument("-out",  help="output map", required=True)
parser.add_argument("-preD", help="prescription dose", type=float, default=2000, required=True)
parser.add_argument("-iso", help="isodoses values as percentage of prescription dose", type=float, default=95, required=True)

args = parser.parse_args()

dose = args.dose
out = args.out
preD=args.preD
iso = args.iso
#xy=args.xy

#if (ct_shape[0]<dose_shape[0] or ct_shape[1]<dose_shape[1] or ct_shape[2]<dose_shape[2]):
#  sys.exit("ERROR: dose size smaller than CT size!")

iso=iso/100
[nn, hs, x0, Map]=unpackVoxels(mhd_read(dose))

isoMap=Map
isoMap[Map < iso*preD]=0
isoMap[Map >= iso*preD]=1
#print(np.argwhere(isoMap[:,:,0]))
#isoidx=np.zeros((isoMap.shape[2],2))
#print(isoidx[0])
print("STRUCTURE:")
print("Z= depth [cm]")
print("X: xmin xmax [cm]")
print("Y: ymin ymax [cm]")
print(" ")
for z in range(0, isoMap.shape[2]):
  slice=isoMap[:,:,z]
  isoidx=np.argwhere(slice)
  if isoidx.size > 0:
    print(np.min(isoidx, axis=0))
    xlim=[x0[0]+hs[0]*np.min(isoidx, axis=0)[0], x0[0]+hs[0]*np.max(isoidx,axis=0)[0]]
    ylim=[x0[1]+hs[1]*np.min(isoidx,axis=0)[1], x0[1]+hs[1]*np.max(isoidx, axis=0)[1]]
    print("Z= "+str(x0[2]+z*hs[2]))
    print("X: "+str(xlim[0])+" "+str(xlim[1]))
    print("Y: "+str(ylim[0])+" "+str(ylim[1]))
    print(" ")
mhd_write(out,packVoxels(np.array(isoMap.shape),hs, x0, isoMap))







