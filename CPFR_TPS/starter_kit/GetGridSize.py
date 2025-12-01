#!/usr/bin/env python3
########################################################################
########################################################################
# Oct 2025 A.Burattini
########################################################################
########################################################################

import argparse

parser = argparse.ArgumentParser(description='get the size of a dose grid keeping the same spacing of the CT')
parser.add_argument("-ct", help="CT map",required=True,nargs=None,metavar='CT')
parser.add_argument("-size", help="grid size [cm]",required=True,nargs=3, type=float)
#parser.add_argument("-o", help="output map",required=True,nargs=None,metavar='out')

args = parser.parse_args()
# exit(0)
########################################################################

import os, sys, shutil, time,re
from math import *
import numpy as np
from mhd_io import *

fin = args.ct
size=args.size
#fout = args.o

[nn,hs,x0,Map] = unpackVoxels(mhd_read(fin))
print(nn, hs, x0)

size_idx=[2*(np.ceil((size[0]/2)/hs[0])),2*(np.ceil((size[1]/2)/hs[1])),(np.ceil(size[2]/hs[2]))]
print("size_idx: "+str(size_idx[0])+" "+str(size_idx[1])+" "+str(size_idx[2]))


coordsZC=[-hs[0]*size_idx[0]*0.5, hs[0]*size_idx[0]*0.5, -hs[1]*size_idx[1]*0.5, hs[1]*size_idx[1]*0.5, 0, hs[2]*size_idx[2]]
print("grid_coord: "+str(coordsZC[0])+" "+str(coordsZC[1])+" "+str(coordsZC[2])+" "+str(coordsZC[3])+" "+str(coordsZC[4])+" "+str(coordsZC[5])+" ")



idxXc=nn[0]/2
idxYc=nn[1]/2
print("centerXY: "+str(idxXc)+" "+str(idxYc))


Xin=x0[0]+(hs[0]*(idxXc-(size_idx[0]/2)))
Xf=x0[0]+(hs[0]*(idxXc+(size_idx[0]/2)))
Yin=x0[1]+(hs[1]*(idxYc-(size_idx[1]/2)))
Yf=x0[1]+(hs[1]*(idxYc+(size_idx[1]/2)))
Zin=0
Zf=hs[2]*size_idx[2]


P=[np.arange(0, nn[0]), np.arange(0, nn[1]), np.arange(0, nn[2])]
Xcoord=x0[0]+hs[0]*np.array(P[0])
Ycoord=x0[1]+hs[1]*np.array(P[1])
Zcoord=x0[2]+hs[2]*np.array(P[2])

Xin=Xcoord[np.argmin(abs(Xcoord-Xin))]
Xf=Xcoord[np.argmin(abs(Xcoord-Xf))]
Yin=Ycoord[np.argmin(abs(Ycoord-Yin))]
Yf=Ycoord[np.argmin(abs(Ycoord-Yf))]
Zin=Zcoord[np.argmin(abs(Zcoord-Zin))]
Zf=Zcoord[np.argmin(abs(Zcoord-Zf))]


print("grid_coord_inREF: "+str(Xin)+" "+str(Xf)+" "+str(Yin)+" "+str(Yf)+" "+str(Zin)+" "+str(Zf)+" ")
print("CHECK: hs[0]: "+str(hs[0])+" "+str((Xf-Xin)/size_idx[0])+" hs[1]: "+str(hs[1])+" "+str((Yf-Yin)/size_idx[1])+"hs[2]: "+str(hs[2])+" "+str((Zf-Zin)/size_idx[2]))

