#!/usr/bin/env python3
########################################################################
########################################################################
# CPFR TPS
# Jan 2026  
########################################################################
########################################################################

import argparse

parser = argparse.ArgumentParser(description='find best field size basing on isodose curves in water')
parser.add_argument("-ptv", help="input  map",required=True,nargs=None,metavar='path')
parser.add_argument("-lutpath", help="LUT PATH",required=True,nargs=None,metavar='path')
parser.add_argument("-min_size", help="X Y [cm]", type=float, required=True,nargs=2)
parser.add_argument("-E", help="energy",required=True)


args = parser.parse_args()
# exit(0)
########################################################################

import os, sys, shutil, time,re
from math import *
import numpy as np
import struct

ptv=args.ptv
LUT=args.lutpath
[minXsize, minYsize]=args.min_size
E=args.E

from mhd_io import *
[nn,hs,x0,Map] = unpackVoxels(mhd_read(ptv))
depth=np.sum(np.max(np.max(np.array(Map), axis=0), axis=0))*hs[2]
print("DEPTH [cm]= "+str(depth))
if (depth > 1.7):
  depth=1.7

from pathlib import Path

lut_dir = Path(LUT)

X_vals = np.array([2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4])
Y_vals = np.array([2, 3, 4])

lut = {}

Xidx=np.argmin(np.abs(X_vals-minXsize))
Yidx=np.argmin(np.abs(Y_vals-minYsize))
flag=False
for y in range(Yidx, len(Y_vals)):
  Y=Y_vals[y]
  for x in range(Xidx, len(X_vals)):
    X=X_vals[x]
    fname = f"LUT_{E}MeV_{X}x{Y}cm2"
    fpath = lut_dir / fname

    if not fpath.exists():
      # file mancante  ^f^r lo saltiamo
      continue

    data = np.loadtxt(fpath)
#    print(data)
    if data.ndim != 2 or data.shape[1] != 5:
      raise ValueError(f"{fname} is not shaped Nx5")

    Zidx=np.argmin(np.abs(data[:,0]-depth))
    for z in range(Zidx, len(data[:,0])):
      if (data[z,2]-data[z,1]>=minXsize) and (data[z,4]-data[z,3]>=minYsize):
        print("Xsize: "+str(X))
        print("Ysize: "+str(Y))
#        print(data)
        flag=True
        break
    if flag==True:
      break
  if flag==True:
    break
if flag==False:
  print("Xsize: "+str(4))
  print("Ysize: "+str(4))
#      lut[(E, X, Y] = data
    
