#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Feb 2023 aschiavi 
########################################################################
########################################################################

import argparse



parser = argparse.ArgumentParser(description='crop a subset of voxels of an mhd map using a mask')
parser.add_argument("-i", help="input  map",required=True,nargs=None,metavar='path')
parser.add_argument("-o", help="output map",required=True,nargs=None,metavar='path')
parser.add_argument("-ixi", help="output map",required=True,nargs=None, type=int, metavar='#')
parser.add_argument("-ixf", help="output map",required=True,nargs=None, type=int, metavar='#')
parser.add_argument("-iyi", help="output map",required=True,nargs=None, type=int, metavar='#')
parser.add_argument("-iyf", help="output map",required=True,nargs=None, type=int, metavar='#')
parser.add_argument("-izi", help="output map",required=True,nargs=None, type=int, metavar='#')
parser.add_argument("-izf", help="output map",required=True,nargs=None, type=int, metavar='#')




args = parser.parse_args()
# exit(0)
########################################################################

import os, sys, shutil, time,re
from math import *
import numpy as np
import struct

from mhd_io import *

fin = args.i
fout = args.o
ixi=args.ixi
ixf=args.ixf
iyi=args.iyi
iyf=args.iyf
izi=args.izi
izf=args.izf

[nn,hs,x0,Map] = unpackVoxels(mhd_read(fin))

newMap=Map[ixi:ixf,iyi:iyf,izi:izf].copy()

x0new = x0+hs*[ixi,iyi,izi]
xFnew = x0+hs*[ixf,iyf,izf]

mhd_write(fout,packVoxels(np.array(newMap.shape),hs,x0new,newMap))
print('cropped map written to file:',fout)

print("new_dims: "+str(ixf-ixi)+" "+str(iyf-iyi)+" "+str(izf-izi))
print("new_offset: "+str(x0new[0])+" "+str(x0new[1])+" "+str(x0new[2]))
print("new_end: "+str(xFnew[0])+" "+str(xFnew[1])+" "+str(xFnew[2]))

