#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Nov 2021 aschiavi 
# 
########################################################################
########################################################################

import argparse
parser = argparse.ArgumentParser(description='clamp voxel values in a mhd map to given range')
parser.add_argument("map", help="path of the map")

parser.add_argument("-o", help="output file path",required=True,nargs=None,metavar='path')

parser.add_argument("-vmin", help="lower value",nargs=None,type=float,metavar='#')
parser.add_argument("-vmax", help="upper value",nargs=None,type=float,metavar='#')


args = parser.parse_args()
# print(args.map)
# exit(0)
########################################################################

import os, sys, shutil, time,re
from math import *
import numpy as np
import struct

from mhd_io import *


fname = args.map
voxels = mhd_read(fname)

print('fname:',fname)

fout = args.o


[nn,hs,x0,Map] = unpackVoxels(voxels)


newMap = Map.copy()

if args.vmin is not None:
	newMap[np.where(newMap<args.vmin)] = args.vmin

if args.vmax is not None:
	newMap[np.where(newMap>args.vmax)] = args.vmax

print('original range = ',np.min(Map),np.max(Map))
print('clamped  range = ',np.min(newMap),np.max(newMap))

voxels['Map']=newMap
mhd_write(fout,voxels)
print('mask written to file:',fout)
