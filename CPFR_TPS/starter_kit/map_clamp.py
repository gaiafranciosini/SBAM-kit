#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Mar 2022 aschiavi 
# Jan 2024 A.Schiavi : moved to SimpleITK IO API
########################################################################
########################################################################
########################################################################
import argparse
parser = argparse.ArgumentParser(description='clamp voxel values to given range')
parser.add_argument("map", help="path of the map")

parser.add_argument("-o", help="output file path",nargs=None,metavar='path')

parser.add_argument("-vmin", help="lower value",nargs=None,type=float,metavar='#')
parser.add_argument("-vmax", help="upper value",nargs=None,type=float,metavar='#')

args = parser.parse_args()
# print(args.map)
# exit(0)
########################################################################
import os, sys, shutil, time,re
from math import *
import numpy as np
from FRED_mhd_io import *

fname = args.map
print('fname:',fname)

[nn,hs,x0,Map,direction] = readMap(fname)

vmin = np.min(Map)
vmax = np.max(Map)
print(f'original range = {vmin} {vmax}')

if args.vmin is None and args.vmax is None:
    print('no range specified -> nothing to do -> map unchanged')
    exit(0)

if args.vmin is not None:
    vmin = np.array([args.vmin],dtype=Map.dtype)[0] # on-the-fly type conversion

if args.vmax is not None:
    vmax = np.array([args.vmax],dtype=Map.dtype)[0]

print(f'clamped  range = {vmin} {vmax}')

newMap = Map.copy()
newMap[np.where(Map<vmin)] = vmin
newMap[np.where(Map>vmax)] = vmax

fout = args.o if args.o != None else fname

writeMap(fout,hs,x0,newMap,direction)
print('clamped map written to file:',fout)

