#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Oct 2016 A.Schiavi
# Sep 2019 adapted from m3d to mhd file format
# Jan 2024 A.Schiavi : moved to SimpleITK IO API
########################################################################
########################################################################
version='1.1'
import sys
if '-V' in sys.argv or '--version' in sys.argv:
    print(f'version: {version}')
    exit(0)
########################################################################
########################################################################
import argparse
parser = argparse.ArgumentParser(description='reports info on 3d map in mhd file format')
parser.add_argument("map",help="list of mhd files",nargs='+',metavar='path')
parser.add_argument("-q", help="short report",action="store_true",default=False)
parser.add_argument("-V","--version",help="version",action="store_true",default=False)

args = parser.parse_args()
# print(args.map)
# exit(0)
########################################################################
########################################################################
import sys
import numpy as np
from FRED_mhd_io import *

for i,f in enumerate(args.map):

    [nn,hs,x0,Map,direction] = readMap(f)

    print('file=',f)
    print('dims=',*nn)
    print('spacing=',*hs)
    print('centroid of first voxel=',*(x0+0.5*hs))
    print('bounding box origin=',*(x0))
    print('transform matrix = ',*direction)
    print('datatype=',Map.dtype)
    print('')
    if args.q:
        continue
    print('x-spatial extent= [',np.float32(x0[0]),',',np.float32(x0[0]+nn[0]*hs[0]),'] => ',np.float32(nn[0]*hs[0]))
    print('y-spatial extent= [',np.float32(x0[1]),',',np.float32(x0[1]+nn[1]*hs[1]),'] => ',np.float32(nn[1]*hs[1]))
    if nn.size == 3:
        print('z-spatial extent= [',np.float32(x0[2]),',',np.float32(x0[2]+nn[2]*hs[2]),'] => ',np.float32(nn[2]*hs[2]))
    print('')
    print('bounding box volume= ',np.float32(np.prod(nn*hs)))
    print('')
    
    print('')
    print('num voxels=',Map.size)
    print('memory occupancy=',np.prod(nn)*Map.dtype.itemsize,'Bytes')
    print('')

    Map=Map.flatten()

    valid = Map[~np.isnan(Map)]

    if len(valid)!=len(Map):
        print('map has NANs')
        print(f'valid values = {len(valid)} = {100.*len(valid)/len(Map):.2f}%')
        print('#')

    print('range=',np.min(valid),np.max(valid))
    print('sum=',np.sum(valid))


    nonzero = (Map!=0).sum()
    print('non-zero voxels=',nonzero,'=',100.0*nonzero/Map.size,' %')
    if i<len(args.map)-1:
        print('\n\n')


