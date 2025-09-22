#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Oct 2016 A.Schiavi
# Sep 2019 adapted from m3d to mhd file format
########################################################################
########################################################################
import sys
import numpy as np
from mhd_io import *

########################################################################
import argparse
parser = argparse.ArgumentParser(description='reports info on 3d map in mhd file format')
parser.add_argument("map",help="list of mhd files",nargs='+',metavar='path')
parser.add_argument("-v", help="verbose report",action="store_true",default=False)
args = parser.parse_args()
# print(args.map)
# exit(0)

files = args.map

for i,f in enumerate(files):

    voxels = mhd_read(f)
    [nn,hs,x0,Map] = unpackVoxels(voxels)

    print('file=',f)
    print('dims=',*nn)
    print('spacing=',*hs)
    print('offset=',*x0)
    print('transform matrix = ',*voxels['TransformMatrix'])
    print('datatype=',Map.dtype)
    print('')
    if not args.v:
        continue
    print("x-spatial_extent: "+str(x0[0])+" "+str(x0[0]+nn[0]*hs[0]))
    print("y-spatial_extent: "+str(x0[1])+" "+str(x0[1]+nn[1]*hs[1])) 
    if nn.size == 3:
        print("z-spatial_extent: "+str(x0[2])+" "+str(x0[2]+nn[2]*hs[2]))
    print('')
    print('bounding box volume= ',np.prod(nn*hs))
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
    if i<len(files)-1:
        print('\n\n')


