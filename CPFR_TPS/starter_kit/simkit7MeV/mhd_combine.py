#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Oct 2016 A.Schiavi
# Oct 2019 adapted from m3d to mhd file format
########################################################################
########################################################################
import sys
import numpy as np
from mhd_io import *
########################################################################

if len(sys.argv)<2:
    print('usage: mapList -sum -avg -stdev')
    print('computes element-wise sum, average and stdev of a list of maps')
    exit(1)

sys.argv.pop(0)

if '-sum' in sys.argv:
    computeSum = True
    del sys.argv[sys.argv.index('-sum')]
else:
    computeSum = False

if '-avg' in sys.argv:
    computeAvg = True
    del sys.argv[sys.argv.index('-avg')]
else:
    computeAvg = False

if '-stdev' in sys.argv:
    computeStdev = True
    del sys.argv[sys.argv.index('-stdev')]
else:
    computeStdev = False


if computeSum==False and computeStdev==False and computeAvg==False:
    print('you must choose at least one of the following operations: -sum, -avg or -stdev')
    exit(1)

files = sys.argv

if len(files)<2:
    print('at least two maps are needed for combine operations')
    exit(1)

voxels = mhd_read(files.pop(0))
[nn,hs,x0,Map] = unpackVoxels(voxels)

s0=1.0
s1=Map.copy()
s2=Map*Map

for i,f in enumerate(files):

    [nnB,hsB,x0B,MapB] = unpackVoxels(mhd_read(f))

    if not np.all(nn==nnB):
        print('error: dimensions do not match')
        exit(1)

    if not np.all(hs==hsB):
        print('error: spacing does not match')
        exit(1)

    if not np.all(x0==x0B):
        print('error: offset does not match')
        exit(1)

    s0+=1.0
    s1+=MapB
    s2+=MapB*MapB

if computeSum:
    voxels['Map']=s1.copy()
    mhd_write('sum',voxels)

if computeAvg:
    voxels['Map']=s1.copy()/s0
    mhd_write('avg',voxels)

if computeStdev:
    Map=(s2-s1*s1/s0)/(s0-1.0);
    Map=np.sqrt(Map)
    voxels['Map']=Map
    mhd_write('stdev',voxels)
