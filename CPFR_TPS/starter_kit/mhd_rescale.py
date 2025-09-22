#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Jul 2020 A.Schiavi
########################################################################
########################################################################
import sys
import numpy as np
from mhd_io import *
########################################################################

if len(sys.argv)<2:
    print('usage: mapList -factor #')
    print('e.g.: map1.mhd map2.mhd map3.mhd -factor 1.602176462E-7')
    print('rescale the values of a list of maps')
    print('rescaling is done in place, i.e. every map is loaded, mulplied by the factor, and saved with the same name')
    print('alternative options:')
    print('\t -multiplier # : map values are muliplied by the multiplies (same as -factor)')
    print('\t -divider # : map values are divided by the divider')
    exit(1)

sys.argv.pop(0)

fac=1.0

nopt=0;
if '-factor' in sys.argv:
    nopt+=1
if '-multiplier' in sys.argv:
    nopt+=1
if '-divider' in sys.argv:
    nopt+=1

if nopt==0:
    print('one options must be present: -factor, -multiplier, -divider')
    exit(1)
if nopt!=1:
    print('options are mutually exclusive: -factor, -multiplier, -divider')
    print('choose one only!')
    exit(1)


if '-factor' in sys.argv:
    try:
        iarg = sys.argv.index('-factor')
        fac = float(sys.argv[iarg+1])
        del sys.argv[iarg:iarg+2]
    except:
        print('error in parsing the rescaling factor')
        exit(1)

if '-multiplier' in sys.argv:
    try:
        iarg = sys.argv.index('-multiplier')
        fac = float(sys.argv[iarg+1])
        del sys.argv[iarg:iarg+2]
    except:
        print('error in parsing the multiplier')
        exit(1)

if '-divider' in sys.argv:
    try:
        iarg = sys.argv.index('-divider')
        div = float(sys.argv[iarg+1])
        if div==0.:
            print('error: divider is zero => this would produce just NANs') 
            exit(1)
        fac  = 1./div
        del sys.argv[iarg:iarg+2]
    except:
        print('error in parsing the divider')
        exit(1)

files = sys.argv

if len(files)<1:
    print('at least one map is needed')
    exit(1)

# print(files)

for i,f in enumerate(files):

    voxels = mhd_read(f)
    [nn,hs,x0,Map] = unpackVoxels(voxels)

    voxels['Map'] = Map*fac

    # mhd_write('test%02d.mhd' % i,voxels)
    mhd_write(f,voxels)
