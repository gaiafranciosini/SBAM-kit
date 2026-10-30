#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Oct 2015 A.Schiavi
# Sep 2019 A.Schiavi ported from map3d to mhd 
########################################################################
########################################################################
import os, sys, shutil, time,re
from math import *
import numpy as np
import struct

from FRED_mhd_io import *

########################################################################
MET_TYPES = ['MET_DOUBLE','MET_FLOAT','MET_INT','MET_LONG','MET_SHORT','MET_CHAR','MET_UCHAR','MET_UINT','MET_ULONG','MET_USHORT']
NUMPY_DTYPES = [np.float64,np.float32,np.int32,np.int64,np.int16,np.int8,np.uint8,np.uint32,np.uint64,np.uint16]
STRING_DTYPES = ['float64','float32','int32','int64','int16','int8','uint8','uint32','uint64','uint16']
########################################################################

if len(sys.argv)<3:
    print('usage: mapFile newDatatype')
    print('convert a whole map to a new datatype')
    print('available datatypes are:',' '.join(STRING_DTYPES))
    exit(1)

newtype = sys.argv[2]
if newtype not in STRING_DTYPES:
    print('required datatype not allowed.')
    print('available datatypes are:',' '.join(STRING_DTYPES))
    exit(1)


nn,hs,x0,Map,direction = readMap(sys.argv[1])

print('fname=',sys.argv[1])

print('dims=',nn)
print('spacing=',hs)
print('offset=',x0)
print('x-spatial extent= [',x0[0],',',x0[0]+nn[0]*hs[0],'] => ',nn[0]*hs[0])
print('y-spatial extent= [',x0[1],',',x0[1]+nn[1]*hs[1],'] => ',nn[1]*hs[1])
if nn.size==3:
    print('z-spatial extent= [',x0[2],',',x0[2]+nn[2]*hs[2],'] => ',nn[2]*hs[2])
print()
print('old datatype=',Map.dtype )
print('range=',np.min(Map),np.max(Map))
print('sum=',np.sum(Map))

# save a copy
writeMap(sys.argv[1].rsplit('.', 1)[0]+'_orig.mhd',hs,x0,Map,direction)

Map = Map.astype(NUMPY_DTYPES[STRING_DTYPES.index(newtype)])

print()
print('new datatype=',Map.dtype )
print('new range=',np.min(Map),np.max(Map))
print('new sum=',np.sum(Map))

writeMap(sys.argv[1],hs,x0,Map,direction)

