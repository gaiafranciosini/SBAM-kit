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
if len(sys.argv)<5+1:
    print('usage: srcMapFile dx dy dz dstMapFile')
    print('shifts map by offset vector [dx,dy,dz]')
    exit(1)

[src,dx,dy,dz,dst] = sys.argv[1:]

# sys.argv.pop(0)

voxels = mhd_read(src)
[nn,hs,x0,Map] = unpackVoxels(voxels)

offset = np.array(list(map(float,[dx,dy,dz])))

print('Original map properties:')
print('dims=',nn)
print('spacing=',hs)
print('offset=',x0)
print('x-spatial extent= [',x0[0],',',x0[0]+nn[0]*hs[0],'] => ',nn[0]*hs[0])
print('y-spatial extent= [',x0[1],',',x0[1]+nn[1]*hs[1],'] => ',nn[1]*hs[1])
print('z-spatial extent= [',x0[2],',',x0[2]+nn[2]*hs[2],'] => ',nn[2]*hs[2])
print('datatype=',Map.dtype )
print('')
print('offset vector = ',offset)
print('')
x0 += offset

mhd_write(dst,packVoxels(nn,hs,x0,Map))

print('New map properties:')
print('dims=',nn)
print('spacing=',hs)
print('offset=',x0)
print('x-spatial extent= [',x0[0],',',x0[0]+nn[0]*hs[0],'] => ',nn[0]*hs[0])
print('y-spatial extent= [',x0[1],',',x0[1]+nn[1]*hs[1],'] => ',nn[1]*hs[1])
print('z-spatial extent= [',x0[2],',',x0[2]+nn[2]*hs[2],'] => ',nn[2]*hs[2])
print('datatype=',Map.dtype )