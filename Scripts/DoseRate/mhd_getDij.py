#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Nov 2016 A.Schiavi
# Sep 2019 A.Schiavi ported from map3d to mhd 
########################################################################
########################################################################
import os, sys, shutil, time,re
from math import *
import numpy as np
import struct

from mhd_io import *

########################################################################

if len(sys.argv)<4:
    print('usage: Dij.bin fieldID pbID')
    print('extracts the dose map for a given pencilbeam from compressed Dij file')
    print('The pencilbeam is identified by its pbID and the referenced fieldID')
    print('e.g. 3 1 identifies pb 1 of field 3')
    print('if pbID is -1 => all pencilbeams of a given field')
    print('if fieldID is -1 => all fields')
    print('therefore you get the whole reconstructed dose map by asking -1 -1')
    exit(1)
    

fname = sys.argv[1]

try:
    FID=int(sys.argv[2])
    print('requested field =',FID)
    PID=int(sys.argv[3])
    print('requested pb =',PID)
except:
    print('Error: cannot parse input parameters')
    exit(-1)

try:
    fin = open(fname,'rb')
except:
    print('IO error: cannot open input file', fname)
    exit(-1)
    

[nx,ny,nz,hx,hy,hz,x0,y0,z0,npb] = struct.unpack('<3i6f1i',fin.read(40)) or die


nn=np.array([nx,ny,nz])
hs=np.around(np.array([hx,hy,hz]),decimals=6)
x0=np.around(np.array([x0,y0,z0]),decimals=6)

print('dims=',nn)
print('spacing=',hs)
# print('offset=',x0)
print('num of pencilbeams=',npb)
found=False

N = np.prod(nn)
M = np.zeros(N,dtype='float32')

while True:
    try:
        [pid,fid,numvxl] = struct.unpack('3i', fin.read(12))
    except:
        break
    print((fid,pid),'nvxl=',numvxl,end='')
    if (pid==PID and fid==FID) or (pid==PID and FID==-1) or (PID==-1 and fid==FID) or (PID==-1 and FID==-1) : 
        found=True
        print(' <-- YES')
        if numvxl>0:
            Ivxl = np.frombuffer(fin.read(numvxl*4),dtype='uint32',count=numvxl)
            Vals = np.frombuffer(fin.read(numvxl*4),dtype='float32',count=numvxl)
            M[Ivxl]+=Vals
    else:
        print(' <-- NO')
        fin.seek(numvxl*8,1) # jump ahead
    if found and not(PID==-1 or FID==-1):
        break

if not found:
    print('Error: could not find pb matching (fieldID,pbID)=',(FID,PID))
    exit(-1)

M = np.reshape(M,nn,order='F')

if PID<0 and FID>=0:
    fname = 'Dose_Field_%d.mhd' % (FID)
    print('Writing dose of al pencilbeams of field %d to %s' % (FID,fname))
if PID>=0 and FID<0:
    fname = 'Dose_AllFields_PB_%04d.mhd' % (PID)    
    print('Writing dose of pb %d for all fields to %s' % (PID,fname))
if PID<0 and FID<0:
    fname = 'Dose_tot'
    print('Writing total reconstructed dose map to %s' % (fname))
if PID>=0 and FID>=0:
    fname = 'Dose_%d_%04d.mhd' % (FID,PID)
    print('Writing dose of pb %d of field %d to %s' % (PID,FID,fname))

mhd_write(fname,packVoxels(nn,hs,x0,M))
