#!/usr/bin/env python3
######################################################################
# Author:       A. Schiavi
# Date:         2020 - 2023
# Purpose:      routines for importing DICOM files in Fred
# Source:       python
######################################################################
import argparse
parser = argparse.ArgumentParser(description='Export ROIs in dicom RTStruct using the phantom grid')
parser.add_argument("-CT",help="CT mhd map",nargs=None,metavar='path',required=True)
parser.add_argument("-rtstruct", help="dicom file(s) containing RTStruct",nargs='+',required=True)
parser.add_argument("-l", help="list ROIs found in DICOM files",action="store_true",default=False)

parser.add_argument("-include", help="substrings of names to be included",nargs='+',metavar='str')
parser.add_argument("-exclude", help="substrings of names to be excluded",nargs='+',metavar='str')

parser.add_argument("-format", help="output file format",choices=['mhd','mha'],default='mhd')
parser.add_argument("-z", help="use compression",action="store_true",default=False)


args = parser.parse_args()
# print(args) ; exit()
######################################################################

import os, sys
scriptDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,scriptDir)

from math import *
import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer
#import pydicom as dicom
import pydicom
import matplotlib
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.cm as cm

matplotlib.use('TkAgg')
import pylab as plt

import SimpleITK as sitk
from mhd_io import *

######################################################################
# global vars
dicoms=[]
cts=[]
cts=[]
rtstruct=[]
rtplan=[]
rtdose=[]
ISO = np.zeros(3)
ROIs = []
nx=ny=nz=0
hx=hy=hz=0
x0=y0=z0=0
Ox=Oy=Oz=0
######################################################################
def banner():
    print(' ')
    print('\tDICOM importer using python libraries for Fred project')
    print(' ')
#################### IMPORT & CHECK DICOM ##########################
#   load fill contour dynamic library
lib = ctypes.cdll.LoadLibrary(os.path.join(scriptDir,"./libFillContours.so"))
print("SCRIPT DIR:: ",scriptDir);
fillContourC = lib.fillContour_f32_f64
fillContourC.restype = None
fillContourC.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                ctypes.c_int , ctypes.c_int ,
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                ctypes.c_int,ctypes.c_int]

######################################################################
def applyFilter(roiName):
    OK = False
    if args.include:
        for s in args.include:
            if s.lower() in roiName.lower():
                OK = True
    else:
        OK = True

    if args.exclude:
        for s in args.exclude:
            if s.lower() in roiName.lower():
                OK = False
    return OK

########################## BUILD ROIs #####################################

def buildROIDriver():
    print('\nLooking for ROIs:')
    if(len(rtstruct)==0):
        print('none found!')
        return

    volumes = []
    for rts in rtstruct:
        for roi in rts.StructureSetROISequence:
            if not applyFilter(roi.ROIName):
                continue
            type = findInterpretedType(rts,roi.ROINumber)
            print(f'    {roi.ROINumber:4d}: {roi.ROIName:30s} {type:20s}')
            volumes.append(roi)

    if args.l: #listMode
        return

    print('\n')
    for roi in volumes:
        buildROI(rts,roi.ROINumber)

######################################################################
def polygonArea(X,Y):
    area = 0
    npt = len(X)
    j = npt-1
    for i in range(npt):
        area += (X[j]+X[i]) * (Y[j]-Y[i])
        j = i
    return area/2;

######################################################################
def buildROI(rts,roiNum):
    global Phantom

    # find ROI
    for obj in rts.StructureSetROISequence:
            if obj.ROINumber == roiNum:
                roi = obj

    # find contour sequence
    for contseq in rts.ROIContourSequence:
        if contseq.ReferencedROINumber==roiNum:
            seq = contseq
    if 'ContourSequence' in seq:
#    if findInterpretedType(rts,roi.ROINumber) == 'GTV':
#        print('GTV ROI num %d, name %s' % (roiNum, roi.ROIName))
#        return;
#    else:
        print('ROI num %d, name %s, num of contours %d' % (roiNum, roi.ROIName, len(seq.ContourSequence)))
    else:
        return 

    ROI = np.zeros_like(Phantom,dtype=np.float32)

    #mask=np.zeros((ny,nx),dtype=np.float32)
    mask=np.zeros((ny,nx),dtype=np.float32)


    myCont = []
    
    for cont in seq.ContourSequence:
        print('.',end='',flush=True)
        xContData = np.array(cont.ContourData[0::3])
        if(xContData[0] < x0): 
            Xcont = (+np.array(cont.ContourData[0::3]) + x0) / hx
            Ycont = (+np.array(cont.ContourData[1::3]) + y0) / hy
        else:
            Xcont = (+np.array(cont.ContourData[0::3]) - x0) / hx
            Ycont = (+np.array(cont.ContourData[1::3]) - y0) / hy
            
#        print("?? ",nx,ny,nz,x0,y0,hx,hy,sep='\t')
#        print("Cont ",np.array(cont.ContourData[0::3]),sep='\t')
#        print("Cont - offX ",Xcont,sep='\t')
#        print("Cont - offY ",Ycont,sep='\t')
        
        lev=5 # oversampling level for voxels on the boundary (1-10)
        fillContourC(mask,ny,nx,Xcont,Ycont,len(Xcont),lev)
        mask = np.abs(mask)

        area = abs(polygonArea(Ycont,Xcont))
        relerr = (np.sum(mask)-area)/area*100.
#        print(area,np.sum(mask),relerr,sep='\t')
#         print(relerr,sep='\t')

        for i in range(nz):
            if abs(Oz+i*hz-cont.ContourData[2])<=hz/2:
                # drawSliceContoursMask(i,mask,Xcont,Ycont)
                # exit(0)
                ROI[i,:,:]+=mask
    
    print('')
    if(np.max(ROI)>2):
        print('error: too many contours overlapping!!!')
        print('check contour definition for ROI ',roi.ROIName)
        exit(1)

    ROI[ROI>1]=2-ROI[ROI>1] # apply rule for subtracting layers (Policlinico)

    print("Numero di voxel diversi da zero nella matrice mask:", np.count_nonzero(ROI))
    name = '_'.join(roi.ROIName.split())
    img_out = sitk.GetImageFromArray(ROI)
    print("Numero di voxel diversi da zero nell'immagine img_out:", sitk.GetArrayViewFromImage(img_out).astype(bool).sum())

    img_out.SetOrigin(imgCT.GetOrigin())
    img_out.SetSpacing(imgCT.GetSpacing())
    img_out.SetDirection(imgCT.GetDirection())
    sitk.WriteImage(img_out,name+'.'+args.format,useCompression = args.z)
    # exit(0)

    roiVolume = np.sum(ROI)*hx*hy*hz
    ROIs.append((name,roiVolume))

##########################################################################################
def drawSliceContoursMask(iz,mask,Xcont,Ycont):
    global nn,hs,x0,Phantom

    plt.ion()
    plt.clf()
    plt.title('Slice %d' % iz)

    Extent = (0,mask.shape[1],0,mask.shape[0])

    slice = Phantom[iz,:,:]

    imCT = plt.imshow(slice, interpolation='none', origin='lower',
                cmap=cm.gray,aspect=1,extent=Extent)        
    im = plt.imshow(mask.squeeze(), interpolation='none', origin='lower',
                cmap=cm.cool,aspect=1,extent=Extent,alpha=0.99)     
    for i in range(mask.shape[0]):
        mask[i,i%2::2]=2
    im2 = plt.imshow(mask.squeeze(), interpolation='nearest', origin='lower',
                cmap=cm.gray,aspect=1,alpha=0.1,extent=Extent)      
            
    plt.fill(Xcont,Ycont,edgecolor='k',label='contour',lw=1,fill=False)
    plt.plot(Xcont[0],Ycont[0],'ro',label='first',lw=1)
    plt.plot(Xcont[1:-1],Ycont[1:-1],'ko',markersize=4)
    plt.plot(Xcont[-1],Ycont[-1],'bo',label='last',lw=1)
    print('start',Xcont[0],Ycont[0])
    Dx = np.max(Xcont)-np.min(Xcont)
    Dy = np.max(Ycont)-np.min(Ycont)
    frac = 0.2
    plt.xlim((np.min(Xcont)-frac*Dx,np.max(Xcont)+frac*Dx))
    plt.ylim((np.min(Ycont)-frac*Dy,np.max(Ycont)+frac*Dy))
    plt.legend()
    plt.show()
    input('Hit return')
##########################################################################################
def findInterpretedType(rts,roiNum):
    for obs in rts.RTROIObservationsSequence:
        if obs.ReferencedROINumber == roiNum:
            return obs.RTROIInterpretedType
    return 'UNKNOWN'
##########################################################################################
##########################################################################################
##########################################################################################

banner()

imgCT   = sitk.ReadImage(args.CT)
Phantom = sitk.GetArrayFromImage(imgCT)
nx,ny,nz = imgCT.GetSize()
hx,hy,hz = imgCT.GetSpacing()
Ox,Oy,Oz = imgCT.GetOrigin()
x0=Ox-hx/2
y0=Oy-hy/2
z0=Oz-hz/2
print('# CT info:')
print('# dims    =',nx,ny,nz)
print('# spacing =',hx,hy,hz)
print('# offset  =',Ox,Oy,Oz)
print('# lowest corner =',x0,y0,z0)
# exit(0)

for fname in args.rtstruct:
    ds = pydicom.dcmread(fname) or die
    if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.3':   # Radiation Therapy Structure Set Storage
        rtstruct.append(ds)
    else:
        print('error: file is not an RTStruct',fname)
        exit(1)

buildROIDriver()

if not args.l:
    print('------------------------------')
    print('%-20s%-10s' %('Name','Volume (cm3)'))
    for r in ROIs:
        print('%-20s%-10f' % (r[0],r[1]/1000.))

