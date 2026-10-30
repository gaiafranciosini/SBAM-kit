#!/usr/bin/env python3
######################################################################
# Author:               A. Schiavi
# Date:                 2020
# Purpose:              routines for importing DICOM files in Fred
# Source:               python
######################################################################

import os, sys
from math import *
import numpy as np
import pydicom as dicom
import SimpleITK as sitk

######################################################################
# global vars
dicoms = []
cts = []
rtstruct = []
rtplan = []
rtdose = []
ISO = np.zeros(3)

x0 = np.zeros(3)
nn = np.ones(3, dtype=np.int32)
hs = np.ones(3)
CT = np.eye(1)
Phantom = np.eye(1)

######################################################################

def banner():
    print(' ')
    print('\tDICOM importer using python libraries for Fred project')
    print(' ')

#################### IMPORT & CHECK DICOM ##########################

def loadDicoms():
    global NZ_orig, ROIs
    print('looking for dicom files:')
    fileList = sys.argv[1:]
    print(fileList)
    for fname in fileList:
        ds = dicom.dcmread(fname)
        dicoms.append(ds)

    for ds in dicoms:
        if ([0x0008, 0x0016] not in ds):  # check if SOPClassUID exists
            continue
        if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.2':  # CT Image Storage
            cts.append(ds)
        if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.3':   # Radiation Therapy Structure Set Storage
            rtstruct.append(ds)
        if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.8': # Radiation Therapy Ion Plan Storage
            rtplan.append(ds)
        if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.2': # Radiation Therapy Dose Storage
            rtdose.append(ds)

    # report
    print('\tnum of CT images     :', len(cts))
    print('\tnum of RTDose files  :', len(rtdose))
    print('\tnum of RTStruct files:', len(rtstruct))
    print('\tnum of RTPlan files  :', len(rtplan))

################################################################################
def consistencyCheck():
    global nn, hs, x0, CT
    if len(cts) == 0:
        print('Error no CT Image Storage found.')
        exit(-1)

    print('Num of CT slices', len(cts))
    print('Frame of Reference UID', cts[0].FrameOfReferenceUID)
    for ct in cts:
        if ct.FrameOfReferenceUID != cts[0].FrameOfReferenceUID:
            print('Error in Frame of Reference UID -> Slices have different FoR !!!')
            exit(-1)

    print('Series Instance UID', cts[0].SeriesInstanceUID)
    for ct in cts:
        if ct.SeriesInstanceUID != cts[0].SeriesInstanceUID:
            print('Error in Series Instance UID -> Slices from different series !!!')
            exit(-1)

    for ct in cts:
        if ct.ImageOrientationPatient[2] != 0 or ct.ImageOrientationPatient[5] != 0:
            print('Error slices not along z-axis !!!')
            exit(-1)

    # sort slices
    cts.sort(key=lambda ct: ct.ImagePositionPatient[2]) # increasing z
    NZ = len(cts)

    # check slice spacing
    hz = cts[1].ImagePositionPatient[2] - cts[0].ImagePositionPatient[2]

    for i in range(NZ - 1):
        hznow = cts[i+1].ImagePositionPatient[2] - cts[i].ImagePositionPatient[2]
        if (abs(hznow - hz) > 1e-2 * hz):
            print('Error: slices are not equally spaced')
            exit(-1)

    NY = cts[0].Rows
    NX = cts[0].Columns
    for ct in cts:
        if ct.Rows != NY or ct.Columns != NX:
            print('Error in Rows and/or Columns')
            exit(-1)

    hx = cts[0].PixelSpacing[0]
    hy = cts[0].PixelSpacing[1]

    hs = np.array([hx, hy, hz])
    for ct in cts:
        if not np.all(cts[0].PixelSpacing == ct.PixelSpacing):
            print('Error in Voxel dimensions (PixelSpacing + SliceThickness)')
            exit(-1)

    [xoff, yoff, zoff] = cts[0].ImagePositionPatient

    for ct in cts:
        if xoff != ct.ImagePositionPatient[0] or yoff != ct.ImagePositionPatient[1]:
            print('Error in ImagePositionPatient offset')
            exit(-1)

    x0 = np.array([xoff, yoff, zoff]) - hs / 2
    nn = np.array((NX, NY, NZ), dtype=np.int32)
    CT = np.zeros(nn, order='F', dtype=np.int16)

##########################################################################################

def loadCT():
    consistencyCheck()

    print('\nLoad CT:')
    for iz in range(CT.shape[2]):
        print('.', iz, end='', flush=True)
        slice = cts[iz].pixel_array.transpose().astype(float)
        slice *= cts[iz].RescaleSlope
        slice += cts[iz].RescaleIntercept
        CT[:, :, iz] = slice

    print('\nwhole CT shape=', CT.shape)

##########################################################################################

def writePhantom():
    print('\nWrite phantom:')
    print('map dims', Phantom.shape)

    # SimpleITK si aspetta l'array nell'ordine (Z, Y, X)
    phantom_zyx = np.transpose(Phantom, (2, 1, 0))
    image = sitk.GetImageFromArray(phantom_zyx)

    # MODIFICA: Salviamo direttamente in mm (senza / 10.0) per essere 1:1 con l'RTStruct
    spacing_mm = hs.tolist()
    origin_mm = x0.tolist()

    image.SetSpacing(spacing_mm)
    image.SetOrigin(origin_mm)

    print('map spacing', spacing_mm, 'mm')
    print('map offset', origin_mm, 'mm')
    print('saving map to Phantom.mhd')

    sitk.WriteImage(image, 'Phantom.mhd')

##########################################################################################

banner()

loadDicoms()
loadCT()

Phantom = CT  # shallow copy
writePhantom()
