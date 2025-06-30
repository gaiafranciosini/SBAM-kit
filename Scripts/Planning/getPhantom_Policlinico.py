#!/usr/bin/env python3
######################################################################
# Author:		A. Schiavi
# Date:			2020
# Purpose:		routines for importing DICOM files in Fred
# Source:		python
######################################################################

import os, sys
from math import *
import numpy as np
import pydicom as dicom
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

x0=np.zeros(3)
nn=np.ones(3,dtype=np.int32)
hs=np.ones(3)
CT=np.eye(1)
Phantom=np.eye(1)

######################################################################

def banner():
	print(' ')
	print('\tDICOM importer using python libraries for Fred project')
	print(' ')
#################### IMPORT & CHECK DICOM ##########################

def loadDicoms():
	global NZ_orig, ROIs
	print('looking for dicom files:')
	# load dicoms
	fileList = sys.argv[1:]
	# print(fileList)
	for fname in fileList:
		try:
			dicoms.append(dicom.read_file(fname))
		except:
			continue

	for ds in dicoms:
		if ([0x0008,0x0016] not in ds):	 # check if SOPClassUID exists
			continue
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.2':  # CT Image Storage
			cts.append(ds)
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.3':	# Radiation Therapy Structure Set Storage
			rtstruct.append(ds)
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.8': # Radiation Therapy Ion Plan Storage
			rtplan.append(ds)
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.2': # Radiation Therapy Dose Storage
			rtdose.append(ds)

	# report
	print('\tnum of CT images     :',len(cts))
	print('\tnum of RTDose files  :',len(rtdose))
	print('\tnum of RTStruct files:',len(rtstruct))
	print('\tnum of RTPlan files  :',len(rtplan))

################################################################################
def consistencyCheck():
	global nn,hs,x0,CT
	# reading dicom error
	if len(cts) == 0:
		print('Error no CT Image Storage found.')
		exit(-1)

	print('Num of CT slices',len(cts))
	# Frame of Reference UID : same for all
	print('Frame of Reference UID',cts[0].FrameOfReferenceUID)
	for ct in cts:
		if ct.FrameOfReferenceUID != cts[0].FrameOfReferenceUID:
			print('Error in Frame of Reference UID -> Slices have different FoR !!!')
			exit(-1)
	# Series Instance UID : same for all
	print('Series Instance UID',cts[0].SeriesInstanceUID)
	for ct in cts:
		if ct.SeriesInstanceUID != cts[0].SeriesInstanceUID:
			print('Error in Series Instance UID -> Slices from different series !!!')
			exit(-1)

	# slices must be along z-axis
	for ct in cts:
		if ct.ImageOrientationPatient[2] != 0 or ct.ImageOrientationPatient[5] != 0: 
			print('Error slices not along z-axis !!!')
			exit(-1)


	# sort slices 
	cts.sort(key = lambda ct : ct.ImagePositionPatient[2]) # increasing z 
	NZ=len(cts)

	# check slice spacing
	hz = cts[1].ImagePositionPatient[2]-cts[0].ImagePositionPatient[2]
	
	for i in range(NZ-1):
		hznow = cts[i+1].ImagePositionPatient[2]-cts[i].ImagePositionPatient[2]
		if (abs(hznow-hz)>1e-2*hz):
			print('Error: slices are not equally spaced')
			exit(-1)

	# Rows and Columns : same for all
	NY= cts[0].Rows
	NX= cts[0].Columns
	for ct in cts:
		if ct.Rows != NY or ct.Columns != NX:
			print('Error in Rows and/or Columns')
			exit(-1)

	# voxel dimensions: same for all
	hx = cts[0].PixelSpacing[0]
	hy = cts[0].PixelSpacing[1]
	

	hs = np.array([hx,hy,hz])
	for ct in cts:
		if not np.all(cts[0].PixelSpacing == ct.PixelSpacing) :
			print('Error in Voxel dimensions (PixelSpacing + SliceThickness)')
			exit(-1)

	# xoffset and yoffset in ImagePositionPatient: same for all 
	[xoff,yoff,zoff] = cts[0].ImagePositionPatient

	
	for ct in cts:
		if xoff != ct.ImagePositionPatient[0] or yoff != ct.ImagePositionPatient[1] :
			print('Error in ImagePositionPatient offset')
			exit(-1)
# 		print 'zslice = ',ct.ImagePositionPatient[2]/10,'cm'
	
	x0 = np.array([xoff,yoff,zoff]) - hs/2
	nn = np.array((NX,NY,NZ),dtype=np.int32)
	CT = np.zeros(nn,order='F',dtype=np.int16) # pixel data are in row-major

##########################################################################################

def loadCT():
	
	consistencyCheck()

	
	print('\nLoad CT:')
	for iz in range(CT.shape[2]):
		print('.',iz,end='',flush=True)
#                print('slice:',iz,'\n')
		slice=cts[iz].pixel_array.transpose().astype(float)
	# HU = (RawPixelValue * RescaleSlope) + RescaleIntercept		
		slice*=cts[iz].RescaleSlope
		slice+=cts[iz].RescaleIntercept
		CT[:,:,iz]=slice

	print('\nwhole CT shape=',CT.shape)
##########################################################################################

def writePhantom():
	print('\nWrite phantom:')
	print('map dims',Phantom.shape)
	print('map spacing',hs/10.,'cm'	)
	print('map offset',x0/10.,'cm')
	print('saving map to Phantom.mhd' )
	voxels = packVoxels(nn,hs/10.,x0/10.,Phantom)
	mhd_write('Phantom',voxels)

##########################################################################################
##########################################################################################
##########################################################################################

banner()

loadDicoms()
loadCT()

# if '-useExternal' in sys.argv:
# 	buildPhantomFromExternal()
# else:
Phantom=CT  # swallow copy 
writePhantom()
