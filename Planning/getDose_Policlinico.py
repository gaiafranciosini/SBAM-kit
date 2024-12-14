#!/usr/bin/env python3
######################################################################
# Author:		A. Schiavi 
# Date:			2020
# Purpose:		routines for importing DICOM files in Fred
# Source:		python
#####################################################################

import os, sys, shutil, re, time
from math import *
import numpy as np
import pydicom as dicom
import struct 

from mhd_io import *

# global vars
dicoms=[]
cts=[]
rtstruct=[]
rtplan=[]
rtdose=[]
isocoordinates = np.zeros(3)
#==================

# VERSION=105

# if '-v' in sys.argv:
# 	print('Fred Dicom importer version %d.%d' % (VERSION/100,VERSION%100))
# 	exit(0)

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
######################### RTDOSE ###################################

def buildDOSE():

	if len(rtdose)==0:
		return 

	rtd = rtdose[0]
		
	print("\n\nLoading RTDOSE ...\n")
	print('Number of frames:',rtd.NumberOfFrames)
	print('Rows:',rtd.Rows)
	print('Columns:',rtd.Columns)
	if rtd.BitsAllocated==16:
		doseArr = np.array(rtd.pixel_array, dtype=np.uint16)
	if rtd.BitsAllocated==32:
		doseArr = np.array(rtd.pixel_array, dtype=np.uint32)
	
	print('shape=',doseArr.shape)
	if (doseArr.shape[0] != rtd.NumberOfFrames):
			print('Error,dose array mismatch!')
			exit(-1)
	if (doseArr.shape[1] != rtd.Rows):
			print('Error,dose array mismatch!')
			exit(-1)
	if (doseArr.shape[2] != rtd.Columns):
			print('Error,dose array mismatch!')
			exit(-1)
	scalingFactor = rtd.DoseGridScaling

	doseArr = doseArr.astype('float32',order='F')
	doseArr *= scalingFactor;
	print('scaling factor=',scalingFactor)
	flat = doseArr.flatten()
	print('min,max=',min(flat),max(flat))
	nvxl = np.array([doseArr.shape[2],doseArr.shape[1],doseArr.shape[0]],dtype=np.int32)
	doseArr = np.reshape(flat,nvxl,order='F')
	
# 	print rtd.GridFrameOffsetVector
	thickness = rtd.GridFrameOffsetVector[1]-rtd.GridFrameOffsetVector[0]
	spacing =  np.array([rtd.PixelSpacing[0],rtd.PixelSpacing[1],thickness])/10.0
	x0 = np.array(rtd.ImagePositionPatient)/10.0
	x0 -= spacing/2
	x0 -= isocoordinates/10
	print('map offset',x0)
	print('saving map to RTDose.mhd' )
	print('nvxl = ',nvxl)
	voxels = packVoxels(nvxl,spacing,x0,doseArr)
	mhd_write('RTDose',voxels)
		

##########################################################################################
##########################################################################################
##########################################################################################

banner()
loadDicoms()
buildDOSE()
