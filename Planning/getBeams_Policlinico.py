#!/usr/bin/env python3
######################################################################
# Author:		A. Schiavi
# Date:			2020
# Purpose:		routines for importing DICOM files in Fred
# Source:		python
######################################################################
import argparse
parser = argparse.ArgumentParser(description='Export field info from dicom RTPLAN')
parser.add_argument("file",help="dicom files",nargs='+',metavar='path')
args = parser.parse_args()
# print(args) ; exit()
######################################################################

import os, sys
from math import *
import numpy as np
import pydicom as dicom
# import matplotlib
# import matplotlib.colors as colors
# import matplotlib.ticker as ticker
# import matplotlib.cm as cm

# matplotlib.use('TkAgg')
# import pylab as plt

# from mhd_io import *

######################################################################
# global vars
dicoms=[]
cts=[]
cts=[]
rtstruct=[]
rtplan=[]
rtdose=[]
ISO = np.zeros(3)
beams = []
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
		# print(ds.SOPClassUID)
		if ([0x0008,0x0016] not in ds):	 # check if SOPClassUID exists
			continue
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.2':  # CT Image Storage
			cts.append(ds)
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.3':	# Radiation Therapy Structure Set Storage
			rtstruct.append(ds)
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.5': # Radiation Therapy Plan Storage
			rtplan.append(ds)
		if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.2': # Radiation Therapy Dose Storage
			rtdose.append(ds)

	# report
	print('\tnum of CT images     :',len(cts))
	print('\tnum of RTDose files  :',len(rtdose))
	print('\tnum of RTStruct files:',len(rtstruct))
	print('\tnum of RTPlan files  :',len(rtplan))
##########################################################################################
def getBeamFrame(psideg,phideg,far):
	global Gf,Gu,Go,Gl

	# far = source distance

	#initial triads for Gantry at zero position
	G0f = np.array([0,1,0])
	G0o = np.array([0,-far,0])
	G0u = np.array([0,0,1])
	
	Rot = np.zeros((3,3))
	Rot2 = np.zeros((3,3))
	
	psi = -pi/180*psideg
	phi = pi/180*phideg

	# CW rotation on Z
	Rot[0,0] = round( cos(psi),3) 
	Rot[1,0] = round(-sin(psi),3)
	Rot[0,1] = round( sin(psi),3)
	Rot[1,1] = round( cos(psi),3)
	Rot[2,2] = 1
	# CW rotation on -Y
	Rot2[0,0] = round( cos(phi),3) 
	Rot2[0,2] = round( sin(phi),3) # CW rotation
	Rot2[2,0] = round(-sin(phi),3)
	Rot2[2,2] = round( cos(phi),3)
	Rot2[1,1] = 1

	# first rotation
	Go = Rot.dot(G0o)
	Gf = Rot.dot(G0f)
	Gu = Rot.dot(G0u)
	# second rotation
	Go = Rot2.dot(Go)
	Gf = Rot2.dot(Gf)
	Gu = Rot2.dot(Gu)
	
	Gl = np.cross(Gu,Gf)

	print('\tBeam Origin Vec = ',Go)
	print('\tBeam Front Vec  = ',Gf)
	print('\tBeam Up Vec     = ',Gu)
	print('\tBeam Left Vec   = ',Gl)

##########################################################################################
def loadPLAN():
	global Gf,Gu,Gl,ISO,Go

	print('\nLoading plan:')

	rtPlan = rtplan[0]

	vt = rtPlan.ReferencedStructureSetSequence[0].ReferencedSOPInstanceUID

	#MU2primary = 2.2e8
	#print ''
	#print 'Using MU to primary conversion: %.3e' % MU2primary
	

	print('\nLoading RTPLAN ...\n')

	nbeam = rtPlan.FractionGroupSequence[0].NumberOfBeams
	print('Number of beams: ',nbeam)

	nfractions = rtPlan.FractionGroupSequence[0].NumberOfFractionsPlanned
	print('Number of fractions: ',nfractions)

	for ibeam in range(nbeam):	
		print('Beam',ibeam+1,'/',nbeam,':')
		beam = rtPlan.BeamSequence[ibeam]
		print('\tBeam Type:',beam.RadiationType)
		masterControlPoint = beam.ControlPointSequence[0]
		GantryAngle = masterControlPoint.GantryAngle
		CouchAngle = masterControlPoint.PatientSupportAngle
		SourceToSurfaceDistance = masterControlPoint.SourceToSurfaceDistance/10.
		print('\tGantryAngle = ',GantryAngle)
		print('\tCouchAngle = ',CouchAngle)
		print('\tSourceToSurfaceDistance = ',SourceToSurfaceDistance,'cm')
		
		getBeamFrame(GantryAngle,CouchAngle,SourceToSurfaceDistance)

		ISO = np.array(masterControlPoint.IsocenterPosition)/10.
		ISO[np.where(abs(ISO)<1e-6)]=0.
		print('\tISO position:',ISO,'cm')
		
		Go += ISO
		# Gf = -Gf 
		Go[np.where(abs(Go)<1e-6)]=+0.0
		Gl[np.where(abs(Gl)<1e-6)]=+0.0
		Gu[np.where(abs(Gu)<1e-6)]=+0.0
		Gf[np.where(abs(Gf)<1e-6)]=+0.0

		mybeam = {'ID':beam.BeamNumber,'O':Go,'f':Gf,'u':Gu,'l':Gl,'distance':SourceToSurfaceDistance}

		beams.append(mybeam)		

##########################################################################################
##########################################################################################

banner()
loadDicoms()
loadPLAN()

fout=open('beams.txt','w')
print('-'*80)
for b in beams:
	line=''
	line += '%-4s' % b['ID']
	frame = (*b['O'],*b['l'],*b['u'],*b['f'],b['distance'])
	# print(frame)
	for x in frame:
		line+='%-8g ' % x
	print(line)
	print(line,file=fout)
print('-'*80)
print('output written to file: beams.txt')
