#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Jan 2021 A.Schiavi
########################################################################
#########################################################################!/usr/bin/env python3
import argparse
parser = argparse.ArgumentParser(description='Convert Fred CT and HU materials into Fluka voxels file')
parser.add_argument("-v", "--verbose", help="increase output verbosity",action="store_true")
parser.add_argument("CT", help="CT map in HU")
parser.add_argument("table", help="dense HU material table")
args = parser.parse_args()

import os,shutil,re
import sys
import numpy as np
import math
import struct

from mhd_io import *

lverbose = args.verbose

def writeFortranRecord(file,record):
    file.write(struct.pack('<1i',len(record)))
    file.write(record)
    file.write(struct.pack('<1i',len(record)))

try:
    if lverbose:
        print('Loading CT map from file ',args.CT)
    [nn,hs,x0,Map] = unpackVoxels(mhd_read(args.CT))
    Map = Map.astype(np.short) # convert to short
    # if Map.dtype!=np.short:
    #   print('CT datatype must be short (integer 16 bit)')
    #   raise 1
    if lverbose:
        print('Loading fred HU table from file ',args.table)
    table = np.loadtxt(args.table,skiprows=1).astype(np.float32)
except:
    print('Fatal error')
    exit(1)

# for i in range(1,table.shape[0]):
#   print(table[i,0],table[i,1]/table[i-1,1])
# exit(0)

# Mapnew=Map.copy()
# HUcut = HUmin+520
# # Mapnew[Mapnew>HUcut]=HUcut
# Mapnew[:] = 1370
# Mapnew[0,0,1] = 1365
# Mapnew[0,0,2] = 1362

# Map=Mapnew.copy()
# HUmin = Map.min()
# HUmax = Map.max()
# print('HU range = ',HUmin,HUmax)
# 

HUmin = Map.min()
HUmax = Map.max()
HUactive = np.unique(Map)
numHUmat = len(HUactive)
print('HU range = ',HUmin,HUmax)
print('HU active = ',HUactive,' = ',len(HUactive),'elements')

HUdense = np.zeros(HUmax-HUmin+1,dtype=np.short)
HUdense[HUactive-HUmin]=1
print('HUdense = ',HUdense,HUdense.sum())
# build HU to VOXEL### link table
KREG=HUdense.copy()

MO = KREG.size-1
NO = KREG.sum()-1
print('MO,NO',MO,NO)
# exit(0)
ireg=-1
for i in range(len(KREG)):
    if KREG[i]!=0:
        ireg+=1
        KREG[i]=ireg

# for i,reg in enumerate(KREG):
#   print(i,reg,HUmin+i)
# exit(0)

if lverbose:
        print('table shape',table.shape)

fin = open(args.table)
line = fin.readline()
cols = line.split()[1:] 
if lverbose:
    print('table header columns',cols)

if len(cols)!= table.shape[1]:
    print('error: columns numbers do not match')
    exit(1)

iHU   =  cols.index('HU')
irho  =  cols.index('rho')
iRSP  =  cols.index('RSP')
iLrad =  cols.index('Lrad')

ielements=list(range(len(cols)))

ielements.remove(iHU)
ielements.remove(irho)
ielements.remove(iRSP)
ielements.remove(iLrad)
if lverbose:
    print(ielements)

fredElements=['H','C','N','O','Na','Mg','P','S','Cl','K','Ca','Fe','I']
flukaElements=['HYDROGEN','CARBON','NITROGEN','OXYGEN','SODIUM','MAGNESIU','PHOSPHO','SULFUR','CHLORINE','POTASSIU','CALCIUM','IRON','IODINE']

iH    =  cols.index('H')
iC    =  cols.index('C')
iN    =  cols.index('N')
iO    =  cols.index('O')
iNa   =  cols.index('Na')
iMg   =  cols.index('Mg')
iP    =  cols.index('P')
iS    =  cols.index('S')
iCl   =  cols.index('Cl')
iK    =  cols.index('K')
iCa   =  cols.index('Ca')
iFe   =  cols.index('Fe')
iI    =  cols.index('I')


table=np.c_[table, np.zeros(table.shape[0])] # add extra column for isActive
iIsActive = table.shape[1]-1

# table[:,iIsActive] = HUactive.T
print(table[:,iIsActive].shape,HUdense.shape)
print(table[:,iIsActive])

for i in range(table.shape[0]):
    hunow = int(table[i,iHU])
    ihu = hunow-HUmin
    isActive = hunow>=HUmin and hunow<=HUmax and HUdense[ihu]>0
    # print(hunow,ihu,isActive)
    table[i,iIsActive] = 1 if isActive else 0
# exit(0)

# build reduced FLUKA materials 
# Dec 2020: fluka non accetta piu' di 510 materiali diversi!!! :-(((
tableHUmin = int(table[ 0][iHU])
tableHUmax = int(table[-1][iHU])

iHUrangemin = table.shape[1]
iHUrangemax = table.shape[1]+1
iCorrFact1 = table.shape[1]+2
iCorrFact2 = table.shape[1]+3

flukaRosettaStone=np.zeros((table.shape[0],4+table.shape[1]),dtype=np.float32)

if numHUmat<200: # build one Fluka material for each HU material
    for i in range(table.shape[0]):
        if(table[i,iIsActive]>0):
            flukaRosettaStone[i,0:iHUrangemin]=table[i,0:iHUrangemin] # copy fred material
            flukaRosettaStone[i,iHUrangemin]=table[i][iHU]
            flukaRosettaStone[i,iHUrangemax]=table[i][iHU]
            flukaRosettaStone[i,iCorrFact1]=table[i,iRSP]/table[i][irho]
            flukaRosettaStone[i,iCorrFact2]=1.0
else: # build HU ranges and use corrfact
    fluHUmat=[0]
    for i in range(1,table.shape[0]-1):
        # if table[i,iHU]<-742 or table[i,iHU]>-740:
        #     continue
        discontinuityFound = False
        for j in ielements:
            # print(i,j,discontinuityFound)
            DerSx = table[i,j]-table[i-1,j]
            DerDx = table[i+1,j]-table[i,j]
            if DerSx==0. and DerDx!=0.:
                discontinuityFound = True
            elif DerSx!=0. and DerDx==0.:
                discontinuityFound = True
            elif DerDx!=0.0 and abs((DerDx-DerSx)/DerDx)>0.20:
                discontinuityFound = True

            # print(i,table[i,iHU],j,DerSx,DerDx,abs((DerSx-DerDx)/DerDx),table[i-1,j],table[i,j],table[i+1,j])
            # print(i,table[i,iHU],DerSx,DerDx,abs((DerSx-DerDx)/DerDx))

            if discontinuityFound:
                # print(i,j,'found')
                break
        if discontinuityFound:
            # print(i,j,table[i,iHU],DerSx,DerDx,abs((DerSx-DerDx)/DerDx))
            fluHUmat.append(i)
    fluHUmat.append(table.shape[0]-1)

    for i in range(table.shape[0]):
        if i in fluHUmat:
            print('fluHUmat',i,table[i,iHU])
    print('not implemented yet!!!! ;)')

    # now build ranged HU materials
    # flukaRosettaStone cols: fluHUrange corrFact1 corrFact2
    ibeg=0
    inext=table.shape[0]
    HUstep=50
    while ibeg<table.shape[0]:
        if ibeg in fluHUmat:
            print('adding ',ibeg)
            if(table[ibeg,iIsActive]>0):
                flukaRosettaStone[ibeg,0:iHUrangemin]=table[ibeg,0:iHUrangemin] # copy fred material
                flukaRosettaStone[ibeg,iHUrangemin]=table[ibeg][iHU]
                flukaRosettaStone[ibeg,iHUrangemax]=table[ibeg][iHU]
                flukaRosettaStone[ibeg,iCorrFact1]=table[ibeg,iRSP]/table[ibeg][irho]
                flukaRosettaStone[ibeg,iCorrFact2]=1.0
            idx = fluHUmat.index(ibeg)+1
            if idx<len(fluHUmat):
                print(ibeg,idx)
                inext = fluHUmat[idx] 
            else:
                inext = table.shape[0]
            ibeg += 1
            continue

        iend = min(ibeg+HUstep,table.shape[0],inext)
        print('try= ',ibeg,iend,'[',inext,']')
        meanHUMat = table[ibeg:iend][:].mean(axis=0)
        while table[ibeg,irho]<2./3.*meanHUMat[irho] or table[iend-1,irho]>3./2.*meanHUMat[irho] or table[ibeg,iRSP]<2./3.*meanHUMat[irho] or table[iend-1,iRSP]>3./2.*meanHUMat[irho]:
          iend-=1
          meanHUMat = table[ibeg:iend][:].mean(axis=0)
        print('got= ',ibeg,iend,'|',table[ibeg,irho]/meanHUMat[irho],table[iend-1,irho]/meanHUMat[irho])
        # print('got= ',ibeg,iend)

        # for v in meanHUMat:
        #   print('%.2e ' % v,end='')
        # name = 'HU<%d' % int(table[iend-1][iHU])
        # print(name)

        # print('%.2e ' % meanHUMat[4:].sum(),end='')
        # print('')
        for i in range(ibeg,iend):
            flukaRosettaStone[i,0:iHUrangemin]=meanHUMat
            flukaRosettaStone[i,iHUrangemin]=table[ibeg][iHU]
            flukaRosettaStone[i,iHUrangemax]=table[iend][iHU]
            flukaRosettaStone[i,iCorrFact1]=table[i,iRSP]/table[i,irho]
            flukaRosettaStone[i,iCorrFact2]=table[i,irho]/meanHUMat[irho]

        ibeg=iend
    # exit(1)

for i in range(table.shape[0]):
    # print(i,end=' ')
    if(table[i,iIsActive]>0):
        print('flukamat = ',i,*flukaRosettaStone[i,:],end='')
        print('')
    # break

# exit(0)

# for i in range(table.shape[0]):
#   print(i,int(table[i,0]),*(flukaRosettaStone[i,iHUrange:]),'|',table[i,irho],
#       flukaRosettaStone[i,iCorrFact2]*flukaRosettaStone[i,irho],'|',
#       table[i,iRSP],
#       flukaRosettaStone[i,iCorrFact1]*flukaRosettaStone[i,iCorrFact2]*flukaRosettaStone[i,irho]
#       )

matLines=[]
matLines.append('MATERIAL         16.    32.066       2.0                              SULFUR    ')
matLines.append('MATERIAL         15. 30.973761       2.2                              PHOSPHO   ')
matLines.append('MATERIAL         17.   35.4527 0.0029947                              CHLORINE  ')
matLines.append('MATERIAL         19.   39.0983     0.862                              POTASSIU  ')
matLines.append('MATERIAL         53. 126.90447      4.93                              IODINE    ')

lFirst=True

for i in range(flukaRosettaStone.shape[0]):
    if flukaRosettaStone[i,iIsActive]==0.0:
        continue # not active => skip it 
    if lFirst:
        lFirst = False
    else:
        if flukaRosettaStone[i,iHUrangemin]==icurr:
            continue # skip it
    icurr = flukaRosettaStone[i,iHUrangemin]
    print(*flukaRosettaStone[i,:])
    
    isHUrange = flukaRosettaStone[i,iHUrangemin]!=flukaRosettaStone[i,iHUrangemax]

    if isHUrange:
        name = 'HU<%d' % int(flukaRosettaStone[i,iHUrangemax])
    else:
        name = 'HU%d' % int(flukaRosettaStone[i,iHUrangemin])


    rho = flukaRosettaStone[i][irho]
    Lrad = flukaRosettaStone[i][iLrad]
    RSP = flukaRosettaStone[i][iRSP]
    composition=[]
    weigths=[]
    for j in ielements:
        if flukaRosettaStone[i][j]>0: 
            composition.append(cols[j])
            weigths.append(flukaRosettaStone[i][j])
            # print(composition[-1],'%e' % weigths[-1])

    if lverbose:
        print('*',i,name,rho,Lrad,composition,weigths)

    # build fluka water for dEdx!!!
    line = 'MATERIAL                      %10f                            0.W%-9s' % (rho,name)
    matLines.append(line)
    line = 'COMPOUND          2.  HYDROGEN        1.    OXYGEN                    W%-9s' % (name)
    matLines.append(line)
    line = 'STERNHEI      3.5017      0.24    2.8004   0.09116    3.4773          W%-9s' % (name)
    matLines.append(line)


    # build fluka material
    line = '%-10s' % ('MATERIAL')
    line += '%10s%10s' % ('','') # what(1) what(2)
    line += '%10f' % (rho) # what(3) density
    line += '%10s' % ('') # what(4)
    line += '%10s' % ('W'+name) # what(5) ref material for dEdx
    line += '%10s' % ('0.') # what(6) 0 = natural isotopic composition
    line+= '%-10s' % (name)
    matLines.append(line)

    line = '%-10s' % ('COMPOUND')
    for iel in range(3):
        if iel < len(composition):
            line += '%10f%10s' % (-weigths[iel]/100.,flukaElements[fredElements.index(composition[iel])])
        else:
            line += '%10s%10s' % ('','')
    line+= '%-10s' % (name)
    matLines.append(line)
    if len(composition)>3:
        line = '%-10s' % ('COMPOUND')
        for iel in range(3,6):
            if iel < len(composition):
                line += '%10f%10s' % (-weigths[iel]/100.,flukaElements[fredElements.index(composition[iel])])
            else:
                line += '%10s%10s' % ('','')
        line+= '%-10s' % (name)
        matLines.append(line)
    if len(composition)>6:
        line = '%-10s' % ('COMPOUND')
        for iel in range(6,9):
            if iel < len(composition):
                line += '%10f%10s' % (-weigths[iel]/100.,flukaElements[fredElements.index(composition[iel])])
            else:
                line += '%10s%10s' % ('','')
        line+= '%-10s' % (name)
        matLines.append(line)
    if len(composition)>9:
        line = '%-10s' % ('COMPOUND')
        for iel in range(9,12):
            if iel < len(composition):
                line += '%10f%10s' % (-weigths[iel]/100.,flukaElements[fredElements.index(composition[iel])])
            else:
                line += '%10s%10s' % ('','')
        line+= '%-10s' % (name)
        matLines.append(line)
    # break

# for line in matLines:
#   print(line)
# exit(0)


assignLines=[]
corrFactLines=[]

assignLines.append('ASSIGNMA      VACUUM     VOXEL                                                  ')
# for i in range(995,1005):
for i in range(table.shape[0]):
    ikreg = int(table[i][iHU])-HUmin
    if ikreg <0 or ikreg >= len(KREG):
        continue
    ireg = KREG[ikreg]
    if ikreg>0 and ireg==0:
        continue # skip this value => not present

    isHUrange = flukaRosettaStone[i,iHUrangemin]!=flukaRosettaStone[i,iHUrangemax]

    if isHUrange:
        name = 'HU<%d' % int(flukaRosettaStone[i,iHUrangemax])
    else:
        name = 'HU%d' % int(flukaRosettaStone[i,iHUrangemin])

# COMPOUND      -0.755  NITROGEN    -0.232    OXYGEN    -0.013     ARGONHU<-1020  
# COMPOUND    0.105909  HYDROGEN  0.339455    CARBON  0.030364  NITROGENHU4       
# COMPOUND    0.519000    OXYGEN  0.000909    SODIUM  0.001000 PHOSPHOHUHU4       
# COMPOUND    0.002000    SULFUR  0.001091  CHLORINE  0.000182  POTASSIUHU4       
# COMPOUND    0.000091      IRON                                        HU4       
# COMPOUND    0.000091      IRON                                        HU4

    line = '%-10s' % ('CORRFACT')
    line += '%10f' % -flukaRosettaStone[i,iCorrFact1] # what(1) relative corr fact dEdx
    line += '%10f' % flukaRosettaStone[i,iCorrFact2] # what(2) corr fact other
    line += '%10s' % ('') # what(3)
    voxelno = ('VOXEL%03d' % (ireg+1)) if (ireg+1) < 1000 else ('VOXE%04d' % (ireg+1))
    line += '%10s' % voxelno
    line += ' '*30
    corrFactLines.append(line)

    line = '%-10s' % ('ASSIGNMA')
    line += '%10s' % name
    line += '%10s' % voxelno
    line += ' '*50
    assignLines.append(line)

# for i,line in enumerate(assignLines):
#   print(line)
# exit(0)

# for line in corrFactLines:
#   print(line)
# exit(0)


# MATERIAL                        1.020910               WATER        0.HU4       
# CORRFACT    1.030090  1.000000            VOXE1004                              
# ASSIGNMA         HU4  VOXE1004                                                  
# COMPOUND    0.105909  HYDROGEN  0.339455    CARBON  0.030364  NITROGENHU4       
# COMPOUND    0.519000    OXYGEN  0.000909    SODIUM  0.001000 PHOSPHOHUHU4       
# COMPOUND    0.002000    SULFUR  0.001091  CHLORINE  0.000182  POTASSIUHU4       
# COMPOUND    0.000091      IRON                                        HU4       


# write voxels file
fout = open('CTHU.vxl','wb')
title = '%-80s' % 'CT from Fred'
record = struct.pack('<80s',str.encode(title))
writeFortranRecord(fout,record)
# print(len(record))
record = struct.pack('<5i',*nn,NO,MO)
writeFortranRecord(fout,record)
record = struct.pack('<3d',*hs)
writeFortranRecord(fout,record)
Map=Map-HUmin
writeFortranRecord(fout,Map.tobytes(order='F'))
writeFortranRecord(fout,KREG[1:].tobytes(order='F'))

for line in matLines:
    writeFortranRecord(fout,str.encode(line))
    print(line)
for i,line in enumerate(assignLines):
    writeFortranRecord(fout,str.encode(line))
    print(line)
for line in corrFactLines:
    writeFortranRecord(fout,str.encode(line))
    print(line)

