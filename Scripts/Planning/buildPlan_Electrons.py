#!/usr/bin/env python3
######################################################################
# Author:       A. Schiavi
# Date:         2020
# Purpose:      build an RTPLAN for Flash electrons project
# Source:       python
######################################################################
import argparse
parser = argparse.ArgumentParser(description='build an RTPLAN using PTV and beams')
parser.add_argument("-ptv",help="PTV mhd map",nargs=None,metavar='path',required=True)
parser.add_argument("-beams",help="file containing beam positions and directions",nargs=None,metavar='path',required=True)
parser.add_argument("-margin",help="radius of internal margin(s) to include around the ptv shadow",nargs='+',metavar='#',type=float,default=[0])
parser.add_argument("-margin_ext",help="radius of external margin(s) to include around the ptv shadow",nargs=None,metavar='#',type=float,default=0)
parser.add_argument("-spotSpacing",help="inter-spot distance in the central core",nargs=None,required=True,metavar='#',type=float)
parser.add_argument("-spotSpacing_ext",help="inter-spot distance in the external core",nargs=None,metavar='#',type=float)

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-paraxial", help="spots are directed parallel to beam axis",action="store_true")
group.add_argument("-pointSource", help="spots start all from beam origin",action="store_true")

outputfiles=parser.add_mutually_exclusive_group(required=True)
outputfiles.add_argument("-fluka_output", help="set fluka output files",action="store_true")
outputfiles.add_argument("-fred_output", help="set fluka output files",action="store_true")
outputfiles.add_argument("-Optimization_output", help="set Optimization output files",action="store_true")

parser.add_argument("-energy",help="beam nominal energy for each field",nargs='+',metavar='#',type=float,default=[0])

parser.add_argument("-FWHM",help="FWHM of a gaussian spot in the central core",nargs=None,metavar='#',default=0.5,type=float)
parser.add_argument("-FWHM_ext",help="FWHM of a gaussian spot in the external core",nargs=None,metavar='#',default=0.5,type=float)

parser.add_argument("-particle",help="beam primary particle",nargs=None,metavar='name',default='electron')
parser.add_argument("-Nprim",help="Number of primaries",nargs=None,required=True,metavar='#',type=int)

parser.add_argument("-Variable_SS",help="RTPLAN with varaible spots spacing",action="store_true")
parser.add_argument("-Variable_Energy",help="RTPLAN with varaible energy",action="store_true")
parser.add_argument("-PreTPS",help="Choice of Fields",action="store_true")

args = parser.parse_args()

# print(args) ; exit()
######################################################################

import os, sys
from math import *
import numpy as np
from numpy.linalg import norm
from beamProjectROI import *

########################################################################
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
########################################################################


class Beam:
    def __init__(self, line):
        tok = line.split()
        self.ID = int(tok[0])
        idx=1
        self.O = np.array(list(map(float,tok[idx:idx+3]))) ; idx +=3
        self.l = np.array(list(map(float,tok[idx:idx+3]))) ; idx +=3
        self.u = np.array(list(map(float,tok[idx:idx+3]))) ; idx +=3
        self.f = np.array(list(map(float,tok[idx:idx+3]))) ; idx +=3
       
        if args.PreTPS :
            self.fromISO = float(tok[idx]); idx +=1
            self.theta=float(tok[idx]); idx += 1
            self.phi=float(tok[idx])
        else:
            self.fromISO = float(tok[idx])   
    def setControlPoints(self,pts):
        self.controlPoints = pts.copy()
    def setFWHM(self,pts):
        self.FWHM = pts.copy()
######################################################################
def pprint(arr):
    return '['+','.join(list(map(str,arr)))+']'

# global vars
beams = []
ptvPoints = []
PTV = None
maskResolution = 0.3 #era 0.7, 0.3 viene molto dettagliata
######################################################################
def loadBeams():
    fin = open(args.beams,'r') or die
    for l in fin.readlines():
        beams.append(Beam(l))
##########################################################################################
def writePlan():
    fout=open('rtplan.txt','w') or die
    fmargins=open('planInfo.txt','w') or die
    print('PLAN:',file=fmargins)
    print('##########INSIDE##########',file=fmargins)
    print('margin_in: {} cm'.format(args.margin),file=fmargins)
    print('spotSpacing_in: {} cm'.format(args.spotSpacing),file=fmargins)
    print('FWHM_in: {} cm'.format(args.FWHM),file=fmargins)
    if args.Variable_SS:
        print('#########OUTSIDE########',file=fmargins)
        print('margin_ext: {} cm'.format(args.margin_ext),file=fmargins)
        print('spotSpacing_ext: {} cm'.format(args.spotSpacing_ext),file=fmargins)
        print('FWHM_ext: {} cm'.format(args.FWHM_ext),file=fmargins)
    
    fbeam=[]
    for b in beams:
        fname ='beam_{}.txt'.format(b.ID) 
        f=open(fname,'w') or die
        fbeam.append(f)
        print('include: {}'.format(fname),file=fout)

    beamE = []
   
    for b in args.energy:
        beamE.append(b)
    Ene = beamE[0]
    for i,b in enumerate(beams):
        if args.paraxial:
            if  '-fred_output' in sys.argv:
                print('RTPLAN:  pb field pb_type      E pb_nprim  FWHM  x0 y0 z0 f_x f_y f_z u_x u_y u_z l_x l_y l_z',file=fbeam[i])
                print('#     pb    field    pb_type      E         pb_nprim     FWHM  x0      y0    z0   f_x  f_y  f_z  u_x  u_y  u_z  l_x  l_y l_z',file=fbeam[i])
              
                
        if args.pointSource:
          if  '-fred_output' in sys.argv:
              print('RTPLAN:  pb field pb_type      E pb_nprim  FWHM  x0 y0 z0 f_x f_y f_z u_x u_y u_z l_x l_y l_z',file=fbeam[i])
              print('#     pb    field    pb_type      E         pb_nprim     FWHM  x0      y0    z0   f_x  f_y  f_z  u_x  u_y  u_z  l_x  l_y l_z',file=fbeam[i])
                

        
        ispot=0
        imargin=0
        uu=[1,0,0]
        Orig = [0,0,0]
        xcord = [0,0,0]
        NEWcord = [0,0,0]

        for spots,fwhm in zip(b.controlPoints,b.FWHM):
            for j in range(spots.shape[0]):
                spot = spots[j]
                ispot+=1
                if args.paraxial:
                    if  '-fred_output' in sys.argv:
    
                        xcord=spot[0]*b.l+spot[1]*b.u+0*b.f
#                        print('xcord',xcord)
#                        print('b.O', b.O)
                        Orig = b.O+xcord
                        
                        
                        print('pb: {:g} {:g} electron {:g} {:d} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}'.format(ispot,b.ID,Ene,args.Nprim,fwhm,Orig[0],Orig[1],Orig[2],b.f[0],b.f[1],b.f[2],b.u[0],b.u[1],b.u[2],b.l[0],b.l[1],b.l[2]),file=fbeam[i]) #fred

                    if  '-fluka_output' in sys.argv:
                        print('{:.6f} {:.6f} {:.6f}'.format(b.f[0],b.f[1],b.f[2]),file=fbeam[i]) #fluka da controllare
                        
                if args.pointSource:
                    module=norm(spot[0]*b.l+spot[1]*b.u+b.fromISO*b.f)
                    uvett=spot[0]*b.l+spot[1]*b.u+b.fromISO*b.f
                    uver=uvett/module

                    
                    #sistema di riferimento Field
                    Ofld=b.O
                    
                    ffld=b.f
                    ufld=b.u
                    lfld=b.l
                    #front vector dello spot
                    fspt=uver
                    
                    #up vector dello spot
                    uspt=ufld-np.dot(ufld,fspt)*fspt
                    uspt/=norm(uspt)
                    
                    #left vector dello spot
                    lspt=np.cross(uspt,fspt)
                    lspt/=norm(lspt)
                    

                    cosBxx=np.dot(lspt,np.array([1,0,0]))
                    cosBxy=np.dot(lspt,np.array([0,1,0]))
                    cosBxz=np.dot(lspt,np.array([0,0,1]))

                    cosBzx=np.dot(fspt,np.array([1,0,0]))
                    cosBzy=np.dot(fspt,np.array([0,1,0]))
                    cosBzz=np.dot(fspt,np.array([0,0,1]))

                   
            
                    if  '-fred_output' in sys.argv:
                        if args.Variable_Energy:
                             print('pb: {:g} {:g} electron {:g} {:d} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}'.format(ispot,b.ID,beamE[i],args.Nprim,fwhm,Ofld[0],Ofld[1],Ofld[2],fspt[0],fspt[1],fspt[2],uspt[0],uspt[1],uspt[2],lspt[0],lspt[1],lspt[2]),file=fbeam[i]) #fred
                        else:
                             print('pb: {:g} {:g} electron {:g} {:d} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}'.format(ispot,b.ID,Ene,args.Nprim,fwhm,Ofld[0],Ofld[1],Ofld[2],fspt[0],fspt[1],fspt[2],uspt[0],uspt[1],uspt[2],lspt[0],lspt[1],lspt[2]),file=fbeam[i]) #fred
                            
                        #if args.PreTPS:
                          #  print('pb: {:g} {:g} electron {:g} {:d} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g} {:g}'.format(ispot,b.ID,beamE[i],args.Nprim,fwhm,Ofld[0],Ofld[1],Ofld[2],fspt[0],fspt[1],fspt[2],uspt[0],uspt[1],uspt[2],lspt[0],lspt[1],lspt[2],b.theta,b.phi),file=fbeam[i]) #fred
                       # else:
                           
                            
                    if  '-Optimization_output' in sys.argv:
                        print('{:g} {:g} {:g}'.format(spot[0],spot[1],b.fromISO),file=fbeam[i]) #Optimization

                    if  '-fluka_output' in sys.argv:
                        print('{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}'.format(fwhm,cosBxx,cosBxy,cosBxz,cosBzx,cosBzy,cosBzz),file=fbeam[i]) #fluka            
        # exit(0)
    
    print('plan written to file: rtplan.txt and beam_#.txt')
##########################################################################################
def removeDuplicates(spotSets,spotSpacing,spotSpacing_ext):
    setPrev = spotSets[0]
    newSets = []
    for i,set in enumerate(spotSets):
        mySet = set.copy()
        if i==0: 
            newSets.append(mySet)
            continue
        setPrev = spotSets[i-1]
        icurr=0
        for j in range(set.shape[0]):
            isDuplicate = False
            for k in range(setPrev.shape[0]):
                if (norm(set[j,:]-setPrev[k,:])<(1e-2*spotSpacing or 1e-2*spotSpacing_ext)):
                    isDuplicate=True
                    break
            if not isDuplicate:
                mySet[icurr,:] = set[j,:]
                icurr+=1
        newSets.append(mySet[0:icurr])

    for i,set in enumerate(spotSets):
        print('{}: {} => {}'.format(i,set.shape[0],newSets[i].shape[0]))

    return newSets



##########################################################################################

loadBeams()

args.margin = sorted(args.margin)
args.energy=sorted(args.energy)

np.set_string_function(pprint,repr=False)

for b in beams:
    print(b.ID,b.O,b.l,b.u,b.f,b.fromISO)

for b in beams:
    print(b.O+b.f*b.fromISO)

PTV = mhd_read(args.ptv) or die
ptvPoints = getROIPoints(PTV)



ncols = int(len(beams)/2) +len(beams)%2 +1

nraws = 2
if len(beams)>7:
    nraws=3
    ncols= int(len(beams)/3) +len(beams)%3

#fig,axs = plt.subplots(3, 2,figsize=(8.3*0.6,11.7*0.6),constrained_layout=True)

fig,axs = plt.subplots(nraws,ncols,figsize=(12,8),constrained_layout=True)

fig.suptitle('Spot spacing = {}  ; margin = {}'.format(args.spotSpacing,args.margin))

axs = axs.flatten()

# plt.ion()
# drawDirections(axs[-1],PTV,beams,args.paraxial)
# plt.show()
# input('Hit return')
# exit(0)

for i,b in enumerate(beams):
    roiPts = beamProjectROI(b,ptvPoints)
    #print(ptvPoints) 
    shadow,extent = buildROImask(roiPts,maskResolution)
    print(extent)
    spotSets=[]
    spotFWHM=[]
    for m in args.margin:
        print('margin=',m)
        if m>0:
            mask = dilateMask(shadow,maskResolution,m)
        elif m<0:
            mask= erodeMask(shadow,maskResolution,-m)
        else:
            mask= shadow
        spotSets.append(layoutSpotGrid(mask,extent,args.spotSpacing))
        spotFWHM.append(args.FWHM)
        if args.Variable_SS :
            mask= donutMask(shadow,maskResolution,-m,args.margin_ext)
            spotSets.append(layoutSpotGrid(mask,extent,args.spotSpacing_ext))
            spotFWHM.append(args.FWHM_ext)
            
    spotSets = removeDuplicates(spotSets,args.spotSpacing,args.spotSpacing_ext)
    # spots=spotSets[0]
    b.setControlPoints(spotSets)
    b.setFWHM(spotFWHM)
    #b.setControlPoints(spotFWHM)
    nspot=0
    for set in spotSets:
        nspot += set.shape[0]
    # print(spots[:][0])
    axs[i].set_title('Beam %s (nspot=%d)'%(b.ID,nspot))
    drawMasks(axs[i],shadow,mask,extent,spotSets)
    # exit(0)
    # break
# exit(0)
#plt.show()



plt.ion()
 
drawDirections(axs[-1],PTV,beams,args.paraxial)
#axs[-1].remove()
# fig.tight_layout()
# axs[-1].axes.get_xaxis().set_visible(False)
# axs[-1].axes.get_yaxis().set_visible(False)
plt.show()
input('Hit return')

plt.savefig('spots.pdf')  
plt.ioff()

writePlan()
