#!/usr/bin/env python
########################################################################
########################################################################
# Fred project
# asch Oct 2020
########################################################################
########################################################################
import os, sys
from math import *
import numpy as np
from numpy import dot,cross
from scipy import ndimage

from mhd_io import *

import matplotlib.path as mpltPath

import matplotlib
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.cm as cm

matplotlib.use('TkAgg')
import pylab as plt

import itertools

########################################################################
########################################################################
########################################################################
def getROIPoints(ROI):
    # carica la ROI
    [nn,hs,x0,Map] = unpackVoxels(ROI)
    # vettori coordinate dei centri dei voxel della ROI
    xc = (np.arange(nn[0])+0.5)*hs[0]+x0[0]
    yc = (np.arange(nn[1])+0.5)*hs[1]+x0[1]
    zc = (np.arange(nn[2])+0.5)*hs[2]+x0[2]

    # matrici "dense" che conterranno la coord di ogni singolo voxel
    XC = np.zeros_like(Map)+xc.reshape((nn[0],1,1))
    YC = np.zeros_like(Map)+yc.reshape((1,nn[1],1))
    ZC = np.zeros_like(Map)+zc.reshape((1,1,nn[2]))

    # identifico gli indici dei voxel che appartengono alla ROI
    I = np.where(Map>0) # active voxels

    # ottengo le coord dei voxel accesi della ROI
    xcroi = XC[I]
    ycroi = YC[I]
    zcroi = ZC[I]
    
    return xcroi,ycroi,zcroi


################################################################################
def beamProjectROI(beam,roiPts):


    # posizione dell'ISOCENTRO
    ISO=beam.O+beam.f*beam.fromISO
    print('ISO=',ISO)


    # distanza del piano ortogonale al beam passante per ISO
    zISO = np.dot(ISO-beam.O,beam.f)

    print('zISO',zISO)

    # # scrivo i punti della ROI (SDR globale)
    # froi = open('roi','w')
    # for i in range(len(xcroi)):
    #   print(xcroi[i],ycroi[i],zcroi[i],file=froi)

    xcroi,ycroi,zcroi = roiPts

    # passo al SDR del beam
    # traslo origine nell'origine del beam 
    xcroi = xcroi - beam.O[0]
    ycroi = ycroi - beam.O[1]
    zcroi = zcroi - beam.O[2]


    # # calcolo le coord dei punti rispetto al beam
    xp = xcroi*beam.l[0]+ycroi*beam.l[1]+zcroi*beam.l[2]
    yp = xcroi*beam.u[0]+ycroi*beam.u[1]+zcroi*beam.u[2]
    zp = xcroi*beam.f[0]+ycroi*beam.f[1]+zcroi*beam.f[2]

    points = np.vstack([xp,yp]).T

    # scrivo i punti proiettati nel SDR del beam
    # (SDR globale)
    # fpoints=open('projected_%s' % beam.ID,'w')
    # for P in points:
    #     print(*P,file=fpoints)
    return points

################################################################################
def buildROImask(pts,res):
    xp,yp = pts[:,0],pts[:,1]
    # calcolo e adatto le dimensioni della griglia
    xmin = np.min(xp)
    xmax = np.max(xp)
    ymin = np.min(yp)
    ymax = np.max(yp)
    print('xmin,xmax,ymin,ymax',xmin,xmax,ymin,ymax)
    # print(xp)
    # print(yp)

    L = min(max(np.abs(np.array([xmin,xmax,ymin,ymax])))+5,10)
    print('L',L)

    nx = int(L/0.2)*2

    mask = np.zeros((nx,nx),dtype=float)

    x0 = -nx/2*res
    y0 = x0

    ix = ((xp-x0)/res+0.5).astype(int)
    iy = ((yp-y0)/res+0.5).astype(int)

    print(xmin,xmax,np.min(ix),np.max(ix))

    print(mask.shape)

    mask[ix,iy] = 1
    extent = (x0,x0+nx*res,y0,y0+nx*res)
    return (mask,extent)
##########################################################################################
def roundKernel(radius,spacing):
    nr=int(radius/spacing)
    if nr<1:
        return np.ones((1,1),dtype=float)
    nk = 2*nr+1
    kernel = np.ones((nk,nk),dtype=float)
    for i in range(nk):
        for j in range(nk):
            kernel[i,j] = 1 if (i-nr)*(i-nr)+(j-nr)*(j-nr)<=nr*nr else 0
    return kernel
##########################################################################################
def dilateMask(mask,spacing,radius):
    # amplia la maschera applicando un kernel di raggio dato
    kernel = roundKernel(radius,spacing)
    grow = ndimage.filters.maximum_filter(mask,footprint=kernel,mode='constant')
    return grow
##########################################################################################
def erodeMask(mask,spacing,radius):
    # amplia la maschera applicando un kernel di raggio dato
    kernel = roundKernel(radius,spacing)
    shrink = ndimage.filters.minimum_filter(mask,footprint=kernel,mode='constant')
    return shrink
##########################################################################################

def donutMask(mask,spacing,rmin,rmax):
    # amplia la maschera applicando un kernel di raggio dato
    kernel_min = roundKernel(rmin,spacing)
    kernel_max = roundKernel(rmax,spacing)
    donut = ndimage.filters.maximum_filter(mask,footprint=kernel_max,mode='constant') - ndimage.filters.minimum_filter(mask,footprint=kernel_min,mode='constant')
    return donut

##########################################################################################

def layoutSpotGrid(mask,extent,spotSpacing):
    # crea la griglia degli spot usando la maschera
    hx = (extent[1]-extent[0])/mask.shape[0]
    hy = (extent[3]-extent[2])/mask.shape[1]

    print('mask.shape', mask.shape[0],mask.shape[1])

    print('hx, hy ',hx,hy)

    I = np.where(mask>0)

    print('Mask range',np.min(mask),np.max(mask))

    xmin = extent[0]+np.min(I[0])*hx
    xmax = extent[0]+np.max(I[0])*hx
    ymin = extent[2]+np.min(I[1])*hy
    ymax = extent[2]+np.max(I[1])*hy    

    print('Grid size',xmin,xmax,ymin,ymax)

    # Lx=xmax-xmin
    # Ly=ymax-ymin
    Lx=2*max(abs(xmax),abs(xmin))
    Ly=2*max(abs(ymax),abs(ymin))
    # print('Lx,Ly = ',Lx,Ly)

    #print('Lx and Ly:',Lx,Ly)

    # default spacing: central spot
    nx = int(Lx/2/spotSpacing)
    xgrid = np.arange(-nx,nx+0.5,1)*spotSpacing
    print('nx,xgrid: ',nx,xgrid)
    
    nx = len(xgrid)

    print('finale nx ',nx)

    

    ny = int(Ly/2/spotSpacing)
    ygrid = np.arange(-ny,ny+0.5,1)*spotSpacing
    ny = len(ygrid)
    
    # print('xgrid',xgrid)
    # print('ygrid',ygrid)

    # # best fit spacing
    # nx = int(Lx/spotSpacing+0.5)+1
    # xgrid = xmin+np.arange(nx)*spotSpacing
    # ny = int(Ly/spotSpacing+0.5)+1
    # ygrid = ymin+np.arange(ny)*spotSpacing

    xp=[]
    yp=[]

    #print('extend[0] and [2]: ', extent[0],extent[2])
    
    for x in xgrid:
        ix = int((x-extent[0])/hx)
        for y in ygrid:
            iy = int((y-extent[2])/hy)
            if mask[ix,iy]>0:
                xp.append(x)
                yp.append(y)
                #print('x and y:',ix,iy)

    spots = np.column_stack([xp,yp])
    return spots.astype(np.float32)
##########################################################################################
def drawMasks(ax,shadow,mask,_extent,spotsSets):

    # plt.ion()
    # plt.clf()



    # Extent = (0,mask.shape[0],0,mask.shape[1])
    Extent = _extent

    im = ax.imshow(mask.T, interpolation='none', origin='lower',
                cmap=cm.coolwarm,aspect=1,extent=Extent,alpha=1.0)

    # for i in range(mask.shape[0]):
    #   mask[i,i%2::2]=2
    shadow = np.ma.masked_where(shadow <=0 , shadow)
    im2 = ax.imshow(shadow.T, interpolation='nearest', origin='lower',
                cmap=cm.gray,aspect=1,alpha=0.5,extent=Extent)      

    # for i in range(mask.shape[0]):
    #   mask[i,i%2::2]=2
    # ax.imshow(mask.T, interpolation='nearest', origin='lower',
    #             cmap=cm.gray,aspect=1,alpha=0.5,extent=Extent)      

    marker = itertools.cycle(('.','+' , '*' , 'x', 'o')) 
    for spots in spotsSets:
        # print('nspots',spots.shape[0])
        ax.scatter(spots[:,0],spots[:,1],marker=next(marker),s=10,color='w')
    # plt.scatter([-8.90],[-9.3],marker='x',s=1,color='g')      

    # print(spots[:])
    # for Xcont,Ycont in Contours:
    #   plt.plot(Xcont,Ycont,'k-',label='contour',lw=0.5)
    # plt.plot(Xcont[0],Ycont[0],'r-o',label='contour',lw=1)
    # print('start',Xcont[0],Ycont[0])
#   print Xcont
#   print Ycont
    plt.legend()
    ax.set_xlabel('cm')
    ax.set_ylabel('cm')

    # exit(0)

########################################################################
def getROICentralSliceZ(ROI):
    [nn,hs,x0,Map] = unpackVoxels(ROI)
    zc = np.arange(nn[2])*hs[2]+x0[2]+hs[2]/2
    I = np.where(Map>0)
    izmin = np.min(I[2])
    izmax = np.max(I[2])
    izmid = int((izmin+izmax)/2)
    zmid = zc[izmid]
    return Map[:,:,izmid],zmid



def drawDirections(ax,ROI,beams,isParaxial):
    ax.set_facecolor('black')

    [nn,hs,x0,Map] = unpackVoxels(ROI)
    
    # slice,z = getROICentralSliceZ(ROI)
    projectedROI = np.sum(Map,axis=2).transpose().squeeze()
    Extent = [x0[0],x0[0]+hs[0]*nn[0],x0[1]+hs[1]*nn[1],x0[1]]
    print("Extent:")
    print(x0[0])
    print(x0[0]+hs[0]*nn[0])
    print(x0[1]+hs[1]*nn[1])
    print(x0[1])

    im = ax.imshow(projectedROI, interpolation='none', origin='upper',
                cmap=cm.seismic,aspect=1,extent=Extent,alpha=1.0)

    uz = np.array([0,0,1]) # z-versor 
    for b in beams:
        fpara = b.f-dot(b.f,uz)*uz
        uperp = cross(uz,fpara)
        print(b.f,fpara,uperp)
        spots = b.controlPoints[-1]
        npts = spots.shape[0]
        # ptsISO = np.zeros((npts,3))
        xperp  = np.zeros(npts)
        for i in range(npts):
            ptISO = b.O + spots[i,0]*b.l + spots[i,1]*b.u + b.fromISO*b.f
            xperp[i] = dot(ptISO,uperp)

        imin = np.argmin(xperp)
        imax = np.argmax(xperp)
        if isParaxial:
            P0 = b.O + spots[imin,0]*b.l + spots[imin,1]*b.u 
            P1 = b.O + spots[imax,0]*b.l + spots[imax,1]*b.u 
        else:
            P0 = b.O
            P1 = b.O

        P2 = b.O + spots[imax,0]*b.l + spots[imax,1]*b.u + b.fromISO*b.f
        P3 = b.O + spots[imin,0]*b.l + spots[imin,1]*b.u + b.fromISO*b.f

        xloop = [P0[0],P1[0],P2[0],P3[0],P0[0]]
        yloop = [P0[1],P1[1],P2[1],P3[1],P0[1]]

        #ax.scatter(xloop,yloop,marker='x',s=1,color='k')
        ax.fill(xloop,yloop,label='Beam_{}'.format(b.ID),alpha=0.5)

    tit = 'Paraxial beam' if isParaxial else 'Point source beam'
    ax.set_title(tit)
    ax.legend(loc='lower left', bbox_to_anchor=(0.8, 0.8),fontsize='small',framealpha=1.0)
    ax.set_xlabel('cm')
    ax.set_ylabel('cm')

    # exit(0)

