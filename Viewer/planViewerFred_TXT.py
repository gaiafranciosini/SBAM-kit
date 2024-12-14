#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# 11 Oct 2020 A.Schiavi 
########################################################################
########################################################################
import argparse
parser = argparse.ArgumentParser(description='Presents CT + Dose + ROIs and DVH for a complete RTplan')
parser.add_argument("-CT",help="CT map",nargs=None,metavar='path',required=True)
parser.add_argument("-roi", help="list of mhd ROIs",nargs='+')
parser.add_argument("-dose",help="dose mhd map",nargs=None,metavar='path')
parser.add_argument("-fileLabel", help="label of the txt file in output", type=str, default="")

# parser.add_argument("-o", help="output file path",required=True,nargs=1,metavar='path')

parser.add_argument("-noDVH", help="do not compute DVHs",action="store_true")

group = parser.add_mutually_exclusive_group()
group.add_argument("--norm-volume", help="express volume in %% => DVH units = [cGy,%%]",action="store_true")
group.add_argument("--norm-dose", help="express dose in %% => DVH units = [%%,cm^3]",action="store_true")
group.add_argument("--norm", help="express dose and volume in %% => DVH units = [%%,%%]",action="store_true")

parser.add_argument("-Dgoal", help="reference planned dose [cGy]",type=float,default=100)


parser.add_argument("-v", "--verbose", help="increase output verbosity",action="store_true")
parser.add_argument("--markers", help="use markers for lineplots",action="store_true")

group = parser.add_mutually_exclusive_group()
group.add_argument("--noreadline", help="do not use built-in python readline",action="store_false",dest='readline',default=False)
group.add_argument("--readline", help="use built-in python readline",action="store_true")

parser.add_argument("-figSize", help="figure canvas size",type=float,default=7)


args = parser.parse_args()
# print(args) ; exit()
########################################################################
import os, sys, shutil, time,re
from math import *
import numpy as np
import struct
import matplotlib
import matplotlib.colors as colors
import matplotlib.ticker as ticker

matplotlib.use('TkAgg')
import warnings
warnings.filterwarnings("ignore")

# matplotlib.use('Qt5Agg')
# print matplotlib.get_backend()
import pylab as plt
from scipy import ndimage
# exit(0)
if args.readline:
    import readline

import itertools

from mhd_io import *

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
def isEyeTransform(transform):
    return np.sum(np.abs(transform.flatten()-np.array([1,0,0,0,1,0,0,0,1])))<1e-5
########################################################################
########################################################################
########################################################################
##########                                                    ##########
##########                  GLOBAL VARIABLES                  ##########
##########                                                    ##########
########################################################################
CT = None
Dose = None
ROIs = []
Maps = []
lverbose = args.verbose
lmarkers = args.markers
marker = itertools.cycle(('.','+' , 'o', '*')) 

binningDVH=100
opacity = 0.5 # dose opacity
roiOpacity = 0.5
ptvContourColor = 'black'
dvhUnits=r'D%'
if args.norm: dvhUnits=r'%%'
if args.norm_volume: dvhUnits=r'D%'
if args.norm_dose: dvhUnits=r'%V'
Dgoal = args.Dgoal
fileLabel = args.fileLabel

########################################################################
########################################################################
########################################################################

def clampToRange(x,xmin,xmax):
    x[x<xmin]=xmin[x<xmin]
    x[x>xmax]=xmax[x>xmax]
    return x
########################################################################
def on_click(event):
    global xnow
    global ax1,ax2,ax3,ax4


    if event.inaxes is not None:
        if event.inaxes==ax1:
            xnow[0]=event.xdata; xnow[1]=event.ydata
        if event.inaxes==ax2:
            xnow[1]=event.xdata; xnow[2]=event.ydata
        if event.inaxes==ax3:
            xnow[0]=event.xdata; xnow[2]=event.ydata
            
        clampToRange(xnow,xmin+hs/2,xmax-hs/2)
        draw()
        # print(event.xdata, event.ydata)
        # print(xnow)
    else:
        pass
#         print('Clicked ouside axes bounds but inside plot window')
########################################################################
def on_key(event):
    global xnow,xmin,xmax,xmin_full,xmax_full,nplot
    global ax1,ax2,ax3
    global lLogMap,lDoseLevels
    global lReverseX,lReverseY,lReverseZ

    sgn = np.ones(3)

    if lReverseX:
        sgn[0]=-1.

    if lReverseY:
        sgn[1]=-1.

    if lReverseZ:
        sgn[2]=-1.

    
    # print(('you pressed', event.key, event.xdata, event.ydata))
    
    if event.key=='q':
        exit(0)

    if event.key=='i': # and event.inaxes is not None:
        l = xmax-xmin
        xmean = (xmax+xmin)*0.5
        l = l * 0.5
        xmin = xmean - l/2
        xmax = xmean + l/2
        draw()
    if event.key=='o': # and event.inaxes is not None:
        l = xmax-xmin
        xmean = (xmax+xmin)*0.5
        l = l / 0.5
        xmin = xmean - l/2
        xmax = xmean + l/2
        draw()

    if event.key=='up': 
        if event.inaxes==ax1:
            xnow[1]+=sgn[1]*hs[1]; draw()
        if event.inaxes==ax2:
            xnow[2]+=sgn[2]*hs[2]; draw()
        if event.inaxes==ax3:
            xnow[2]+=sgn[2]*hs[2]; draw()
    if event.key=='down': 
        if event.inaxes==ax1:
            xnow[1]-=sgn[1]*hs[1]; draw()
        if event.inaxes==ax2:
            xnow[2]-=sgn[2]*hs[2]; draw()
        if event.inaxes==ax3:
            xnow[2]-=sgn[2]*hs[2]; draw()
    if event.key=='left': 
        if event.inaxes==ax1:
            xnow[0]-=sgn[0]*hs[0]; draw()
        if event.inaxes==ax2:
            xnow[1]+=sgn[1]*hs[1]; draw()
        if event.inaxes==ax3:
            xnow[0]-=sgn[0]*hs[0]; draw()
    if event.key=='right': 
        if event.inaxes==ax1:
            xnow[0]+=sgn[0]*hs[0]; draw()
        if event.inaxes==ax2:
            xnow[1]-=sgn[1]*hs[1]; draw()
        if event.inaxes==ax3:
            xnow[0]+=sgn[0]*hs[0]; draw()

    if event.key=='pageup':
        l = (xmax-xmin)/4
        if event.inaxes==ax1:
            xmin[1]+=sgn[1]*l[1]; xmax[1]+=sgn[1]*l[1]; draw()
        if event.inaxes==ax2:
            xmin[2]+=sgn[2]*l[2]; xmax[2]+=sgn[2]*l[2]; draw()
        if event.inaxes==ax3:
            xmin[2]+=sgn[2]*l[2]; xmax[2]+=sgn[2]*l[2]; draw()
    if event.key=='pagedown': 
        l = (xmax-xmin)/4
        if event.inaxes==ax1:
            xmin[1]-=sgn[1]*l[1]; xmax[1]-=sgn[1]*l[1]; draw()
        if event.inaxes==ax2:
            xmin[2]-=sgn[2]*l[2]; xmax[2]-=sgn[2]*l[2]; draw()
        if event.inaxes==ax3:
            xmin[2]-=sgn[2]*l[2]; xmax[2]-=sgn[2]*l[2]; draw()
    if event.key=='home': 
        l = (xmax-xmin)/4
        if event.inaxes==ax1:
            xmin[0]-=sgn[0]*l[0]; xmax[0]-=sgn[0]*l[0]; draw()
        if event.inaxes==ax2:
            xmin[1]+=sgn[1]*l[1]; xmax[1]+=sgn[1]*l[1]; draw()
        if event.inaxes==ax3:
            xmin[0]-=sgn[0]*l[0]; xmax[0]-=sgn[0]*l[0]; draw()
    if event.key=='end': 
        l = (xmax-xmin)/4
        if event.inaxes==ax1:
            xmin[0]+=sgn[0]*l[0]; xmax[0]+=sgn[0]*l[0]; draw()
        if event.inaxes==ax2:
            xmin[1]-=sgn[1]*l[1]; xmax[1]-=sgn[1]*l[1]; draw()
        if event.inaxes==ax3:
            xmin[0]+=sgn[0]*l[0]; xmax[0]+=sgn[0]*l[0]; draw()

                    
    if event.key=='u':
        resetView()
        # xmin=xmin_full
        # xmax=xmax_full
        draw()
    if event.key=='c':
        l = xmax-xmin
        if event.inaxes is not None:
            if event.inaxes==ax1:
                xmin[0]=event.xdata-l[0]/2; xmax[0]=event.xdata+l[0]/2; 
                xmin[1]=event.ydata-l[1]/2; xmax[1]=event.ydata+l[1]/2; 
            if event.inaxes==ax2:
                xmin[1]=event.xdata-l[1]/2; xmax[1]=event.xdata+l[1]/2; 
                xmin[2]=event.ydata-l[2]/2; xmax[2]=event.ydata+l[2]/2; 
            if event.inaxes==ax3:
                xmin[0]=event.xdata-l[0]/2; xmax[0]=event.xdata+l[0]/2; 
                xmin[2]=event.ydata-l[2]/2; xmax[2]=event.ydata+l[2]/2; 
            draw()

    if event.key=='y':
        lLogMap = not lLogMap
        draw()

    if event.key=='l':
        lDoseLevels = not lDoseLevels
        draw()

    if event.key=='0':
        xnow = xnow*0
        draw()

    if event.key=='h':
        help()

    if event.key=='w':
        restoreEnv()
        draw()

    if event.key=='W':
        saveEnv()

    if event.key=='r':
        if event.inaxes is None:
            return
        ax = getCloserAxis(event.xdata,event.ydata,event.inaxes)
        # print('closer axis',ax)
        if ax==None:
            return
        if ax == 'x':
            if event.inaxes==ax1:
                lReverseX = not lReverseX
            if event.inaxes==ax2:
                lReverseY = not lReverseY
            if event.inaxes==ax3:
                lReverseX = not lReverseX
        if ax == 'y':
            if event.inaxes==ax1:
                lReverseY = not lReverseY
            if event.inaxes==ax2:
                lReverseZ = not lReverseZ
            if event.inaxes==ax3:
                lReverseZ = not lReverseZ
        draw()

    if event.key=='p':
        nplot +=1
        fname = 'plot%02d.pdf' % nplot
        plt.savefig(fname)
        print('Saved plot to file ',fname)

    if event.key=='P':
        nplot +=1
        fname = 'plot%02d.png' % nplot
        plt.savefig(fname)
        print('Saved plot to file ',fname)

########################################################################
def help():
    print("#############################################################################")
    print("input line commands: ")
    print("h: this help")
    print("q: quit")
    print("return: reload from files")
    print("p, pdf: save plot as pdf")
    print("P, PNG: save plot as png")
    print("iso: move position to isocenter, i.e. (0,0,0)")
    print("x pos: move current position to x = pos")
    print("y pos: move current position to y = pos")
    print("z pos: move current position to z = pos")
    print("cm name: change colormap to name [e.g. bone, jet, hot, ...]")
    print("l , levels: show dose using TPS-like level colormap [on/off]")
    print("y: toggle lin/log color map")
    
    print('')
    print("interactive commands (mouse and keyboard events): ")
    print("h: this help")
    print("q: quit")
    # print("e: reload from files")
    print("i: zoom in")
    print("o: zoom out")
    print("arrows up,down,left,right: move current position by one voxel in the given direction")
    print("u: show all")
    print("c: center on cursor position")
    print("0: move position to ISO, i.e. (0,0,0)")
    print("r: reverse closer axis")
    print("l: show dose using TPS-like level colormap [on/off]")
    
    print("y: toggle lin/log color map")
    print("#############################################################################")
    
########################################################################
def calcDVH():
    global Dose,ROIs,hs
    global lverbose
    global Dgoal,voxelVolume
    global fileLabel

    for roi in ROIs:

        D = Dose['Map'][np.where(roi['Map']>0)].flatten() # dose values in the ROI

        D *= 100. # from Gy (standard fred mhd) to cGy (standard DVH visualization)
       

        myLabel=roi['name']
        maxDose = np.max(D)
        nroi = len(D)
        roiVolume = 1.0*nroi*voxelVolume
        roi['volume']=roiVolume
        print('Calcuting DVH for',roi['name'])
        if lverbose:
            print('planned dose',Dgoal,'cGy')            
            print('maxDose',maxDose,'cGy')
            print('ROI num voxels',nroi)
            print('ROI volume %.1f cm^3' % (roiVolume))
        binSize = Dgoal/binningDVH
        nbin = max(binningDVH,int(maxDose/binSize)+1)
        if lverbose:
            print('nbin',nbin)

        row = [roi['name'],np.min(D),np.max(D),np.mean(D),np.std(D),roiVolume]
        roi['stats']=row

        bins = (np.arange(nbin)+0.5)*binSize
        hist = np.zeros(nbin)

        for d in D:
            imax = int(d/binSize)
            hist[0:imax+1]+=1

        hist*=voxelVolume

        roi['DVH']=(bins,hist)
        roiLabel=myLabel.split(r'/')[-1]
        print(roiLabel)  
        #data = np.column_stack([bins/Dgoal*100, hist/roiVolume*100])    
        np.savetxt(str(roiLabel)+fileLabel+".txt", np.c_[bins, hist/roiVolume*100], fmt='%.3e', delimiter="  ")
        
        # print(bins)

def roiStats():
    global Dose,ROIs,hs
    global lverbose
    global Dgoal,voxelVolume

    table=[['ROI','Min','Max','Mean','Stdev','Volume (cm^3)']]
    for roi in ROIs:
        table.append(roi['stats'])

    print('----------------------------------------------------------------------------------------------------')
    print('ROI Statistics\n')
    for i,row in enumerate(table):
        line='%-20s' % str(row[0])
        for v in row[1:]:
            if i==0:
                line+='%-20s' % v
            else:
                line+='%-20.1f' % v
        print(line)
        if i==0: print('')
    print('----------------------------------------------------------------------------------------------------')
    print('')    
########################################################################
def loadMaps():
    global CT,Dose,ROIs,Maps
    global nn,hs,x0
    global voxelVolume

    # remove possible duplicates from ROI list preserving order
    if args.roi is not None:
        for i in range(len(args.roi)):
            for j in reversed(range(i+1,len(args.roi))):
                if args.roi[i]==args.roi[j]:
                    del args.roi[j]
        if args.CT in args.roi: args.roi.remove(args.CT)

        if args.dose is not None and args.dose in args.roi: args.roi.remove(args.dose)


    print('Loading CT map')
    CT = mhd_read(args.CT)
    if not isEyeTransform(CT['TransformMatrix']):
        print('error: TransformMatrix != Identity not yet implemented'); exit(1)
    Maps.append({'name':'CT','map':CT,'cmap':'gray'})
    Maps[-1]['normColors'] = colors.Normalize(vmin=np.min(Maps[-1]['map']['Map']), vmax=np.max(Maps[-1]['map']['Map']))

    nn,hs,x0=(CT['nn'],CT['hs'],CT['x0'])
    voxelVolume = np.prod(hs)

    if args.dose is not None:
        print('Loading dose map')
    Dose = mhd_read(args.dose) if args.dose else None
    if Dose is not None:
        if not isEyeTransform(Dose['TransformMatrix']):
            print('error: TransformMatrix != Identity not yet implemented'); exit(1)
        if np.any(nn!=Dose['nn']):
            print('error: CT and Dose dimensions do not match')
            exit(1)
        if np.any(abs(hs-Dose['hs'])>1e-4*hs):
            print('warning: CT and Dose spacing are not the same')
        if np.any(abs(x0-Dose['x0'])>1e-4):
            print('warning: CT and Dose origin are not the same')

        Maps.append({'name':'Dose','map':Dose,'cmap':'jet'})
        try:
            Maps[-1]['normColors'] = colors.Normalize(vmin=np.min(Maps[-1]['map']['Map']), vmax=np.max(Maps[-1]['map']['Map']))
            Map = Dose['Map']
            arr = Map[np.where(Map>0.0)].flatten()
            mapGlobalLogMin = log10(np.min(arr))
            mapGlobalLogMax = log10(np.max(arr))
            # print(mapGlobalLogMin,mapGlobalLogMax)
            Maps[-1]['normColorsLog'] = colors.Normalize(vmin=mapGlobalLogMin, vmax=mapGlobalLogMax)
        except:
            pass
    # print(Maps[-1])
    # exit(0)

    ROIs = []
    if args.roi is not None:
        for fname in args.roi:
            print('Loading roi map',fname[:-4])

            roi = mhd_read(fname)
            if np.any(nn!=roi['nn']):
                print('error: CT and ROI dimensions do not match')
                exit(1)
            if np.any(abs(hs-roi['hs'])>1e-4*hs):
                print('warning: CT and ROI spacing are not the same')
            if np.any(abs(x0-roi['x0'])>1e-4):
                print('warning: CT and ROI origin are not the same')

            icol = len(ROIs)%9
            cmap = matplotlib.cm.get_cmap('Set1')
            roi['color']=cmap(icol/8.0)
            roi['Map']=roi['Map']*(1+icol)
            roi['name']=fname[:-4]
            ROIs.append(roi)
            if not isEyeTransform(ROIs[-1]['TransformMatrix']):
                print('error: TransformMatrix != Identity not yet implemented'); exit(1)
            Maps.append({'name':fname[:-4],'map':roi,'cmap':'Set1'})
            Maps[-1]['normColors'] = colors.Normalize(vmin=1, vmax=9)


    if Dose is not None and not args.noDVH:
        calcDVH()
        roiStats()
    # exit(0)

    resetView()
########################################################################
def getCentroid(roi):
    global nn,hs,x0
    xc = np.arange(nn[0])*hs[0]+x0[0]+hs[0]/2
    yc = np.arange(nn[1])*hs[1]+x0[1]+hs[1]/2
    zc = np.arange(nn[2])*hs[2]+x0[2]+hs[2]/2
    Map = roi['Map']
    I = np.where(Map>0)
    ixmin = np.min(I[0])
    ixmax = np.max(I[0])
    iymin = np.min(I[1])
    iymax = np.max(I[1])
    izmin = np.min(I[2])
    izmax = np.max(I[2])
    Px = (xc[ixmin]+xc[ixmax])/2
    Py = (yc[iymin]+yc[iymax])/2
    Pz = (zc[izmin]+zc[izmax])/2
    # print(ixmin,ixmax,)
    # print(Px,Py,Pz)
    return np.array([Px,Py,Pz])

########################################################################
def getCloserAxis(xdata,ydata,axes):
    fracx = abs((xdata-axes.viewLim.x0)/(axes.viewLim.x1-axes.viewLim.x0)-0.5)
    fracy = abs((ydata-axes.viewLim.y0)/(axes.viewLim.y1-axes.viewLim.y0)-0.5)
    threshold = 0.2
    if fracx<0.5-threshold and  fracy<0.5-threshold:
        return  None
    if fracx>fracy:
        return 'y'
    else:
        return 'x'

########################################################################
class Formatter(object):
    def __init__(self, im, islice):
        self.im = im
        self.islice = islice
    def __call__(self, x, y):
#         [xl,xr,yb,yt] = self.im.get_extent()
#         # print('extent = ',self.im.get_extent())
#         # print('shape = ',self.im.get_array().shape)
#         [ny,nx] = self.im.get_array().shape
#         [hx,hy] = [abs(xr-xl)/nx,abs(yt-yb)/ny]
# #       print('spacing = ',[hx,hy])
# #       print('x,y',x,y)
#         [i,j] = [int(abs(y-yb)/hy), int(abs(x-xl)/hx)]
#         if xl>xr:
#             j=nx-1-j
#         if yb>yt:
#             i=ny-1-i
#         # print('i,j',i,j)
#         val = self.im.get_array()[i,j]
#         # val = 0
        if self.islice==1:
            return 'x=%.01f, y=%.01f' % (x, y)
        if self.islice==2:
            return 'y=%.01f, z=%.01f' % (x, y)
        if self.islice==3:
            return 'x=%.01f, z=%.01f' % (x, y)
#       
#       
#       return 'x y = %f %f' % (x,y)
########################################################################
def resetView():
    global nn,hs,x0
    global xmin,xmax,xnow
    global xmin_full,xmax_full
    
    xnow = x0+nn*hs/2.
    L = np.ones(3)*np.max(hs*nn)

    xmin =xnow - L/2
    xmax =xnow + L/2

    xmin_full=x0*1
    xmax_full=x0+hs*nn
    # print(hs*nn,L)
    # print(xmin,xmax)
    # print(xmin_full,xmax_full)
    # print(x0)
    # print(xnow)
    # exit(0)

########################################################################
def myLog_formatter(x, pos):
    return "%.2g" % pow(10,x)
########################################################################
def draw():
    global lfirst,lLineout,lLogMap,cmapname
    global lDoseLevels,clevels_cmap,clevels_norm,clevels_bounds
    global lReverseX,lReverseY,lReverseZ
    global nn,hs,x0,Maps
    global xmin,xmax,xnow
    global xmin_full,xmax_full
    global ax1,ax2,ax3,cbaxes
    global opacity,marker
    global Dgoal,voxelVolume
    
    xc = np.arange(nn[0])*hs[0]+x0[0]+hs[0]/2
    yc = np.arange(nn[1])*hs[1]+x0[1]+hs[1]/2
    zc = np.arange(nn[2])*hs[2]+x0[2]+hs[2]/2
#   print(zc)



    # if lfirst:
    #     resetView()

    idx = ((xnow-x0)/hs).astype('int32')
    nx=nn[0];ny=nn[1];nz=nn[2]
    idx[np.where(nn==1)]=0
    xnow = x0+idx*hs+hs/2

    # print('\nCurrent position=',xnow,'index=',idx,'value=',Map[idx[0],idx[1],idx[2]])
    
    plt.figure(1)
    plt.clf()
    plt.style.use('dark_background')    
    fig = plt.gcf()
    rect = fig.patch
    rect.set_facecolor('black')

    fig.canvas.mpl_connect('button_press_event', on_click)
    fig.canvas.mpl_connect('key_press_event', on_key)
    
    plt.subplot(221)
    ax1=plt.gca()
    ax1.set_facecolor('black')

    for m in Maps:
        Map = m['map']['Map']
        slice_xy=np.copy(Map[:,:,idx[2]]).transpose().squeeze()
        if m['name'] != 'CT':
            slice_xy = np.ma.masked_where(slice_xy <=0 , slice_xy)

        if lLogMap and m['name']=='Dose':
            I=np.where(slice_xy>0)
            if len(I[0])>0:         
                minval = np.min(slice_xy[np.where(slice_xy>0)])
                slice_xy[np.where(slice_xy<=0)] = minval
                slice_xy=np.log10(slice_xy)   

        Extent = [xmin_full[0],xmax_full[0],xmax_full[1],xmin_full[1]]
        mycmap = m['cmap']

        myopacity = roiOpacity
        if m['name']=='Dose':
            myopacity = opacity
        if m['name']=='CT':
            myopacity = 1.0

        if m['name'][0:3].upper()=='PTV' and (Dose is not None and opacity>0):
            maxval = np.max(slice_xy)
            if maxval>0:
                img = slice_xy/(Dgoal/100.)
                smoothed = ndimage.uniform_filter(img,size=1)
                levels = [0.5,1.0,1.5]
                ax1.contour(smoothed, levels, linewidths=1, linestyles='dotted',colors=ptvContourColor, origin='upper', extent=Extent)
        else:
            if lLogMap and 'normColorsLog' in m.keys():
                im1 = ax1.imshow(slice_xy,interpolation='nearest', norm = m['normColorsLog'],extent=Extent,aspect=1,cmap=mycmap,alpha=myopacity)
            else:
                if lDoseLevels and m['name']=='Dose':
                    img = slice_xy/(Dgoal/100.)
                    smoothed = ndimage.uniform_filter(img,size=4)
                    smoothed = np.ma.masked_where(smoothed <=0 , smoothed)
                    ax1.contourf(smoothed,clevels_bounds, linestyles='dotted',
                        colors=clevels_cmap(clevels_bounds), origin='upper', extent=Extent,alpha=myopacity)
                    # im1 = ax1.imshow(slice_xy/(Dgoal/100.),interpolation='nearest',norm = clevels_norm,extent=Extent,aspect=1,cmap=clevels_cmap,alpha=myopacity)
                else:
                    im1 = ax1.imshow(slice_xy,interpolation='nearest',norm = m['normColors'],extent=Extent,aspect=1,cmap=mycmap,alpha=myopacity)


        # im1B = ax1.imshow(sliceCT_xy,interpolation='nearest', norm = mapNormB,extent=ExtentB ,aspect=1,cmap='bone')
        # im1  = ax1.imshow(slice_xy,  interpolation='nearest', norm = mapNorm ,extent=Extent  ,aspect=1,cmap=topcmap,alpha=opacity)

        ax1.set_xlim(xmin[0],xmax[0]) if not lReverseX else ax1.set_xlim(xmax[0],xmin[0]) 
        ax1.set_ylim(xmin[1],xmax[1]) if not lReverseY else ax1.set_ylim(xmax[1],xmin[1]) 

    # print('full range',Extent)
    # print('view range',xmin[0],xmax[0],xmin[1],xmax[1])
    # plt.title('XY slice at z=%.2f' % xnow[2])
    # plt.xlabel('x (cm)')
    # plt.ylabel('y (cm)')
    ax1.axvline(xnow[0], color='blue', lw=1.0, alpha=1.0,dashes=(1,1),antialiased=False)
    ax1.axhline(xnow[1], color='blue', lw=1.0, alpha=1.0,dashes=(1,1),antialiased=False)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    ax1.format_coord = Formatter(im1,1)

    plt.text(0.03, 0.5, 'R',transform=ax1.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.97, 0.5, 'L',transform=ax1.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.5, 0.03, 'P',transform=ax1.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.5, 0.97, 'A',transform=ax1.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    
    plt.subplot(222)
    ax2=plt.gca()
    ax2.set_facecolor('black')
    for m in Maps:
        Map = m['map']['Map']
        slice_yz=np.copy(Map[idx[0],:,:]).transpose().squeeze()
        slice_yz = slice_yz[::-1,:]
        if m['name'] != 'CT':
            slice_yz = np.ma.masked_where(slice_yz <=0 , slice_yz)

        if lLogMap and m['name']=='Dose':
            I=np.where(slice_yz>0)
            if len(I[0])>0:         
                minval = np.min(slice_yz[np.where(slice_yz>0)])
                slice_yz[np.where(slice_yz<=0)] = minval
                slice_yz=np.log10(slice_yz)   

        Extent = [xmin_full[1],xmax_full[1],xmin_full[2],xmax_full[2]]
        mycmap = m['cmap']

        myopacity = roiOpacity
        if m['name']=='Dose':
            myopacity = opacity
        if m['name']=='CT':
            myopacity = 1.0

        if m['name'][0:3].upper()=='PTV' and (Dose is not None and opacity>0):
            maxval = np.max(slice_yz)
            if maxval>0:
                img = slice_yz/(Dgoal/100.)
                smoothed = ndimage.uniform_filter(img,size=1)
                levels = [0.5,1.0,1.5]
                ax2.contour(smoothed, levels, linewidths=1, linestyles='dotted',colors=ptvContourColor, origin='upper', extent=Extent)
        else:
            if lLogMap and 'normColorsLog' in m.keys():
                im2 = ax2.imshow(slice_yz,interpolation='nearest', norm = m['normColorsLog'],extent=Extent,aspect=1,cmap=mycmap,alpha=myopacity)
            else:
                if lDoseLevels and m['name']=='Dose':
                    img = slice_yz/(Dgoal/100.)
                    smoothed = ndimage.uniform_filter(img,size=4)
                    smoothed = np.ma.masked_where(smoothed <=0 , smoothed)
                    ax2.contourf(smoothed,clevels_bounds, linestyles='dotted',
                        colors=clevels_cmap(clevels_bounds), origin='upper', extent=Extent,alpha=myopacity)
                    # im2 = ax2.imshow(slice_yz/(Dgoal/100.),interpolation='nearest',norm = clevels_norm,extent=Extent,aspect=1,cmap=clevels_cmap,alpha=myopacity)
                else:
                    im2 = ax2.imshow(slice_yz,interpolation='nearest',norm = m['normColors'],extent=Extent,aspect=1,cmap=mycmap,alpha=myopacity)


    # plt.title('ZY slice at x=%.2f' % xnow[0])
    # plt.xlabel('z (cm)')
    # plt.ylabel('y (cm)')
    
    ax2.axvline(xnow[1], color='blue', lw=1.0, alpha=1.0,dashes=(1,1),antialiased=False)
    ax2.axhline(xnow[2], color='blue', lw=1.0, alpha=1.0,dashes=(1,1),antialiased=False)

    ax2.set_ylim(xmin[2],xmax[2]) if not lReverseZ else ax2.set_ylim(xmax[2],xmin[2]) 
    ax2.set_xlim(xmax[1],xmin[1]) if not lReverseY else ax2.set_xlim(xmin[1],xmax[1]) 
    
    ax2.format_coord = Formatter(im2,2)
    ax2.axes.get_xaxis().set_visible(False)
    ax2.axes.get_yaxis().set_visible(False)

    plt.text(0.03, 0.5, 'A',transform=ax2.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.97, 0.5, 'P',transform=ax2.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.5, 0.03, 'I',transform=ax2.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.5, 0.97, 'S',transform=ax2.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)


    plt.subplot(224)
    ax3=plt.gca()
    ax3.set_facecolor('black')

    # Now adding the colorbar
    cbaxes = fig.add_axes([0.93, 0.55, 0.02, 0.40])     

    for m in Maps:
        Map = m['map']['Map']
        slice_xz=np.copy(Map[:,idx[1],:]).transpose().squeeze()
        if m['name'] != 'CT':
            slice_xz = np.ma.masked_where(slice_xz <=0 , slice_xz)

        if lLogMap and m['name']=='Dose':
            I=np.where(slice_xz>0)
            if len(I[0])>0:         
                minval = np.min(slice_xz[np.where(slice_xz>0)])
                slice_xz[np.where(slice_xz<=0)] = minval
                slice_xz=np.log10(slice_xz)   

        Extent = [xmin_full[0],xmax_full[0],xmax_full[2],xmin_full[2]]
        mycmap = m['cmap']

        myopacity = roiOpacity
        if m['name']=='Dose':
            myopacity = opacity
        if m['name']=='CT':
            myopacity = 1.0

        if m['name'][0:3].upper()=='PTV' and (Dose is not None and opacity>0):
            maxval = np.max(slice_xz)
            if maxval>0:
                img = slice_xz/(Dgoal/100.)
                smoothed = ndimage.uniform_filter(img,size=1)
                levels = [0.5,1.0,1.5]
                ax3.contour(smoothed, levels, linewidths=1, linestyles='dotted',colors=ptvContourColor, origin='upper', extent=Extent)
        else:
            if lLogMap and 'normColorsLog' in m.keys():
                im3 = ax3.imshow(slice_xz,interpolation='nearest', norm = m['normColorsLog'],extent=Extent,aspect=1,cmap=mycmap,alpha=myopacity)
                cbar = fig.colorbar(im3, cax = cbaxes,format=matplotlib.ticker.FuncFormatter(myLog_formatter), ticklocation='left')
            else:
                if lDoseLevels and m['name']=='Dose':
                    img = slice_xz/(Dgoal/100.)
                    smoothed = ndimage.uniform_filter(img,size=4)
                    smoothed = np.ma.masked_where(smoothed <=0 , smoothed)
                    CS = ax3.contourf(smoothed,clevels_bounds, linestyles='dotted',
                        colors=clevels_cmap(clevels_bounds), origin='upper', extent=Extent,alpha=myopacity)
                    cbar = fig.colorbar(CS, cax = cbaxes, ticklocation='left',alpha=1.0, ticks=clevels_bounds)
                        # cmap = clevels_cmap, values = [clevels_bounds[0],clevels_bounds[-1]])
                    # im3 = ax3.imshow(slice_xz/(Dgoal/100.),interpolation='nearest',norm = clevels_norm,extent=Extent,aspect=1,cmap=clevels_cmap,alpha=myopacity)
                else:
                    im3 = ax3.imshow(slice_xz,interpolation='nearest',norm = m['normColors'],extent=Extent,aspect=1,cmap=mycmap,alpha=myopacity)
                    if m['name']=='Dose':
                        cbar = fig.colorbar(im3, cax = cbaxes, ticklocation='left', drawedges=False)

    # plt.title('XZ slice at y=%.2f' % xnow[1])
    # plt.xlabel('x (cm)')
    # plt.ylabel('z (cm)')
    ax3.axvline(xnow[0], color='blue', lw=1.0, alpha=1.0,dashes=(1,1),antialiased=False)
    ax3.axhline(xnow[2], color='blue', lw=1.0, alpha=1.0,dashes=(1,1),antialiased=False)

    ax3.set_xlim(xmin[0],xmax[0]) if not lReverseX else ax3.set_xlim(xmax[0],xmin[0]) 
    ax3.set_ylim(xmin[2],xmax[2]) if not lReverseZ else ax3.set_ylim(xmax[2],xmin[2]) 

    ax3.format_coord = Formatter(im3,3)

    ax3.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)

    plt.text(0.03, 0.5, 'R',transform=ax3.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.97, 0.5, 'L',transform=ax3.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.5, 0.03, 'I',transform=ax3.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)
    plt.text(0.5, 0.97, 'S',transform=ax3.transAxes,color='orange',horizontalalignment='center',verticalalignment='center',fontsize=10)

    if not args.noDVH and Dose is not None and len(ROIs):    
        plt.subplot(223)
        ax4=plt.gca()
        # ax4.set_facecolor('blue')
        ax4.set_facecolor('black')
        ax4.axes.get_xaxis().set_visible(True)
        ax4.axes.get_yaxis().set_visible(True)
        for roi in ROIs:
            if 'DVH' not in roi.keys():
                continue 

            myLabel = roi['name']
            myColor = roi['color'] 

            origx,origy = roi['DVH']
            x=origx.copy() # deep copy
            y=origy.copy()

            if dvhUnits[0]=='%':
                x*=100./Dgoal
            if dvhUnits[1]=='%':
                y*=100/roi['volume']

                      
            if lmarkers:
                plt.plot(x,y,label=myLabel,marker = next(marker),lw=0.4,color=myColor,markevery=10)
            else:
                plt.plot(x,y,label=myLabel,lw=2.0,color=myColor)

        if dvhUnits==r'%%':
            plt.xlabel('Dose [%]')
            plt.ylabel('Volume [%]')
            plt.plot([95.],[95.],'r*',lw=3,label=r'95% 95%',markerfacecolor='r',
                 markeredgewidth=0.7, markeredgecolor='yellow',markersize=8)
            plt.xticks(np.arange(0,120+1, 10))
            plt.yticks(np.arange(0,100+1, 10))
        elif dvhUnits==r'D%':
            plt.xlabel('Dose [cGy]')
            plt.ylabel('Volume [%]')
            plt.plot([0.95*Dgoal],[95.],'r*',lw=3,label=r'95% 95%',markerfacecolor='r',
                 markeredgewidth=0.7, markeredgecolor='yellow',markersize=8)
            plt.yticks(np.arange(0,100+1, 10))
        elif dvhUnits==r'%V':
            plt.xlabel('Dose [%]')
            plt.ylabel('Volume [cm^3]')
            plt.xticks(np.arange(0,120+1, 10))
        else:
            plt.xlabel('Dose [cGy]')
            plt.ylabel('Volume [cm^3]')

        # plt.legend(loc='upper right',bbox_to_anchor=(0.9, 1.55))
        fig.legend(bbox_to_anchor=(0.70, 0.6))
        plt.grid(color='lightgray', linestyle='-.', linewidth=0.4)
        ax4.set_axisbelow(True)
        

    if args.noDVH or Dose is None:
        plt.subplot(223)
        ax4=plt.gca()
        ax4.axes.get_xaxis().set_visible(False)
        ax4.axes.get_yaxis().set_visible(False)

        for roi in ROIs:
            myLabel = roi['name']
            myColor = roi['color']
            plt.plot([0,1],[0,1],label=myLabel,lw=2.0,color=myColor)
        if len(ROIs):
            fig.legend(loc='center',bbox_to_anchor=(0.25, 0.25))
            plt.cla()
        

    
    # hsp=0.2
    plt.subplots_adjust(hspace=0,wspace=0,left=0,right=1.0,top=1.0,bottom=0.0)
    fig.tight_layout()
    plt.show()


def saveEnv():
    global xmin,xmax,xnow
    try:
        fout = open('.env','w')
        print(' '.join(list(map(str,xmin))),file=fout)
        print(' '.join(list(map(str,xmax))),file=fout)
        print(' '.join(list(map(str,xnow))),file=fout)
        fout.close()
        print('environment saved')
    except:
        pass

def restoreEnv():
    global xmin,xmax,xnow
    try:
        fin = open('.env')
        lines = fin.readlines()
        xmin = np.array(list(map(float,lines[0].split())))
        xmax = np.array(list(map(float,lines[1].split())))
        xnow = np.array(list(map(float,lines[2].split())))
        fin.close()
        print('environment restored')
    except:
        pass

########################################################################
########################################################################
########################################################################


xnow=np.zeros([3,1])
xmin=np.zeros([3,1])
xmax=np.zeros([3,1])
xmin_full=np.zeros([3,1])
xmax_full=np.zeros([3,1])

if len(sys.argv)<2:
    print('usage: mapFile')
    print('interactive three-slices display of a map')
    exit(1)

lfirst=True
lLineout=False
lLogMap=False
lDoseLevels=False
lReverseX=lReverseY=lReverseZ=False
lReverseY=True
if '-revX' in sys.argv: lReverseX=True
if '-revY' in sys.argv: lReverseY=True
if '-revZ' in sys.argv: lReverseZ=True
nplot=-1


fig = plt.figure(1,figsize=(args.figSize,args.figSize))

#fig.canvas.set_window_title("plan viewer for Fred")
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_key)
# cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8]) 


cmapname='jet'
plt.rcParams['image.cmap'] = cmapname

# make a color map of fixed colors: for TPS-like contour levels
clevels_cmap = colors.ListedColormap(
[[180./255.,51./255.,32./255.],
[152./255.,203./255.,0./255.],
[6./255.,55./255.,198./255.],
[243./255.,124./255.,191./255.],
[52./255.,179./255.,252./255.],
[103./255.,243./255.,79./255.],
[252./255.,246./255.,32./255.],
[250./255.,198./255.,39./255.],
[245./255.,124./255.,44./255.],
[242./255.,65./255.,99./255.]]
)

clevels_bounds=np.array([0.0,0.01,0.05,0.10,0.20,0.30,0.60,0.70,0.80,0.90,1.10])
clevels_norm = colors.BoundaryNorm(clevels_bounds, clevels_cmap.N)

lredraw=True

nplot=0

loadMaps()

for roi in ROIs:
    if (roi['name'][0:3]).upper()=='PTV':
        xnow = getCentroid(roi)

# exit(0)

while True:
    
    plt.ion()

    if lfirst:
        print('# dims=',nn)
        print('# spacing=',hs)
        print('# offset=',x0)
        print('# x-spatial extent= [',x0[0],',',x0[0]+nn[0]*hs[0],'] => ',nn[0]*hs[0])
        print('# y-spatial extent= [',x0[1],',',x0[1]+nn[1]*hs[1],'] => ',nn[1]*hs[1])
        print('# z-spatial extent= [',x0[2],',',x0[2]+nn[2]*hs[2],'] => ',nn[2]*hs[2])
    
    if lredraw:
        draw()
    lredraw = True
        
    if lfirst:
        lfirst=False
    
    # choice = input('Press return to redraw;  input "h" for help or "q" to exit: ')
    choice = input('>>> ')
#   print(choice)
    
    # get new position
    m = re.match(r'\s*([xyz])\s*=\s*(\S+)',choice)
    if m==None:
        m = re.match(r'\s*([xyz])\s+(\S+)',choice)
    if m:
        try:
            if m.groups()[0]=='x':
                xnow[0] = float(m.groups()[1])
            if m.groups()[0]=='y':
                xnow[1] = float(m.groups()[1])
            if m.groups()[0]=='z':
                xnow[2] = float(m.groups()[1])
        except:
            pass

    # get new index
    m = re.match(r'\s*i([xyz])\s*=\s*(\S+)',choice)
    if m==None:
        m = re.match(r'\s*i([xyz])\s+(\S+)',choice)
    if m:
        try:
            if m.groups()[0]=='x':
                xnow[0] = float(x0[0]+(int(m.groups()[1])+0.5)*hs[0])
            if m.groups()[0]=='y':
                xnow[1] = float(x0[1]+(int(m.groups()[1])+0.5)*hs[1])
            if m.groups()[0]=='z':
                xnow[2] = float(x0[2]+(int(m.groups()[1])+0.5)*hs[2])
        except:
            pass

    m = re.match(r'\s*(cm)\s*',choice)
    if m:
        print('available colormaps: ',','.join(plt.colormaps()))
            
    m = re.match(r'\s*(cm)\s+(\S+)',choice)
    if m:
        try:
            cmapname = m.groups()[1]
            if cmapname not in plt.colormaps():
                print('available colormaps: ',','.join(plt.colormaps()))
            else:
                plt.rcParams['image.cmap'] = cmapname
        except:
            pass

    if choice == "iso":
        xnow = xnow*0

    if choice == "q":
        break

    if choice == "W":
        saveEnv()
        break

    if choice == "w":
        restoreEnv()
        break

    if choice == "h":
        help()
        lredraw = False
        
    if choice in ["p","pdf"]:
        nplot +=1
        fname = 'plot%02d.pdf' % nplot
        plt.savefig(fname)
        print("Saved plot to file ",fname)
        # if lLineout:
        #     plt.figure(1)
        #     fname = 'plot%02dL.pdf' % nplot
        #     plt.savefig(fname)          
        #     plt.figure(2)
        lredraw = False

    if choice in ["P","PNG"]:
        nplot +=1
        fname = 'plot%02d.png' % nplot
        plt.savefig(fname)  
        # if lLineout:
        #     plt.figure(1)
        #     fname = 'plot%02dL.png' % nplot
        #     plt.savefig(fname)          
        #     plt.figure(2)
        lredraw = False
        print("Saved plot to file ",fname)

    if choice in ["l","levels"]:
        lDoseLevels = not lDoseLevels

    if choice in ["y"]:
        lLogMap = not lLogMap

    m = re.match(r'\s*opacity\s+(\S+)',choice) or re.match(r'\s*op\s+(\S+)',choice,re.IGNORECASE)
    if m:
        try:
            opacity = float(m.groups()[0])
        except:
            pass
            
    m = re.match(r'\s*roiOpacity\s+(\S+)',choice) or re.match(r'\s*roiOp\s+(\S+)',choice,re.IGNORECASE)
    if m:
        try:
            roiOpacity = float(m.groups()[0])
        except:
            pass
        
    m = re.match(r'\s*Dgoal\s+(\S+)',choice,re.IGNORECASE)
    if m:
        try:
            Dgoal = float(m.groups()[0])
            calcDVH()
        except:
            pass

    m = re.match(r'\s*dvhUnits\s+(\S+)',choice,re.IGNORECASE)
    if m:
        try:
            if (m.groups()[0]).upper() in [r'DV',r'D%',r'%V',r'%%']:
                dvhUnits = (m.groups()[0]).upper()
        except:
            pass

    m = re.match(r'^\s*stats\s*$',choice,re.IGNORECASE)
    if m:
        roiStats()
        lredraw=False

    clampToRange(xnow,xmin+hs/2,xmax-hs/2)
            
plt.ioff()
