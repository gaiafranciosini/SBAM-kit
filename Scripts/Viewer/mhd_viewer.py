#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Oct 2016 A.Schiavi 
# Sep 2019 adapted from m3d to mhd file format
# Jan 2022 A.Schiavi : multiple maps ; new GUI interaction
########################################################################
########################################################################
import argparse
parser = argparse.ArgumentParser(description='mhd viewer for one or more maps with 2D slices and lineouts')
parser.add_argument("map", help="list of mhd maps",nargs='+')
parser.add_argument("-CT",help="CT map",nargs=None,metavar='path')

# parser.add_argument("-p", "--plot", help="interactive visualization",action="store_true")
parser.add_argument("-pdf", help="save image to PDF",action="store_true")
parser.add_argument("-png", help="save image to PNG",action="store_true")

parser.add_argument("-v", "--verbose", help="increase output verbosity",action="store_true")

parser.add_argument("--title",help="plot title",nargs=None,metavar='string')

parser.add_argument("-revX", help="reverse x axis",action="store_true",default=False)
parser.add_argument("-revY", help="reverse y axis",action="store_true",default=False)
parser.add_argument("-revZ", help="reverse z axis",action="store_true",default=False)

parser.add_argument("-DCO", help="dose cutoff (or just value cutoff)",type=float)
parser.add_argument("-opacity", help="map opacity (0-1)",type=float)

parser.add_argument("-legend", help="list of CSV labels for the maps, e.g. 'aa,bb,cc'")

parser.add_argument("-lineout", help="plot profile(s)",action="store_true",default=False)
parser.add_argument("-gridLineout", help="plot grid in profile(s)",action="store_true",default=False)

parser.add_argument("-GI", help="use adapted colormap for gamma-index visualization",action="store_true",default=False)

parser.add_argument("-ISO", help="position of ISO center",nargs=3,type=float,metavar='#',default=[0,0,0])

args = parser.parse_args()
# print(args) ; exit(0)
########################################################################
########################################################################
########################################################################

import os, sys, shutil, time,re
from math import *
import numpy as np
import struct
import matplotlib
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.use('TkAgg')
# matplotlib.use('Qt5Agg')
# print matplotlib.get_backend()
import pylab as plt
# import readline

from mhd_io import *

class Map( object ):
    def __init__( self , fname):
        self.path = fname
        [self.nn,self.hs,self.xmin,self.map] = unpackVoxels(mhd_read(fname))
        self.vmin = np.min(self.map)
        self.vmax = np.max(self.map)
        self.L = self.hs*self.nn
        self.xmax = self.xmin + self.L

        try:
            self.logVmin = log10(np.min(self.map[np.where(self.map>0.0)]))
            self.logVmax = log10(np.max(self.map[np.where(self.map>0.0)]))
        except:
            self.logVmin = None
            self.logVmax = None

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
def clampToRange(x,xmin,xmax):
    x[x<xmin]=xmin[x<xmin]
    x[x>xmax]=xmax[x>xmax]
    return x
########################################################################
def on_click(event):
    global xnow
    global axXY,axYZ,axZX,axLineoutX,axLineoutY,axLineoutZ
    global maps

    # print(event.xdata, event.ydata)

    if event.inaxes is not None and event.dblclick:
        if event.inaxes in axXY:
            xnow[0]=event.xdata; xnow[1]=event.ydata
        if event.inaxes in axYZ:
            xnow[1]=event.xdata; xnow[2]=event.ydata
        if event.inaxes in axZX:
            xnow[2]=event.xdata; xnow[0]=event.ydata
        if event.inaxes==axLineoutX:
            xnow[0]=event.xdata;
        if event.inaxes==axLineoutY:
            xnow[1]=event.xdata;
        if event.inaxes==axLineoutZ:
            xnow[2]=event.xdata;

        draw()
        # print(event.xdata, event.ydata)
        # print(xnow)
    else:
        pass
#         print('Clicked ouside axes bounds but inside plot window')
########################################################################
def on_key(event):
    global xnow,xmin,xmax,xmin_full,xmax_full,nplot
    global axXY,axYZ,axZX,axLineoutX,axLineoutY,axLineoutZ
    global lLogMap,lContLevels,isEmptyMap
    global lReverseX,lReverseY,lReverseZ,lLineout,lRestoreFullView
    global mycmap
    global opacity

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
    if event.key=='e':
        loadMaps()
        draw()

    if event.key=='1':
        opacity = 0.1; draw(); return
    if event.key=='2':
        opacity = 0.2; draw(); return
    if event.key=='3':
        opacity = 0.3; draw(); return
    if event.key=='4':
        opacity = 0.4; draw(); return
    if event.key=='5':
        opacity = 0.5; draw(); return
    if event.key=='6':
        opacity = 0.6; draw(); return
    if event.key=='7':
        opacity = 0.7; draw(); return
    if event.key=='8':
        opacity = 0.8; draw(); return
    if event.key=='9':
        opacity = 0.9; draw(); return
    if event.key=='0':
        opacity = 1.0; draw(); return

    if event.key=='up': 
        if event.inaxes in axXY:
            xnow[1]+=sgn[1]*hs[1];
        if event.inaxes in axYZ:
            xnow[2]+=sgn[2]*hs[2];
        if event.inaxes in axZX:
            xnow[0]+=sgn[0]*hs[0];
        draw()

    if event.key=='down': 
        if event.inaxes in axXY:
            xnow[1]-=sgn[1]*hs[1];
        if event.inaxes in axYZ:
            xnow[2]-=sgn[2]*hs[2];
        if event.inaxes in axZX:
            xnow[0]-=sgn[0]*hs[0];
        draw()

    if event.key=='left': 
        if event.inaxes in axXY:
            xnow[0]-=sgn[0]*hs[0];
        if event.inaxes in axYZ:
            xnow[1]-=sgn[1]*hs[1];
        if event.inaxes in axZX:
            xnow[2]-=sgn[2]*hs[2];
        draw()

    if event.key=='right': 
        if event.inaxes in axXY:
            xnow[0]+=sgn[0]*hs[0];
        if event.inaxes in axYZ:
            xnow[1]+=sgn[1]*hs[1];
        if event.inaxes in axZX:
            xnow[2]+=sgn[2]*hs[2];
        draw()

                    
    if event.key=='u':
        lRestoreFullView = True
        draw()

    if event.key=='y':
        lLogMap = not lLogMap
        draw()

    if event.key=='h':
        help()

    if event.key=='w':
        restoreEnv()
        draw()

    if event.key=='j':
        lLineout = not lLineout
        draw()

    if event.key=='W':
        saveEnv()

    if event.key=='b':
        plt.rcParams['image.cmap'] = mycmap if plt.rcParams['image.cmap']=='gray' else 'gray'
        draw()


    if event.key=='r':
        print(axYZ)
        if event.inaxes is None:
            return
        ax = getCloserAxis(event.xdata,event.ydata,event.inaxes)
        print('closer axis',ax)
        if ax==None:
            return
        if ax == 'x':
            if event.inaxes in axXY:
                lReverseX = not lReverseX
            if event.inaxes in axYZ:
                lReverseY = not lReverseY
            if event.inaxes in axZX:
                lReverseZ = not lReverseZ
        if ax == 'y':
            if event.inaxes in axXY:
                lReverseY = not lReverseY
            if event.inaxes in axYZ:
                lReverseZ = not lReverseZ
            if event.inaxes in axZX:
                lReverseX = not lReverseX
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

# Declare and register callbacks
def on_xlims_change(event_ax):
    global axXY,axYZ,axZX,axLineoutX,axLineoutY,axLineoutZ
    # print("updated xlims: ", event_ax.get_xlim())
    if event_ax in [*axXY,axLineoutX]:
        axZX[0].set_ylim(event_ax.get_xlim())
    if event_ax in [*axYZ,axLineoutY]:
        axXY[0].set_ylim(event_ax.get_xlim())
    if event_ax in [*axZX,axLineoutZ]:
        axYZ[0].set_ylim(event_ax.get_xlim())

def on_ylims_change(event_ax):
    global axXY,axYZ,axZX,axLineoutX,axLineoutY,axLineoutZ
    # print("updated ylims: ", event_ax.get_ylim())
    if event_ax in [*axXY,axLineoutY]:
        axYZ[0].set_xlim(event_ax.get_ylim())
    if event_ax in [*axYZ,axLineoutZ]:
        axZX[0].set_xlim(event_ax.get_ylim())
    if event_ax in [*axZX,axLineoutX]:
        axXY[0].set_xlim(event_ax.get_ylim())

########################################################################
def help():
    print("#############################################################################")
    print("input line commands: ")
    print("h: this help")
    print("q: quit")
    # print("return: reload from files")
    print("p, pdf: save plot as pdf")
    print("P, PNG: save plot as png")
    print("iso: move position to isocenter (if set) or to the origin (0,0,0)")
    print("x pos: move current position to x = pos")
    print("y pos: move current position to y = pos")
    print("z pos: move current position to z = pos")
    print("cm name: change colormap to name [e.g. bone, jet, hot, ...]")
    print("j , lineout: show lineouts through current position")
    print("y: toggle lin/log color map")
    print('')
    print("GI: toggle on/off special color map for gamma-index")
    print("vmin: set minimum value for map visualization")
    print("vmax: set maximum value for map visualization")
    print("DCO, vcutoff: set value cutoff below which maps are not shown")
    

    
    print('')
    print("interactive commands (mouse and keyboard events): ")
    print("h: this help")
    print("q: quit")
    print("e: reload from files")
    print("j: show/hide lineouts")
    print("arrows up,down,left,right: move current position by one voxel in the given direction")
    print("u: show all")
    print("double click: update current position with mouse")
    print("r: reverse closer axis")
    print("y: toggle lin/log color map")
    print("1,2,...,9,0: set opacity to 10%,20%,...,90%,100%")
    print("#############################################################################")
    
    

########################################################################
def loadCT():
    global mapCT
    mapCT= Map(args.CT)
    
########################################################################
def loadMaps():
    global maps
    global cbNorm
    global mapGlobalMin,mapGlobalMax,mapGlobalLogMin,mapGlobalLogMax

    maps = []
    cbNorm = []

    for m in args.map:
        maps.append(Map(m))
        cbNorm.append(colors.Normalize(vmin=maps[-1].vmin, vmax=maps[-1].vmax))

    mapGlobalMin = np.min([m.vmin for m in maps])
    mapGlobalMax = np.max([m.vmax for m in maps])

    try:
        mapGlobalLogMin = np.min([m.logVmin for m in maps])
        mapGlobalLogMax = np.max([m.logVmax for m in maps])
    except:
        mapGlobalLogMin = -1
        mapGlobalLogMax = -1

    # print(mapGlobalMin,mapGlobalMax)
    # print(mapGlobalLogMin,mapGlobalLogMax)
    # exit(0)

########################################################################
def getColorNormalization():
    global cbNorm
    global Vmin,Vmax
    global lGammaIndexMode

    # print('Vmin,Vmax = ',Vmin,Vmax)

    if Vmin is not None:
        vmin = Vmin
        if lLogMap and Vmin>0:
            vmin = log10(Vmin)
        if lGammaIndexMode and Vmin<0:
            vmin = 0
    else:
        vmin = mapGlobalMin
        if lLogMap:
            vmin = mapGlobalLogMin
        if lGammaIndexMode:
            vmin = 0

    if Vmax is not None:
        vmax = Vmax
        if lLogMap and Vmin>0:
            vmax = log10(Vmax)
    else:
        vmax = mapGlobalMax
        if lLogMap:
            vmax = mapGlobalLogMax

    if lGammaIndexMode:
        return GINormalize(vmin=vmin,vmax=vmax)

    return colors.Normalize(vmin=vmin, vmax=vmax)



########################################################################
class GINormalize(colors.Normalize):
    def __init__(self, vmin=0, vmax=None, clip=False):
        self.midpoint = 1.0
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # result, is_scalar = self.process_value(value)
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.array(np.interp(value, x, y), mask=np.ma.masked_less(value,0).mask, copy=False)        

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
        if self.islice==1:
            return 'x=%.01f, y=%.01f' % (x, y)
        if self.islice==2:
            return 'y=%.01f, z=%.01f' % (x, y)
        if self.islice==3:
            return 'z=%.01f, x=%.01f' % (x, y)
        if self.islice==11:
            return 'x=%.01f, v=%.01f' % (x, y)
        if self.islice==12:
            return 'y=%.01f, v=%.01f' % (x, y)
        if self.islice==13:
            return 'z=%.01f, v=%.01f' % (x, y)
########################################################################
def resetView():
    global maps,mapCT
    global xmin_full,xmax_full,xmin,xmax,xnow
    
    xmin_full=maps[0].xmin.copy()
    xmax_full=maps[0].xmax.copy()
    for m in maps:
        xmin_full = np.minimum(xmin_full,m.xmin)
        xmax_full = np.maximum(xmax_full,m.xmax)

    xmin = xmin_full.copy()
    xmax = xmax_full.copy()
    xnow = (xmin+xmax)/2.

    if args.CT != None:
        xmin_full = np.minimum(xmin_full,mapCT.xmin)
        xmax_full = np.maximum(xmax_full,mapCT.xmax)
    xmin = xmin_full.copy()
    xmax = xmax_full.copy()

########################################################################
def myLog_formatter(x, pos):
    return "%.2g" % pow(10,x)

def getSlice(array,slice,idx):
    if slice == 'xy':
        vxls = np.copy(array[:,:,idx[2]]).transpose().squeeze()
    if slice == 'yz':
        vxls = np.copy(array[idx[0],:,:]).transpose().squeeze()
    if slice == 'zx':
        vxls = np.copy(array[:,idx[1],:]).squeeze()
    return vxls

def getLineout(array,dir,idx):
    global xc,yc,zc
    if dir == 'x':
        vxls = np.copy(array[:,idx[1],idx[2]]).squeeze()
        binc = xc.copy()
    if dir == 'y':
        vxls = np.copy(array[idx[0],:,idx[2]]).squeeze()
        binc = yc.copy()
    if dir == 'z':
        vxls = np.copy(array[idx[0],idx[1],:]).squeeze()
        binc = zc.copy()

    return binc,vxls


def getExtent(myMap,slice):
    global maps
    if slice == 'xy':
        ext = [myMap.xmin[0],myMap.xmax[0],myMap.xmax[1],myMap.xmin[1]]
    if slice == 'yz':
        ext = [myMap.xmin[1],myMap.xmax[1],myMap.xmax[2],myMap.xmin[2]]
    if slice == 'zx':
        ext = [myMap.xmin[2],myMap.xmax[2],myMap.xmax[0],myMap.xmin[0]]
    return ext

def getCurrIdx(myMap):
    global xnow
    idx = ((xnow-myMap.xmin)/myMap.hs).astype('int32')
    idx[np.where(myMap.nn==1)]=0
    if np.all(idx>=0) and np.all(idx<myMap.nn): # is in range
        return idx
    else:
        return None



def drawSlice(axes,myMap,slice,Cmap,Norm,masked):
    global mapCT
    global lLogMap

    if args.CT!=None:
        idx = getCurrIdx(mapCT)
        if idx is not None:
            img = getSlice(mapCT.map,slice,idx)
            axes.imshow(img,interpolation='nearest',aspect=1,extent=getExtent(mapCT,slice),cmap='bone')

    idx = getCurrIdx(myMap)
    if idx is None:
        return None

    #get slice
    img = getSlice(myMap.map,slice,idx)
    if masked and args.DCO!=None:
        img = np.ma.masked_where(img<args.DCO,img)

    if lLogMap:
        img = np.log10(img)
   # if lLogMap:
   #          slice_xy[np.where(slice_xy<=0)] = mapGlobalNonZeroMin
   #          slice_xy=np.log10(slice_xy)   
    
    # print('extent',getExtent(myMap,slice))

    im1 = axes.imshow(img,interpolation='nearest',aspect=1,extent=getExtent(myMap,slice),cmap=Cmap,norm=Norm,alpha=opacity)

    return im1

def drawLineout(axes,myMap,dir,ls,label=None):
    idx = getCurrIdx(myMap)
    if idx is not None:
        pc,vxls = getLineout(myMap.map,dir,idx)
        axes.set_facecolor('white')
        im1 = axes.plot(pc,vxls,ls,label=label)
        return im1
    else:
        return None



def getCurrView():
    global axXY,axYZ,axZX,axLineoutX,axLineoutY,axLineoutZ
    global xmin,xmax
    # print('XY',axXY[0].get_xlim(),axXY[0].get_ylim())
    # print('YZ',axYZ[0].get_xlim(),axYZ[0].get_ylim())
    # print('ZX',axZX[0].get_xlim(),axZX[0].get_ylim())
    xmin[0] = min(axXY[0].get_xlim())
    xmax[0] = max(axXY[0].get_xlim())
    xmin[1] = min(axYZ[0].get_xlim())
    xmax[1] = max(axYZ[0].get_xlim())
    xmin[2] = min(axZX[0].get_xlim())
    xmax[2] = max(axZX[0].get_xlim())
    # print('xmin',xmin)
    # print('xmax',xmax)



########################################################################
def draw():
    global lfirst,lLineout,lLogMap,lRestoreFullView
    global lContLevels,clevels_cmap,clevels_norm,clevels_bounds
    global lReverseX,lReverseY,lReverseZ
    global nn,hs,x0,Map
    global xmin,xmax,xnow
    global xmin_full,xmax_full
    global axXY,axYZ,axZX,axLineoutX,axLineoutY,axLineoutZ
    global maps,mapCT
    global idx
    global xc,yc,zc

    mapA = maps[0]  # a local view!

    nn = mapA.nn
    hs = mapA.hs
    x0 = mapA.xmin
    xc = np.arange(nn[0])*hs[0]+x0[0]+hs[0]/2
    yc = np.arange(nn[1])*hs[1]+x0[1]+hs[1]/2
    zc = np.arange(nn[2])*hs[2]+x0[2]+hs[2]/2

    if lfirst or lRestoreFullView:
        resetView()
        lRestoreFullView = False
    else:
        getCurrView()
    # print(xnow)

    if lfirst and np.any(ISO!=np.zeros(3)) and np.all(xmin<=ISO) and np.all(ISO<=xmax):
        xnow = ISO.copy()

    print('Current position=',xnow)
    
    plt.clf()   
    fig = plt.gcf()

    Cmap = mycmap
    if lGammaIndexMode:
        Cmap='seismic'
        Norm = GINormalize(vmin=0,vmax=mapGlobalMax)

    Norm = getColorNormalization()

    nmaps = len(maps)
    nrows = nmaps

    if lLineout:
        nrows += 1

    ncols = 3
    
    axXY = []
    axYZ = []
    axZX = []

    axXY.append(plt.subplot(nrows,ncols,1))
    axYZ.append(plt.subplot(nrows,ncols,2))
    axZX.append(plt.subplot(nrows,ncols,3))

    axXY_A = axXY[0]
    axYZ_A = axYZ[0]
    axZX_A = axZX[0]

    if args.legend is None: # set default to basenames without extension
        args.legend='basenames'

    if args.legend is not None:
        tokens = args.legend.split(',')
        # print(tokens)
        if len(tokens)>nmaps:
            print('Error: too many labels in legend for the given maps')
            exit(1)
        if len(tokens)==1 and 'basenames'.startswith(tokens[0]):
            # print('basename!',tokens[0]);exit(0)
            labels = [os.path.basename(m.path)[:-4] for m in maps]
        elif len(tokens)==1 and 'path'.startswith(tokens[0]):
            # print('basename!',tokens[0]);exit(0)
            labels = [m.path for m in maps]
        elif len(tokens)==1 and tokens[0]=='none':
            labels = [None]*nmaps
            args.legend = None
        else:
            labels = [None]*nmaps
            for i,lab in enumerate(tokens):
                labels[i] = lab.strip()

    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # cycle = ['r','b','g']

    for i in range(1,nmaps):
        axXY.append(plt.subplot(nrows,ncols,ncols*i+1,sharex=axXY_A,sharey=axXY_A))
        axYZ.append(plt.subplot(nrows,ncols,ncols*i+2,sharex=axYZ_A,sharey=axYZ_A))
        axZX.append(plt.subplot(nrows,ncols,ncols*i+3,sharex=axZX_A,sharey=axZX_A))

    if lLineout:
        axLineoutX = plt.subplot(nrows,ncols,nmaps*ncols+1,sharex=axXY_A)
        axLineoutY = plt.subplot(nrows,ncols,nmaps*ncols+2,sharex=axYZ_A)
        axLineoutZ = plt.subplot(nrows,ncols,nmaps*ncols+3,sharex=axZX_A)

    for i in range(nmaps):
        im1 = drawSlice(axXY[i],maps[i],'xy',Cmap,Norm,True)
        axXY[i].axvline(xnow[0], color='white', lw=1, alpha=0.5,dashes=(3,3))
        axXY[i].axhline(xnow[1], color='white', lw=1, alpha=0.5,dashes=(3,3))
        axXY[i].format_coord = Formatter(im1,1)
        if args.legend is not None:
            axXY[i].set_ylabel(labels[i],fontsize=16,color=cycle[i%len(cycle)])


        im2 = drawSlice(axYZ[i],maps[i],'yz',Cmap,Norm,True)
        axYZ[i].axvline(xnow[1], color='white', lw=1, alpha=0.5,dashes=(3,3))
        axYZ[i].axhline(xnow[2], color='white', lw=1, alpha=0.5,dashes=(3,3))
        axYZ[i].format_coord = Formatter(im2,2)
    
        im3 = drawSlice(axZX[i],maps[i],'zx',Cmap,Norm,True)
        axZX[i].axvline(xnow[2], color='white', lw=1, alpha=0.5,dashes=(3,3))
        axZX[i].axhline(xnow[0], color='white', lw=1, alpha=0.5,dashes=(3,3))
        axZX[i].format_coord = Formatter(im3,3)
        
        # colorbar 
        divider = make_axes_locatable(axZX[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        if lLogMap:
            cbar = plt.colorbar(im3, cax=cax,format=matplotlib.ticker.FuncFormatter(myLog_formatter))
        else:
            cbar = plt.colorbar(im3, cax=cax)



    axXY_A.set_xlim(xmin[0],xmax[0]) if not lReverseX else axXY_A.set_xlim(xmax[0],xmin[0]) 
    axXY_A.set_ylim(xmin[1],xmax[1]) if not lReverseY else axXY_A.set_ylim(xmax[1],xmin[1]) 
    axXY_A.set_title('XY slice at z=%.2f' % xnow[2])

    axYZ_A.set_xlim(xmin[1],xmax[1]) if not lReverseY else axYZ_A.set_xlim(xmax[1],xmin[1]) 
    axYZ_A.set_ylim(xmin[2],xmax[2]) if not lReverseZ else axYZ_A.set_ylim(xmax[2],xmin[2]) 
    axYZ_A.set_title('YZ slice at x=%.2f' % xnow[0])

    axZX_A.set_xlim(xmin[2],xmax[2]) if not lReverseZ else axZX_A.set_xlim(xmax[2],xmin[2]) 
    axZX_A.set_ylim(xmin[0],xmax[0]) if not lReverseX else axZX_A.set_ylim(xmax[0],xmin[0]) 
    axZX_A.set_title('ZX slice at y=%.2f' % xnow[1])

    if lLineout:
        lineStyle = ['r-','b-','g-','c-','m-','y-','k-']
        for i in range(nmaps):
            drawLineout(axLineoutX,maps[i],'x',lineStyle[i],label=labels[i])
            drawLineout(axLineoutY,maps[i],'y',lineStyle[i],label=labels[i])
            drawLineout(axLineoutZ,maps[i],'z',lineStyle[i],label=labels[i])
        axLineoutX.set_xlabel('x')
        axLineoutY.set_xlabel('y')
        axLineoutZ.set_xlabel('z')
        axLineoutX.format_coord = Formatter(None,11)
        axLineoutY.format_coord = Formatter(None,12)
        axLineoutZ.format_coord = Formatter(None,13)
        if args.legend is not None:
            axLineoutX.legend()
            axLineoutY.legend()
            axLineoutZ.legend()
        if args.gridLineout:
            axLineoutX.grid()
            axLineoutY.grid()
            axLineoutZ.grid()


    else:
        axLineoutX = None
        axLineoutY = None
        axLineoutZ = None

    for ax in [*axXY,*axYZ,*axZX]:
        ax.callbacks.connect('xlim_changed', on_xlims_change)
        ax.callbacks.connect('ylim_changed', on_ylims_change)

    for ax in [axLineoutX,axLineoutY,axLineoutZ]:
        if ax is not None:
            ax.callbacks.connect('xlim_changed', on_xlims_change)

    hsp=0.1
    plt.subplots_adjust(hspace=hsp)

    plt.show()

    if args.pdf:
        plt.savefig('plot.pdf')

    if args.png:
        plt.savefig('plot.png')

    if args.png or args.pdf:
        exit(0)

    return

########################################################################
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
########################################################################
def restoreEnv():
    global xmin,xmax,xnow
    try:
        fin = open('.env')
        lines = fin.readlines()
        xmin = np.array(map(float,lines[0].split()))
        xmax = np.array(map(float,lines[1].split()))
        xnow = np.array(map(float,lines[2].split()))
        fin.close()
        print('environment restored')
    except:
        pass

########################################################################
########################################################################
########################################################################
opacity = 1.0

if args.CT!=None:
    lCT=True
    loadCT()
    print('CT map loaded with shape',mapCT.map.shape)
    opacity = 0.5
    # exit(0)

if args.opacity is not None:
    if args.opacity<=0 or args.opacity>1:
        print('Opacity range error: must be in (0,1]')
        exit(1)
    opacity = args.opacity

if args.ISO!=None:
    ISO = np.array(args.ISO)
else:
    ISO = np.zeros(3)

# exit(0)

xnow=np.zeros(3)
xmin=np.zeros(3)
xmax=np.zeros(3)
xmin_full=np.zeros(3)
xmax_full=np.zeros(3)

Vmin=None
Vmax=None

lInteractive = not (args.pdf or args.png)
lfirst=True
lLineout= len(args.map) > 1 or args.lineout
lLogMap=False
lContLevels=False
lGammaIndexMode=args.GI
lReverseX=args.revX
lReverseY=args.revY
lReverseZ=args.revZ
nplot=-1

plt.rcParams["keymap.back"] = ['ctrl+left']
plt.rcParams["keymap.forward"] = ['ctrl+right']

from cycler import cycler
plt.rcParams['axes.prop_cycle'] = cycler(color='rbgcmyk')

fig = plt.figure(1,figsize=(12,10))

# fig.canvas.set_window_title(sys.argv[1])
fig.canvas.mpl_connect('button_press_event', on_click)
fig.canvas.mpl_connect('key_press_event', on_key)

plt.rcParams['image.cmap'] = 'jet'

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
mycmap = 'jet'


lredraw=True
lReloadMap=True
lRestoreFullView=True

nplot=0

while True:
    
    if lReloadMap:
        loadMaps()
        lReloadMap=False        
    
    plt.ion()
    
    if lredraw:
        draw()
    lredraw = True
        
    if lfirst:
        lfirst=False
    
    choice = input('Press return to redraw;  input "h" for help or "q" to exit: ')
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
                mycmap = cmapname
                plt.rcParams['image.cmap'] = cmapname
        except:
            pass

    if choice.lower() == "iso":
        xnow = ISO

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
        lReloadMap = False
        
    if choice in ["p","pdf"]:
        nplot +=1
        fname = 'plot%02d.pdf' % nplot
        plt.savefig(fname)
        print("Saved plot to file ",fname)
        lredraw = False
        lReloadMap = False

    if choice in ["P","PNG"]:
        nplot +=1
        fname = 'plot%02d.png' % nplot
        plt.savefig(fname)  
        lredraw = False
        lReloadMap = False
        print("Saved plot to file ",fname)

    if choice in ["j","lineout"]:
        lLineout = not lLineout

    if choice in ["y"]:
        lLogMap = not lLogMap

    m = re.match(r'\s*dco\s+(\S+)',choice,re.IGNORECASE) or re.match(r'\s*vcutoff\s+(\S+)',choice,re.IGNORECASE)
    if m:
        try:
            args.DCO = float(m.groups()[0])
        except:
            pass

    if choice.upper() in ["GI"]:
        lGammaIndexMode = not lGammaIndexMode
        if lGammaIndexMode:
            plt.rcParams['axes.facecolor'] = 'lightgray'
        else:
            plt.rcParams['axes.facecolor'] = 'white'

    m = re.match(r'\s*(vmin)\s+(\S+)',choice,re.IGNORECASE)
    if m:
        try:
            Vmin = float(m.groups()[1])
        except:
            Vmin = None

    m = re.match(r'\s*(vmax)\s+(\S+)',choice,re.IGNORECASE)
    if m:
        try:
            Vmax = float(m.groups()[1])
        except:
            Vmax = None

            
plt.ioff()
