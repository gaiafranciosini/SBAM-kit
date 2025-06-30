#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
   # Apr 2020 
########################################################################
########################################################################
#import sys
#import numpy as np
#import argparse
#import struct
#from mhd_io import *
#from math import *



#def usage():
#	print('usage: -A mapA -o mapC -mm #')
#	print('get a new mask whihch have a distance 1/2/3/4 cm from PTV')
#
#
#
#
#if len(sys.argv)==1:
#	usage()
#	exit(1)
#sys.argv.pop(0)
#
#
#fnameC = 'newmap.mhd'
#margin = mm
#
#
#if '-o' in sys.argv:
#	try:
#		fnameC = sys.argv[sys.argv.index('-o')+1]
#		del sys.argv[sys.argv.index('-o')+1]
#		del sys.argv[sys.argv.index('-o')]
#	except:
#		print('error in parsing output filename')
#		exit(1)
#
#operator = 'none'
#
#if '-A' in sys.argv:
#	try:
#		fnameA = sys.argv[sys.argv.index('-A')+1]
#		del sys.argv[sys.argv.index('-A')+1]
#		del sys.argv[sys.argv.index('-A')]
#	except:
#		print('error in parsing filename for map A')
#		exit(1)
#
#else:
#	print('error: map A not defined')
#	usage()
#	exit(1)
#
#if '-mm' in sys.argv:
#	try:
#		margin = sys.argv[sys.argv.index('-mm')+1]
#		del sys.argv[sys.argv.index('-mm')+1]
#		del sys.argv[sys.argv.index('-mm')]
#	except:
#		print('error in parsing margin for the map')
#		exit(1)
#	else:
#		print('error: margin is not defined')
#		usage()
#		exit(1)
#
#
#if '-v' in sys.argv:
#	verbose = True
#	del sys.argv[sys.argv.index('-v')]
#else:
#	verbose = False


import argparse
parser = argparse.ArgumentParser(description='get a new mask which have a distance 1/2/3/4 cm from PTV')
parser.add_argument("map", help="path of the reference map")
parser.add_argument("-o", help="output file path",required=True,nargs=None,metavar='path')
parser.add_argument("-mm", help="margin of the mask",nargs=None,type=float,metavar='#')

args = parser.parse_args()

import sys
import numpy as np
import argparse
import struct
from mhd_io import *
from math import *

fnameA = args.map

fnameC = args.o 

mm=args.mm







#diamo in input il file mhd da studiare, in questo caso è il PTV
voxelsA = mhd_read(fnameA)
[nnA,hsA,x0A,MapA] = unpackVoxels(voxelsA)
#print(MapA)

#vado a trovare gli indici dove la mappa A e' =1 
I=np.where(MapA!=0)
PTV = MapA
nvxlPTV=len(I[0])
#print(nvxlPTV)
#costruisco l'opportuna maschera
mask=np.zeros_like(MapA)
Imask = np.where(mask==1)



#qui vado a definire il valore minimo e massimo del PTV

ixminptv=np.min(I[0])
#print(ixminptv)
ixmaxptv=np.max(I[0])
#print(ixmaxptv)
iyminptv=np.min(I[1])
#print(iyminptv)
iymaxptv=np.max(I[1])
#print(iymaxptv)
izminptv=np.min(I[2])
#print(izminptv)
izmaxptv=np.max(I[2])
#print(izmaxptv)



# exit(1)

#qui costruisco il margine della maschera che nel seguente caso volevamo di 4 cm 

 # maximum margin (cm)

ixminmm = int(max(0,ixminptv-mm/hsA[0]))
ixmaxmm = int(min(nnA[0]-1,ixmaxptv+mm/hsA[0]))
iyminmm = int(max(0,iyminptv-mm/hsA[1]))
iymaxmm = int(min(nnA[1]-1,iymaxptv+mm/hsA[1]))
izminmm = int(max(0,izminptv-mm/hsA[2]))
izmaxmm = int(min(nnA[2]-1,izmaxptv+mm/hsA[2]))


mask[ixminmm:ixmaxmm,iyminmm:iymaxmm,izminmm:izmaxmm]=1


mask[I]=2

Imask = np.where(mask==1)
#print(Imask)
Iptv = I

nvxlmask=len(Imask[0])




#vado a definire il numero di voxel sul bordo del PTV
onBoundary =np.zeros([nvxlPTV])

for idx in range(nvxlPTV):
	i = Iptv[0][idx]
	j = Iptv[1][idx]
	k = Iptv[2][idx]
	somma_primi_vicini=PTV[i-1,j,k]+PTV[i+1,j,k]+PTV[i,j-1,k]+PTV[i,j+1,k]+PTV[i,j,k-1]+PTV[i,j,k+1]
	if somma_primi_vicini!=6:
		onBoundary[idx]=1
print('num voxels inside PTV',nvxlPTV)
print('num voxels on PTV boundary',np.sum(onBoundary))
print('num voxels inside mask',nvxlmask)

Iboundary = (Iptv[0][np.where(onBoundary==1)],Iptv[1][np.where(onBoundary==1)],Iptv[2][np.where(onBoundary==1)])


# mask2[:]=PTV[:]
mask2 = np.zeros_like(PTV,dtype=np.float32)

nvxlBoundary=len(Iboundary[0])


I_zeros=np.where(MapA==0)

mask2[I]=0
mask2[I_zeros]=0


#loop sui voxel del PTV o sul numero di voxel sul bordo del PTV 
#i,j,k=indici del voxel del PTV
nstep = nvxlmask/100

print('looping...')
for i in range(nvxlmask):
	# mask2[Imask[0][i],Imask[1][i],Imask[2][i]]=i
# 	# print('.')
# 	distance_min=10000000
	if i>=nstep:
		print(int(100.*i/nvxlmask),'%')
		nstep += nvxlmask/100
	 #dsq = ((Imask[0][i]-Iboundary[0])*(Imask[0][i]-Iboundary[0]))*hsA[0]*hsA[0]
	dsqx = (Imask[0][i]-Iboundary[0])*hsA[0]
	dsqy = (Imask[1][i]-Iboundary[1])*hsA[1]
	dsqz = (Imask[2][i]-Iboundary[2])*hsA[2]
	dsq = dsqx*dsqx + dsqy*dsqy + dsqz*dsqz
# 	# for j in range(nvxlBoundary):
# 	# 	deltax2=((Imask[0][i]-Iboundary[0][j])*(Imask[0][i]-Iboundary[0][j]))*hsA[0]*hsA[0]
# 	# 	deltay2=((Imask[1][i]-Iboundary[1][j])*(Imask[1][i]-Iboundary[1][j]))*hsA[1]*hsA[1]
# 	# 	deltaz2=((Imask[2][i]-Iboundary[2][j])*(Imask[2][i]-Iboundary[2][j]))*hsA[2]*hsA[2]
# 	# 	distance_now=sqrt(deltax2+deltay2+deltaz2)
# 	# 	if (distance_min>distance_now): 
# 	# 		distance_min=distance_now
# 	# mask2[Imask[0][i],Imask[1][i],Imask[2][i]]=distance_min
	mask2[Imask[0][i],Imask[1][i],Imask[2][i]]=sqrt(np.min(dsq))


	if mask2[Imask[0][i],Imask[1][i],Imask[2][i]]<=mm:
		mask2[Imask[0][i],Imask[1][i],Imask[2][i]]=1
	else:
		mask2[Imask[0][i],Imask[1][i],Imask[2][i]]=0
		














mhd_write(fnameC,packVoxels(nnA,hsA,x0A,mask2))

				



		








		
			






		



			






				
			
				

				



						

	
		
			






















			




				












			




				




			
				

			






















			









	








			
	











			







































  




 




























	

