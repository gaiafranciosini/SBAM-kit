#!/usr/bin/env python3
########################################################################
########################################################################
# Nov 2025 A.Burattini
########################################################################
########################################################################

import argparse
import numpy as np
from mhd_io import *

parser = argparse.ArgumentParser(
    description='Evaluates isodoses and saves them on a mhd map'
)
parser.add_argument("-dose", help="Dose map", required=True)
parser.add_argument("-out",  help="output map", required=True)
parser.add_argument("-preD", help="prescription dose", type=float, default=2000, required=True)
parser.add_argument("-iso", help="isodoses values as percentage of prescription dose", nargs="+", type=float, default=[10, 50, 80, 95, 100, 105, 110], required=True)

args = parser.parse_args()

dose = args.dose
out = args.out
preD=args.preD
iso = args.iso


#if (ct_shape[0]<dose_shape[0] or ct_shape[1]<dose_shape[1] or ct_shape[2]<dose_shape[2]):
#  sys.exit("ERROR: dose size smaller than CT size!")


[nn, hs, x0, Map]=unpackVoxels(mhd_read(dose))

iso_sort=np.sort(np.array(iso))/100
print(iso_sort)
isoMap=Map
isoMap[Map < iso_sort[0]*preD]=0

if (len(iso_sort)>1):
  for d in range(1, len(iso_sort)):
    isoMap[(Map >= iso_sort[d-1]*preD) & (Map < iso_sort[d]*preD)]=iso_sort[d-1]*preD
isoMap[Map >= iso_sort[-1]*preD]=iso_sort[-1]*preD

mhd_write(out,packVoxels(np.array(isoMap.shape),hs, x0, isoMap))







