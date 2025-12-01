#!/usr/bin/env python3
########################################################################
########################################################################
# Oct 2025 A.Burattini
########################################################################
########################################################################

import argparse
import numpy as np
from mhd_io import *

parser = argparse.ArgumentParser(
    description='Get the size of a dose grid keeping the same spacing of the CT'
)
parser.add_argument("-ct",   help="CT map",   required=True)
parser.add_argument("-dose", help="Dose map", required=True)
parser.add_argument("-out",  help="output map", required=True)

args = parser.parse_args()

ct   = args.ct
dose = args.dose
fout = args.out

[nn_ct, hs_ct, x0_ct, Map_ct]=unpackVoxels(mhd_read(ct))
[nn_dose, hs_dose, x0_dose, Map_dose]=unpackVoxels(mhd_read(dose))

ct_val=np.array(Map_ct)
dose_val=np.array(Map_dose)

print(np.shape(dose_val))

ct_shape=np.shape(Map_ct)
dose_shape=np.shape(Map_dose)

if (ct_shape[0]<dose_shape[0] or ct_shape[1]<dose_shape[1] or ct_shape[2]<dose_shape[2]):
  sys.exit("ERROR: dose size smaller than CT size!")

ct_coords=[hs_ct[0]*np.arange(0, ct_shape[0])+x0_ct[0],hs_ct[1]*np.arange(0, ct_shape[1])+x0_ct[1],hs_ct[2]*np.arange(0, ct_shape[2])+x0_ct[2]]

#print(ct_coords[0])
zero_idx_ct=np.array([np.argmin(np.abs(ct_coords[0])), np.argmin(np.abs(ct_coords[1])), np.argmin(np.abs(ct_coords[2]))])
print("ZERO: ", zero_idx_ct)

ct_center=np.array([ct_shape[0] // 2, ct_shape[1] // 2])
dose_center=np.array([dose_shape[0] // 2, dose_shape[1] // 2])

print(ct_center, dose_center)

new_dose=np.zeros((ct_shape[0], ct_shape[1], ct_shape[2]))
print(np.shape(new_dose))

if zero_idx_ct[0]-dose_center[0]<0 or zero_idx_ct[1]-dose_center[1]<0:
  sys.exit("ERROR: dose size and CT size not compatible!")

start_idx = [zero_idx_ct[0] - dose_center[0],zero_idx_ct[1] - dose_center[1],0]

end_idx = [zero_idx_ct[0] - dose_center[0] + dose_shape[0], zero_idx_ct[1] - dose_center[1] + dose_shape[1], dose_shape[2]]
new_dose[start_idx[0]:end_idx[0], start_idx[1]:end_idx[1], start_idx[2]:end_idx[2]]=dose_val

print(np.shape(new_dose))
print(ct_shape)


mhd_write(fout,packVoxels(np.array(new_dose.shape),hs_ct, x0_ct, new_dose))







