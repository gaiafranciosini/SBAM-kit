import os, sys, shutil, time,re
import SimpleITK as sitk
import matplotlib.pyplot as plt
from scipy.ndimage import binary_fill_holes
from skimage import measure, morphology
from math import *
import numpy as np
import struct
from mhd_io import *
import argparse
import textwrap

def main():
    parser = argparse.ArgumentParser(description='Get your own CT')
    parser.add_argument("map", help="path of the CT")
    args = parser.parse_args()
    CT = args.map
    voxels = mhd_read(CT)
    [nn, hs, xmin, Map] = unpackVoxels(voxels)
    print(type(Map))
    ct_image = sitk.ReadImage(CT)
    ct_values = sitk.GetArrayFromImage(ct_image)
    ct_spacing = ct_image.GetSpacing()
    ct_origin = ct_image.GetOrigin()
    rep_index = np.unravel_index(np.argmax(ct_values), ct_values.shape)
    print('repery indexes: ', rep_index)
    x = ct_origin[0] + rep_index[2] * ct_spacing[0]  
    y = ct_origin[1] + rep_index[1] * ct_spacing[1]  
    z = ct_origin[2] + rep_index[0] * ct_spacing[2]
    rep_coords = (x*0.1,y*0.1,z*0.1)
    print('repery coordinates [cm]: ', rep_coords)


if __name__ == "__main__":
    main()

