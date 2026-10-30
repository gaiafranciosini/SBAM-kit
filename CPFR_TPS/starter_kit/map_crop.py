#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Feb 2023 aschiavi
########################################################################
########################################################################

import argparse
import sys
import numpy as np
import SimpleITK as sitk

from FRED_mhd_io import readMap


def write_map_custom(dst_path, Map, hs, x0, direction=None):
    """Scrive la mappa MHD convertendo la matrice da ordinamento Fortran (x,y,z)

    all'ordinamento C (z,y,x) richiesto da SimpleITK.
    """
    array_zyx = np.transpose(Map, (2, 1, 0))
    img = sitk.GetImageFromArray(array_zyx)

    # Conversione dello spacing da cm a mm per SimpleITK
    spacing_mm = [float(s * 10.0) for s in hs]
    img.SetSpacing(spacing_mm)

    # Conversione dell'origine da cm a mm (centro del primo voxel)
    if direction is not None:
        dir_flat = np.array(direction).flatten().astype(float).tolist()
        img.SetDirection(dir_flat)
        R = np.array(direction).reshape(3, 3)
        origin_mm = (x0 + 0.5 * R.dot(hs)) * 10.0
    else:
        origin_mm = (x0 + 0.5 * hs) * 10.0

    img.SetOrigin([float(o) for o in origin_mm])
    sitk.WriteImage(img, dst_path)


def main():
    parser = argparse.ArgumentParser(description='crop a subset of voxels of an mhd map using a mask')
    parser.add_argument("-i", help="input map", required=True, metavar='path')
    parser.add_argument("-o", help="output map", required=True, metavar='path')
    parser.add_argument("-ixi", help="start index x", required=True, type=int, metavar='#')
    parser.add_argument("-ixf", help="end index x", required=True, type=int, metavar='#')
    parser.add_argument("-iyi", help="start index y", required=True, type=int, metavar='#')
    parser.add_argument("-iyf", help="end index y", required=True, type=int, metavar='#')
    parser.add_argument("-izi", help="start index z", required=True, type=int, metavar='#')
    parser.add_argument("-izf", help="end index z", required=True, type=int, metavar='#')

    args = parser.parse_args()

    fin = args.i
    fout = args.o
    ixi, ixf = args.ixi, args.ixf
    iyi, iyf = args.iyi, args.iyf
    izi, izf = args.izi, args.izf

    # Lettura via readMap (restituisce hs e x0 in cm)
    nn, hs, x0, Map, direction = readMap(fin)

    newMap = Map[ixi:ixf, iyi:iyf, izi:izf].copy()

    idx_i = np.array([ixi, iyi, izi])
    idx_f = np.array([ixf, iyf, izf])

    x0new = x0 + hs * idx_i
    xFnew = x0 + hs * idx_f

    # Scrittura della mappa ritagliata
    write_map_custom(fout, newMap, hs, x0new, direction=direction)

    print('cropped map written to file:', fout)
    print("new_dims: " + str(ixf - ixi) + " " + str(iyf - iyi) + " " + str(izf - izi))
    print("new_offset: " + str(x0new[0]) + " " + str(x0new[1]) + " " + str(x0new[2]))
    print("new_end: " + str(xFnew[0]) + " " + str(xFnew[1]) + " " + str(xFnew[2]))
    print("new_L: " + str(xFnew[0] - x0new[0]) + " " + str(xFnew[1] - x0new[1]) + " " + str(xFnew[2] - x0new[2]))


if __name__ == '__main__':
    main()
