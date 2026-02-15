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

def getQuaternion(axis, theta):
    q0=np.cos(theta/2)
    qq=axis*np.sin(theta/2)
    return np.array([q0, qq[0], qq[1], qq[2]])

def getRotMatrix(quaternion):
    w,x,y,z=quaternion
    R = np.array([
        [1 - 2*(y**2 + z**2),   2*(x*y - z*w),     2*(x*z + y*w)],
        [2*(x*y + z*w),         1 - 2*(x**2 + z**2), 2*(y*z - x*w)],
        [2*(x*z - y*w),         2*(y*z + x*w),     1 - 2*(x**2 + y**2)]
    ])
    return R

def rotate_ct_sitk(ct_image, R, fill_value=0):

    # Imposta la trasformazione affine
    transform = sitk.AffineTransform(3)
    transform.SetMatrix(R.flatten())

    # Centro di rotazione nel centro fisico del volume
    size = ct_image.GetSize()
    center_index = np.array(size) / 2
#    print(center_index)
    center_physical = ct_image.TransformContinuousIndexToPhysicalPoint(center_index)
    transform.SetCenter(center_physical)
#    print(center_physical)
    # Resampling
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(ct_image)  # mantiene dimensioni, spacing, origin
    resampler.SetTransform(transform)
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    resampler.SetDefaultPixelValue(fill_value)  # HU dell’aria

    rotated_image = resampler.Execute(ct_image)
#    filled_rot_image = binary_fill_holes(rotated_image).astype(np.uint8)

    # 4. Riempimento buchi con scipy.ndimage.binary_fill_holes
    array = sitk.GetArrayFromImage(rotated_image)  # (z, y, x)
    binary = (array > 0).astype(np.uint8)  # binarizza se servisse
    filled = binary_fill_holes(binary).astype(np.uint8)

    # 5. Ricostruisci immagine SimpleITK con metadata originali
    filled_image = sitk.GetImageFromArray(filled)
    filled_image.CopyInformation(rotated_image)

    return filled_image

def xyz_to_sitk(volume_xyz, spacing, origin, direction=(1,0,0,0,1,0,0,0,1)):
    image = sitk.GetImageFromArray(np.transpose(volume_xyz, (2,1,0)))  # da (x,y,z) → (z,y,x)
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    image.SetDirection(direction)
    return image

def calculate_isocenter(ptv_array, spacing, origin):
    coords = np.array(np.nonzero(ptv_array)).T
#    print(coords.shape)
    physical_coords = coords * spacing + origin
#    for i in range(physical_coords.shape[0]):
#    print(physical_coords)
    isocenter = physical_coords.mean(axis=0)
    return isocenter[::-1]

def main():
    parser = argparse.ArgumentParser(description='Rotates ROIs')
    parser.add_argument("map", help="path of the CT")
    parser.add_argument("normal", nargs=3, help="shooting direction")
    parser.add_argument("out", help="rotated map mhd name")
    args = parser.parse_args()
    CT = args.map
    sitkCT=sitk.ReadImage(CT)
    voxels = mhd_read(CT)
    [nn, hs, xmin, Map] = unpackVoxels(voxels)
#    print(nn)
#    print(type(hs), type(xmin))
    sitkCT=xyz_to_sitk(Map, 10*hs, 10*xmin)
    normal=np.array(args.normal, dtype=float)
    normal=normal/np.linalg.norm(normal)
    beam=np.array([0,0,-1])
    theta=-np.arccos(np.dot(normal,beam)) #a ruotare e' il sistema di riferimento quindi va invertito il segno

    if np.isclose(theta, 0):
        R = np.eye(3)
        axis = np.cross(normal, beam)
    else:
        axis = np.cross(normal, beam)
        axis = axis / np.linalg.norm(axis)
        q = getQuaternion(axis, theta)
        R = getRotMatrix(q)
        U, _, Vt = np.linalg.svd(R)
        R = U @ Vt  # forza ortogonalità numerica
    print("ROTATION ANGLE: "+str(-np.degrees(theta))+"° around AXIS:", axis) 
    print("")
    print("ROTATION MATRIX")
    print(R)
    print("")
    RotMap=rotate_ct_sitk(sitkCT, R, fill_value=0)
    sitk.WriteImage(RotMap,args.out)
    print("Rotated image saved as", args.out)
    print("")
    ptv_array_zyx=sitk.GetArrayFromImage(RotMap)
#    print(ptv_array.shape)
    iso_ptv=calculate_isocenter(ptv_array_zyx, hs[::-1], xmin[::-1])
    print("ISO_PTV[cm]:"+" "+str(iso_ptv[0])+" "+str(iso_ptv[1])+" "+str(iso_ptv[2]))
    isoptv_idx = np.floor((iso_ptv - xmin) / (hs + 1e-9)).astype(np.int64)
    print("ISO_PTV[voxel]:"+" "+str(isoptv_idx[0])+" "+str(isoptv_idx[1])+" "+str(isoptv_idx[2]))
    print("OFFSET: ", xmin)
    print("SPACING: ", hs)
    print("")
#    ptv_array=np.transpose(ptv_array_zyx,(2,1,0))
#    print(np.argwhere(ptv_array[ptv_array>1]))

if __name__ == "__main__":
    main()
