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

def rotate_ct_sitk(ct_image, R, fill_value=-1000):

    # Imposta la trasformazione affine
    transform = sitk.AffineTransform(3)
    transform.SetMatrix(R.flatten())

    # Centro di rotazione nel centro fisico del volume
    size = ct_image.GetSize()
    center_index = np.array(size) / 2
    print(center_index)
    center_physical = ct_image.TransformContinuousIndexToPhysicalPoint(center_index)
    transform.SetCenter(center_physical)
    print(center_physical)
    # Resampling
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(ct_image)  # mantiene dimensioni, spacing, origin
    resampler.SetTransform(transform)
    resampler.SetInterpolator(sitk.sitkBSpline)
    resampler.SetDefaultPixelValue(fill_value)  # HU dell’aria

    rotated_image = resampler.Execute(ct_image)
    return rotated_image

def xyz_to_sitk(volume_xyz, spacing, origin, direction=(1,0,0,0,1,0,0,0,1)):
    image = sitk.GetImageFromArray(np.transpose(volume_xyz, (2,1,0)))  # da (x,y,z) → (z,y,x)
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    image.SetDirection(direction)
    return image

def main():
    parser = argparse.ArgumentParser(description='Rotation a la Hamilton')
    parser.add_argument("map", help="path of the CT")
    parser.add_argument("normal", nargs=3, help="shooting direction")
    parser.add_argument("out", help="rotated map mhd name")
    args = parser.parse_args()
    CT = args.map
    sitkCT=sitk.ReadImage(CT)
    voxels = mhd_read(CT)
    [nn, hs, xmin, Map] = unpackVoxels(voxels)
    print(nn)
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
    print("rotation angle: "+str(-np.degrees(theta))+"° around axis:", axis) 
    print(R)
    RotMap=rotate_ct_sitk(sitkCT, R, fill_value=-1000)
    sitk.WriteImage(RotMap,args.out)
    print("Rotated image saved as", args.out)
if __name__ == "__main__":
    main()
