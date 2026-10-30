import argparse
import os
import re
import sys
import time
import numpy as np
import SimpleITK as sitk
from scipy.ndimage import binary_fill_holes

from FRED_mhd_io import readMap, writeMap


def getQuaternion(axis, theta):
    q0 = np.cos(theta / 2)
    qq = axis * np.sin(theta / 2)
    return np.array([q0, qq[0], qq[1], qq[2]])


def getRotMatrix(quaternion):
    w, x, y, z = quaternion
    R = np.array([
        [1 - 2 * (y**2 + z**2), 2 * (x * y - z * w),     2 * (x * z + y * w)],
        [2 * (x * y + z * w),     1 - 2 * (x**2 + z**2), 2 * (y * z - x * w)],
        [2 * (x * z - y * w),     2 * (y * z + x * w),     1 - 2 * (x**2 + y**2)]
    ])
    return R


def rotate_ct_sitk(ct_image, R, fill_value=0):
    transform = sitk.AffineTransform(3)
    transform.SetMatrix(R.flatten())

    size = ct_image.GetSize()
    center_index = np.array(size) / 2.0
    center_physical = ct_image.TransformContinuousIndexToPhysicalPoint(center_index)
    transform.SetCenter(center_physical)

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(ct_image)
    resampler.SetTransform(transform)
    resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    resampler.SetDefaultPixelValue(fill_value)

    rotated_image = resampler.Execute(ct_image)

    array = sitk.GetArrayFromImage(rotated_image)
    binary = (array > 0).astype(np.uint8)
    filled = binary_fill_holes(binary).astype(np.uint8)

    filled_image = sitk.GetImageFromArray(filled)
    filled_image.CopyInformation(rotated_image)

    return filled_image


def xyz_to_sitk(volume_xyz, spacing, origin, direction=(1, 0, 0, 0, 1, 0, 0, 0, 1)):
    image = sitk.GetImageFromArray(np.transpose(volume_xyz, (2, 1, 0)))
    image.SetSpacing([float(s) for s in spacing])
    image.SetOrigin([float(o) for o in origin])
    
    dir_flat = np.array(direction).flatten().astype(float).tolist()
    image.SetDirection(dir_flat)
    
    return image


def calculate_isocenter(ptv_array, spacing, origin):
    coords = np.array(np.nonzero(ptv_array)).T
    physical_coords = coords * spacing + origin
    isocenter = physical_coords.mean(axis=0)
    return isocenter[::-1]


def main():
    parser = argparse.ArgumentParser(description='Rotates ROIs')
    parser.add_argument("map", help="path of the CT")
    parser.add_argument("normal", nargs=3, help="shooting direction")
    parser.add_argument("out", help="rotated map mhd name")
    args = parser.parse_args()

    CT = args.map

    nn, hs, xmin, Map, direction = readMap(CT)

    sitkCT = xyz_to_sitk(Map, 10 * hs, 10 * xmin, direction)

    normal = np.array(args.normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    beam = np.array([0, 0, -1])
    theta = -np.arccos(np.dot(normal, beam))

    if np.isclose(theta, 0):
        R = np.eye(3)
        axis = np.cross(normal, beam)
    else:
        axis = np.cross(normal, beam)
        axis = axis / np.linalg.norm(axis)
        q = getQuaternion(axis, theta)
        R = getRotMatrix(q)
        U, _, Vt = np.linalg.svd(R)
        R = U @ Vt

    print("ROTATION ANGLE: " + str(-np.degrees(theta)) + "° around AXIS:", axis)
    print("\nROTATION MATRIX")
    print(R)
    print("")

    RotMap = rotate_ct_sitk(sitkCT, R, fill_value=0)
    sitk.WriteImage(RotMap, args.out)
    print("Rotated image saved as", args.out, "\n")

    ptv_array_zyx = sitk.GetArrayFromImage(RotMap)
    iso_ptv = calculate_isocenter(ptv_array_zyx, hs[::-1], xmin[::-1])

    print("ISO_PTV[cm]: " + f"{iso_ptv[0]} {iso_ptv[1]} {iso_ptv[2]}")
    isoptv_idx = np.floor((iso_ptv - xmin) / (hs + 1e-9)).astype(np.int64)
    print("ISO_PTV[voxel]: " + f"{isoptv_idx[0]} {isoptv_idx[1]} {isoptv_idx[2]}")
    print("OFFSET: ", xmin)
    print("SPACING: ", hs)
    print("")


if __name__ == "__main__":
    main()
