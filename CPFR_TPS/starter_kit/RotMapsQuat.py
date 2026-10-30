import argparse
import numpy as np
import SimpleITK as sitk

from FRED_mhd_io import readMap


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


def rotate_ct_sitk(ct_image, R, fill_value=-1000):
    # Imposta la trasformazione affine
    transform = sitk.AffineTransform(3)
    transform.SetMatrix(R.flatten())

    # Centro di rotazione nel centro fisico del volume
    size = ct_image.GetSize()
    center_index = np.array(size) / 2.0
    center_physical = ct_image.TransformContinuousIndexToPhysicalPoint(center_index)
    transform.SetCenter(center_physical)

    # Resampling con interpolatore B-Spline
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(ct_image)  # Mantiene dimensioni, spacing, origin
    resampler.SetTransform(transform)
    resampler.SetInterpolator(sitk.sitkBSpline)
    resampler.SetDefaultPixelValue(fill_value)  # HU dell'aria per CT

    rotated_image = resampler.Execute(ct_image)
    return rotated_image


def xyz_to_sitk(volume_xyz, spacing, origin, direction=(1, 0, 0, 0, 1, 0, 0, 0, 1)):
    image = sitk.GetImageFromArray(np.transpose(volume_xyz, (2, 1, 0)))  # da (x,y,z) → (z,y,x)

    # Conversione in liste native di float per evitare TypeError con SimpleITK
    image.SetSpacing([float(s) for s in spacing])
    image.SetOrigin([float(o) for o in origin])

    dir_flat = np.array(direction).flatten().astype(float).tolist()
    image.SetDirection(dir_flat)

    return image


def main():
    parser = argparse.ArgumentParser(description='Rotation a la Hamilton')
    parser.add_argument("map", help="path of the CT")
    parser.add_argument("normal", nargs=3, help="shooting direction")
    parser.add_argument("out", help="rotated map mhd name")
    args = parser.parse_args()

    CT = args.map

    # Caricamento tramite readMap
    nn, hs, xmin, Map, direction = readMap(CT)

    # Conversione cm -> mm per SimpleITK
    sitkCT = xyz_to_sitk(Map, 10 * hs, 10 * xmin, direction)

    normal = np.array(args.normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    beam = np.array([0, 0, -1])
    theta = -np.arccos(np.dot(normal, beam))  # Segno invertito per ruotare il sistema di riferimento

    if np.isclose(theta, 0):
        R = np.eye(3)
        axis = np.cross(normal, beam)
    else:
        axis = np.cross(normal, beam)
        axis = axis / np.linalg.norm(axis)
        q = getQuaternion(axis, theta)
        R = getRotMatrix(q)
        U, _, Vt = np.linalg.svd(R)
        R = U @ Vt  # Forza ortogonalità numerica

    print("ROTATION ANGLE: " + str(-np.degrees(theta)) + "° around AXIS:", axis)
    print("\nROTATION MATRIX")
    print(R)
    print("")

    RotMap = rotate_ct_sitk(sitkCT, R, fill_value=-1000)
    sitk.WriteImage(RotMap, args.out)
    print("Rotated image saved as", args.out)


if __name__ == "__main__":
    main()
