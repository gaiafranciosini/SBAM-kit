#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Jul 2020 A.Schiavi
########################################################################
########################################################################
import sys
import numpy as np
import SimpleITK as sitk

from FRED_mhd_io import readMap


def write_map_custom(dst_path, Map, hs, x0, direction=None):
    """Scrive la mappa MHD convertendo la matrice da ordinamento Fortran (x,y,z)

    all'ordinamento C (z,y,x) richiesto da SimpleITK.
    """
    # Trasposta da (nx, ny, nz) a (nz, ny, nx)
    array_zyx = np.transpose(Map, (2, 1, 0))
    img = sitk.GetImageFromArray(array_zyx)

    # Conversione dello spacing da cm a mm
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
    if len(sys.argv) < 5 + 1:
        print('usage: srcMapFile dx dy dz dstMapFile')
        print('shifts map by offset vector [dx,dy,dz]')
        sys.exit(1)

    src, dx, dy, dz, dst = sys.argv[1:6]

    # Lettura via FRED_mhd_io
    nn, hs, x0, Map, direction = readMap(src)

    offset = np.array(list(map(float, [dx, dy, dz])))

    print('Original map properties:')
    print('dims=', nn)
    print('spacing=', hs)
    print('offset=', x0)
    print('x-spatial extent= [', x0[0], ',', x0[0] + nn[0] * hs[0], '] => ', nn[0] * hs[0])
    print('y-spatial extent= [', x0[1], ',', x0[1] + nn[1] * hs[1], '] => ', nn[1] * hs[1])
    print('z-spatial extent= [', x0[2], ',', x0[2] + nn[2] * hs[2], '] => ', nn[2] * hs[2])
    print('datatype=', Map.dtype)
    print('')
    print('offset vector = ', offset)
    print('')

    # Applicazione dello spostamento
    x0 += offset

    # Scrittura della mappa shiftata
    write_map_custom(dst, Map, hs, x0, direction=direction)

    print('New map properties:')
    print('dims=', nn)
    print('spacing=', hs)
    print('offset=', x0)
    print('x-spatial extent= [', x0[0], ',', x0[0] + nn[0] * hs[0], '] => ', nn[0] * hs[0])
    print('y-spatial extent= [', x0[1], ',', x0[1] + nn[1] * hs[1], '] => ', nn[1] * hs[1])
    print('z-spatial extent= [', x0[2], ',', x0[2] + nn[2] * hs[2], '] => ', nn[2] * hs[2])


if __name__ == '__main__':
    main()
