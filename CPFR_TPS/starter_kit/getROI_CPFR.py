import argparse
import os
import sys
from math import *
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import pydicom
import matplotlib

matplotlib.use('TkAgg')
import pylab as plt
import SimpleITK as sitk

######################################################################
# CLI Arguments
######################################################################
parser = argparse.ArgumentParser(description='Export ROIs in dicom RTStruct using the phantom grid (solid ROIs only)')
parser.add_argument("-CT", help="CT mhd map", nargs=None, metavar='path', required=True)
parser.add_argument("-rtstruct", help="dicom file(s) containing RTStruct", nargs='+', required=True)
parser.add_argument("-l", help="list ROIs found in DICOM files", action="store_true", default=False)

parser.add_argument("-include", help="substrings of names to be included", nargs='+', metavar='str')
parser.add_argument("-exclude", help="substrings of names to be excluded", nargs='+', metavar='str')

parser.add_argument("-format", help="output file format", choices=['mhd', 'mha'], default='mhd')
parser.add_argument("-z", help="use compression", action="store_true", default=False)
parser.add_argument("-binary", help="binarize output mask (threshold 0.5) instead of keeping partial-volume values",
                    action="store_true", default=False)

args = parser.parse_args()

scriptDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, scriptDir)

######################################################################
# Global variables
######################################################################
dicoms = []
cts = []
rtstruct = []
rtplan = []
rtdose = []
ISO = np.zeros(3)
ROIs = []
nx = ny = nz = 0
hx = hy = hz = 0
x0 = y0 = z0 = 0
Ox = Oy = Oz = 0
Phantom = None
imgCT = None

def banner():
    print(' ')
    print('\tDICOM importer using python libraries for Fred project')
    print(' ')

######################################################################
# Dynamic Library Setup (C/C++ contours loader)
######################################################################
lib_path = os.path.join(scriptDir, "./libFillContours.so")
if not os.path.exists(lib_path):
    print(f"[ERROR] Dynamic library not found at: {lib_path}")
    sys.exit(1)

lib = ctypes.cdll.LoadLibrary(lib_path)
print("SCRIPT DIR:: ", scriptDir)
fillContourC = lib.fillContour_f32_f64
fillContourC.restype = None
fillContourC.argtypes = [
    ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
    ctypes.c_int, ctypes.c_int,
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ctypes.c_int, ctypes.c_int
]

######################################################################
# Helper Functions
######################################################################
def applyFilter(roiName):
    OK = False
    if args.include:
        for s in args.include:
            if s.lower() in roiName.lower():
                OK = True
    else:
        OK = True

    if args.exclude:
        for s in args.exclude:
            if s.lower() in roiName.lower():
                OK = False
    return OK

def findInterpretedType(rts, roiNum):
    if hasattr(rts, 'RTROIObservationsSequence'):
        for obs in rts.RTROIObservationsSequence:
            if obs.ReferencedROINumber == roiNum:
                return getattr(obs, 'RTROIInterpretedType', 'UNKNOWN')
    return 'UNKNOWN'

######################################################################
# ROI Processing Logic
######################################################################
def buildROIDriver():
    print('\nLooking for ROIs:')
    if len(rtstruct) == 0:
        print('none found!')
        return

    volumes = []
    for rts in rtstruct:
        for roi in rts.StructureSetROISequence:
            if not applyFilter(roi.ROIName):
                continue
            type_str = findInterpretedType(rts, roi.ROINumber)
            print(f'    {roi.ROINumber:4d}: {roi.ROIName:30s} {type_str:20s}')
            volumes.append((rts, roi))

    if args.l:
        return

    print('\n')
    for rts, roi in volumes:
        buildROI(rts, roi.ROINumber)

def buildROI(rts, roiNum):
    global Phantom, imgCT, nx, ny, nz, hx, hy, hz, Ox, Oy, Oz, x0, y0, z0

    roi = None
    for obj in rts.StructureSetROISequence:
        if obj.ROINumber == roiNum:
            roi = obj
            break

    if roi is None:
        return

    seq = None
    for contseq in rts.ROIContourSequence:
        if contseq.ReferencedROINumber == roiNum:
            seq = contseq
            break

    if seq is not None and 'ContourSequence' in seq:
        print('ROI num %d, name %s, num of contours %d' % (roiNum, roi.ROIName, len(seq.ContourSequence)))
    else:
        return

    ROI = np.zeros_like(Phantom, dtype=np.float32)

    # 1. Raccogli e ordina le quote Z dei contorni
    z_coords = sorted(list(set(cont.ContourData[2] for cont in seq.ContourSequence)))
    
    # 2. Calcola lo spessore Z medio associato a ciascun contorno DICOM
    if len(z_coords) > 1:
        z_spacings = np.diff(z_coords)
        default_dz = float(np.median(z_spacings))
    else:
        default_dz = hz

    for cont in seq.ContourSequence:
        print('.', end='', flush=True)

        mask = np.zeros((ny, nx), dtype=np.float32)

        # Mappatura coordinate X, Y
        Xcont = (np.array(cont.ContourData[0::3], dtype=np.float64) - x0) / hx
        Ycont = (np.array(cont.ContourData[1::3], dtype=np.float64) - y0) / hy

        lev = 5
        fillContourC(mask, ny, nx, Xcont, Ycont, len(Xcont), lev)
        mask = np.abs(mask)

        # 3. Determina l'intervallo Z da coprire (estrusione 3D per riempire i vuoti)
        zc = cont.ContourData[2]
        
        # Intervallo Z coperto da questo contorno
        z_min = zc - default_dz / 2.0
        z_max = zc + default_dz / 2.0

        # Mappa l'intervallo reale sulle fette della TAC
        iz_min = max(0, int(round((z_min - Oz) / hz)))
        iz_max = min(nz, int(round((z_max - Oz) / hz)) + 1)

        # Se l'intervallo calcolato e' vuoto, assegna almeno la fetta piu' vicina
        if iz_min >= iz_max:
            iz_closest = int(round((zc - Oz) / hz))
            if 0 <= iz_closest < nz:
                ROI[iz_closest] = np.maximum(ROI[iz_closest], mask)
        else:
            for iz in range(iz_min, iz_max):
                ROI[iz] = np.maximum(ROI[iz], mask)

    print('')

    # Binarizzazione opzionale
    if args.binary:
        ROI = (ROI >= 0.5).astype(np.float32)

    nz_voxels = np.count_nonzero(ROI)
    print("Numero di voxel diversi da zero nella matrice ROI:", nz_voxels)
    if nz_voxels == 0:
        print(f"[WARNING] ROI '{roi.ROIName}' resulted completely empty!")

    name = '_'.join(roi.ROIName.split())
    img_out = sitk.GetImageFromArray(ROI)

    img_out.SetOrigin(imgCT.GetOrigin())
    img_out.SetSpacing(imgCT.GetSpacing())
    img_out.SetDirection(imgCT.GetDirection())

    sitk.WriteImage(img_out, name + '.' + args.format, useCompression=args.z)

    roiVolume = np.sum(ROI) * hx * hy * hz
    ROIs.append((name, roiVolume))

######################################################################
# Main Execution
######################################################################
if __name__ == '__main__':
    banner()

    imgCT = sitk.ReadImage(args.CT)
    Phantom = sitk.GetArrayFromImage(imgCT)
    nx, ny, nz = imgCT.GetSize()
    hx, hy, hz = imgCT.GetSpacing()
    Ox, Oy, Oz = imgCT.GetOrigin()
    x0 = Ox - hx / 2
    y0 = Oy - hy / 2
    z0 = Oz - hz / 2

    print('# CT info:')
    print('# dims    =', nx, ny, nz)
    print('# spacing =', hx, hy, hz)
    print('# offset  =', Ox, Oy, Oz)

    for fname in args.rtstruct:
        try:
            ds = pydicom.dcmread(fname)
            if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.481.3':
                rtstruct.append(ds)
            else:
                print('error: file is not an RTStruct', fname)
                sys.exit(1)
        except Exception as e:
            print(f'error reading DICOM file {fname}: {e}')
            sys.exit(1)

    buildROIDriver()

    if not args.l:
        print('------------------------------')
        print('%-20s%-10s' % ('Name', 'Volume (cm3)'))
        for r in ROIs:
            print('%-20s%-10f' % (r[0], r[1] / 1000.))
