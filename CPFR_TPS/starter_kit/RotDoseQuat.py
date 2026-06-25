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

def rotate_dose_sitk(dose_image, Inv, fill_value=0.0):

    # Imposta la trasformazione affine
    transform = sitk.AffineTransform(3)
    # Usa la matrice Inv (es. R.T calcolata in precedenza)
    transform.SetMatrix(Inv.flatten())

    # Centro di rotazione nel centro fisico del volume
    size = dose_image.GetSize()
    center_index = np.array(size) / 2
    center_physical = dose_image.TransformContinuousIndexToPhysicalPoint(center_index)
    transform.SetCenter(center_physical)

    # Resampling
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(dose_image)  
    resampler.SetTransform(transform)
    
    # Interpolazione Lineare per evitare artefatti sui gradienti di dose
    resampler.SetInterpolator(sitk.sitkLinear) 
    
    # Il background deve essere 0 Gy (assenza di dose)
    resampler.SetDefaultPixelValue(fill_value)
    rotated_image = resampler.Execute(dose_image)
    return rotated_image
def xyz_to_sitk(volume_xyz, spacing, origin, direction=(1,0,0,0,1,0,0,0,1)):
    image = sitk.GetImageFromArray(np.transpose(volume_xyz, (2,1,0)))  # da (x,y,z) → (z,y,x)
    image.SetSpacing(spacing)
    image.SetOrigin(origin)
    image.SetDirection(direction)
    return image
if __name__ == "__main__":
    # 1. Controllo degli argomenti passati da terminale (ora passiamo la matrice invece dei vettori)
    if len(sys.argv) < 4:
        print("Uso: python3 RotDoseQuat.py <input_dose.mhd> <matrice_rotazione.txt> <output_dose.mhd>")
        sys.exit(1)

    input_file = sys.argv[1]
    matrix_file = sys.argv[2]
    output_file = sys.argv[3]

    print(f"Caricamento mappa di dose: {input_file}...")

    # 2. Lettura dell'immagine con SimpleITK
    try:
        dose_image = sitk.ReadImage(input_file)
    except Exception as e:
        print(f"Errore durante la lettura del file {input_file}: {e}")
        sys.exit(1)

    # 3. Caricamento della matrice di rotazione originale
    try:
        R_forward = np.loadtxt(matrix_file)
        print("Matrice di rotazione originale caricata con successo:\n", R_forward)
    except Exception as e:
        print(f"Errore durante la lettura del file matrice {matrix_file}: {e}")
        sys.exit(1)

    # 4. Calcolo della CONTRO-ROTAZIONE
    # Per una matrice di rotazione ortogonale, l'inversa è uguale alla trasposta.
    # Usiamo l'inversa per annullare esattamente la trasformazione fatta sulla CT.
    R_inv = np.linalg.inv(R_forward)

    # 5. Applichiamo la trasformazione
    print("Contro-rotazione della dose in corso...")
    rotated_image = rotate_dose_sitk(dose_image, R_inv)

    # 6. Salvataggio del risultato
    try:
        sitk.WriteImage(rotated_image, output_file)
        print(f"Salvataggio completato: {output_file}")
    except Exception as e:
        print(f"Errore durante il salvataggio: {e}")
        sys.exit(1)



