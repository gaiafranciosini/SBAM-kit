#!/usr/bin/env python3
"""
ptv_isocenter.py

Stampa SOLO il baricentro (centro di massa) di una ROI .mhd, in cm,
come tre numeri separati da spazio su una sola riga: "x y z".

Uso:
  python3 ptv_isocenter.py imgs/PTV_plan.mhd
  ISO=$(python3 ptv_isocenter.py imgs/PTV_plan.mhd)
  python3 map_viewer_isodose.py ... -ISO $ISO
"""

import sys
import argparse
import numpy as np
import SimpleITK as sitk


def main():
    ap = argparse.ArgumentParser(description="Baricentro (cm) di una ROI .mhd")
    ap.add_argument("roi", help="ROI .mhd (maschera binaria)")
    ap.add_argument("--threshold", type=float, default=0.0,
                    help="i voxel con valore > soglia fanno parte della ROI (default 0)")
    args = ap.parse_args()

    img = sitk.ReadImage(args.roi)
    arr = sitk.GetArrayFromImage(img)            # (Z, Y, X)

    mask = arr > args.threshold
    if not mask.any():
        sys.exit("ERRORE: la ROI e' vuota: nessun voxel sopra soglia.")

    sp = np.array(img.GetSpacing(), dtype=float)             # (X, Y, Z), mm
    org = np.array(img.GetOrigin(), dtype=float)
    D = np.array(img.GetDirection(), dtype=float).reshape(3, 3)

    idx = np.argwhere(mask)[:, ::-1].astype(float)           # (X, Y, Z)
    pts = org + (idx * sp) @ D.T                             # mm, coordinate fisiche
    centroid_cm = pts.mean(axis=0) * 0.1                     # mm -> cm

    print(" ".join(f"{x:.4f}" for x in centroid_cm))


if __name__ == "__main__":
    main()
