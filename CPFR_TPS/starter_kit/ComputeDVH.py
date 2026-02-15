#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
import subprocess
from pathlib import Path

import numpy as np

try:
    import SimpleITK as sitk
except Exception as e:
    print("Errore: SimpleITK non disponibile. Installa con: pip install SimpleITK", file=sys.stderr)
    raise

BINNING_DVH_DEFAULT = 100

def read_mhd_to_numpy(path):
    """Legge un .mhd/.raw con SimpleITK e restituisce (array numpy, spacing_mm, origin_mm, direction)."""
    img = sitk.ReadImage(str(path))
    arr = sitk.GetArrayFromImage(img)  # z,y,x
    spacing = img.GetSpacing()         # (sx, sy, sz) in mm
    origin = img.GetOrigin()
    direction = img.GetDirection()
    # Converte in ordine x,y,z coerente per shape [z,y,x]
    # Per i calcoli volumetrici basta lo spacing, indipendente dall'ordine
    return arr, spacing, origin, direction

def roi_name_from_path(p: Path) -> str:
    """Estrae il nome ROI come in C++: basename senza estensione."""
    return p.stem

def compute_dvh_cumulative_leq(D_cGy, roi_values, spacing_mm, Dgoal, binning=BINNING_DVH_DEFAULT):
    """
    Calcola DVH cumulativo '≤ d' (stessa logica del C++):
      bins[j] = (j+0.5)*binSize
      hist[j] = volume frazionario (in %) con dose <= bins[j]
    D_cGy: np.ndarray delle dosi in cGy
    roi_values: np.ndarray della ROI (int o float); attivi dove >0
    """
    mask = roi_values > 0
    if not np.any(mask):
        return None  # ROI vuota

    sx, sy, sz = spacing_mm
    voxel_cm3 = (sx/10.0) * (sy/10.0) * (sz/10.0)

    # Pesi come nel C++: sommo roi_value * voxelVolume (int -> 1*voxel, float -> frazioni)
    weights = roi_values[mask].astype(np.float64) * voxel_cm3
    roi_volume_cm3 = weights.sum()
    n_vox = int(mask.sum())

    D_vals = D_cGy[mask].astype(np.float64)
    maxDose = float(D_vals.max())

    binSize = float(Dgoal) / float(binning)
    nbin = max(binning, int(maxDose / binSize) + 1)

    # Bins (centri) e edges per istogramma
    bins = (np.arange(nbin) + 0.5) * binSize
    edges = np.arange(nbin + 1) * binSize

    # Istogramma pesato
    hist_counts, _ = np.histogram(D_vals, bins=edges, weights=weights)  # per intervalli [edge_i, edge_{i+1})
    cdf_leq = np.cumsum(hist_counts)  # volume con dose <= edge_{i+1}
    # Riporta al centro-bin come nel C++ (che accumulava fino a imax inclusivo)
    hist_percent = (cdf_leq / roi_volume_cm3) * 100.0

    info = {
        "voxel_cm3": voxel_cm3,
        "roi_volume_cm3": roi_volume_cm3,
        "n_vox": n_vox,
        "maxDose_cGy": maxDose,
        "binSize": binSize,
        "nbin": nbin,
    }
    return bins, hist_percent, info

def main():
    parser = argparse.ArgumentParser(
        description="Calcolo DVH da dose + ROI (traduzione Python del tool C++).",
        add_help=True
    )
    parser.add_argument("-Dgoal", type=float, required=True, help="Dose di riferimento (stesso uso del C++).")
    parser.add_argument("-roi", nargs="+", required=True, help="Lista di file ROI (.mhd).")
    parser.add_argument("-dose", required=True, help="File dose (.mhd).")
    parser.add_argument("-type", choices=["int", "float"], required=True, help="Tipo ROI (int o float).")
    parser.add_argument("-fileLabel", required=True, help="Etichetta da appendere al nome file DVH.")
    parser.add_argument("-dir", required=True, help="Directory di output per i file DVH.")
    parser.add_argument("-plot", action="store_true", help="Se presente, lancia plotDVH.py come nel C++.")
    parser.add_argument("--binningDVH", type=int, default=BINNING_DVH_DEFAULT, help="Numero di bin base (default 100).")

    args = parser.parse_args()

    dose_path = Path(args.dose)
    roi_paths = [Path(p) for p in args.roi]
    out_dir = Path(args.dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Leggi dose
    print(f"Dose file: {dose_path}")
    dose_arr, spacing_mm, origin, direction = read_mhd_to_numpy(dose_path)
    # SimpleITK -> array [z,y,x]; manteniamo la stessa forma per ROI
    dose_arr = dose_arr.astype(np.float32)
    # Come nel C++: D[n] = dose[n] * 100 (Gy -> cGy usualmente)
    D_cGy = dose_arr * 100.0

    all_roi_names_for_plot = []

    for roi_path in roi_paths:
        # ROI
        roi_arr, spacing_mm_roi, _, _ = read_mhd_to_numpy(roi_path)

        # Check dimensioni/spacing coerenti
        if roi_arr.shape != D_cGy.shape:
            print(f"Errore: shape diversa tra dose {D_cGy.shape} e ROI {roi_arr.shape} per {roi_path}", file=sys.stderr)
            sys.exit(1)
        if spacing_mm_roi != spacing_mm:
            # Non blocco, ma avviso
            print(f"Warning: spacing differente tra dose {spacing_mm} e ROI {spacing_mm_roi} per {roi_path}", file=sys.stderr)

        if args.type == "int":
            roi_vals = roi_arr.astype(np.int16)
        else:
            roi_vals = roi_arr.astype(np.float32)

        res = compute_dvh_cumulative_leq(
            D_cGy=D_cGy,
            roi_values=roi_vals,
            spacing_mm=spacing_mm,
            Dgoal=args.Dgoal,
            binning=args.binningDVH
        )

        roi_name = roi_name_from_path(roi_path)
        all_roi_names_for_plot.append(roi_name)

        if res is None:
            print(f"ROI vuota (nessun voxel >0): {roi_name}")
            continue

        bins, hist_percent, info = res

        print(f"ROI: {roi_name} | vox attivi: {info['n_vox']} | Vol ROI: {info['roi_volume_cm3']:.3f} cm^3")
        print(f"  VoxelVolume: {info['voxel_cm3']:.6f} cm^3  | MaxDose: {info['maxDose_cGy']:.3f} cGy")
        print(f"  binSize: {info['binSize']:.4f}  | nbin: {info['nbin']}  | Dgoal: {args.Dgoal}")

        out_path = out_dir / f"{roi_name}{args.fileLabel}.txt"
        with open(out_path, "w") as f:
            for b, h in zip(bins, hist_percent):
                f.write(f"{b} {h}\n")
        print(f"  -> scritto: {out_path}")

    # Comando plot come in C++
    if args.plot and all_roi_names_for_plot:
        # Costruisce: python3 plotDVH.py -label1 <fileLabel> -dir1 <dir> -roi "<roi1 roi2 ...>" -Dgoal <Dgoal>
        roi_str = " ".join(all_roi_names_for_plot)
        cmd = [
            sys.executable, "plotDVH.py",
            "-label1", args.fileLabel,
            "-dir1", str(out_dir),
            "-roi", roi_str,
            "-Dgoal", str(args.Dgoal)
        ]
        print("Eseguo:", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"plotDVH.py ha restituito errore: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()
