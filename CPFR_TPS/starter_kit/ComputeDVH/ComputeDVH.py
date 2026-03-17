#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import subprocess
from pathlib import Path

import numpy as np

try:
    import SimpleITK as sitk
except Exception as e:
    print("Error: SimpleITK not available. Install it with: pip install SimpleITK", file=sys.stderr)
    raise


BINNING_DVH_DEFAULT = 2000


def read_mhd_to_numpy(path: Path):
    """
    Read a .mhd/.raw file using SimpleITK and return:
      - arr: numpy array in [z, y, x] order
      - spacing_mm: tuple (sx, sy, sz) in mm
      - origin, direction: image metadata
    """
    img = sitk.ReadImage(str(path))
    arr = sitk.GetArrayFromImage(img)      # shape: (z, y, x)
    spacing_mm = img.GetSpacing()          # (sx, sy, sz) in mm
    origin = img.GetOrigin()
    direction = img.GetDirection()
    return arr, spacing_mm, origin, direction


def roi_name_from_path(p: Path) -> str:
    """Return the filename without extension."""
    return p.stem


def compute_dvh_cumulative_geq(
    D_cGy: np.ndarray,
    roi_values: np.ndarray,
    spacing_mm,
    Dgoal: float,
    binningDVH: int = BINNING_DVH_DEFAULT,
):
    """
    Compute cumulative DVH in the form:
        V(d) = percentage volume receiving dose >= d

    Steps:
      - select voxels where roi > 0
      - voxel weight = roi_value * voxelVolume
      - compute bin index as floor(D/binSize)
      - accumulate weights per bin
      - compute reverse cumulative sum (>= threshold)
    """
    mask = roi_values > 0
    if not np.any(mask):
        return None

    # Compute voxel volume in cm^3 (convert mm -> cm)
    sx, sy, sz = spacing_mm
    voxel_cm3 = (sx / 10.0) * (sy / 10.0) * (sz / 10.0)

    roi_vals = roi_values[mask].astype(np.float64)
    weights = roi_vals * voxel_cm3

    roi_volume_cm3 = float(weights.sum())
    n_vox = int(mask.sum())

    D_vals = D_cGy[mask].astype(np.float64)
    maxDose = float(D_vals.max())

    binSize = float(Dgoal) / float(binningDVH)
    if binSize <= 0:
        raise ValueError("binSize <= 0: check Dgoal and binningDVH.")

    nbin = max(int(binningDVH), int(maxDose / binSize) + 1)

    # Bin centers
    bins = (np.arange(nbin, dtype=np.float64) + 0.5) * binSize

    # Compute bin index
    idx = np.floor(D_vals / binSize).astype(np.int64)
    idx = np.clip(idx, 0, nbin - 1)

    # Weighted count per bin
    counts = np.bincount(idx, weights=weights, minlength=nbin).astype(np.float64)

    # Reverse cumulative sum (dose >= threshold)
    hist = np.cumsum(counts[::-1])[::-1]

    hist_percent = (hist / roi_volume_cm3) * 100.0

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
        description="Compute DVH from dose and ROI volumes.",
        add_help=True,
    )
    parser.add_argument("-Dgoal", type=float, required=True, help="Reference dose value.")
    parser.add_argument("-roi", nargs="+", required=True, help="List of ROI files (.mhd).")
    parser.add_argument("-dose", required=True, help="Dose file (.mhd).")
    parser.add_argument("-type", choices=["int", "float"], required=True, help="ROI type (int or float).")
    parser.add_argument("-fileLabel", required=True, help="Label appended to output DVH filename.")
    parser.add_argument("-dir", required=True, help="Output directory for DVH files.")
    parser.add_argument("-plot", action="store_true", help="If set, execute plotDVH.py.")
    parser.add_argument(
        "--binningDVH",
        type=int,
        default=BINNING_DVH_DEFAULT,
        help=f"Number of bins (default {BINNING_DVH_DEFAULT}).",
    )

    args = parser.parse_args()

    dose_path = Path(args.dose)
    roi_paths = [Path(p) for p in args.roi]
    out_dir = Path(args.dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read dose
    print(f"Dose file: {dose_path}")
    dose_arr, spacing_mm, origin, direction = read_mhd_to_numpy(dose_path)
    dose_arr = dose_arr.astype(np.float32)

    # Convert dose to cGy
    D_cGy = dose_arr * 100.0

    all_roi_names_for_plot = []

    for roi_path in roi_paths:
        print(f"roitype: {args.type}")

        roi_arr, spacing_mm_roi, _, _ = read_mhd_to_numpy(roi_path)

        # Check shape consistency
        if roi_arr.shape != D_cGy.shape:
            print(
                f"Error: mismatched shape between dose {D_cGy.shape} and ROI {roi_arr.shape} for {roi_path}",
                file=sys.stderr,
            )
            sys.exit(1)

        # Warn if spacing differs
        if spacing_mm_roi != spacing_mm:
            print(
                f"Warning: spacing differs between dose {spacing_mm} and ROI {spacing_mm_roi} for {roi_path}",
                file=sys.stderr,
            )

        if args.type == "int":
            roi_vals = roi_arr.astype(np.int16)
        else:
            roi_vals = roi_arr.astype(np.float32)

        roi_name = roi_name_from_path(roi_path)
        all_roi_names_for_plot.append(roi_name)

        res = compute_dvh_cumulative_geq(
            D_cGy=D_cGy,
            roi_values=roi_vals,
            spacing_mm=spacing_mm,
            Dgoal=args.Dgoal,
            binningDVH=args.binningDVH,
        )

        if res is None:
            print(f"Empty ROI (no voxel > 0): {roi_name}")
            continue

        bins, hist_percent, info = res

        print(f"VoxelVolume: {info['voxel_cm3']:.6f}")
        print(f"ROI: {roi_name} Non-zero voxels: {info['n_vox']} VolumeRoi: {info['roi_volume_cm3']}")
        print(f"MaxDose: {info['maxDose_cGy']:.6f}")
        print(f"{info['maxDose_cGy']} {info['binSize']} {info['nbin']} {args.Dgoal}")

        out_path = out_dir / f"{roi_name}{args.fileLabel}.txt"
        with open(out_path, "w") as f:
            for b, h in zip(bins, hist_percent):
                f.write(f"{b:.6f} {h:.6f}\n")

    print(" ".join(all_roi_names_for_plot))

    if args.plot and all_roi_names_for_plot:
        cmd = [
            sys.executable,
            "plotDVH.py",
            "-label1",
            args.fileLabel,
            "-dir1",
            str(out_dir),
            "-roi",
            *all_roi_names_for_plot,
            "-Dgoal",
            str(args.Dgoal),
        ]
        try:
            subprocess.run(cmd, check=True)
        except FileNotFoundError:
            print("Error: plotDVH.py not found in current directory.", file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"plotDVH.py returned an error: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
