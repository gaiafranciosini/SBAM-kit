#!/usr/bin/env python3
"""
PTV coverage check against the 95% / 90% isodose (3D + slice-by-slice).
Writes a compact English report to a text file (and one summary line to stdout).

numpy arrays from SimpleITK are [Z, Y, X]:
  axial    = XY plane, scrolling along Z (numpy axis 0)
  coronal  = XZ plane, scrolling along Y (numpy axis 1)
  sagittal = YZ plane, scrolling along X (numpy axis 2)
"""
import argparse
import os
import sys

import numpy as np
import SimpleITK as sitk

# plane -> (numpy axis, axis label, plane name)
PLANES = {
    "axial": (0, "Z", "XY"),
    "coronal": (1, "Y", "XZ"),
    "sagittal": (2, "X", "YZ"),
}


def same_geometry(a, b):
    return (
        a.GetSize() == b.GetSize()
        and np.allclose(a.GetSpacing(), b.GetSpacing(), atol=1e-6)
        and np.allclose(a.GetOrigin(), b.GetOrigin(), atol=1e-6)
        and np.allclose(a.GetDirection(), b.GetDirection(), atol=1e-6)
    )


def overlap_stats(ptv, iso, thr):
    """3D coverage and Dice of the PTV against the region iso >= thr."""
    reached = np.isfinite(iso) & (iso >= thr)
    n_ptv = np.count_nonzero(ptv)
    n_iso = np.count_nonzero(reached)
    n_in = np.count_nonzero(ptv & reached)
    dice = 2.0 * n_in / (n_ptv + n_iso) if (n_ptv + n_iso) else np.nan
    return n_in, 100.0 * n_in / n_ptv, dice


def slice_rows(ptv, iso, img, plane, thrs):
    """One row per slice containing at least 1 PTV voxel:
    (index, position_mm, n_ptv_voxels, [coverage% for each threshold])."""
    axis = PLANES[plane][0]
    p3 = np.moveaxis(ptv, axis, 0)
    i3 = np.moveaxis(iso, axis, 0)
    rows = []
    for k in range(p3.shape[0]):
        p = p3[k]
        n = np.count_nonzero(p)
        if n == 0:
            continue
        s = i3[k]
        ok = np.isfinite(s)
        covs = [100.0 * np.count_nonzero(p & ok & (s >= t)) / n for t in thrs]
        idx = [0, 0, 0]              # SimpleITK index is (x, y, z)
        idx[2 - axis] = k
        pos = img.TransformIndexToPhysicalPoint(tuple(idx))[2 - axis]
        rows.append((k, pos, n, covs))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--isodose", required=True, help="isodose map (.mhd)")
    ap.add_argument("-p", "--ptv", required=True, help="PTV map (.mhd), same geometry")
    ap.add_argument("-o", "--output", required=True, help="output text report")
    ap.add_argument("--label", default=None, help="run label (default: isodose file name)")
    ap.add_argument("--summary", default=None,
                    help="tab-separated file to which ONE summary line is appended")
    ap.add_argument("--slice-threshold", type=float, default=100.0,
                    help="min %% of PTV covered in a slice to call it covered (default 100)")
    ap.add_argument("--table", action="store_true",
                    help="also write the full per-slice table for each plane")
    args = ap.parse_args()

    label = args.label or os.path.basename(args.isodose)

    iso_img = sitk.ReadImage(args.isodose)
    ptv_img = sitk.ReadImage(args.ptv)
    if not same_geometry(iso_img, ptv_img):
        sys.exit("ERROR: PTV and isodose map have different geometry.")

    iso = sitk.GetArrayFromImage(iso_img).astype(np.float64)
    ptv = sitk.GetArrayFromImage(ptv_img) > 0
    n_ptv = int(np.count_nonzero(ptv))
    if n_ptv == 0:
        sys.exit("ERROR: no PTV voxels (PTV > 0).")

    # isodose scale: 0-1 or 0-100
    if np.nanmax(iso) <= 1.5:
        thrs, scale = (0.95, 0.90), "0-1"
    else:
        thrs, scale = (95.0, 90.0), "0-100"

    n95, cov95, dice95 = overlap_stats(ptv, iso, thrs[0])
    n90, cov90, dice90 = overlap_stats(ptv, iso, thrs[1])

    L = []
    add = L.append
    add(f"PTV COVERAGE REPORT - {label}")
    add(f"Isodose map : {args.isodose}")
    add(f"PTV map     : {args.ptv}")
    add(f"Isodose scale {scale}; a PTV voxel is covered if isodose >= threshold")
    add(f"PTV voxels  : {n_ptv}")
    add("")
    add("3D VOLUMETRIC")
    add("coverage = covered PTV voxels / PTV voxels")
    add("Dice     = 2*|PTV and isodose region| / (|PTV| + |isodose region|)")
    add(f"{'iso':>5} {'covered':>9} {'below':>9} {'coverage%':>10} {'Dice':>8}")
    add(f"{'95%':>5} {n95:9d} {n_ptv - n95:9d} {cov95:10.2f} {dice95:8.4f}")
    add(f"{'90%':>5} {n90:9d} {n_ptv - n90:9d} {cov90:10.2f} {dice90:8.4f}")
    add("")
    if n95 == n_ptv:
        add("RESULT: the whole PTV is inside the 95% isodose.")
    elif n90 == n_ptv:
        add(f"RESULT: PTV NOT fully inside the 95% isodose ({cov95:.2f}% covered), "
            "but fully inside the 90% isodose.")
    else:
        add(f"RESULT: PTV NOT fully inside the 95% isodose ({cov95:.2f}%) "
            f"nor the 90% isodose ({cov90:.2f}%).")
    add("")

    add(f"SLICE-BY-SLICE (only slices with >=1 PTV voxel; slice OK if coverage >= "
        f"{args.slice_threshold:.1f}%)")
    add("depth = distance in mm from the first slice containing PTV "
        "(NOT radiological depth, NOT distance from the beam)")
    add(f"{'plane':<9} {'axis':<4} {'idx range':<11} {'mm range':<19} {'n':>3} "
        f"{'OK@95':>8} {'OK@90':>8}  worst cov95% (idx, depth)")

    per_plane = {}
    for plane, (_, ax, _) in PLANES.items():
        rows = slice_rows(ptv, iso, iso_img, plane, thrs)
        per_plane[plane] = rows
        n = len(rows)
        ok95 = sum(r[3][0] >= args.slice_threshold for r in rows)
        ok90 = sum(r[3][1] >= args.slice_threshold for r in rows)
        first_mm = rows[0][1]
        worst = min(rows, key=lambda r: r[3][0])
        add(f"{plane:<9} {ax:<4} {rows[0][0]:>4}-{rows[-1][0]:<6} "
            f"{first_mm:8.2f}->{rows[-1][1]:<9.2f} {n:3d} "
            f"{ok95:>4}/{n:<3} {ok90:>4}/{n:<3}  "
            f"{worst[3][0]:.2f} (idx {worst[0]}, {abs(worst[1] - first_mm):.2f} mm)")
        per_plane[plane] = (rows, ok95, ok90)

    if args.table:
        for plane, (_, ax, pname) in PLANES.items():
            rows = per_plane[plane][0]
            first_mm = rows[0][1]
            add("")
            add(f"{plane.upper()} ({pname}, along {ax})")
            add(f"{'idx':>5} {'mm':>9} {'depth_mm':>9} {'PTV_vox':>8} {'cov95%':>8} {'cov90%':>8}")
            for k, pos, n, covs in rows:
                add(f"{k:5d} {pos:9.2f} {abs(pos - first_mm):9.2f} {n:8d} "
                    f"{covs[0]:8.2f} {covs[1]:8.2f}")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")

    if args.summary:
        new = not os.path.exists(args.summary) or os.path.getsize(args.summary) == 0
        with open(args.summary, "a", encoding="utf-8") as f:
            if new:
                f.write("label\tPTV_voxels\tcov95_%\tcov90_%\tdice95\tdice90\t"
                        "axial_OK95\tcoronal_OK95\tsagittal_OK95\n")
            ok = [f"{per_plane[p][1]}/{len(per_plane[p][0])}" for p in PLANES]
            f.write(f"{label}\t{n_ptv}\t{cov95:.2f}\t{cov90:.2f}\t{dice95:.4f}\t"
                    f"{dice90:.4f}\t{ok[0]}\t{ok[1]}\t{ok[2]}\n")

    print(f"[{label}] cov95={cov95:.2f}% cov90={cov90:.2f}% "
          f"dice95={dice95:.3f} dice90={dice90:.3f} -> {args.output}")


if __name__ == "__main__":
    main()
