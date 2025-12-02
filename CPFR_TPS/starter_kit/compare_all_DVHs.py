#!/usr/bin/env python3
"""
Single-figure comparative DVH plotter
-------------------------------------

Scans directories named like:  sim{E}MeV_{slit}mm   (e.g., sim7MeV_5mm)
and, within each, reads DVH text files from:  DVH{E}MeV_{slit}mm
Each DVH .txt is expected to have two whitespace-separated columns:
    1) dose in cGy
    2) volume percentage (%)

Creates ONE single figure containing ALL DVH curves from all ROIs/energies/slits.

Usage:
    python3 plot_dvh_all.py --base-dir . --out "DVH_ALL.png" [--legend right|inside|none]

Color/style scheme:
- Same ROI -> same color
- Different slit (mm) -> different linestyle/marker
"""

import argparse
import os
import re
import sys
from glob import glob
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt

SIM_DIR_RE = re.compile(r"^sim(?P<E>\d+)MeV_(?P<slit>\d+)p$")

def read_two_column_txt(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read two-column DVH text file (dose cGy, volume %)."""
    try:
        arr = np.loadtxt(path, comments="#", ndmin=2)
        if arr.size == 0 or arr.shape[1] < 2:
            raise ValueError("Not enough columns")
        return arr[:, 0], arr[:, 1]
    except Exception:
        xs, ys = [], []
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                line = line.replace(",", " ")
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    xs.append(float(parts[0])); ys.append(float(parts[1]))
                except Exception:
                    continue
        if not xs:
            raise ValueError(f"No numeric pairs in '{path}'")
        return np.asarray(xs, float), np.asarray(ys, float)

def pretty_roi_name(stem: str) -> str:
    s = stem
    if s.endswith("_plan"):
        s = s[:-5]
    return s.replace("_", " ")

def scan_sim_dirs(base_dir: str):
    for name in sorted(os.listdir(base_dir)):
        full = os.path.join(base_dir, name)
        if not os.path.isdir(full):
            continue
        m = SIM_DIR_RE.match(name)
        if not m:
            continue
        yield full, int(m.group("E")), int(m.group("slit"))

def collect_all_curves(base_dir: str):
    """
    Return list of dicts with keys:
      roi, E, slit, dose (np.ndarray), vol (np.ndarray), label (str), path (str)
    """
    curves = []
    found_any = False
    for sim_path, E, slit in scan_sim_dirs(base_dir):
        dvh_dir = os.path.join(sim_path, f"DVH{E}MeV_{slit}p")
        if not os.path.isdir(dvh_dir):
            print(f"[WARN] DVH dir not found: {dvh_dir}", file=sys.stderr)
            continue
        found_any = True
        for txt in sorted(glob(os.path.join(dvh_dir, "*.txt"))):
            stem = os.path.splitext(os.path.basename(txt))[0]
            roi = pretty_roi_name(stem)
            try:
                x, y = read_two_column_txt(txt)
            except Exception as e:
                print(f"[WARN] Skip '{txt}': {e}", file=sys.stderr)
                continue
            label = f"{roi} | {E}MeV {slit}p"
            curves.append(dict(roi=roi, E=E, slit=slit, dose=x, vol=y, label=label, path=txt))

    if not found_any:
        print(f"[ERROR] No sim{{E}}MeV_{{slit}}p/DVH folders found under '{base_dir}'", file=sys.stderr)
    return sorted(curves, key=lambda c: (c["roi"].lower(), c["E"], c["slit"]))

def make_single_figure(curves, out_path: str, legend_mode: str = "right", dpi: int = 300):
    if not curves:
        print("[ERROR] No curves to plot.", file=sys.stderr)
        sys.exit(1)

    # === Palette: stesso ROI -> stesso colore ===
    rois = sorted({c["roi"] for c in curves}, key=str.lower)
    # usa tab20 se >10 ROI, altrimenti tab10
    cmap = plt.colormaps["tab20"] if len(rois) > 10 else plt.colormaps["tab10"]
    roi_color = {roi: cmap(i % cmap.N) for i, roi in enumerate(rois)}

    # === Stili: slit diverso -> stile diverso (linestyle+marker) ===
    slits = sorted({c["slit"] for c in curves})
    linestyles = ["-", "--", "-.", ":"]
    markers = ["", "o", "s", "^", "d", "v", "x", "P", "*"]
    def style_for_slit(s):
        i = slits.index(s)
        return linestyles[i % len(linestyles)], markers[i % len(markers)]

    plt.figure(figsize=(16, 8))
    for c in curves:
        order = np.argsort(c["dose"])
        xd = c["dose"][order]; yd = c["vol"][order]
        color = roi_color[c["roi"]]
        ls, mk = style_for_slit(c["slit"])
        # linewidth un po' più spesso per leggibilità; marker poco fitto
        plt.plot(xd, yd, label=c["label"], color=color, linestyle=ls, marker=mk, markevery=max(1, len(xd)//25), linewidth=1.8, alpha=0.95)

    plt.title("DVH - All ROIs, Energies, Slits")
    plt.xlabel("Dose [cGy]"); plt.ylabel("Volume [%]")
    plt.grid(True, alpha=0.3)

    if legend_mode == "right":
        n = len(curves)
        ncol = 1 if n <= 16 else 2 if n <= 32 else 3
        plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5),
                   ncol=ncol, fontsize="x-small", frameon=False)
    elif legend_mode == "inside":
        plt.legend(fontsize="x-small", frameon=False)
    elif legend_mode == "none":
        pass
    else:
        print(f"[WARN] Unknown legend mode '{legend_mode}', defaulting to 'right'", file=sys.stderr)
        n = len(curves)
        ncol = 1 if n <= 16 else 2 if n <= 32 else 3
        plt.legend(loc="center left", bbox_to_anchor=(1.0, 0.5),
                   ncol=ncol, fontsize="small", frameon=False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"[INFO] Saved: {out_path}")

def main():
    ap = argparse.ArgumentParser(description="Plot ALL DVHs into a single figure.")
    ap.add_argument("--base-dir", default=".",
                    help="Directory containing sim{E}MeV_{slit}p folders (default: .)")
    ap.add_argument("--out", default="DVH_ALL.png",
                    help="Output PNG path (default: DVH_ALL.png)")
    ap.add_argument("--legend", choices=["right", "inside", "none"], default="right",
                    help="Legend placement (default: right)")
    args = ap.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    out_path = os.path.abspath(args.out)

    print(f"[INFO] Scanning: {base_dir}")
    curves = collect_all_curves(base_dir)
    print(f"[INFO] Total curves: {len(curves)}")
    if curves:
        preview = ", ".join(c["label"] for c in curves[:5])
        print(f"[INFO] Example labels: {preview} ...")

    make_single_figure(curves, out_path, legend_mode=args.legend)

if __name__ == "__main__":
    main()
