#!/usr/bin/env python3
import argparse
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def main():
    parser = argparse.ArgumentParser(
        description="Tripanello XY/YZ/ZX con slider e colorbar condivisa: tutti i piani usano la stessa dinamica [min,max] globale."
    )
    parser.add_argument("file", type=str, help="Percorso al file .mhd")
    parser.add_argument("--vmin", type=float, default=None, help="Forza il minimo display (default: min globale del volume)")
    parser.add_argument("--vmax", type=float, default=None, help="Forza il massimo display (default: max globale del volume)")
    args = parser.parse_args()

    # Leggi immagine: vol shape (Z, Y, X), spacing (sx, sy, sz)
    image = sitk.ReadImage(args.file)
    sx, sy, sz = image.GetSpacing()
    vol = sitk.GetArrayFromImage(image)
    Z, Y, X = vol.shape

    # Dinamica globale comune a tutte le viste
    gmin = float(np.min(vol)) if args.vmin is None else args.vmin
    gmax = float(np.max(vol)) if args.vmax is None else args.vmax
    if gmax <= gmin:
        raise ValueError("vmax deve essere > vmin.")

    # Slice iniziali
    iz, iy, ix = Z // 2, Y // 2, X // 2

    # Figura + griglia
    fig = plt.figure(figsize=(13, 6))
    gs = fig.add_gridspec(3, 3, height_ratios=[20, 1, 1.2], hspace=0.4)

    ax_xy = fig.add_subplot(gs[0, 0])
    ax_yz = fig.add_subplot(gs[0, 1])
    ax_zx = fig.add_subplot(gs[0, 2])

    # Imshow con stessa colormap e stessi limiti per TUTTE le viste
    common_kwargs = dict(cmap="gray", vmin=gmin, vmax=gmax, origin="lower")

    # XY: fissato Z
    im_xy = ax_xy.imshow(
        vol[iz, :, :],
        extent=[0, sx * X, 0, sy * Y],
        aspect="equal",
        **common_kwargs
    )
    ax_xy.set_title(f"XY (Z={iz+1}/{Z})")
    ax_xy.set_xlabel("X (mm)"); ax_xy.set_ylabel("Y (mm)")

    # YZ: fissato X (trasposta per avere assi (Z, Y) → (col, row))
    im_yz = ax_yz.imshow(
        vol[:, :, ix].T,
        extent=[0, sz * Z, 0, sy * Y],
        aspect="equal",
        **common_kwargs
    )
    ax_yz.set_title(f"YZ (X={ix+1}/{X})")
    ax_yz.set_xlabel("Z (mm)"); ax_yz.set_ylabel("Y (mm)")

    # ZX: fissato Y
    im_zx = ax_zx.imshow(
        vol[:, iy, :],
        extent=[0, sx * X, 0, sz * Z],
        aspect="equal",
        **common_kwargs
    )
    ax_zx.set_title(f"ZX (Y={iy+1}/{Y})")
    ax_zx.set_xlabel("X (mm)"); ax_zx.set_ylabel("Z (mm)")

    # Slider
    ax_s_xy = fig.add_subplot(gs[1, 0])
    ax_s_yz = fig.add_subplot(gs[1, 1])
    ax_s_zx = fig.add_subplot(gs[1, 2])

    s_xy = Slider(ax_s_xy, "Z-slice", 0, Z - 1, valinit=iz, valfmt="%d")
    s_yz = Slider(ax_s_yz, "X-slice", 0, X - 1, valinit=ix, valfmt="%d")
    s_zx = Slider(ax_s_zx, "Y-slice", 0, Y - 1, valinit=iy, valfmt="%d")

    def up_xy(_):
        z = int(s_xy.val)
        im_xy.set_data(vol[z, :, :])
        ax_xy.set_title(f"XY (Z={z+1}/{Z})")
        fig.canvas.draw_idle()

    def up_yz(_):
        x = int(s_yz.val)
        im_yz.set_data(vol[:, :, x].T)
        ax_yz.set_title(f"YZ (X={x+1}/{X})")
        fig.canvas.draw_idle()

    def up_zx(_):
        y = int(s_zx.val)
        im_zx.set_data(vol[:, y, :])
        ax_zx.set_title(f"ZX (Y={y+1}/{Y})")
        fig.canvas.draw_idle()

    s_xy.on_changed(up_xy)
    s_yz.on_changed(up_yz)
    s_zx.on_changed(up_zx)

    # Colorbar condivisa con etichette min/max globali
    cax = fig.add_subplot(gs[2, :])
    cbar = fig.colorbar(im_xy, cax=cax, orientation="horizontal")
    cbar.set_label("Intensità (valori voxel)")
    cax.text(0.0, -1.3, f"Min: {gmin:.3f}", ha="left", va="center", transform=cax.transAxes)
    cax.text(1.0, -1.3, f"Max: {gmax:.3f}", ha="right", va="center", transform=cax.transAxes)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
