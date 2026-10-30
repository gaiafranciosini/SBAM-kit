#!/usr/bin/env python3
"""
normalize_dose_to_isodose.py

Legge una mappa di dose in formato .mhd, calcola la dose massima e
normalizza la mappa rispetto ad essa per ottenere una mappa di isodose
(valori relativi, frazione 0-1 oppure percentuale), salvando il
risultato in un nuovo file .mhd.

Uso da riga di comando (bash):
    python normalize_dose_to_isodose.py -i dose.mhd -o isodose.mhd
    python normalize_dose_to_isodose.py -i dose.mhd -o isodose.mhd --percent
    python normalize_dose_to_isodose.py -i dose.mhd -o isodose.mhd --preview isodose.png -v
"""

import argparse
import numpy as np
import SimpleITK as sitk


########################################################################
def readMap(fname):
    """Legge una mappa .mhd e restituisce dimensioni, spacing, origine,
    array (ordinato x,y,z) e matrice di direzione."""
    img = sitk.ReadImage(fname)
    arr = sitk.GetArrayFromImage(img)
    nn = np.array(img.GetSize(), dtype=np.int32)
    hs = np.array(img.GetSpacing(), dtype=np.float64)
    x0 = np.array(img.GetOrigin(), dtype=np.float64)
    direction = np.array(img.GetDirection(), dtype=np.float64)
    # riordina l'array in modo che gli indici siano [x,y,z] (Fortran order),
    # senza toccare la memoria sottostante
    arr = arr.reshape(-1, order='A').reshape(nn, order='F')
    return nn, hs, x0, arr, direction


def writeMap(fname, hs, x0, arr, direction=(1, 0, 0, 0, 1, 0, 0, 0, 1)):
    """Salva un array [x,y,z] come mappa .mhd con spacing, origine e
    direzione date."""
    nn = arr.shape
    if arr.flags.f_contiguous:
        img = sitk.GetImageFromArray(
            arr.reshape(-1, order='A').reshape(nn[2], nn[1], nn[0], order='C')
        )
    else:
        img = sitk.GetImageFromArray(
            np.asfortranarray(arr).reshape(-1, order='A').reshape(nn[2], nn[1], nn[0], order='C')
        )
    img.SetSpacing(tuple(np.array(hs, dtype=np.double)))
    img.SetDirection(tuple(np.array(direction, dtype=np.double)))
    img.SetOrigin(tuple(np.array(x0, dtype=np.double)))
    sitk.WriteImage(img, fname)


########################################################################
def make_preview(isodose_map, out_png, percent=False):
    """Genera un PNG della slice centrale (piano XY) usando una colormap
    'TPS-like' a livelli fissi di isodose (ispirata a M. Schwarz)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    clevels_cmap = colors.ListedColormap([
        [102/255, 255/255, 255/255],
        [0/255, 255/255, 0/255],
        [197/255, 248/255, 190/255],
        [0/255, 0/255, 255/255],
        [122/255, 0/255, 255/255],
        [255/255, 125/255, 0/255],
        [255/255, 255/255, 0/255],
        [255/255, 0/255, 0/255],
        [255/255, 184/255, 205/255],
        [197/255, 0/255, 190/255],
        [100/255, 144/255, 255/255],
    ])

    if percent:
        bounds = np.array([25, 50, 60, 70, 80, 90, 95, 100, 105, 107, 110])
    else:
        bounds = np.array([0.25, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00, 1.05, 1.07, 1.10])

    norm = colors.BoundaryNorm(bounds, clevels_cmap.N, clip=True)

    zc = isodose_map.shape[2] // 2
    slice_xy = isodose_map[:, :, zc]
    slice_xy = np.ma.masked_where(slice_xy < bounds[0], slice_xy)

    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(slice_xy.T, origin="upper", cmap=clevels_cmap, norm=norm)
    ax.set_title("Isodose - slice centrale (XY)")
    plt.colorbar(im, ax=ax, boundaries=bounds, ticks=bounds)
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Anteprima salvata in: {out_png}")


########################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Normalizza una mappa di dose .mhd rispetto alla dose massima per ottenere una mappa di isodose"
    )
    parser.add_argument("-i", "--input", required=True, metavar="path",
                         help="File .mhd della mappa di dose in ingresso")
    parser.add_argument("-o", "--output", required=True, metavar="path",
                         help="File .mhd della mappa di isodose in uscita")
    parser.add_argument("--percent", action="store_true",
                         help="Esprime l'isodose in percentuale (0-100) invece che come frazione (0-1)")
    parser.add_argument("--preview", metavar="path",
                         help="Salva anche un PNG di anteprima (slice centrale) con colormap TPS-like")
    parser.add_argument("-v", "--verbose", action="store_true",
                         help="Stampa informazioni aggiuntive")
    args = parser.parse_args()

    if args.verbose:
        print(f"Carico la mappa di dose: {args.input}")

    nn, hs, x0, dose_map, direction = readMap(args.input)

    print(f"Dimensioni mappa: {nn}")
    print(f"Spacing (unita' del file, tipicamente mm): {hs}")

    max_dose = np.nanmax(dose_map)
    print(f"Dose massima: {max_dose}")

    if max_dose == 0 or np.isnan(max_dose):
        raise ValueError("La dose massima e' zero o NaN: impossibile normalizzare.")

    isodose_map = dose_map.astype(np.float64) / max_dose
    if args.percent:
        isodose_map *= 100.0

    writeMap(args.output, hs, x0, isodose_map.astype(np.float32), direction)
    print(f"Mappa di isodose salvata in: {args.output}")

    if args.preview:
        make_preview(isodose_map, args.preview, percent=args.percent)


if __name__ == "__main__":
    main()
