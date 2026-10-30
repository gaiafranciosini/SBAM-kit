#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verifica la copertura del PTV con le isodose 95/90/80/50/10% della dose
prescritta: copertura e Dice volumetrici 3D, analisi slice per slice nei tre
piani (solo slice con PTV) e verifica del criterio V(Rx) >= 95%.

La dose prescritta NON ha un valore di default: va sempre indicata con
--dose-prescription-gy.

Esempio:
    python3 copertura_ptv_isodose.py -i isodose.mhd -p PTV_plan.mhd \
        --dose-prescription-gy <Gy> --normalization prescription
"""

import argparse

import numpy as np
import SimpleITK as sitk

W = 112
LEVELS = (100, 95, 90, 80, 50, 10)          # % della dose prescritta
MAIN_LEVELS = (95, 90, 80, 50, 10)          # livelli delle analisi principali
VOLUME_CRITERION = 95.0                     # criterio clinico V(Rx) >= 95%
# (nome, asse del piano, piano, asse numpy). SimpleITK e' (X,Y,Z), numpy [Z,Y,X].
PLANES = (("ASSIALE", "Z", "XY", 0), ("CORONALE", "Y", "XZ", 1), ("SAGITTALE", "X", "YZ", 2))


def same_geometry(a, b):
    return (a.GetSize() == b.GetSize()
            and np.allclose(a.GetSpacing(), b.GetSpacing(), atol=1e-6)
            and np.allclose(a.GetOrigin(), b.GetOrigin(), atol=1e-6)
            and np.allclose(a.GetDirection(), b.GetDirection(), atol=1e-6))


def cov_dice(ptv, iso, thr):
    """Copertura 3D e Dice 3D tra PTV e regione isodose >= thr.
    Coverage = V(P∩I)/V(P);  Dice = 2 V(P∩I) / (V(P) + V(I))."""
    region = np.isfinite(iso) & (iso >= thr)
    inter = int(np.count_nonzero(ptv & region))
    n_ptv = int(np.count_nonzero(ptv))
    n_reg = int(np.count_nonzero(region))
    return {
        "covered": inter, "vol_ptv": n_ptv, "vol_iso": n_reg,
        "coverage": 100.0 * inter / n_ptv if n_ptv else np.nan,
        "dice": 2.0 * inter / (n_ptv + n_reg) if n_ptv + n_reg else np.nan,
    }


def slice_rows(ptv, iso, img, axis, columns):
    """Per ogni slice con almeno un voxel di PTV: (idx, posizione mm, voxel PTV,
    [voxel coperti per ciascuna soglia in columns])."""
    sitk_axis = 2 - axis
    rows = []
    for i in range(ptv.shape[axis]):
        p = np.take(ptv, i, axis=axis)
        total = int(np.count_nonzero(p))
        if total == 0:
            continue
        d = np.take(iso, i, axis=axis)
        ok = p & np.isfinite(d)
        idx = [0, 0, 0]
        idx[sitk_axis] = i
        pos = img.TransformIndexToPhysicalPoint(tuple(idx))[sitk_axis]
        rows.append((i, pos, total, [int(np.count_nonzero(ok & (d >= thr))) for _, thr in columns]))
    return rows


def print_slices(title, ax, rows, columns, slice_thr):
    print(f"\n{'=' * W}\n{title}\n{'=' * W}")
    if not rows:
        print("Nessuna slice contenente PTV.")
        return
    a = ax.lower()
    first = rows[0][1]
    print(f"PTV presente lungo {ax} dalla slice {rows[0][0]} alla {rows[-1][0]} "
          f"({first:.2f} -> {rows[-1][1]:.2f} mm), {len(rows)} slice con PTV.")
    print(f"{ax}_depth_PTV_mm = distanza dalla PRIMA slice con PTV lungo {ax} "
          "(NON e' profondita' radiologica ne' distanza dal fascio).\n")
    print(f"{a + '_idx':>6} {a + '_mm':>9} {a + '_depth_PTV_mm':>17} {'voxel_PTV':>10}"
          + "".join(f" {'>=' + n:>10} {'cov_' + n:>10}" for n, _ in columns))
    print("-" * W)
    ok_slices = [0] * len(columns)
    for idx, pos, total, counts in rows:
        line = f"{idx:6d} {pos:9.2f} {abs(pos - first):17.2f} {total:10d}"
        for j, c in enumerate(counts):
            cov = 100.0 * c / total
            ok_slices[j] += cov >= slice_thr
            line += f" {c:10d} {cov:10.2f}"
        print(line)
    print(f"\n{'-' * W}\nRIASSUNTO (slice con copertura >= {slice_thr:.2f}%)")
    for (n, _), k in zip(columns, ok_slices):
        print(f"  {n:>8}: {k}/{len(rows)}")


def main():
    ap = argparse.ArgumentParser(
        description="Copertura del PTV con isodose 95/90/80/50/10% della dose prescritta "
                    "(3D e slice per slice nei piani assiale, coronale, sagittale).")
    ap.add_argument("-i", "--isodose", required=True, help="Mappa di isodose .mhd")
    ap.add_argument("-p", "--ptv", required=True, help="Mappa PTV .mhd")
    ap.add_argument("--dose-prescription-gy", type=float, required=True,
                    help="Dose prescritta in Gy (obbligatoria, nessun default)")
    ap.add_argument("--normalization", choices=("dmax", "prescription"), default="dmax",
                    help="Riferimento della mappa: dmax (vecchie mappe) o prescription (nuove)")
    ap.add_argument("--dmax-gy", type=float,
                    help="Dmax effettiva (Gy) corrispondente a 1.0/100 nella mappa; "
                         "obbligatoria con --normalization dmax")
    ap.add_argument("--isodose-scale", choices=("auto", "fraction", "percent"), default=None,
                    help="fraction: 1=100%%; percent: 100=100%%. Default: fraction con "
                         "prescription, auto con dmax")
    ap.add_argument("--tolerance-gy", type=float, default=0.05,
                    help="Tolleranza (Gy) sottratta alle soglie (solo limite inferiore); 0 = nessuna")
    ap.add_argument("--slice-threshold", type=float, default=100.0,
                    help="%% minima di PTV coperta in una slice per considerarla coperta")
    args = ap.parse_args()

    rx, tol = args.dose_prescription_gy, args.tolerance_gy
    if not np.isfinite(rx) or rx <= 0:
        ap.error("--dose-prescription-gy deve essere un valore finito > 0.")
    if args.normalization == "dmax" and (args.dmax_gy is None or not np.isfinite(args.dmax_gy)
                                         or args.dmax_gy <= 0):
        ap.error("Con --normalization dmax indicare --dmax-gy > 0.")
    if not np.isfinite(tol) or tol < 0:
        ap.error("--tolerance-gy deve essere finita e non negativa (Gy).")
    if not (np.isfinite(args.slice_threshold) and 0 <= args.slice_threshold <= 100):
        ap.error("--slice-threshold deve essere tra 0 e 100 (%).")

    # --- lettura (SimpleITK (X,Y,Z) -> numpy [Z,Y,X])
    iso_img, ptv_img = sitk.ReadImage(args.isodose), sitk.ReadImage(args.ptv)
    if not same_geometry(iso_img, ptv_img):
        raise RuntimeError("ERRORE: PTV e mappa di isodose non hanno la stessa geometria.")
    iso = sitk.GetArrayFromImage(iso_img).astype(np.float64)
    mask = sitk.GetArrayFromImage(ptv_img) > 0
    total = int(np.count_nonzero(mask))
    if total == 0:
        raise ValueError("ERRORE: nessun voxel con PTV > 0.")

    # --- scala e soglie
    req = args.isodose_scale or ("fraction" if args.normalization == "prescription" else "auto")
    if req == "fraction" or (req == "auto" and np.nanmax(iso) <= 1.5):
        scale, scale_label = 1.0, "frazione (1.0 = 100% del riferimento)"
    else:
        scale, scale_label = 100.0, "percentuale (100 = 100% del riferimento)"
    ref = rx if args.normalization == "prescription" else args.dmax_gy
    nom_gy = {L: L / 100.0 * rx for L in LEVELS}
    tol_gy = {L: max(0.0, g - tol) for L, g in nom_gy.items()}   # tolleranza: solo limite inferiore
    nom = {L: scale * g / ref for L, g in nom_gy.items()}
    tolr = {L: scale * g / ref for L, g in tol_gy.items()}
    main_cols = [(f"{L}%", tolr[L]) for L in MAIN_LEVELS]

    print(f"\n{'=' * W}\nCONTROLLO COPERTURA PTV - ISODOSE "
          f"{'/'.join(str(L) + '%' for L in MAIN_LEVELS)} (SOGLIE CON TOLLERANZA)\n{'=' * W}")
    print(f"Scala isodose in ingresso : {scale_label}")
    print(f"Normalizzazione mappa     : {args.normalization} "
          f"(dose di riferimento {ref:.3f} Gy; la mappa NON viene modificata)")
    print(f"Dose prescritta           : {rx:.3f} Gy")
    print(f"Tolleranza                : {tol:.3f} Gy (abbassa SOLO la soglia minima)")
    print(f"Voxel totali PTV          : {total}")
    print("Soglie (livello di DOSE, % di prescrizione; non la quota del volume PTV):")
    for L in LEVELS:
        print(f"  {L:3d}% = {nom_gy[L]:7.3f} Gy (isodose >= {nom[L]:.6f}) -> con tolleranza "
              f">= {tol_gy[L]:.3f} Gy (isodose >= {tolr[L]:.6f})")
    print("Un voxel del PTV e' COPERTO se isodose >= soglia (eventualmente ridotta della tolleranza).")

    # --- copertura e Dice 3D
    print(f"\n{'-' * W}\nCOPERTURA E DICE VOLUMETRICI 3D (soglie con tolleranza)\n{'-' * W}")
    print(f"{'Livello':>8} {'iso >=':>10} {'PTV coperti':>12} {'PTV sotto':>10} "
          f"{'Coverage%':>10} {'Vol iso':>10} {'Dice':>8} {'Dice%':>7}")
    res = {}
    for L in MAIN_LEVELS:
        r = res[L] = cov_dice(mask, iso, tolr[L])
        print(f"{L:7d}% {tolr[L]:10.6f} {r['covered']:12d} {total - r['covered']:10d} "
              f"{r['coverage']:10.2f} {r['vol_iso']:10d} {r['dice']:8.4f} {100 * r['dice']:7.2f}")
    print("Dice = 2 V(PTV∩isodose) / (V(PTV) + V(isodose)): 1 = coincidenza perfetta, 0 = nessuna "
          "sovrapposizione.\nLa coverage dice QUANTO PTV e' coperto; il Dice tiene conto anche "
          "dell'isodose fuori dal PTV.")

    # --- confronto nominale / tollerato
    print(f"\n{'=' * W}\nCONFRONTO NOMINALE / CON TOLLERANZA - VOLUME 3D E DICE\n{'=' * W}")
    print(f"{'Dose % Rx':>10} {'Soglia nom. Gy':>15} {'Soglia tol. Gy':>15} {'Cov. nom.':>11} "
          f"{'Cov. tol.':>11} {'Dice nom.':>10} {'Dice tol.':>10}")
    for L in LEVELS:
        n, t = cov_dice(mask, iso, nom[L]), cov_dice(mask, iso, tolr[L])
        print(f"{L:10d} {nom_gy[L]:15.3f} {tol_gy[L]:15.3f} {n['coverage']:10.2f}% "
              f"{t['coverage']:10.2f}% {n['dice']:10.4f} {t['dice']:10.4f}")

    # --- slice per slice: livelli principali (con tolleranza), solo slice con PTV
    for name, ax, plane, axis in PLANES:
        rows = slice_rows(mask, iso, iso_img, axis, main_cols)
        print_slices(f"ANALISI {name}: PIANO {plane}, SCORRENDO LUNGO {ax}", ax, rows,
                     main_cols, args.slice_threshold)

    # --- criterio V(Rx) >= 95%, nominale e con tolleranza
    n, t = cov_dice(mask, iso, nom[100]), cov_dice(mask, iso, tolr[100])
    print(f"\n{'=' * W}\nVERIFICA DEL CRITERIO V({rx:.2f} Gy) >= {VOLUME_CRITERION:.0f}% "
          f"- NOMINALE E CON TOLLERANZA\n{'=' * W}")
    print(f"Dose prescritta                       : {rx:.3f} Gy")
    print(f"Soglia nominale / con tolleranza      : {nom_gy[100]:.3f} / {tol_gy[100]:.3f} Gy")
    print(f"Voxel PTV totali                      : {total}")
    print(f"Voxel coperti nominalmente            : {n['covered']} ({n['coverage']:.2f}%)")
    print(f"Voxel coperti con tolleranza          : {t['covered']} ({t['coverage']:.2f}%)")
    print(f"Voxel aggiunti grazie alla tolleranza : {t['covered'] - n['covered']}")
    print(f"Dice nominale / con tolleranza        : {n['dice']:.4f} / {t['dice']:.4f}")
    print(f"Criterio clinico V({rx:.2f} Gy) >= {VOLUME_CRITERION:.0f}%  : "
          f"{'RAGGIUNTO' if n['coverage'] >= VOLUME_CRITERION else 'NON RAGGIUNTO'}")
    print(f"Criterio con tolleranza (esplorativo) : "
          f"{'RAGGIUNTO' if t['coverage'] >= VOLUME_CRITERION else 'NON RAGGIUNTO'}")
    print("NOTA: il risultato clinico si basa sulla soglia NOMINALE; la tolleranza non cambia la "
          "prescrizione.\nI voxel sopra la soglia restano sempre coperti.")
    rx_cols = [("Rx", nom[100]), ("Rx-tol", tolr[100])]
    for name, ax, plane, axis in PLANES:
        rows = slice_rows(mask, iso, iso_img, axis, rx_cols)
        print_slices(f"ANALISI V({rx:.2f} Gy) {name} (lungo {ax}) - NOMINALE E CON TOLLERANZA",
                     ax, rows, rx_cols, args.slice_threshold)
    print("Il test sulle singole slice non sostituisce il criterio volumetrico 3D.")

    # --- risultato finale
    print(f"\n{'=' * W}\nRISULTATO VOLUMETRICO 3D\n{'=' * W}")
    if res[95]["covered"] == total:
        print(">>> L'INTERO PTV 3D e' contenuto nella regione con isodose >= 95%.")
    else:
        print(f">>> L'isodose del 95% NON contiene completamente il PTV 3D "
              f"(copertura {res[95]['coverage']:.2f}%).")
        if res[90]["covered"] == total:
            print(">>> L'INTERO PTV 3D e' contenuto nella regione con isodose >= 90%.")
        else:
            print(f">>> Nemmeno l'isodose del 90% contiene completamente il PTV 3D "
                  f"(copertura {res[90]['coverage']:.2f}%).")
    print()
    for L in MAIN_LEVELS:
        print(f">>> {L:2d}%: copertura 3D {res[L]['coverage']:.2f}% | "
              f"Dice 3D {res[L]['dice']:.4f} ({100 * res[L]['dice']:.2f}%)")
    print(f"\n{'-' * W}\nNOTA: analizzate solo le slice con almeno 1 voxel di PTV "
          "(assiale XY lungo Z, coronale XZ lungo Y, sagittale YZ lungo X).\n"
          "Le soglie nominali sono (% prescrizione x dose prescritta) / riferimento; quelle tollerate "
          "sono max(0, soglia_Gy - tolleranza_Gy) / riferimento, moltiplicate per la scala della mappa.")


if __name__ == "__main__":
    main()
