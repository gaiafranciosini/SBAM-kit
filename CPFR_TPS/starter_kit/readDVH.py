#!/usr/bin/env python3
import argparse
import os
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Evaluates DVH metrics")
    parser.add_argument(
        "-DVH", type=str, nargs="+", required=True, help="DVH file(s)"
    )
    parser.add_argument(
        "-preD", type=float, required=True, help="Prescription dose (in Gy)"
    )
    args = parser.parse_args()

    vol_thresholds = np.array([95, 90, 80, 70, 50, 40, 30, 25, 20, 10])

    for dvh_file in args.DVH:
        roi_name = os.path.splitext(os.path.basename(dvh_file))[0]
        print(roi_name)

        # Carica dose (in cGy) e volume (in %)
        dose, volume = np.loadtxt(dvh_file, unpack=True)

        # Calcolo delle dosi ai vari volumi (convertite in Gy)
        for v in vol_thresholds:
            d_vol = 0.01 * dose[np.argmin(np.abs(volume - v))]
            print(f"     D(V{v}) = {d_vol:.3f} Gy")

        # 1. Dose di prescrizione dall'argomento in ingresso (Gy)
        dpre = args.preD
        print(f"        D_pre = {dpre:.3f} Gy")

        # 2. Dose massima (al punto V < 0.01%) in Gy
        dmax = 0.01 * dose[np.argmin(np.abs(volume - 0.01))]
        print(f"       D_max = {dmax:.3f} Gy (V<0.01%)")

        # 3. Rapporto percentuale Dmax / Dpre
        perc = 100.0 * (dmax / dpre)
        print(f"       Dmax is {perc:.3f}% of D_pre")

        print("")
        print("###########")
        print("")


if __name__ == "__main__":
    main()
