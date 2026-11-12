import SimpleITK as sitk
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="De-crop (pad) a dose map with zeros and align its geometry to the reference CT"
    )
    parser.add_argument("-i", "--input", required=True, help="Path of the cropped dose map (.mhd)")
    parser.add_argument("-o", "--output", required=True, help="Path to save the de-cropped dose map")
    parser.add_argument("-cz", "--cropZ", type=int, required=True,
                        help="Number of slices that were cut at the start of Z (cropZ)")
    parser.add_argument("-r", "--reference", default=None,
                        help="Reference CT (e.g. Phantom.mhd): origin, spacing and direction are copied from it")
    parser.add_argument("--tol", type=float, default=1e-4,
                        help="Tolerance on spacing mismatch between dose and CT (mm)")

    args = parser.parse_args()

    # 1. Leggi la dose ritagliata
    dose_img = sitk.ReadImage(args.input)

    # 2. Padding di cropZ fette di zeri all'inizio dell'asse Z (nessun taglio su X, Y e nessuno in coda)
    pad_lower = [0, 0, args.cropZ]
    pad_upper = [0, 0, 0]
    padded = sitk.ConstantPad(dose_img, pad_lower, pad_upper, 0.0)

    print(f"Dose in ingresso : size={dose_img.GetSize()} origin={dose_img.GetOrigin()}")
    print(f"Dose dopo padding: size={padded.GetSize()} origin={padded.GetOrigin()}")

    # 3. Allineamento geometrico alla CT (se fornita)
    if args.reference is not None:
        ref = sitk.ReadImage(args.reference)

        if padded.GetSize() != ref.GetSize():
            raise SystemExit(
                f"ERRORE: dimensioni diverse. Dose={padded.GetSize()} CT={ref.GetSize()}. "
                f"Controlla cropZ ({args.cropZ})."
            )

        if not np.allclose(padded.GetSpacing(), ref.GetSpacing(), atol=args.tol):
            raise SystemExit(
                f"ERRORE: spacing diverso. Dose={padded.GetSpacing()} CT={ref.GetSpacing()}"
            )

        # Differenza di origine prima della correzione (attesa: ~ mezzo voxel in x, y, z)
        diff = np.array(ref.GetOrigin()) - np.array(padded.GetOrigin())
        half = np.array(ref.GetSpacing()) / 2.0
        print(f"Differenza origine CT - dose: {tuple(diff)}  (mezzo voxel = {tuple(half)})")
        if not np.allclose(diff, half, atol=1e-3):
            print("ATTENZIONE: lo shift NON e' esattamente mezzo voxel. "
                  "Verifica cropZ e l'origine della CT ritagliata usata in FLUKA.")

        # Copia origin, spacing e direction dalla CT (i valori dei voxel non cambiano)
        padded.CopyInformation(ref)
        print(f"Dose finale: origin={padded.GetOrigin()} (allineata alla CT)")
    else:
        print("Nessuna CT di riferimento: applicato solo il padding, origine non corretta.")

    # 4. Salva
    sitk.WriteImage(padded, args.output)
    print(f"De-cropping completato! Aggiunte {args.cropZ} fette di zeri all'inizio dell'asse Z.")


if __name__ == "__main__":
    main()
