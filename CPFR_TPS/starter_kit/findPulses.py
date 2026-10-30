import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Calcola impulsi per copertura DVH")
    parser.add_argument("input_file", help="file path")
    parser.add_argument("preD", type=float, help="Prescription dose (Gy)")
    parser.add_argument("preV", type=float, help="Prescription volume (%)")
    # Nuovo parametro facoltativo (default 100%)
    parser.add_argument("--preD_percent", type=float, default=100.0, help="Perc. della dose prescritta al preV (es. 95 per D95%=95%D_presc)")
    args = parser.parse_args()

    data = np.loadtxt(args.input_file, comments="#")
    dose_cgy = np.array(data[:, 0])
    volume = np.array(data[:, 1])

    # NOTA: Se usi ComputeDVH.py con DVH <= D, il volume da cercare è (100 - preV)
    # Se hai già corretto ComputeDVH.py per il DVH >= D standard, usa direttamente args.preV
    target_volume_index = args.preV
    targetDose = dose_cgy[np.argmin(np.abs(volume - target_volume_index))]

    # Calcola la dose richiesta in cGy applicando la percentuale
    required_dose_cgy = 100.0 * args.preD * (args.preD_percent / 100.0)

    # Calcolo impulsi arrotondato all'intero più vicino
    pulses = int(0.5 + required_dose_cgy / targetDose)
    print(pulses)

if __name__ == "__main__":
    main()
