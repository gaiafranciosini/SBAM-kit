import argparse
import numpy as np

def main():
    # Parser degli argomenti
    parser = argparse.ArgumentParser(description="Leggi file txt con due colonne: dose_cgy e volume")
    parser.add_argument("input_file", help="file path")
    parser.add_argument("preD", type=float, help="Prescrisption dose")
    parser.add_argument("preV", type=float, help="Prescrisption volume")
    args = parser.parse_args()

    # Carica i dati con NumPy, saltando le righe che iniziano con '#'
    data = np.loadtxt(args.input_file, comments="#")

    # Prima colonna → dose_cgy, seconda → volume
    dose_cgy = np.array(data[:, 0])
    volume = np.array(data[:, 1])

#    print("Dose (cGy):", dose_cgy)
#    print("Volume:", volume)

    rescaling_factor=100*np.array([args.preD])/dose_cgy[np.argmin(np.abs(volume-np.array([args.preV])))]
    print(rescaling_factor[0])

if __name__ == "__main__":
    main()
