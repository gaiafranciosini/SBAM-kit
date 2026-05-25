import SimpleITK as sitk
import argparse

def main():
    parser = argparse.ArgumentParser(description='De-crop (Pad) a dose map with zeros')
    parser.add_argument("-i", "--input", required=True, help="Path of the cropped dose map (.mhd)")
    parser.add_argument("-o", "--output", required=True, help="Path to save the de-cropped dose map")
    parser.add_argument("-cz", "--cropZ", type=int, required=True, help="Number of slices that were cut (your cropZ variable)")
    
    args = parser.parse_args()

    # 1. Leggi l'immagine della dose ritagliata
    dose_img = sitk.ReadImage(args.input)

    # 2. Definisci quanto "spazio" aggiungere
    # SimpleITK richiede le dimensioni come lista: [Pad_X, Pad_Y, Pad_Z]
    # Aggiungiamo 'cropZ' fette all'inizio (lower) dell'asse Z. Nessun taglio su X e Y.
    pad_lower = [0, 0, args.cropZ] 
    
    # Dal tuo script Bash precedente, non hai tagliato fette alla fine (izf era dimZ),
    # quindi non dobbiamo riaggiungere nulla alla fine (upper)
    pad_upper = [0, 0, 0] 

    # 3. Applica il padding
    # Il valore di riempimento è 0.0, che per la dose significa "zero Gray"
    decrapped_img = sitk.ConstantPad(dose_img, pad_lower, pad_upper, 0.0)

    # 4. Salva il risultato
    sitk.WriteImage(decrapped_img, args.output)
    print(f"De-cropping completato! Aggiunte {args.cropZ} fette di zeri all'inizio dell'asse Z.")

if __name__ == "__main__":
    main()
