#!/usr/bin/env python3

import struct
import sys
import numpy as np
# Assicurati che mhd_io sia nella stessa cartella o nel PYTHONPATH
from mhd_io import *


def read_all_pencilbeams(fname, FID):
    """Legge il file binario Dij ed estrae un dizionario {PID: Matrice_Dose}

    per tutti i pencil beam appartenenti al campo FID.
    Restituisce anche i metadati dell'header originale.
    """
    try:
        fin = open(fname, "rb")
    except:
        print(f"IO error: cannot open input file {fname}")
        return None, None

    # Leggi l'header da 40 byte
    header_bytes = fin.read(40)
    [nx, ny, nz, hx, hy, hz, x0, y0, z0, npb] = struct.unpack(
        "<3i6f1i", header_bytes
    )
    nn = np.array([nx, ny, nz])

    pb_doses = {}

    while True:
        try:
            [pid, fid, numvxl] = struct.unpack("3i", fin.read(12))
        except:
            break

        if fid == FID:
            if pid not in pb_doses:
                pb_doses[pid] = np.zeros(np.prod(nn), dtype="float32")

            if numvxl > 0:
                Ivxl = np.frombuffer(
                    fin.read(numvxl * 4), dtype="uint32", count=numvxl
                )
                Vals = np.frombuffer(
                    fin.read(numvxl * 4), dtype="float32", count=numvxl
                )
                pb_doses[pid][Ivxl] += Vals
        else:
            fin.seek(numvxl * 8, 1)

    fin.close()

    # Rimodella le matrici in 3D (Fortran order)
    for pid in pb_doses:
        pb_doses[pid] = np.reshape(pb_doses[pid], nn, order="F")

    return pb_doses, header_bytes


# ==============================================================================
# MAIN PROGRAM: RISCALAMENTO DIJ BASATO SUI MASSIMI ROBUSTI
# ==============================================================================
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py Dijopt.bin Dijfluka.bin fieldID")
        sys.exit(1)

    file_opt = sys.argv[1]
    file_fluka = sys.argv[2]
    FID = int(sys.argv[3])

    mu = 0.2  # Coefficiente di attenuazione del Tungsteno (W)
    top_voxels = 10  # Numero di voxel più caldi su cui fare la media del picco

    print(f"--- Caricamento dei file Dij per il Campo {FID} ---")
    print("Leggo OPT (Target)...")
    pb_opt, _ = read_all_pencilbeams(file_opt, FID)
    print("Leggo FLUKA (Monte Carlo)...")
    pb_fluka, header_originale = read_all_pencilbeams(file_fluka, FID)

    if not pb_opt or not pb_fluka:
        print("Errore nel caricamento dei file o campo non trovato.")
        sys.exit(-1)

    # Definiamo il nome del nuovo file Dij binario in uscita
    out_binary_fname = f"Dijfluka_riscalata_Field_{FID}.bin"

    print(
        f"\n-> Apertura file binario di uscita in scrittura: {out_binary_fname}"
    )
    try:
        fout = open(out_binary_fname, "wb")
    except:
        print(f"Errore: Impossibile creare il file di uscita {out_binary_fname}")
        sys.exit(-1)

    # 1. Scrivi l'header identico a quello originale (40 byte)
    fout.write(header_originale)

    print(f"\n--- Calcolo spessori (indipendente da shift spaziali) ---")
    print(f"{'PB (j)':<10}{'Max OPT (vxl)':<15}{'Max FLUKA (vxl)':<15}{'Delta_x (W)':<15}")
    print("-" * 55)

    pids_comuni = set(pb_opt.keys()).intersection(set(pb_fluka.keys()))

    for j in sorted(pids_comuni):
        M_opt_j = pb_opt[j]
        M_fluka_j = pb_fluka[j]

        # Estraiamo i 10 voxel più alti di OPT e di FLUKA indipendentemente da dove si trovano
        # Questo elimina l'errore dovuto allo shift geometrico visibile nei profili
        max_opt_robusto = np.mean(
            np.sort(M_opt_j.flatten())[-top_voxels:]
        )
        max_fluka_robusto = np.mean(
            np.sort(M_fluka_j.flatten())[-top_voxels:]
        )

        delta_x_j = 0.0

        # Calcoliamo lo spessore solo se il raggio ha effettivamente depositato dose
        if max_fluka_robusto > 1e-6 and max_opt_robusto > 1e-6:
            rapporto_picchi = max_opt_robusto / max_fluka_robusto
            delta_x_j = -(1.0 / mu) * np.log(rapporto_picchi)

        print(
            f"{j:<10}{max_opt_robusto:<15.4e}{max_fluka_robusto:<15.4e}{delta_x_j:<15.4f}"
        )

        # Applicazione dell'attenuazione esponenziale alla matrice 3D del fascio j
        M_fluka_j_riscalata = M_fluka_j * np.exp(-mu * delta_x_j)

        # Ripristiniamo la matrice in un array monodimensionale (F order) per estrarre gli indici
        M_flat = np.reshape(M_fluka_j_riscalata, -1, order="F")

        # Selezioniamo solo i voxel attivi (dose > 0) per preservare la compressione Dij
        indici_attivi = np.where(M_flat > 0)[0].astype("uint32")
        valori_attivi = M_flat[indici_attivi].astype("float32")
        numvxl = len(indici_attivi)

        # 2. Scrittura del blocco compresso nel nuovo file binario (struttura Fred)
        fout.write(struct.pack("3i", j, FID, numvxl))

        if numvxl > 0:
            # Array degli indici (Ivxl)
            fout.write(indici_attivi.tobytes())
            # Array dei valori riscalati (Vals)
            fout.write(valori_attivi.tobytes())

    fout.close()
    print("\n" + "=" * 55)
    print(f"Nuovo file Dij binario compresso generato con successo!")
    print(f"File di uscita: {out_binary_fname}")
