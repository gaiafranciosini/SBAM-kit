#!/usr/bin/env python3
import struct
import sys
import numpy as np
from mhd_io import mhd_read

def read_header(fin):
    header_bytes = fin.read(40)
    if len(header_bytes) < 40:
        return None, None
    [nx, ny, nz, hx, hy, hz, x0, y0, z0, npb] = struct.unpack("<3i6f1i", header_bytes)
    return np.array([nx, ny, nz]), header_bytes

def read_next_pb_raw(fin):
    block_header = fin.read(12)
    if len(block_header) < 12:
        return None
    [pid, fid, numvxl] = struct.unpack("3i", block_header)
    
    if numvxl > 0:
        Ivxl_bytes = fin.read(numvxl * 4)
        Vals_bytes = fin.read(numvxl * 4)
    else:
        Ivxl_bytes = b""
        Vals_bytes = b""
        
    return pid, fid, numvxl, Ivxl_bytes, Vals_bytes


if __name__ == "__main__":
    if len(sys.argv) < 7:
        print("Usage: python script.py Dijopt.bin Dijfluka.bin fieldID metodo mask.mhd spessori_in.txt")
        print("NOTA: Il parametro 'metodo' viene ignorato in questa versione ottimizzata per il CORE del fascio.")
        sys.exit(1)

    file_opt = sys.argv[1]
    file_fluka = sys.argv[2]
    FID = int(sys.argv[3])
    # metodo = int(sys.argv[4]) # Sovrascritto dalla logica del Core
    file_mask = sys.argv[5]
    file_spessori_in = sys.argv[6]

    # --- CONFIGURAZIONE FISICA E GEOMETRICA (TUTTO IN CM) ---
    mu = 0.20*19.3             # Coefficiente di attenuazione lineare in cm^-1
    soglia_core = 0.70   # Soglia adattiva: considera solo i voxel > 70% del rispettivo picco massimo
    top_n_voxels = 20    # Numero di voxel di picco per stabilizzare la media della cresta

    # --- 1. LETTURA FILE SPESSORI INIZIALI ---
    print("Lettura file spessori iniziali (in cm)...")
    spessori_iniziali = {}
    righe_originali_meta = {}
    with open(file_spessori_in, "r") as f_in:
        for riga in f_in:
            riga_strip = riga.strip()
            if not riga_strip or riga_strip.startswith("#"): continue
            parti = riga_strip.split()
            if len(parti) < 6: continue
            j_val = int(parti[0])
            spessori_iniziali[j_val] = float(parti[5])
            righe_originali_meta[j_val] = " ".join(parti[0:5])

    # --- 2. CARICAMENTO MASCHERA ED ESTRAZIONE INDICI ---
    print(f"Caricamento maschera {file_mask}...")
    voxels_mask = mhd_read(file_mask)
    mask_data = voxels_mask['Map']
    
    # Appiattimento in ordine Fortran ('F') per raccordarsi alla struttura dei file binari Dij
    mask_flat = np.reshape(mask_data, -1, order="F")
    indici_maschera_set = set(np.where(mask_flat > 0)[0])
    del mask_data, mask_flat # Liberazione immediata della memoria RAM
    print("-> Maschera convertita correttamente in set di ricerca rapida.")

    # --- 3. APERTURA FILE DIJ IN STREAMING (ZONA RAM SICURA) ---
    f_opt = open(file_opt, "rb")
    f_fluka = open(file_fluka, "rb")

    nn_opt, _ = read_header(f_opt)
    nn_fluka, header_originale = read_header(f_fluka)

    # Indicizzazione rapida del file Target (legge solo intestazioni di 12 byte, salta i vettori pesanti)
    print("Indicizzazione rapida del file OPT (Target)...")
    opt_offsets = {}
    while True:
        pos = f_opt.tell()
        block = f_opt.read(12)
        if len(block) < 12: break
        [pid, fid, numvxl] = struct.unpack("3i", block)
        if fid == FID:
            opt_offsets[pid] = (pos, numvxl)
        f_opt.seek(numvxl * 8, 1)

    # Configurazione file di output
    out_binary_fname = f"Dijfluka_riscalata_Core_Field_{FID}.bin"
    out_text_fname = f"Spessori_Aggiornati_Core_Field_{FID}.txt"

    fout = open(out_binary_fname, "wb")
    ftxt = open(out_text_fname, "w")
    fout.write(header_originale)
    
    print(f"\n--- AVVIO CALCOLO OTTIMIZZAZIONE DEL CORE (Scala: cm) ---")
    print(f"{'j':<6} | {'X_init (cm)':<11} | {'Fattore scalo k':<16} | {'X_fine (cm)':<11} | {'Delta_X (cm)':<12}")
    print("-" * 70)

    while True:
        res_fluka = read_next_pb_raw(f_fluka)
        if res_fluka is None: break
        
        j, fid, numvxl_fluka, Ivxl_f_bytes, Vals_f_bytes = res_fluka
        if numvxl_fluka == 0 and len(Ivxl_f_bytes) == 0 and j == 0: break 
        if fid != FID: continue

        fattore_scalo_k = 1.0
        delta_x_j = 0.0

        if j in opt_offsets:
            pos_opt, numvxl_opt = opt_offsets[j]
            f_opt.seek(pos_opt + 12)
            Ivxl_o_bytes = f_opt.read(numvxl_opt * 4)
            Vals_o_bytes = f_opt.read(numvxl_opt * 4)

            # Conversione dei flussi binari in array leggeri (Nessuna matrice 3D densa allocata)
            Ivxl_o = np.frombuffer(Ivxl_o_bytes, dtype="uint32", count=numvxl_opt)
            Vals_o = np.frombuffer(Vals_o_bytes, dtype="float32", count=numvxl_opt)
            
            Ivxl_f = np.frombuffer(Ivxl_f_bytes, dtype="uint32", count=numvxl_fluka)
            Vals_f = np.frombuffer(Vals_f_bytes, dtype="float32", count=numvxl_fluka)

            # Estrazione posizioni interne alla maschera geometrica
            vxl_f_in_mask = [i for i, vxl in enumerate(Ivxl_f) if vxl in indici_maschera_set]
            vxl_o_in_mask = [i for i, vxl in enumerate(Ivxl_o) if vxl in indici_maschera_set]

            if len(vxl_f_in_mask) > 0 and len(vxl_o_in_mask) > 0:
                vals_f_mask = Vals_f[vxl_f_in_mask]
                vals_o_mask = Vals_o[vxl_o_in_mask]

                # Individuazione indipendente delle vette del profilo di dose
                max_f_locale = np.max(vals_f_mask)
                max_o_locale = np.max(vals_o_mask)

                if max_f_locale > 1e-5 and max_o_locale > 1e-5:
                    # Filtro adattivo: isoliamo solo le creste dei due profili
                    core_f = vals_f_mask[vals_f_mask > (soglia_core * max_f_locale)]
                    core_o = vals_o_mask[vals_o_mask > (soglia_core * max_o_locale)]

                    core_f_sorted = np.sort(core_f)
                    core_o_sorted = np.sort(core_o)
                    
                    n_vxl_f = min(top_n_voxels, len(core_f_sorted))
                    n_vxl_o = min(top_n_voxels, len(core_o_sorted))

                    media_core_f = np.mean(core_f_sorted[-n_vxl_f:])
                    media_core_o = np.mean(core_o_sorted[-n_vxl_o:])

                    if media_core_f > 1e-6:
                        # Fattore k pulito: es. 2.5 / 13.0 -> ~0.19
                        fattore_scalo_k = media_core_o / media_core_f

                    if fattore_scalo_k > 1e-6:
                        # Calcolo dello spessore necessario a ritroso dal fattore k
                        delta_x_j = - (1.0 / mu) * np.log(fattore_scalo_k)

        # --- BILANCIO E AGGIORNAMENTO GEOMETRICO ---
        x_iniziale = spessori_iniziali.get(j, 0.0)
        x_finale = x_iniziale + delta_x_j
        
        # Vincolo fisico inferiore: lo spessore del compensatore non può scendere sotto lo zero
        if x_finale < 0: 
            x_finale = 0.0
            delta_x_j = -x_iniziale

        # Stampa del log in tempo reale per monitorare l'evoluzione di ogni pencil beam
        print(f"{j:<6} | {x_iniziale:<11.4f} | {fattore_scalo_k:<16.4f} | {x_finale:<11.4f} | {delta_x_j:+.4f}")

        # --- 4. EMISSIONE REPORT TXT ---
        if j in righe_originali_meta:
            ftxt.write(f"{righe_originali_meta[j]} {x_finale:.6f}\n")
        else:
            ftxt.write(f"{j} 0 0 0 0 {x_finale:.6f}\n")

        # --- 5. RISCALAMENTO DIRETTO E TRASMISSIONE SU DISCO DEL FILE DIJ ---
        fout.write(struct.pack("3i", j, FID, numvxl_fluka))
        if numvxl_fluka > 0:
            Vals_f_numpy = np.frombuffer(Vals_f_bytes, dtype="float32", count=numvxl_fluka)
            # Moltiplicazione diretta sul vettore compresso nativo di FLUKA
            Vals_f_riscalati = Vals_f_numpy * fattore_scalo_k
            
            fout.write(Ivxl_f_bytes)                 
            fout.write(Vals_f_riscalati.tobytes())   

    print("\n")
    f_opt.close(); f_fluka.close(); fout.close(); ftxt.close()
    print("=" * 70)
    print("Elaborazione e riscalamento completati con successo!")
    print(f"Mappa Dij riscalata: {out_binary_fname}")
    print(f"Report spessori:     {out_text_fname}")