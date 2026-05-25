#!/usr/bin/env python3
########################################################################
# Calcola la Dij ottimizzata (moltiplicata per le fluenze) mantenendo
# lo stesso formato binario sparso originale.
#
# Usage:
#    python dij_optiplan_scaled.py <Dij.bin> <optiPlan.txt> <fieldID>
########################################################################
import sys
import struct
import numpy as np

# -----------------------------------------------------------------------
def load_optiplan(fname, fieldID):
    """
    Legge optiPlan.txt e restituisce un dizionario {pbID: NumParticles}
    filtrato per il fieldID richiesto.
    """
    weights = {}
    with open(fname, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            fid = int(parts[1])
            pid = int(parts[2])
            npart = float(parts[3])
            if fid == fieldID:
                weights[pid] = npart
    if not weights:
        print(f'ATTENZIONE: nessun pencil beam trovato per fieldID={fieldID} in {fname}')
    else:
        print(f'Trovati {len(weights)} pencil beam per fieldID={fieldID}')
    return weights

# -----------------------------------------------------------------------
def scale_dij(dij_fname, out_fname, weights, fieldID):
    """
    Legge il file Dij binario, riscala i valori dei pencil beam appartenenti
    al fieldID usando i pesi di optiPlan, e scrive il risultato nel nuovo file.
    I pencil beam non appartenenti al fieldID o non ottimizzati vengono saltati
    o scritti con peso 0 (azzerati).
    """
    with open(dij_fname, 'rb') as fin, open(out_fname, 'wb') as fout:

        # --- header ---
        header = fin.read(40)
        if len(header) < 40:
            print('Errore: file Dij troppo corto o corrotto')
            sys.exit(-1)

        # Copia l'header originale direttamente nel file di output
        fout.write(header)

        nx, ny, nz, hx, hy, hz, x0, y0, z0, npb = struct.unpack('<3i6f1i', header)
        nn = np.array([nx, ny, nz])
        
        print(f'Dims     : {nn}')
        print(f'Num PB nel file Dij originale: {npb}')

        pb_scaled = 0
        pb_zeroed = 0
        pb_skipped = 0

        while True:
            chunk = fin.read(12)
            if len(chunk) < 12:
                break  # EOF

            pid, fid, numvxl = struct.unpack('3i', chunk)

            # Condizione: fa parte del campo e abbiamo un peso valido in optiPlan
            if fid == fieldID:
                weight = weights.get(pid, 0.0)
                
                # Scriviamo comunque il chunk di intestazione del PB (id, field, nvxl)
                fout.write(chunk)
                
                if numvxl > 0:
                    # Leggi indici e valori originali
                    ivxl_bytes = fin.read(numvxl * 4)
                    vals_bytes = fin.read(numvxl * 4)
                    
                    if weight > 0:
                        # Converti i valori, riscalali per la fluenza e riconverti in bytes
                        vals = np.frombuffer(vals_bytes, dtype='float32').copy()
                        vals_scaled = vals * weight
                        vals_bytes_scaled = vals_scaled.tobytes()
                        
                        # Scrivi indici immutati e valori riscaldati
                        fout.write(ivxl_bytes)
                        fout.write(vals_bytes_scaled)
                        pb_scaled += 1
                    else:
                        # Se il peso è zero (o assente), azzeriamo i valori di dose
                        fout.write(ivxl_bytes)
                        fout.write(np.zeros(numvxl, dtype='float32').tobytes())
                        pb_zeroed += 1
            else:
                # Se appartiene a un altro fieldID, saltiamo il blocco nel file sorgente
                # senza scriverlo nell'output (stiamo filtrando solo il fieldID richiesto)
                fin.seek(numvxl * 8, 1)  # salta indici e valori
                pb_skipped += 1

    print(f'\n=== Risultato elaborazione ===')
    print(f'PB riscalati con successo : {pb_scaled}')
    print(f'PB azzerati (peso 0)     : {pb_zeroed}')
    print(f'PB di altri campi ignorati: {pb_skipped}')

# -----------------------------------------------------------------------
def main():
    if len(sys.argv) < 4:
        print('Usage: python dij_optiplan_scaled.py <Dij.bin> <optiPlan.txt> <fieldID>')
        sys.exit(1)

    dij_fname  = sys.argv[1]
    opti_fname = sys.argv[2]
    fieldID    = int(sys.argv[3])

    out_fname = f'Dij_optimized_field{fieldID}.bin'

    print(f'\n=== Dij file sorgente : {dij_fname}')
    print(f'=== optiPlan          : {opti_fname}')
    print(f'=== fieldID filtrato  : {fieldID}')
    print(f'=== Dij ottimizzata   : {out_fname}\n')

    # 1) Carica i pesi dal piano ottimizzato
    weights = load_optiplan(opti_fname, fieldID)

    # 2) Riscala la Dij e genera il nuovo file binario
    scale_dij(dij_fname, out_fname, weights, fieldID)
    
    print(f'\nNuovo file Dij ottimizzato generato con successo: {out_fname}')

# -----------------------------------------------------------------------
if __name__ == '__main__':
    main()