import pandas as pd
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def calcola_intensita_max(region_file, optiplan_file, field_id_richiesto):
    try:
        # 1. Caricamento del Regionfile
        df_regions = pd.read_csv(region_file, sep=r'\s+')
        df_regions.columns = df_regions.columns.str.strip()
        
        # 2. Caricamento dell'Optiplan
        df_optiplan = pd.read_csv(
            optiplan_file, 
            sep=r'\s+', 
            comment='#', 
            names=['isource', 'fieldID', 'pbID', 'NumParticles']
        )
        
        field_id_richiesto = int(field_id_richiesto)

        # 3. Filtraggio per fieldID
        df_field = df_optiplan[df_optiplan['fieldID'] == field_id_richiesto].copy()
        
        if df_field.empty:
            print(f"[-] Nessun dato trovato per il fieldID: {field_id_richiesto}")
            return

        # 4. Ricerca massimo PIN per intensita
        indice_max = df_field['NumParticles'].idxmax()
        riga_max = df_field.loc[indice_max]
        id_pin_max = int(riga_max['pbID'])
        intensita_i_max = riga_max['NumParticles']

        # 5. Recupero dettagli dal Regionfile
        colonna_id = 'SPB_ID'
        if colonna_id not in df_regions.columns:
            print(f"[-] Errore: Colonna '{colonna_id}' non trovata. Colonne: {list(df_regions.columns)}")
            return

        info_geo_max = df_regions[df_regions[colonna_id] == id_pin_max]
        if info_geo_max.empty:
            print(f"[-] Errore: PIN {id_pin_max} non trovato nel Regionfile.")
            return

        io_simu_max = info_geo_max.iloc[0]['N_Part']

        # 6. Output dati PIN Max
        print("\n" + "="*50)
        print(f" ANALISI FIELD ID: {field_id_richiesto}")
        print("="*50)
        print(f"PIN con Intensita Max:  {id_pin_max}")
        print(f"Intensita (I_max):      {intensita_i_max:.6e}")
        print(f"Io_simu (N_Part):       {int(io_simu_max)}")
        print("-" * 50)

        # 7. Calcolo Npulse basato sul PIN massimo
        # Formula: Num_ele_max = intensita_i * 1e9 / io_simu
        num_ele_max = intensita_i_max * 1e9 / io_simu_max
        num_ele_per_pulse = 2e-6 / 1.6e-19
        npulse_float = num_ele_max / num_ele_per_pulse
        npulse = math.ceil(npulse_float)
        
        print(f"N. pulse richiesto: {npulse_float:.2f}  Arrotondato: {npulse}")

        # 8. Calcolo spessori
        mu = 0.051944 * 19.3  # Valore di mu fornito

        df_merged = pd.merge(
            df_field, 
            df_regions[[colonna_id, 'N_Part', 'Xmin_cell', 'Xmax_cell', 'Ymin_cell', 'Ymax_cell']], 
            left_on='pbID', 
            right_on=colonna_id
        )

        # Calcolo Io_reale per ogni PIN
        df_merged['Io_reale'] = (df_merged['N_Part'] * npulse * num_ele_per_pulse) / 1e9

        # Funzione per calcolare lo spessore con i tuoi vincoli
        def calcola_x(row):
            I = row['NumParticles']
            Io = row['Io_reale']
            if I <= 0 or Io <= 0:
                return 2.0
            val_x = -np.log(I / Io) / mu
            if val_x < 0: return 0.0
            if val_x > 2.0: return 8.0
            return val_x

        df_merged['x_thickness'] = df_merged.apply(calcola_x, axis=1)

        # 9. Preparazione e salvataggio file output
        df_output = df_merged[['pbID', 'Xmin_cell', 'Xmax_cell', 'Ymin_cell', 'Ymax_cell', 'x_thickness']].copy()
        df_output.rename(columns={'pbID': 'PIN_ID'}, inplace=True)
        
        output_name = f"risultati_field_{field_id_richiesto}.txt"
        df_output.to_csv(output_name, sep='\t', index=False, float_format="%.6f")
        print(f"[+] File generato: {output_name}")

        # 10. Visualizzazione 
        print("[*] Generazione grafico 3D...")
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')

        x_centers = ((df_output['Xmin_cell'] + df_output['Xmax_cell']) / 2).values
        y_centers = ((df_output['Ymin_cell'] + df_output['Ymax_cell']) / 2).values
        z_base = np.zeros_like(x_centers)
        dz = df_output['x_thickness'].values

        # Larghezza cubetti (0.2 se la cella = 0.2)
        dx = dy = 0.2 

        # Colormap
        colors = plt.cm.viridis(dz / 2.0)

        ax.bar3d(x_centers, y_centers, z_base, dx, dy, dz, color=colors, alpha=0.8, shade=True)

        ax.set_title(f'Volume Attenuatore 3D - Field {field_id_richiesto}')
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('Y [cm]')
        ax.set_zlabel('Spessore [cm]')
        ax.set_zlim(0, 2.1) # Un po piu di 2 per vederlo bene

        plt.savefig(f"volume_3d_field_{field_id_richiesto}.png")
        print(f"[+] Grafico salvato: volume_3d_field_{field_id_richiesto}.png")
        plt.show()

    except Exception as e:
        print(f"[-] Si  verificato un errore: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Utilizzo: python3 MakeAttenutator.py regionfile.txt optiplan.txt 1")
    else:
        calcola_intensita_max(sys.argv[1], sys.argv[2], sys.argv[3])