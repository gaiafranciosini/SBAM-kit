#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import argparse

parser = argparse.ArgumentParser(description='Presents label of the TXT files')
parser.add_argument("-label1", help="label of TXTs files 1", type=str, default="")
parser.add_argument("-dir1", help="folder of TXTs files 2", type=str, default="")
parser.add_argument("-roi", help="list of ROIs name", type=str, nargs='+')
args = parser.parse_args()

fileLabel_1 = args.label1
fileDir_1 = args.dir1

my_Color = ['r','g','b','k','y','c','m','tab:orange','tab:purple','tab:pink','tab:brown','yellowgreen','slategray','hotpink', 'tan','lightcoral', 'tomato']

# Configurazione font e figura
plt.rcParams.update({'font.size': 15})
fig = plt.figure(figsize=(8, 8))
ax1 = plt.subplot(111)

plt.xlabel('Dose [Gy/s]', fontsize=15)
plt.ylabel('Volume [%]', fontsize=15)

# Loop per il plot delle ROI
c = 0
for ROIname in args.roi:
    try:
        filepath = fileDir_1 + '/' + ROIname + fileLabel_1 + '.txt'
        x1, y1 = np.loadtxt(filepath, comments=['#', 'idx'], usecols=(0, 1), unpack=True)
        plt.plot(x1, y1, color=my_Color[c % len(my_Color)], label=ROIname, linewidth=1.8)
        c += 1
    except Exception as e:
        print(f"Errore nel caricamento di {ROIname}: {e}")

# Aggiunta linee verticali tratteggiate
plt.axvline(x=10, color='grey', linestyle='--', linewidth=1.5, label='10 Gy/s')
plt.axvline(x=40, color='darkred', linestyle='--', linewidth=1.5, label='40 Gy/s')

# Regolazione posizione per far stare la legenda sotto
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.15, box.width, box.height * 0.85])

# Creazione della legenda unica con ROI e linee verticali
handles, labels = ax1.get_legend_handles_labels()
legend1 = plt.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                     fancybox=True, shadow=True, ncol=3, fontsize=12)

# Titolo o etichetta superiore (opzionale, basato sul tuo fileLabel)
plt.title(f"Label: {fileLabel_1}", fontsize=15)

plt.grid(True, which='both', linestyle=':', alpha=0.6)
plt.show()
