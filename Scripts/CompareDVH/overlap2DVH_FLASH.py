#!/usr/bin/env python3

import os
import matplotlib
matplotlib.use(os.environ.get('FRED_MATPLOTLIB_BACKEND', 'Agg'))
import pylab as plt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
import argparse

# Parser degli argomenti
parser = argparse.ArgumentParser(description='Presents label of the TXT files')
parser.add_argument("-label1", help="label of TXTs files 1", type=str, default="")
parser.add_argument("-label2", help="label of TXTs files 2", type=str, default="")
parser.add_argument("-label3", help="label of TXTs files 3", type=str, default="")  # Nuovo parametro label3
parser.add_argument("-dir1", help="folder of TXTs files 1", type=str, default="")
parser.add_argument("-dir2", help="folder of TXTs files 2", type=str, default="")
parser.add_argument("-dir3", help="folder of TXTs files 3", type=str, default="")  # Nuovo parametro dir3
parser.add_argument("-Dgoal", help="reference planned dose [cGy]", type=float, default=100)
parser.add_argument("-roi", help="list of ROIs name", type=str, nargs='+')
args = parser.parse_args()

# Definizione variabili globali e colori
fileLabel_1 = args.label1
fileLabel_2 = args.label2
fileLabel_3 = args.label3  # Nuova variabile per label3
fileDir_1 = args.dir1
fileDir_2 = args.dir2
fileDir_3 = args.dir3  # Nuova variabile per dir3
my_Color = ['r', 'g', 'b', 'k', 'y', 'c', 'm', 'tab:orange', 'tab:purple', 'tab:pink', 
            'tab:brown', 'yellowgreen', 'slategray', 'hotpink', 'tan', 'lightcoral', 'tomato']

# Crea figure con due subplot affiancati
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))

# Configurazione del primo grafico (scala lineare)
ax1.set_xlabel('Dose [cGy]', fontsize=11)
ax1.set_ylabel('Volume [%]', fontsize=11)

# Configurazione del secondo grafico (zoom su Y da 0 a 20)
ax2.set_xlabel('Dose [cGy]', fontsize=11)
ax2.set_ylabel('Volume [%]', fontsize=11)
ax2.set_ylim(0, 10)  # Limiti asse Y da 0 a 20
ax2.set_xlim(0, 6000)  # Limiti asse Y da 0 a 20

for c, ROIname in enumerate(args.roi):
    # Caricamento dei dati per i primi due file
    x1, y1 = np.loadtxt(f"{fileDir_1}/{ROIname}{fileLabel_1}.txt", comments=['#', 'idx'], usecols=(0, 1), unpack=True)
    x2, y2 = np.loadtxt(f"{fileDir_2}/{ROIname}{fileLabel_2}.txt", comments=['#', 'idx'], usecols=(0, 1), unpack=True)
    
    # Caricamento dei dati per il terzo file (nuovo)
    x3, y3 = np.loadtxt(f"{fileDir_3}/{ROIname}{fileLabel_3}.txt", comments=['#', 'idx'], usecols=(0, 1), unpack=True)

    # Interpolazione per allineare x2 e x3 a x1
    f2 = interp1d(x2, y2, bounds_error=False, fill_value=0)
    y2_interpolated = f2(x1)
    
    f3 = interp1d(x3, y3, bounds_error=False, fill_value=0)
    y3_interpolated = f3(x1)

    # Linea continua per il primo file
    ax1.plot(x1, y1, color=my_Color[c], linestyle = 'solid', label=ROIname, linewidth=1.8)
    # Linea tratteggiata per il secondo file
    ax1.plot(x1, y2_interpolated, color=my_Color[c], linestyle='solid', linewidth=1.8,alpha=0.)
    # Banda trasparente tra label1 e label2
    ax1.fill_between(x1, y1, y2_interpolated, color=my_Color[c], alpha=0.2)  
    
    # Linea tratteggiata separata per il terzo file (label3)
    ax1.plot(x1, y3_interpolated, color=my_Color[c], linestyle='dashed', linewidth=1.8)


    # Linea continua per il primo file
    ax2.plot(x1, y1, color=my_Color[c], linestyle = 'solid', label=ROIname, linewidth=1.8)
    # Linea tratteggiata per il secondo file
    ax2.plot(x1, y2_interpolated, color=my_Color[c], linestyle='solid', linewidth=1.8, alpha=0.)
    # Banda trasparente tra label1 e label2
    ax2.fill_between(x1, y1, y2_interpolated, color=my_Color[c], alpha=0.2)  
    
    # Linea tratteggiata separata per il terzo file (label3)
    ax2.plot(x1, y3_interpolated, color=my_Color[c], linestyle='dashed', linewidth=1.8)

# Aggiungi stella su entrambi i grafici
for ax in [ax1, ax2]:
    ax.plot([args.Dgoal], [100], 'r*', markersize=17, markerfacecolor='r', markeredgewidth=1, markeredgecolor='yellow')
    star = Line2D([], [], marker='*', linewidth=0, color='r', markeredgewidth=1, markeredgecolor='yellow', markersize=17, label='V100% 100%')

    # Legenda in alto con fileLabel
    line1 = Line2D([0], [0], label=fileLabel_1, color='k', linestyle='solid')
#    line2 = Line2D([0], [0], label=fileLabel_2, color='k', linestyle='dotted')
    line3 = Line2D([0], [0], label=fileLabel_3, color='k', linestyle='dashed')  # Nuova legenda per il file 3
    legend_top = ax.legend(handles=[star, line1, line3], loc='lower center', bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=4)
    ax.add_artist(legend_top)

# Legenda ROIname centrata sotto entrambi i grafici
handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles=handles, labels=labels, loc='lower center', bbox_to_anchor=(0.5, -0.001), ncol=3, fancybox=True, shadow=True)

# Aggiungi griglia e layout
ax1.grid()
ax2.grid()
plt.tight_layout()
plt.subplots_adjust(top=0.85, bottom=0.2)  # Spazio per le legende sopra e sotto
plt.show()


