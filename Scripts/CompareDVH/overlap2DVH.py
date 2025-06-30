




#!/usr/bin/env python3

import os
import matplotlib
matplotlib.use(os.environ['FRED_MATPLOTLIB_BACKEND'])
import pylab as plt
import numpy as np
from matplotlib.lines import Line2D
import argparse

# Parser degli argomenti
parser = argparse.ArgumentParser(description='Presents label of the TXT files')
parser.add_argument("-label1", help="label of TXTs files 1", type=str, default="")
parser.add_argument("-label2", help="label of TXTs files 2", type=str, default="")
parser.add_argument("-label3", help="label of TXTs files 3 (optional)", type=str, default="")
parser.add_argument("-dir1", help="folder of TXTs files 1", type=str, default="")
parser.add_argument("-dir2", help="folder of TXTs files 2", type=str, default="")
parser.add_argument("-dir3", help="folder of TXTs files 3 (optional)", type=str, default="")
parser.add_argument("-Dgoal", help="reference planned dose [cGy]", type=float, default=100)
parser.add_argument("-roi", help="list of ROIs name", type=str, nargs='+')
args = parser.parse_args()

# Definizione variabili globali e colori
fileLabel_1 = args.label1
fileLabel_2 = args.label2
fileLabel_3 = args.label3  # Nuovo parametro per il terzo file
fileDir_1 = args.dir1
fileDir_2 = args.dir2
fileDir_3 = args.dir3  # Nuovo parametro per la cartella del terzo file
my_Color = ['r', 'g', 'b', 'k', 'y', 'c', 'm', 'tab:orange', 'tab:purple', 'tab:pink', 'tab:brown', 'yellowgreen', 'slategray', 'hotpink', 'tan', 'lightcoral', 'tomato']

# Crea figure con due subplot affiancati
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))

# Configurazione del primo grafico (scala lineare)
ax1.set_xlabel('Dose [cGy]', fontsize=11)
ax1.set_ylabel('Volume [%]', fontsize=11)

# Aggiungi il primo e secondo file
for c, ROIname in enumerate(args.roi):
    x1, y1 = np.loadtxt(fileDir_1 + '/' + ROIname + fileLabel_1 + '.txt', comments=['#', 'idx'], usecols=(0, 1), unpack=True)
    x2, y2 = np.loadtxt(fileDir_2 + '/' + ROIname + fileLabel_2 + '.txt', comments=['#', 'idx'], usecols=(0, 1), unpack=True)
    ax1.plot(x1, y1, color=my_Color[c], label=ROIname, linewidth=1.8)
    ax1.plot(x2, y2, color=my_Color[c], linestyle='dashed', linewidth=1.8)

    # Se esiste un terzo file, aggiungi il grafico
    if fileLabel_3 and fileDir_3:
        x3, y3 = np.loadtxt(fileDir_3 + '/' + ROIname + fileLabel_3 + '.txt', comments=['#', 'idx'], usecols=(0, 1), unpack=True)
        ax1.plot(x3, y3, color=my_Color[c], linestyle='dotted', linewidth=1.8)

# Configurazione del secondo grafico (scala logaritmica)
ax2.set_xlabel('Dose [cGy]', fontsize=11)
ax2.set_ylabel('Volume [%]', fontsize=11)
ax2.set_yscale('log')  # Imposta l'asse Y in scala logaritmica

# Aggiungi il primo e secondo file
for c, ROIname in enumerate(args.roi):
    x1, y1 = np.loadtxt(fileDir_1 + '/' + ROIname + fileLabel_1 + '.txt', comments=['#', 'idx'], usecols=(0, 1), unpack=True)
    x2, y2 = np.loadtxt(fileDir_2 + '/' + ROIname + fileLabel_2 + '.txt', comments=['#', 'idx'], usecols=(0, 1), unpack=True)
    ax2.plot(x1, y1, color=my_Color[c], label=ROIname, linewidth=1.8)
    ax2.plot(x2, y2, color=my_Color[c], linestyle='dashed', linewidth=1.8)

    # Se esiste un terzo file, aggiungi il grafico
    if fileLabel_3 and fileDir_3:
        x3, y3 = np.loadtxt(fileDir_3 + '/' + ROIname + fileLabel_3 + '.txt', comments=['#', 'idx'], usecols=(0, 1), unpack=True)
        ax2.plot(x3, y3, color=my_Color[c], linestyle='dotted', linewidth=1.8)

# Aggiungi stella su entrambi i grafici
for ax in [ax1, ax2]:
    ax.plot([args.Dgoal], [100], 'r*', markersize=17, markerfacecolor='r', markeredgewidth=1, markeredgecolor='yellow')
    star = Line2D([], [], marker='*', linewidth=0, color='r', markeredgewidth=1, markeredgecolor='yellow', markersize=17, label='V100% 100%')

    # Legenda in alto con fileLabel
    line1 = Line2D([0], [0], label=fileLabel_1, color='k', linestyle='solid')
    line2 = Line2D([0], [0], label=fileLabel_2, color='k', linestyle='dashed')
    if fileLabel_3:  # Se il terzo file è presente, aggiungi la legenda
        line3 = Line2D([0], [0], label=fileLabel_3, color='k', linestyle='dotted')
        legend_top = ax.legend(handles=[star, line1, line2, line3], loc='lower center', bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=4)
    else:
        legend_top = ax.legend(handles=[star, line1, line2], loc='lower center', bbox_to_anchor=(0, 1.02, 1, 0.2), ncol=3)
    
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
