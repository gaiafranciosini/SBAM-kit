#!/usr/bin/env python3

import os
import matplotlib
if 'FRED_MATPLOTLIB_BACKEND' in os.environ:
    matplotlib.use(os.environ['FRED_MATPLOTLIB_BACKEND'])
import matplotlib.pyplot as plt
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
parser.add_argument("-log", help="show logarithmic plot side by side", action="store_true")
parser.add_argument("-star", help="star position: multiplier_for_Dgoal Y_value (e.g., -star 0.95 95)", type=float, nargs=2, default=[1.0, 100.0])
args = parser.parse_args()

# Font size globale
FONT_SIZE = 15
plt.rcParams.update({'font.size': FONT_SIZE})

# Calcolo coordinate stella
star_x = args.Dgoal * args.star[0]
star_y = args.star[1]

# Variabili
fileLabel_1, fileLabel_2, fileLabel_3 = args.label1, args.label2, args.label3
fileDir_1, fileDir_2, fileDir_3 = args.dir1, args.dir2, args.dir3
my_Color = ['r', 'g', 'b', 'k', 'y', 'c', 'm', 'tab:orange', 'tab:purple', 'tab:pink', 'tab:brown', 'yellowgreen', 'slategray', 'hotpink', 'tan', 'lightcoral', 'tomato']

# Configurazione Figure
if args.log:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8))
    axes = [ax1, ax2]
else:
    fig, ax1 = plt.subplots(1, 1, figsize=(9, 8))
    axes = [ax1]

for ax in axes:
    ax.set_xlabel('Dose [cGy]', fontsize=FONT_SIZE)
    ax.set_ylabel('Volume [%]', fontsize=FONT_SIZE)
    if ax == axes[-1] and args.log:
        ax.set_yscale('log')
    
    for c, ROIname in enumerate(args.roi):
        color = my_Color[c % len(my_Color)]
        # File 1
        x1, y1 = np.loadtxt(f"{fileDir_1}/{ROIname}{fileLabel_1}.txt", comments=['#', 'idx'], usecols=(0, 1), unpack=True)
        ax.plot(x1, y1, color=color, label=ROIname, linewidth=2.0)
        # File 2
        x2, y2 = np.loadtxt(f"{fileDir_2}/{ROIname}{fileLabel_2}.txt", comments=['#', 'idx'], usecols=(0, 1), unpack=True)
        ax.plot(x2, y2, color=color, linestyle='dashed', linewidth=2.0)
        # File 3
        if fileLabel_3 and fileDir_3:
            x3, y3 = np.loadtxt(f"{fileDir_3}/{ROIname}{fileLabel_3}.txt", comments=['#', 'idx'], usecols=(0, 1), unpack=True)
            ax.plot(x3, y3, color=color, linestyle='dotted', linewidth=2.0)

    # Posizionamento stella personalizzata
    ax.plot([star_x], [star_y], 'r*', markersize=17, markerfacecolor='r', markeredgewidth=1, markeredgecolor='yellow')
    ax.grid(True, which='both', linestyle='--', alpha=0.5)

# Elementi legende
star_label = f"V{int(args.star[0]*100)}% {int(star_y)}%"
star_handle = Line2D([], [], marker='*', linewidth=0, color='r', markeredgewidth=1, markeredgecolor='yellow', markersize=17, label=star_label)
l1 = Line2D([0], [0], label=fileLabel_1, color='k', linestyle='solid')
l2 = Line2D([0], [0], label=fileLabel_2, color='k', linestyle='dashed')
handles_top = [star_handle, l1, l2]
if fileLabel_3:
    handles_top.append(Line2D([0], [0], label=fileLabel_3, color='k', linestyle='dotted'))

# Legenda Superiore
for ax in axes:
    ax.legend(handles=handles_top, loc='lower center', bbox_to_anchor=(0.5, 1.02), ncol=len(handles_top), fontsize=FONT_SIZE-2, frameon=False)

# Legenda Inferiore (ROIs)
handles_roi, labels_roi = ax1.get_legend_handles_labels()
if args.log:
    fig.legend(handles=handles_roi, labels=labels_roi, loc='lower center', bbox_to_anchor=(0.5, 0.05), ncol=4, fancybox=True, shadow=True, fontsize=FONT_SIZE)
    plt.subplots_adjust(top=0.85, bottom=0.22, wspace=0.3)
else:
    ax1.legend(handles=handles_roi, labels=labels_roi, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fancybox=True, shadow=True, fontsize=FONT_SIZE)
    plt.subplots_adjust(top=0.88, bottom=0.25)

plt.show()