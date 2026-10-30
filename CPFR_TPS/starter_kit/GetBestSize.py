#!/usr/bin/env python3
########################################################################
########################################################################
# CPFR TPS
# Jan 2026
########################################################################
########################################################################

import argparse

parser = argparse.ArgumentParser(description='find best field size basing on isodose curves in water')
parser.add_argument("-ptv", help="input map",required=True,nargs=None,metavar='path')
parser.add_argument("-lutpath", help="LUT PATH",required=True,nargs=None,metavar='path')
parser.add_argument("-min_size", help="X Y [cm]", type=float, required=True,nargs=2)
parser.add_argument("-E", help="energy",required=True)


args = parser.parse_args()
# exit(0)
########################################################################

import os, sys, shutil, time, re
from math import *
import numpy as np
import struct

ptv=args.ptv
LUT=args.lutpath
[minXsize, minYsize]=args.min_size
E=args.E

# --- MODIFICA 1: Utilizzo di FRED_mhd_io ---
from FRED_mhd_io import readMap
nn, hs, x0, Map, direction = readMap(ptv)
# -------------------------------------------

depth=np.sum(np.max(np.max(np.array(Map), axis=0), axis=0))*hs[2]
#depth = depth + 0.05  # Aggiunge un margine di sicurezza di 0.5 mm
if depth > 1.7:
    depth = 1.7
print("DEPTH [cm]= "+str(depth))
if (depth > 1.7):
  depth=1.7

from pathlib import Path

lut_dir = Path(LUT)

# --- MODIFICA 2: Estensione dei valori X e Y ---
# X_vals fino a 5 con passo 0.2. Uso round(1) per evitare problemi di floating point nei nomi dei file
X_vals = np.round(np.arange(2.0, 5.1, 0.2), 1) 
# Y_vals fino a 5 con passo 1
Y_vals = np.array([2, 3, 4, 5])
# -----------------------------------------------

lut = {}

Xidx=np.argmin(np.abs(X_vals-minXsize))
Yidx=np.argmin(np.abs(Y_vals-minYsize))
flag=False
for y in range(Yidx, len(Y_vals)):
  Y=int(Y_vals[y])
  for x in range(Xidx, len(X_vals)):
    X=X_vals[x]
    
    # Formattazione per i nomi dei file: se X è intero (es 2.0), lo formatto senza virgola
    X_str = int(X) if X.is_integer() else X
    
    fname = f"LUT_{E}MeV_{X_str}x{Y}cm2"
    fpath = lut_dir / fname

    if not fpath.exists():
      # file mancante lo saltiamo
      continue

    data = np.loadtxt(fpath)
#   print(data)
    if data.ndim != 2 or data.shape[1] != 5:
      raise ValueError(f"{fname} is not shaped Nx5")

    Zidx=np.argmin(np.abs(data[:,0]-depth))
    #print(f"DEPTH calcolata: {depth:.8f}")
    #print(f"Riga LUT selezionata (Zidx): {Zidx}")
    for z in range(Zidx, len(data[:,0])):
      if (data[z,2]-data[z,1]>=minXsize) and (data[z,4]-data[z,3]>=minYsize):
        print("Xsize: "+str(X_str))
        print("Ysize: "+str(Y))
#       print(data)
        flag=True
        break
    if flag==True:
      break
  if flag==True:
    break
if flag==False:
  print("Xsize: "+str(5))
  print("Ysize: "+str(5))
#      lut[(E, X, Y)] = data
