import argparse 
import numpy as np
import os

parser = argparse.ArgumentParser(description="Evaluates DVH metrics")
parser.add_argument("-DVH", type=str, nargs='+', required=True, help="DVH file")
#parser.add_argument("-", type=float, required=True, help="Valore B")
args = parser.parse_args()

vol=np.array([95, 90, 80, 70, 50, 40, 30, 25, 20, 10])

files=args.DVH
for dvh in files:
  roi = os.path.splitext(os.path.basename(dvh))[0]
  print(roi)
  dose, volume= np.loadtxt(dvh, unpack=True)
  for v in vol:
    d_vol=0.01*dose[np.argmin(np.abs(volume-v))]
    print("	D(V"+str(v)+") = "+str(f"{d_vol:.3f}")+" Gy")  
  dpre=(1/95)*dose[np.argmin(np.abs(volume-95))]
  print("        D_pre = "+str(f"{dpre:.3f}")+" Gy")
  dmax=0.01*dose[np.argmin(np.abs(volume-0.01))]
  print("	D_max = "+str(f"{dmax:.3f}")+" Gy (V<0.01%)")
  perc=100*(dmax/dpre)
  print("	Dmax is "+str(f"{perc:.3f}")+"% of D_pre")
  print("")
  print("###########")
  print("")

