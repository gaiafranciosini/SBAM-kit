import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-cpu", type=int, help="Number of CPU available", required=True)
#parser.add_argument("-A", type=int, help="Number of angles", required=True)
parser.add_argument("-E", type=int, help="Number of energies", required=True)
parser.add_argument("-P", type=float, help="Total number of primaries", required=True)

args = parser.parse_args()

cpu = args.cpu
E = args.E
P = args.P

inps_per_setup = cpu // E
primaries_per_inp = P / inps_per_setup

# Stampa in formato scientifico con massimo 3 decimali, E maiuscola
print("INPs: "+str(inps_per_setup))
print(f"PRIMARIES: {primaries_per_inp:.3E}")
print(f"CPUs: "+str(inps_per_setup*E))
