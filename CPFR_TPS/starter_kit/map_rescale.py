#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Jul 2020 - 2025 A.Schiavi
########################################################################
########################################################################
import sys
import numpy as np
from FRED_mhd_io import *
########################################################################
import argparse
parser = argparse.ArgumentParser(description='rescales the values of a list of maps\ne.g.: map1.mhd map2.mhd map3.mhd -factor 1.602176462E-7')
parser.add_argument("map", help="list of maps", nargs='+')

group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-factor","-multiplier", help="multiplies voxel values by a factor",type=float,metavar='#')
group.add_argument("-divider", help="divides voxel values by the divider",type=float,metavar='#')

parser.add_argument("-v", "--verbose", help="increase output verbosity",action="store_true",dest='verbose')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)  # Exit with a non-zero status for no arguments

args = parser.parse_args()
# print(args)
# exit(0)
########################################################################


if args.divider!=None and args.divider==0:
    print('error: divider is zero => this would produce just NANs') 
    exit(1)

for i,f in enumerate(args.map):

    nn,hs,x0,Map,direction = readMap(f)
    if args.factor != None:
        Map*=args.factor
    if args.divider != None:
        Map/=args.divider

    writeMap(f,hs,x0,Map,direction)
