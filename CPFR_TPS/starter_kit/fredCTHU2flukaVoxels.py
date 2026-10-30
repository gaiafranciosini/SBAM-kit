#!/usr/bin/env python3
########################################################################
########################################################################
# Fred project
# Jan 2021 A.Schiavi
########################################################################
########################################################################
import argparse

parser = argparse.ArgumentParser(
    description="Convert Fred CT and HU materials into Fluka voxels file"
)
parser.add_argument(
    "-v", "--verbose", help="increase output verbosity", action="store_true"
)
parser.add_argument("CT", help="CT map in HU")
parser.add_argument("table", help="dense HU material table")
args = parser.parse_args()

import math
import os
import re
import shutil
import struct
import sys
import numpy as np

# Sostituita la vecchia libreria mhd_io con FRED_mhd_io
from FRED_mhd_io import readMap

lverbose = args.verbose


def writeFortranRecord(file, record):
    file.write(struct.pack("<1i", len(record)))
    file.write(record)
    file.write(struct.pack("<1i", len(record)))


try:
    if lverbose:
        print("Loading CT map from file ", args.CT)
    # Utilizzo di readMap da FRED_mhd_io
    nn, hs, x0, Map, direction = readMap(args.CT)
    Map = Map.astype(np.short)  # convert to short

    if lverbose:
        print("Loading fred HU table from file ", args.table)
    table = np.loadtxt(args.table, skiprows=1).astype(np.float32)
except Exception as e:
    print(f"Fatal error: {e}")
    exit(1)

HUmin = Map.min()
HUmax = Map.max()
HUactive = np.unique(Map)
numHUmat = len(HUactive)
print("HU range = ", HUmin, HUmax)
print("HU active = ", HUactive, " = ", len(HUactive), "elements")

HUdense = np.zeros(HUmax - HUmin + 1, dtype=np.short)
HUdense[HUactive - HUmin] = 1
print("HUdense = ", HUdense, HUdense.sum())

# build HU to VOXEL### link table
KREG = HUdense.copy()

MO = KREG.size - 1
NO = KREG.sum() - 1
print("MO,NO", MO, NO)

ireg = -1
for i in range(len(KREG)):
    if KREG[i] != 0:
        ireg += 1
        KREG[i] = ireg

if lverbose:
    print("table shape", table.shape)

fin = open(args.table)
line = fin.readline()
cols = line.split()[1:]
if lverbose:
    print("table header columns", cols)

if len(cols) != table.shape[1]:
    print("error: columns numbers do not match")
    exit(1)

iHU = cols.index("HU")
irho = cols.index("rho")
iRSP = cols.index("RSP")
iLrad = cols.index("Lrad")

ielements = list(range(len(cols)))

ielements.remove(iHU)
ielements.remove(irho)
ielements.remove(iRSP)
ielements.remove(iLrad)
if lverbose:
    print(ielements)

fredElements = [
    "H",
    "C",
    "N",
    "O",
    "Na",
    "Mg",
    "P",
    "S",
    "Cl",
    "K",
    "Ca",
    "Fe",
    "I",
]
flukaElements = [
    "HYDROGEN",
    "CARBON",
    "NITROGEN",
    "OXYGEN",
    "SODIUM",
    "MAGNESIU",
    "PHOSPHO",
    "SULFUR",
    "CHLORINE",
    "POTASSIU",
    "CALCIUM",
    "IRON",
    "IODINE",
]

iH = cols.index("H")
iC = cols.index("C")
iN = cols.index("N")
iO = cols.index("O")
iNa = cols.index("Na")
iMg = cols.index("Mg")
iP = cols.index("P")
iS = cols.index("S")
iCl = cols.index("Cl")
iK = cols.index("K")
iCa = cols.index("Ca")
iFe = cols.index("Fe")
iI = cols.index("I")

table = np.c_[
    table, np.zeros(table.shape[0])
]  # add extra column for isActive
iIsActive = table.shape[1] - 1

print(table[:, iIsActive].shape, HUdense.shape)
print(table[:, iIsActive])

for i in range(table.shape[0]):
    hunow = int(table[i, iHU])
    ihu = hunow - HUmin
    isActive = hunow >= HUmin and hunow <= HUmax and HUdense[ihu] > 0
    table[i, iIsActive] = 1 if isActive else 0

# build reduced FLUKA materials
tableHUmin = int(table[0][iHU])
tableHUmax = int(table[-1][iHU])

iHUrangemin = table.shape[1]
iHUrangemax = table.shape[1] + 1
iCorrFact1 = table.shape[1] + 2
iCorrFact2 = table.shape[1] + 3

flukaRosettaStone = np.zeros(
    (table.shape[0], 4 + table.shape[1]), dtype=np.float32
)

if numHUmat < 200:  # build one Fluka material for each HU material
    for i in range(table.shape[0]):
        if table[i, iIsActive] > 0:
            flukaRosettaStone[i, 0:iHUrangemin] = table[
                i, 0:iHUrangemin
            ]  # copy fred material
            flukaRosettaStone[i, iHUrangemin] = table[i][iHU]
            flukaRosettaStone[i, iHUrangemax] = table[i][iHU]
            flukaRosettaStone[i, iCorrFact1] = (
                table[i, iRSP] / table[i][irho]
            )
            flukaRosettaStone[i, iCorrFact2] = 1.0
else:  # build HU ranges and use corrfact
    fluHUmat = [0]
    for i in range(1, table.shape[0] - 1):
        discontinuityFound = False
        for j in ielements:
            DerSx = table[i, j] - table[i - 1, j]
            DerDx = table[i + 1, j] - table[i, j]
            if DerSx == 0.0 and DerDx != 0.0:
                discontinuityFound = True
            elif DerSx != 0.0 and DerDx == 0.0:
                discontinuityFound = True
            elif DerDx != 0.0 and abs((DerDx - DerSx) / DerDx) > 0.20:
                discontinuityFound = True

            if discontinuityFound:
                break
        if discontinuityFound:
            fluHUmat.append(i)
    fluHUmat.append(table.shape[0] - 1)

    for i in range(table.shape[0]):
        if i in fluHUmat:
            print("fluHUmat", i, table[i, iHU])
    print("not implemented yet!!!! ;)")

    ibeg = 0
    inext = table.shape[0]
    HUstep = 50
    while ibeg < table.shape[0]:
        if ibeg in fluHUmat:
            print("adding ", ibeg)
            if table[ibeg, iIsActive] > 0:
                flukaRosettaStone[ibeg, 0:iHUrangemin] = table[
                    ibeg, 0:iHUrangemin
                ]
                flukaRosettaStone[ibeg, iHUrangemin] = table[ibeg][iHU]
                flukaRosettaStone[ibeg, iHUrangemax] = table[ibeg][iHU]
                flukaRosettaStone[ibeg, iCorrFact1] = (
                    table[ibeg, iRSP] / table[ibeg][irho]
                )
                flukaRosettaStone[ibeg, iCorrFact2] = 1.0
            idx = fluHUmat.index(ibeg) + 1
            if idx < len(fluHUmat):
                print(ibeg, idx)
                inext = fluHUmat[idx]
            else:
                inext = table.shape[0]
            ibeg += 1
            continue

        iend = min(ibeg + HUstep, table.shape[0], inext)
        print("try= ", ibeg, iend, "[", inext, "]")
        meanHUMat = table[ibeg:iend][:].mean(axis=0)
        while (
            table[ibeg, irho] < 2.0 / 3.0 * meanHUMat[irho]
            or table[iend - 1, irho] > 3.0 / 2.0 * meanHUMat[irho]
            or table[ibeg, iRSP] < 2.0 / 3.0 * meanHUMat[irho]
            or table[iend - 1, iRSP] > 3.0 / 2.0 * meanHUMat[irho]
        ):
            iend -= 1
            meanHUMat = table[ibeg:iend][:].mean(axis=0)
        print(
            "got= ",
            ibeg,
            iend,
            "|",
            table[ibeg, irho] / meanHUMat[irho],
            table[iend - 1, irho] / meanHUMat[irho],
        )

        for i in range(ibeg, iend):
            flukaRosettaStone[i, 0:iHUrangemin] = meanHUMat
            flukaRosettaStone[i, iHUrangemin] = table[ibeg][iHU]
            flukaRosettaStone[i, iHUrangemax] = table[iend][iHU]
            flukaRosettaStone[i, iCorrFact1] = (
                table[i, iRSP] / table[i, irho]
            )
            flukaRosettaStone[i, iCorrFact2] = (
                table[i, irho] / meanHUMat[irho]
            )

        ibeg = iend

for i in range(table.shape[0]):
    if table[i, iIsActive] > 0:
        print("flukamat = ", i, *flukaRosettaStone[i, :], end="")
        print("")

matLines = []
matLines.append(
    "MATERIAL         16.    32.066       2.0                              SULFUR    "
)
matLines.append(
    "MATERIAL         15. 30.973761       2.2                              PHOSPHO   "
)
matLines.append(
    "MATERIAL         17.   35.4527 0.0029947                              CHLORINE  "
)
matLines.append(
    "MATERIAL         19.   39.0983     0.862                              POTASSIU  "
)
matLines.append(
    "MATERIAL         53. 126.90447      4.93                              IODINE    "
)

lFirst = True

for i in range(flukaRosettaStone.shape[0]):
    if flukaRosettaStone[i, iIsActive] == 0.0:
        continue  # not active => skip it
    if lFirst:
        lFirst = False
    else:
        if flukaRosettaStone[i, iHUrangemin] == icurr:
            continue  # skip it
    icurr = flukaRosettaStone[i, iHUrangemin]
    print(*flukaRosettaStone[i, :])

    isHUrange = (
        flukaRosettaStone[i, iHUrangemin]
        != flukaRosettaStone[i, iHUrangemax]
    )

    if isHUrange:
        name = "HU<%d" % int(flukaRosettaStone[i, iHUrangemax])
    else:
        name = "HU%d" % int(flukaRosettaStone[i, iHUrangemin])

    rho = flukaRosettaStone[i][irho]
    Lrad = flukaRosettaStone[i][iLrad]
    RSP = flukaRosettaStone[i][iRSP]
    composition = []
    weigths = []
    for j in ielements:
        if flukaRosettaStone[i][j] > 0:
            composition.append(cols[j])
            weigths.append(flukaRosettaStone[i][j])

    if lverbose:
        print("*", i, name, rho, Lrad, composition, weigths)

    # build fluka water for dEdx!!!
    
    # MATERIAL command
    line = "%-10s" % "MATERIAL"
    line += "%10s%10s" % ("", "")
    line += "%10f" % rho
    line += "%10s%10s" % ("", "")
    line += "%10s" % "0."
    line += "%-10s" % ("W" + name)
    matLines.append(line)

    # COMPOUND command
    line = "%-10s" % "COMPOUND"
    line += "%10s" % "2."
    line += "%10s" % "HYDROGEN"
    line += "%10s" % "1."
    line += "%10s" % "OXYGEN"
    line += "%10s%10s" % ("", "")
    line += "%-10s" % ("W" + name)
    matLines.append(line)

    # STERNHEI command
    line = "%-10s" % "STERNHEI"
    line += "%10s" % "3.5017"
    line += "%10s" % "0.24"
    line += "%10s" % "2.8004"
    line += "%10s" % "0.09116"
    line += "%10s" % "3.4773"
    line += "%10s" % ""
    line += "%-10s" % ("W" + name)
    matLines.append(line)

    # build fluka material
    line = "%-10s" % ("MATERIAL")
    line += "%10s%10s" % ("", "")  # what(1) what(2)
    line += "%10f" % (rho)  # what(3) density
    line += "%10s" % ("")  # what(4)
    line += "%10s" % ("W" + name)  # what(5) ref material for dEdx
    line += "%10s" % ("0.")  # what(6) 0 = natural isotopic composition
    line += "%-10s" % (name)
    matLines.append(line)

    line = "%-10s" % ("COMPOUND")
    for iel in range(3):
        if iel < len(composition):
            line += "%10f%10s" % (
                -weigths[iel] / 100.0,
                flukaElements[fredElements.index(composition[iel])],
            )
        else:
            line += "%10s%10s" % ("", "")
    line += "%-10s" % (name)
    matLines.append(line)
    if len(composition) > 3:
        line = "%-10s" % ("COMPOUND")
        for iel in range(3, 6):
            if iel < len(composition):
                line += "%10f%10s" % (
                    -weigths[iel] / 100.0,
                    flukaElements[fredElements.index(composition[iel])],
                )
            else:
                line += "%10s%10s" % ("", "")
        line += "%-10s" % (name)
        matLines.append(line)
    if len(composition) > 6:
        line = "%-10s" % ("COMPOUND")
        for iel in range(6, 9):
            if iel < len(composition):
                line += "%10f%10s" % (
                    -weigths[iel] / 100.0,
                    flukaElements[fredElements.index(composition[iel])],
                )
            else:
                line += "%10s%10s" % ("", "")
        line += "%-10s" % (name)
        matLines.append(line)
    if len(composition) > 9:
        line = "%-10s" % ("COMPOUND")
        for iel in range(9, 12):
            if iel < len(composition):
                line += "%10f%10s" % (
                    -weigths[iel] / 100.0,
                    flukaElements[fredElements.index(composition[iel])],
                )
            else:
                line += "%10s%10s" % ("", "")
        line += "%-10s" % (name)
        matLines.append(line)

assignLines = []
corrFactLines = []

assignLines.append(
    "ASSIGNMA     VACUUM      VOXEL                                                                  "
)
for i in range(table.shape[0]):
    ikreg = int(table[i][iHU]) - HUmin
    if ikreg < 0 or ikreg >= len(KREG):
        continue
    ireg = KREG[ikreg]
    if ikreg > 0 and ireg == 0:
        continue  # skip this value => not present

    isHUrange = (
        flukaRosettaStone[i, iHUrangemin]
        != flukaRosettaStone[i, iHUrangemax]
    )

    if isHUrange:
        name = "HU<%d" % int(flukaRosettaStone[i, iHUrangemax])
    else:
        name = "HU%d" % int(flukaRosettaStone[i, iHUrangemin])

    line = "%-10s" % ("CORRFACT")
    line += "%10f" % -flukaRosettaStone[i, iCorrFact1]
    line += "%10f" % flukaRosettaStone[i, iCorrFact2]
    line += "%10s" % ("")
    voxelno = (
        ("VOXEL%03d" % (ireg + 1))
        if (ireg + 1) < 1000
        else ("VOXE%04d" % (ireg + 1))
    )
    line += "%10s" % voxelno
    line += " " * 30
    corrFactLines.append(line)

    line = "%-10s" % ("ASSIGNMA")
    line += "%10s" % name
    line += "%10s" % voxelno
    line += " " * 50
    assignLines.append(line)

# write voxels file
fout = open("CTHU.vxl", "wb")
title = "%-80s" % "CT from Fred"
record = struct.pack("<80s", str.encode(title))
writeFortranRecord(fout, record)

record = struct.pack("<5i", *nn, NO, MO)
writeFortranRecord(fout, record)
record = struct.pack("<3d", *hs)
writeFortranRecord(fout, record)
Map = Map - HUmin
writeFortranRecord(fout, Map.tobytes(order="F"))
writeFortranRecord(fout, KREG[1:].tobytes(order="F"))

for line in matLines:
    writeFortranRecord(fout, str.encode(line))
    print(line)
for i, line in enumerate(assignLines):
    writeFortranRecord(fout, str.encode(line))
    print(line)
for line in corrFactLines:
    writeFortranRecord(fout, str.encode(line))
    print(line)
