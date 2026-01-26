#!/bin/bash
echo "Energies to plot (space separated):"
read -a energies
echo "Number of pulses to plot (space separated):"
read -a pulses
echo "Generating new plot"

dirs=()
plotname="DVH_"
for E in "${energies[@]}"; do
plotname+="${E}_"
done
plotname+="MeV_"
for P in "${pulses[@]}"; do
plotname+="${P}_"
done
plotname+="pulses.png"

for E in "${energies[@]}"; do
for P in "${pulses[@]}"; do
dirs+=("sim${E}MeV/DVH${E}MeV_${P}pulses")
done
done

python3 starter_kit/GetDVHPlot.py "${dirs[@]}" --out "./${plotname}"
