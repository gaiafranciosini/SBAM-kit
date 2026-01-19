#!/bin/bash

echo "Beam Energy (MeV):"
read -a energies
CHOICE="${energies[*]}"
if [[ "$CHOICE" != "7" && "$CHOICE" != "9" && "$CHOICE" != "7 9" && "$CHOICE" != "9 7" ]]; then
    echo 'Energy not available, choose "7" or "9" or "7 9"'
    exit 1
fi
echo "ENERGY: ${energies[@]}" > setup.out

echo "Number of pulses:"
read pulses
echo "PULSES: ${pulses}" >> setup.out

#echo "Prescription dose [Gy] and volume [%]:"
#"read -a preDV
#
#preD=${preDV[0]}
#preV=${preDV[1]}
#echo "PERCENTAGE PRESCRIPTION DOSE: ${preD}" >> setup.out
#echo "PERCENTAGE PRESCRIPTION VOLUME: ${preV}" >> setup.out

echo "How many primaries do you want to generate?"
read primaries
echo "PRIMARIES: ${primaries}" >> setup.out

echo "How many CPUs are available?"
read available_CPUs
echo "AVAILABLE CPUs: ${available_CPUs}" >> setup.out

num_energies=${#energies[@]}

if (( available_CPUs < num_energies )); then
    echo "ERROR: You selected $num_energies energies but only $available_CPUs CPUs are available."
    echo "Please provide at least as many CPUs as the number of energies."
    exit 1
fi

python3 starter_kit/eval_cpu.py -cpu "${available_CPUs}" -E "${#energies[@]}"  -P "${primaries}" > cpu_setup.out

INPs=$(awk -F': ' '/INPs/{print $2}' cpu_setup.out)
echo "INP FILES: ${INPs}" >> setup.out
primaries_per_CPU=$(awk -F': ' '/PRIMARIES/{print $2}' cpu_setup.out)
echo "PRIMARIES PER CPU: ${primaries_per_CPU}" >> setup.out
CPUs=$(awk -F': ' '/CPUs/{print $2}' cpu_setup.out)
echo "USED CPUs: ${CPUs}" >> setup.out

echo "You are using ${CPUs} CPUs, running ${primaries_per_CPU} per CPU"
echo "Each setup (beam energy+slits aperture) simulation is distributed over ${INPs} CPUs"
echo " "
echo "Press any key to proceed or CMD+Z to exit"
echo " "

bash starter_kit/card_modifier.sh starter_kit/EF70mm_start.inp starter_kit/tmp.inp "START" 1 "${primaries_per_CPU}" 1
mv starter_kit/tmp.inp starter_kit/EF70mm_start.inp

read wait

echo "Copy and customize: "
echo "nohup bash FULL_CHAIN_70mm.sh setup.out imgs/<<Phantom.mhd>> imgs/<<PTV.mhd>> imgs/<<MARKER.mhd>> imgs/<<OAR1.mhd>> imgs/<<OAR2.mhd>> ... >full_log.txt 2>&1 & tail -f full_log.txt"
