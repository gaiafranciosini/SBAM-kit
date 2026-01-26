#!/bin/bash
#++++++++++++++++++++++++++++++++++
echo "New number of pulses:"
read new_pulses
input="input_log.out"
read -r file CT PTV MARKER rest < "$input"
read -r _ roi_line < <(grep '^ROI:' "$input")
ROIs=(${roi_line#ROI: })

# Inizializza variabili
#ape=0
energies=()
pulses=""
#preD=""
#preV=""
primaries=""
available_CPUs=""
INPs=""
primaries_per_CPU=""
CPUs=""

# Legge il file riga per riga
while IFS= read -r line; do
    case "$line" in
        "ENERGY:"*)
            energies=(${line#ENERGY: })
            ;;
        "PULSES:"*)
            pulses=(${line#PULSES: })
            ;;
#        "PERCENTAGE PRESCRIPTION DOSE:"*)
#            preD=${line#PERCENTAGE PRESCRIPTION DOSE: }
#            ;;
#        "PERCENTAGE PRESCRIPTION VOLUME:"*)
#            preV=${line#PERCENTAGE PRESCRIPTION VOLUME: }
#            ;;
        "PRIMARIES:"*)
            primaries=${line#PRIMARIES: }
            ;;
        "AVAILABLE CPUs:"*)
            available_CPUs=${line#AVAILABLE CPUs: }
            ;;
        "INP FILES:"*)
            INPs=${line#INP FILES: }
            ;;
        "PRIMARIES PER CPU:"*)
            primaries_per_CPU=${line#PRIMARIES PER CPU: }
            ;;
        "USED CPUs:"*)
            CPUs=${line#USED CPUs: }
            ;;
    esac
done < "$file"


#CALIBRATION FACTORS (primaries/pulse)
kFLASH_7MeV=9.382238805934787e+11
kCONV_7MeV=1.42155133423254e+10
kFLASH_9MeV=7.106140035174675e+11
kCONV_9MeV=1.07668788411738e+10


all_rois_path=()
ROIs_plan=()
for roi in "${ROIs[@]}"; do
all_rois_path+=("../imgs/${roi}_plan.mhd")
ROIs_plan+=("${roi}_plan")
done 

echo "ORIGINAL # of PULSES: ${pulses}"
echo "NEW # of PULSES: ${new_pulses}"

pids=()
for E in "${energies[@]}"; do
  echo "Rescaling ${E}MeV test"
  cd sim${E}MeV
  [[ -d DVH${E}MeV ]] && mv DVH${E}MeV DVH${E}MeV_${pulses}pulses
  [[ -f ${E}MeV.png ]] && mv ${E}MeV.png ${E}MeV_${pulses}pulses.png
  [[ -f DOSE_${E}MeV.mhd ]] && cp DOSE_${E}MeV.mhd DOSE_${E}MeV_${pulses}pulses.mhd
  mkdir "DVH${E}MeV_${new_pulses}pulses"
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -divider "${pulses}"
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -multiplier "${new_pulses}"
  echo "Map rescaled"
  ../starter_kit/ComputeDVH/ComputeDVH.x -Dgoal 2000 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV.mhd -type float -fileLabel ${E}MeV -dir DVH${E}MeV_${new_pulses}pulses > trash.out
  echo "New DVH curves created"
  python3 ../starter_kit/plotDVH.py -label1 ${E}MeV -dir1 DVH${E}MeV_${new_pulses}pulses -roi PTV_plan "${ROIs_plan[@]}" > trash.out
  python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV_${new_pulses}pulses/PTV_plan${E}MeV.txt > DVH${E}MeV_${new_pulses}pulses/DVH_PTV.out
#  mv ${E}MeV.png ${E}MeV_${new_pulses}pulses.png
  echo "DVH plot saved"
  pids+=("$!")
  cd ..
  echo " "
done 

ex_pulses=${pulses}
pulses=(${ex_pulses} ${new_pulses})

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
echo "${plotname} generated"
echo "To compare DVHs with other energies or pulses values use:" 
echo "bash newDVHplot.sh"










