#!/bin/bash

#read -p "CT path:  " CT
#read -p "PTV path:  " PTV
#read -p "MARKER path:  " MARKER
#
#echo "ROIs paths:  " 
#read -a ROIpaths
#
#for path in "${ROIpaths[@]}"; do
#  # Rimuove "imgs/" e ".mhd"
#  name="${path##*/}"         # rimuove il percorso → osso.mhd
#  name="${name%.mhd}"        # rimuove l'estensione → osso
#  ROIs+=("$name")     # aggiunge all'array
#done
#
#+++++++++++++++++++++++++++++++++


wait_with_spinner_and_report() {
  local pid="$1"
  local spin='|/-\'
  local i=0
  local msg="${SPINNER_MSG:-Working}"

  # Spinner mentre il processo è vivo
  while kill -0 "$pid" 2>/dev/null; do
    i=$(( (i+1) % 4 ))
    printf "\r%s %s" "$msg" "${spin:$i:1}"
    sleep 0.2
  done

  # Linea finale "Done!" come nel tuo snippet
  printf "\rDone!   \n"

  # Report esito come nel tuo codice
  local fail=0
  if ! wait "$pid"; then
    ((fail++))
    echo "Process PID $pid terminated with an error"
  fi

  if (( fail == 0 )); then
    echo "All simulations terminated successfully"
  else
    echo "$fail process(es) terminated with an error"
  fi

  # Ritorna 0/1 in base all'esito
  (( fail == 0 ))

#USAGE:
#( sleep 5 ) &
#pid=$!
#SPINNER_MSG="Simulazione in corso"
#wait_with_spinner_and_report "$pid"



}

#++++++++++++++++++++++++++++++++++


CT="$1"
PTV="$2"
MARKER="$3"

if [ "$#" -lt 3 ]; then
  echo "Use: CT PTV MARKER [ROI1 ROI2 ROI3 ROI4 ...]"
  exit 1
fi

shift 3

# Array che conterrà i nomi puliti
ROIs=()

# Ciclo sugli argomenti rimanenti
for path in "$@"; do
  # Rimuove "imgs/" e ".mhd"
  name="${path##*/}"         # rimuove il percorso → osso.mhd
  name="${name%.mhd}"        # rimuove l'estensione → osso
  ROIs+=("$name")     # aggiunge all'array
done
#++++++++++++++++++++++++++++++++++

echo "Shaper orientation angles (degrees):"
read -a angles

echo "Beam Energy (MeV):"
read -a energies
CHOICE="${energies[*]}"
if [[ "$CHOICE" != "7" && "$CHOICE" != "9" && "$CHOICE" != "7 9" && "$CHOICE" != "9 7" ]]; then
    echo 'Energy not available, choose "7" or "9" or "7 9"'
    exit 1
fi

echo "Prescription dose [Gy] and volume [%]:"
read -a preDV

preD=${preDV[0]}
preV=${preDV[1]}

echo "How many primaries do you want to generate?"
read primaries

echo "How many CPUs are available?"
read available_CPUs

python3 starter_kit/eval_cpu.py -cpu "${available_CPUs}"  -A "${#angles[@]}"  -E "${#energies[@]}"  -P "${primaries}" > cpu_setup.out

INPs=$(awk -F': ' '/INPs/{print $2}' cpu_setup.out)
primaries_per_CPU=$(awk -F': ' '/PRIMARIES/{print $2}' cpu_setup.out)
CPUs=$(awk -F': ' '/CPUs/{print $2}' cpu_setup.out)

echo "You are using ${CPUs} CPUs, running ${primaries_per_CPU} per CPU"
echo "Each setup (beam energy+shaper angle) simulation is distributed over ${INPs} CPUs"
echo " "
echo "Press any key to proceed or CMD+Z to exit"
echo " "

bash starter_kit/card_modifier.sh starter_kit/EF70mm_start.inp starter_kit/tmp.inp "START" 1 "${primaries_per_CPU}" 1
mv starter_kit/tmp.inp starter_kit/EF70mm_start.inp

read wait

echo "Good luck!"


#1. Get BEAM DIRECTION

python3 starter_kit/GetDirection.py "$CT" -PTV "$PTV" -marker "$MARKER"  > out.out &
echo " " &
pid=$!
SPINNER_MSG="Computing beam direction"
wait_with_spinner_and_report "$pid"


output=$(grep "BEAM DIRECTION:" out.out | head -n 1)

read -r Vx Vy Vz < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')

echo "Beam direction: Vx:${Vx} Vy:${Vy} Vz:${Vz}"

echo "Rotating PTV, ROIs and CT"
python3 starter_kit/RotROIsQuat.py "$PTV" "$Vx" "$Vy" "$Vz" imgs/PTV_ROT.mhd > ptv.out
echo "PTV rotated and saved as imgs/PTV_ROT.mhd"
for roi in "${ROIs[@]}"; do
  python3 starter_kit/RotROIsQuat.py "imgs/${roi}.mhd" "$Vx" "$Vy" "$Vz" imgs/${roi}_ROT.mhd > trash.out
  echo "${roi} rotated and saved as imgs/${roi}_ROT.mhd"
done
python3 starter_kit/RotMapsQuat.py "$CT" "$Vx" "$Vy" "$Vz" imgs/CT_ROT.mhd > trash.out
echo "CT rotated and saved as imgs/CT_ROT.mhd"

echo "Shifting ROIs and CT to center the PTV in [0,0,z]"

output=$(grep -F "ISO_PTV[voxel]:" ptv.out | head -n 1)
read -r idxPTVx idxPTVy idxPTVz < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')
echo "X_ptv[idx]: $idxPTVx | Y_ptv[idx]: $idxPTVy | Z_ptv[idx]: $idxPTVz"

output=$(grep -F "ISO_PTV[cm]:" ptv.out | head -n 1)
read -r PTVx PTVy PTVz < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')
echo "X_ptv: $PTVx | Y_ptv: $PTVy | Z_ptv: $PTVz"

SHx=$(echo "$PTVx * -1" | bc -l)
SHy=$(echo "$PTVy * -1" | bc -l)

python3 starter_kit/mhd_shift.py "imgs/CT_ROT.mhd" "${SHx}" "${SHy}" 0 imgs/CT_SH.mhd > trash.out
echo "CT shifted  and saved as imgs/CT_SH.mhd"
python3 starter_kit/mhd_shift.py "imgs/PTV_ROT.mhd" "${SHx}" "${SHy}" 0 imgs/PTV_SH.mhd > trash.out
echo "PTV_shifted  and saved as imgs/PTV_SH.mhd"

for roi in "${ROIs[@]}"; do
  python3 starter_kit/mhd_shift.py "imgs/${roi}_ROT.mhd" "${SHx}" "${SHy}" 0  imgs/${roi}_SH.mhd > trash.out
  echo "${roi} shifted and saved as imgs/${roi}_SH.mhd"
done
echo

# 3. FIELD SIZE dimensions (modify slit size with applicator) !!!!!!!!!!!!!!!!!!!!!!!
echo "Optimising field size for each shaper orientation"

cp "starter_kit/EF70mm_start.inp" "starter_kit/EF70mm.inp"
echo "inp copied"
width=8.0
height=4.0
echo "Slit size for 70mm applicator: ${width}x${height}cm^2"

python3 starter_kit/GetFieldSize.py imgs/PTV_SH.mhd -rot "${angles[@]}" > rectangles.out &

pid=$!
SPINNER_MSG="Optimising geometry"
wait_with_spinner_and_report "$pid"

echo "Geometry optimization details reported in rectangles.out"
echo " "
echo "Cropping PTV, CT and ROIs"
python3 starter_kit/mhd_info.py imgs/CT_SH.mhd > info.out
output=$(grep "dims=" info.out | head -n 1)

read -r dimX dimY dimZ < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')

cropZ=$(( idxPTVz - 1 )) #to make sure that the first slice of the PTV is not cut out
python3 starter_kit/mhd_crop.py \-i imgs/CT_SH.mhd \-o imgs/CT_CROP.mhd \-ixi 0  \-ixf "$dimX" \-iyi 0 \-iyf "$dimY" \-izi "$cropZ" \-izf "$dimZ" > crop.out
python3 starter_kit/mhd_crop.py \-i imgs/PTV_SH.mhd \-o imgs/PTV_CROP.mhd \-ixi 0  \-ixf "$dimX" \-iyi 0 \-iyf "$dimY" \-izi "$cropZ" \-izf "$dimZ" > trash.out
for roi in "${ROIs[@]}"; do
  python3 starter_kit/mhd_crop.py \-i imgs/${roi}_SH.mhd \-o imgs/${roi}_CROP.mhd \-ixi 0  \-ixf "$dimX" \-iyi 0 \-iyf "$dimY" \-izi "$cropZ" \-izf "$dimZ" > trash.out
done

echo "CT, PTV and ROIs cropped"

output=$(grep "new_offset:" crop.out | head -n 1)
read -r X0crop Y0crop Z0crop < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')
echo "Images cropped at [cm]: $X0crop $Y0crop $Z0crop"
#***for USRBIN***

output=$(grep "new_end:" crop.out | head -n 1)
read -r Xfcrop Yfcrop Zfcrop < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')

output=$(grep "new_dims:" crop.out | head -n 1)
read -r dimXcrop dimYcrop dimZcrop < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')
#*** ***

echo "Cropped images dimensions: $dimXcrop $dimYcrop $dimZcrop"

bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 4 "${Xfcrop}" 1
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 5 "${Yfcrop}" 1
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 6 "${Zfcrop}" 1
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 1 "${X0crop}" 2
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 2 "${Y0crop}" 2
#bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 3 "${Z0crop}" 2 #dose map evaluated starting from 0 by default
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 4 "${dimXcrop}" 2
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 5 "${dimYcrop}" 2
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 6 "${dimZcrop}" 2

X0voxels=$(echo "${X0crop} * -1" | bc -l)
Y0voxels=$(echo "${Y0crop} * -1" | bc -l)

bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "VOXELS" 1 "${X0crop}"
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "VOXELS" 2 "${Y0crop}"

#mv starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp

echo "FLUKA inp modified"
SHz=$(echo "${Z0crop} * -1" | bc -l)
python3 starter_kit/mhd_shift.py "imgs/CT_CROP.mhd" 0 0 ${SHz} imgs/CT_plan.mhd > trash.out
python3 starter_kit/mhd_shift.py "imgs/PTV_CROP.mhd" 0 0 ${SHz} imgs/PTV_plan.mhd > trash.out
for roi in ${ROIs[@]}; do
  python3 starter_kit/mhd_shift.py "imgs/${roi}_CROP.mhd" 0 0 ${SHz} imgs/${roi}_plan.mhd > trash.out
done

echo "z-offset removed from CT; PTV and ROIs"

mv starter_kit/fredCTHU2flukaVoxels.py imgs/
mv starter_kit/hu2materials_fredEMGPU.txt imgs/
mv starter_kit/mhd_clamp.py imgs/

cd imgs/

file=CT_plan.mhd
filename="${file%%.*}"
filenameclamped="${filename}_clamped.mhd"

python3 mhd_clamp.py ${file} -o ${filenameclamped} -vmin -1000 -vmax 1376 > trash.out \
&& python3 fredCTHU2flukaVoxels.py ${filenameclamped} hu2materials_fredEMGPU.txt > trash.out
mv CTHU.vxl ${filename}.vxl
rm ${filenameclamped}

mv ${filename}.vxl CT.vxl
echo "CT .mhd file converted to .vxl as required by FLUKA"
cd ..

mv imgs/fredCTHU2flukaVoxels.py starter_kit/
mv imgs/hu2materials_fredEMGPU.txt starter_kit/
mv imgs/mhd_clamp.py starter_kit/

python3 starter_kit/mhd_astype.py imgs/CT_plan.mhd float32 > trash.out
python3 starter_kit/mhd_astype.py imgs/PTV_plan.mhd float32 > trash.out
for roi in ${ROIs[@]}; do
python3 starter_kit/mhd_astype.py imgs/${roi}_plan.mhd float32 > trash.out
done

echo "CT, PTV and ROIs values converted to float32" 

Ws=()
Hs=()
for E in "${energies[@]}"; do
for deg in "${angles[@]}"; do
  read -r W < <(grep "Width_${deg}:" rectangles.out | awk '{
    for(i=1; i <= NF; i++) {
      if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
        print $i; exit;
      }
    }
  }')
  Ws+=("$W")

  read -r H < <(grep "Height_${deg}:" rectangles.out | awk '{
    for(i=1; i <= NF; i++) {
      if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
        print $i; exit;
      }
    }
  }')
  Hs+=("$H")

# 4. MODIFY FLUKA file - field size
  mkdir -p "sim${E}MeV_${deg}deg"
  cp starter_kit/EF70mm.inp sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp 

  Xi_right=$(echo "scale=3; $W / 2" | bc -l)
  Xf_right=$(echo "scale=3; ${Xi_right} + ${height}" | bc -l) 
  Xf_left=$(echo "scale=3; -1 * $W / 2" | bc -l)
  Xi_left=$(echo "scale=3; ${Xf_left} - ${height}" | bc -l)
  Yf_down=$(echo "scale=3; -1 * $H / 2" | bc -l)
  Yi_down=$(echo "scale=3; ${Yf_down} - ${height}" | bc -l)
  Yi_up=$(echo "scale=3; $H / 2" | bc -l)
  Yf_up=$(echo "scale=3; ${Yi_up} + ${height}" | bc -l)
#  echo "$Xi_right $Xf_right $Xi_left $Xf_left $Yi_down $Yf_down $Yi_up $Yf_up"
#  sed -i -E "s/^((RPP[[:space:]]+lam1[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Yi_up} ${Yf_up} /" sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp 
#  sed -i -E "s/^((RPP[[:space:]]+lam2[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Yi_down} ${Yf_down} /" sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp 
#  sed -i -E "s/^((RPP[[:space:]]+lam3[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Xi_left} ${Xf_left} /" sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp 
#  sed -i -E "s/^((RPP[[:space:]]+lam4[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Xi_right} ${Xf_right} /" sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp 

  # lam1: aggiorna Ymin (k=3) e Ymax (k=4)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    RPP lam1 3 "${Yi_up}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    RPP lam1 4 "${Yf_up}"

  # lam2: aggiorna Ymin (k=3) e Ymax (k=4)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    RPP lam2 3 "${Yi_down}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    RPP lam2 4 "${Yf_down}"

  # lam3: aggiorna Xmin (k=1) e Xmax (k=2)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    RPP lam3 1 "${Xi_left}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    RPP lam3 2 "${Xf_left}"

  # lam4: aggiorna Xmin (k=1) e Xmax (k=2)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    RPP lam4 1 "${Xi_right}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV_${deg}deg/tmp.inp" \
    "sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp" \
    RPP lam4 2 "${Xf_right}"



degRot=$(echo "${deg} * -1" | bc -l)

  bash starter_kit/card_modifier.sh sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp sim${E}MeV_${deg}deg/tmp.inp "ROT-DEFI" 3 "${degRot}" 2
  bash starter_kit/card_modifier.sh sim${E}MeV_${deg}deg/tmp.inp sim${E}MeV_${deg}deg/EF70mm${E}MeV_${deg}deg.inp "SOURCE" 7 "${E}MeV" 1
  rm sim${E}MeV_${deg}deg/tmp.inp
  cp imgs/CT.vxl sim${E}MeV_${deg}deg/CT.vxl
  cp starter_kit/simkit${E}MeV/* sim${E}MeV_${deg}deg/

  cd sim${E}MeV_${deg}deg
  bash crea_input_homemade.sh EF70mm${E}MeV_${deg}deg.inp ${INPs} > trash.out
  echo "FLUKA input files created for ${E}MeV - ${deg}° treatment dose evaluation"
  pids=()
  what_pid=()
  mkdir -p logs
  for run in $(seq 1 "${INPs}"); do
    nohup ${FLUPRO}/flutil/rfluka run_${run}R.inp -e fluka_EF_${E}MeV.exe -N0 -M1 > "logs/run_${run}.log" 2>&1 &
    pids+=($!)
    what_pid+=("${E}MeV_${deg}°_run${run}")
  done
  printf "%s\n" "${pids[@]}" > processall.pid
  cd ..

done
done

echo "Waiting all the simulations to be completed"
fail=0
i=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then 
    ((fail++))
    echo "Process PID ${pid} (${what_pid[i]}) terminated with an error"
  fi
  (( i++ ))
done

if ((fail==0)); then
  echo "All simulations terminated succesfully"
else
  echo "$fail process(es) terminated with an error"
fi

all_rois_path=()
ROIs_plan=()
for roi in "${ROIs[@]}"; do
all_rois_path+=("../imgs/${roi}_plan.mhd")
ROIs_plan+=("${roi}_plan")
done 



pids=()
for E in "${energies[@]}"; do
for deg in "${angles[@]}"; do
  cd sim${E}MeV_${deg}deg
  for run in $(seq 1 "${INPs}")
  do
    nohup ./bnn2mhd run_${run}R001_fort.23 dose_tot_run_${run}.mhd -Gy > trash.out 2>&1 & 
    pids+=("$!")
  done
  cd ..
done 
done

echo "Waiting all .bnn files to be converted to .mhd maps"
fail=0
i=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then 
    ((fail++))
    echo "Process PID $pid (${what_pid[i]}) terminated with an error"
  fi
  (( i++ ))
done

if ((fail==0)); then
  echo "All conversions terminated successfully"
else
  echo "$fail process(es) terminated with an error"
fi


pids=()
for E in "${energies[@]}"; do
for deg in "${angles[@]}"; do
  cd sim${E}MeV_${deg}deg
  nohup mhd_combine.py -avg dose_tot_run* > trash.out 2>&1 &
  pids+=("$!")
  cd ..
done 
done

echo "Waiting maps to be combined"
fail=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then 
    ((fail++))
    echo "Process PID $pid terminated with an error"
  fi
done

if ((fail==0)); then
  echo "All maps combined - dose maps created successfully"
else
  echo "$fail process(es) terminated with an error"
fi
pids=()
for E in "${energies[@]}"; do
for deg in "${angles[@]}"; do
  cd sim${E}MeV_${deg}deg
  cp ../starter_kit/mhd_smooth.x ./
  ./mhd_smooth.x avg.mhd -o avg_smooth.mhd > trash.out
  mv avg_smooth.mhd DOSE_${E}MeV_${deg}deg.mhd

  python3 ../starter_kit/mhd_info.py -v DOSE_${E}MeV_${deg}deg.mhd > dose_info.out
  output=$(grep "range=" dose_info.out | head -n 1)

  read -r min max < <(echo "$output" | awk '{
    count = 0;
    for (i=1; i<=NF; i++) {
      if ($i ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) {
        printf "%s ", $i;
        count++;
        if (count == 2) break;
      }
    }
  }')

  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV_${deg}deg.mhd -divider ${max}
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV_${deg}deg.mhd -multiplier 50

  echo "Creating dose map and DVH for ${E}MeV ${deg}°"
  mkdir -p DVH${E}MeV_${deg}deg
  ../starter_kit/ComputeDVH/ComputeDVH.x -Dgoal 2000 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_${deg}deg.mhd -type float -fileLabel ${E}MeV${deg}deg -dir DVH${E}MeV_${deg}deg > trash.out
  cp ../starter_kit/find_prescription_point.py ./
  rescale_factor=$(python3 find_prescription_point.py "DVH${E}MeV_${deg}deg/PTV_plan${E}MeV${deg}deg.txt" ${preD} ${preV})
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV_${deg}deg.mhd -multiplier "${rescale_factor}"
  ../starter_kit/ComputeDVH/ComputeDVH.x -Dgoal 2000 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_${deg}deg.mhd -type float -fileLabel ${E}MeV${deg}deg -dir DVH${E}MeV_${deg}deg > trash.out
  python3 ../starter_kit/plotDVH.py -label1 ${E}MeV${deg}deg -dir1 DVH${E}MeV_${deg}deg -roi PTV_plan "${ROIs_plan[@]}" > trash.out
  nohup  python3 ../starter_kit/mhd_viewer_RayS.py DOSE_${E}MeV_${deg}deg.mhd -CT ../imgs/CT_plan.mhd -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -png > trash.out &
  pids+=("$!")
  cd ..
done 
done

echo "Waiting png images to be created and saved"
fail=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then 
    ((fail++))
    echo "Process PID $pid terminated with an error"
  fi
done

if ((fail==0)); then
  echo "All dose maps and DVHs created!"
else
  echo "$fail process(es) terminated with an error"
fi


mv starter_kit/compare_all_DVHs.py ./
python3 compare_all_DVHs.py --base-dir . --out DVH_ALL.png --legend right > trash.out
pid="$!"

SPINNER_MSG="Generating final plot"
if ! wait_with_spinner_and_report "$pid"; then
  continue
fi
mv compare_all_DVHs.py starter_kit/compare_all_DVHs.py
echo "Wait just another moment and enjoy!"
xdg-open DVH_ALL.png &






