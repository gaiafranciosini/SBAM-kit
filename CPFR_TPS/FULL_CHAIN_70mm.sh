#!/bin/bash
#++++++++++++++++++++++++++++++++++

file="$1"
CT="$2"
PTV="$3"
MARKER="$4"

if [ "$#" -lt 4 ]; then
  echo "Use: CT PTV MARKER [ROI1 ROI2 ROI3 ROI4 ...]"
  exit 1
fi
echo "$1 $2 $3 $4" > input_log.out
shift 4
# Array che conterr   i nomi puliti
ROIs=()

# Ciclo sugli argomenti rimanenti
for path in "$@"; do
  # Rimuove "imgs/" e ".mhd"
  name="${path##*/}"         # rimuove il percorso  ^f^r osso.mhd
  name="${name%.mhd}"        # rimuove l'estensione  ^f^r osso
  ROIs+=("$name")     # aggiunge all'array
done
#++++++++++++++++++++++++++++++++++
echo "ROI: ${ROIs[@]}" >> input_log.out

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


#1. Get BEAM DIRECTION

nohup python3 starter_kit/GetDirection.py "$CT" -PTV "$PTV" -marker "$MARKER"  > out.out &
pid=$!
echo "Computing beam direction"

  fail=0
  if ! wait "$pid"; then
    ((fail++))
    echo "Process PID $pid terminated with an error"
  fi

  if (( fail == 0 )); then
    echo "Process terminated successfully"
  else
    echo "$fail process terminated with an error"
  fi

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

angles=(0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85)

python3 starter_kit/GetFieldSize.py imgs/PTV_SH.mhd -rot "${angles[@]}" > rectangles.out &
pid=$!
echo "Optimising geometry"

  fail=0
  if ! wait "$pid"; then
    ((fail++))
    echo "Process PID $pid terminated with an error"
  fi

  if (( fail == 0 )); then
    echo "Process terminated successfully"
  else
    echo "$fail process terminated with an error"
  fi

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
if ((cropZ<0)); then
    cropZ=0
fi

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


output=$(grep "new_L:" crop.out | head -n 1)
read -r newLx newLy newLz < <(echo "$output" | awk '{
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

Xgrid=8 #cm     -Xgrid/2 |------0------| +Xgrid/2
Ygrid=8 #cm     -Ygrid/2 |------0------| +Ygrid/2
Zgrid=6 #cm            0 |-------------| +Zgrid

python3 starter_kit/GetGridSize.py -ct imgs/CT_CROP.mhd -size ${Xgrid} ${Ygrid} ${Zgrid} > grid_size.out

output=$(grep "grid_coord: " grid_size.out | head -n 1)
read -r gridXin gridXf gridYin gridYf gridZin gridZf < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 6) break;
    }
  }
}')


output=$(grep "size_idx: " grid_size.out | head -n 1)
read -r gridXidx gridYidx gridZidx < <(echo "$output" | awk '{
  count = 0;
  for (i=1; i<=NF; i++) {
    if ($i ~ /^[-+]?[0-9]/) {
      printf "%s ", $i;
      count++;
      if (count == 3) break;
    }
  }
}')

echo "Cropped images dimensions: $dimXcrop $dimYcrop $dimZcrop"
echo "New images length: $newLx $newLy $newLz"
echo "DOSE GRID SIZE: ${gridXidx} ${gridYidx} $gridZidx}" 
echo "X: ${Xgrid} (${gridXin} ${gridXf}) [cm]" 
echo "Y: ${Ygrid} (${gridYin} ${gridYf}) [cm]" 
echo "Z: ${Zgrid} (${gridZin} ${gridZf}) [cm]" 

bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 4 "${gridXf}" 1
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 5 "${gridYf}" 1
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 6 "${gridZf}" 1
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 1 "${gridXin}" 2
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 2 "${gridYin}" 2
#bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 3 "${Z0crop}" 2 #dose map evaluated starting from 0 by default
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 4 "${gridXidx}" 2
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 5 "${gridYidx}" 2
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 6 "${gridZidx}" 2

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

read -r deg < <(grep "BEST_ANGLE:" rectangles.out | awk '{
  for(i=1; i <= NF; i++) {
    if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
      print $i; exit;
    }
  }
}')


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


for E in "${energies[@]}"; do
mkdir -p "sim${E}MeV"
rm -f sim${E}MeV/*
python3 starter_kit/GetBestSize.py -ptv imgs/PTV_plan.mhd -lutpath starter_kit/LUTs -min_size ${W} ${H} -E ${E} > "sim${E}MeV/field_size.out"
#usage: GetBestSize.py [-h] -ptv path -lutpath path -min_size MIN_SIZE MIN_SIZE -E E

read -r bestXsize < <(grep "Xsize:" sim${E}MeV/field_size.out | awk '{
  for(i=1; i <= NF; i++) {
    if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
      print $i; exit;
    }
  }
}')

read -r bestYsize < <(grep "Ysize:" sim${E}MeV/field_size.out | awk '{
  for(i=1; i <= NF; i++) {
    if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
      print $i; exit;
    }
  }
}')

# 4. MODIFY FLUKA file - field size
  cp starter_kit/EF70mm.inp "sim${E}MeV/EF70mm${E}MeV.inp"

  python3 starter_kit/GetSlitMargins.py ${width} ${height} ${bestXsize} ${bestYsize} > "sim${E}MeV/aperture.out"
  read -r Xi_right Xf_right Xi_left Xf_left Yi_down Yf_down Yi_up Yf_up < <(
  awk '/^APERTURE:/{ 
    for(i=2;i<=NF;i++) 
      if($i ~ /^[+-]?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?$/) 
        printf "%s ", $i; 
    print ""
  }' "sim${E}MeV/aperture.out"
)

#  echo "$Xi_right $Xf_right $Xi_left $Xf_left $Yi_down $Yf_down $Yi_up $Yf_up"
#  sed -i -E "s/^((RPP[[:space:]]+lam1[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Yi_up} ${Yf_up} /" sim${E}MeV/EF70mm${E}MeV.inp 
#  sed -i -E "s/^((RPP[[:space:]]+lam2[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Yi_down} ${Yf_down} /" sim${E}MeV/EF70mm${E}MeV.inp 
#  sed -i -E "s/^((RPP[[:space:]]+lam3[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Xi_left} ${Xf_left} /" sim${E}MeV/EF70mm${E}MeV.inp 
#  sed -i -E "s/^((RPP[[:space:]]+lam4[[:space:]])+-?[0-9.]+[[:space:]]+-?[0-9.]+)/\1 ${Xi_right} ${Xf_right} /" sim${E}MeV/EF70mm${E}MeV.inp 

  # lam1: aggiorna Ymin (k=3) e Ymax (k=4)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    "sim${E}MeV/tmp.inp" \
    RPP lam1 3 "${Yi_up}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/tmp.inp" \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    RPP lam1 4 "${Yf_up}"

  # lam2: aggiorna Ymin (k=3) e Ymax (k=4)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    "sim${E}MeV/tmp.inp" \
    RPP lam2 3 "${Yi_down}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/tmp.inp" \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    RPP lam2 4 "${Yf_down}"

  # lam3: aggiorna Xmin (k=1) e Xmax (k=2)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    "sim${E}MeV/tmp.inp" \
    RPP lam3 1 "${Xi_left}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/tmp.inp" \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    RPP lam3 2 "${Xf_left}"

  # lam4: aggiorna Xmin (k=1) e Xmax (k=2)
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    "sim${E}MeV/tmp.inp" \
    RPP lam4 1 "${Xi_right}"
  bash starter_kit/object_card_modifier.sh \
    "sim${E}MeV/tmp.inp" \
    "sim${E}MeV/EF70mm${E}MeV.inp" \
    RPP lam4 2 "${Xf_right}"

  read -r Wfield Hfield < <(
  awk '/^FIELD_SIZE:/{
    for(i=2;i<=NF;i++)
      if($i ~ /^[+-]?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?$/)
        printf "%s ", $i;
    print ""
  }' "sim${E}MeV/aperture.out"
)

python3 starter_kit/GetBEV.py -ptv imgs/PTV_SH.mhd -slitsize "${width}" "${height}" -W "${bestXsize}" -H "${bestYsize}" -angle "${deg}" -showsave 1

echo "---------------------------------------------------------------"
echo "                      FIELD SIZE ${E}MeV                       "
echo "---------------------------------------------------------------"
echo " "
echo "                    WIDTH: ${bestXsize} cm                     "
echo "                   HEIGHT: ${bestYsize} cm                     "
echo " "
echo "---------------------------------------------------------------"




degRot=$(echo "${deg} * -1" | bc -l)

  bash starter_kit/card_modifier.sh "sim${E}MeV/EF70mm${E}MeV.inp" "sim${E}MeV/tmp.inp" "ROT-DEFI" 3 "${degRot}" 2
  bash starter_kit/card_modifier.sh "sim${E}MeV/tmp.inp" "sim${E}MeV/EF70mm${E}MeV.inp" "SOURCE" 7 "${E}MeV" 1
  rm sim${E}MeV/tmp.inp
  cp imgs/CT.vxl "sim${E}MeV/CT.vxl"
  cp starter_kit/simkit${E}MeV/* "sim${E}MeV/"



  cd sim${E}MeV
  bash crea_input_homemade.sh EF70mm${E}MeV.inp ${INPs} > trash.out
  echo "FLUKA input files created for ${E}MeV treatment dose evaluation"
  pids=()
  what_pid=()
  mkdir -p logs

    for run in $(seq 1 "${INPs}"); do
        nohup ${FLUPRO}/flutil/rfluka run_${run}R.inp -e fluka_EF_${E}MeV.exe -N0 -M1 > "logs/run_${run}.log" 2>&1 &
        pids+=($!)
        what_pid+=("${E}MeV_run${run}")
    done

  printf "%s\n" "${pids[@]}" > processall.pid
  cd ..

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
  cd sim${E}MeV
  for run in $(seq 1 "${INPs}")
  do
    nohup ../starter_kit/bnn2mhd run_${run}R001_fort.23 dose_tot_run_${run}.mhd -Gy > trash.out 2>&1 & 
    pids+=("$!")
  done
  cd ..
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
  cd sim${E}MeV
  check=(dose_tot_run*)
  if (( ${#check[@]} == 1 )); then
    cp "${check[0]}" dose_tot_run_copy.mhd
  fi
  nohup ../starter_kit/mhd_combine.py -avg dose_tot_run* > trash.out 2>&1 &
  pids+=("$!")
  cd ..
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
  cd sim${E}MeV
  cp ../starter_kit/mhd_smooth.x ./
  ./mhd_smooth.x avg.mhd -o avg_smooth.mhd > trash.out
  mv avg_smooth.mhd DOSE_${E}MeV_GRID.mhd

  python3 ../starter_kit/mhd_refill.py -ct ../imgs/CT_plan.mhd -dose DOSE_${E}MeV_GRID.mhd -out DOSE_${E}MeV.mhd
  python3 ../starter_kit/mhd_astype.py DOSE_${E}MeV.mhd float32
  python3 ../starter_kit/mhd_info.py -v DOSE_${E}MeV.mhd > dose_info.out
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

#  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -divider ${max}
#  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -multiplier 50
#
  echo "Creating dose map and DVH for ${E}MeV"
  mkdir -p DVH${E}MeV
#  echo "RESCALE FACTOR= ${rescale_factor}" > rescaling.out
#  echo "TOTAL RESCALE= ${rescale_factor} x 50 x (1/${max})" >> rescaling.out
  factor="kFLASH_${E}MeV"
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -multiplier "${!factor}"
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -multiplier "${pulses}"
  ../starter_kit/ComputeDVH/ComputeDVH.x -Dgoal 2000 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV.mhd -type float -fileLabel ${E}MeV -dir DVH${E}MeV > trash.out
  python3 ../starter_kit/plotDVH.py -label1 ${E}MeV -dir1 DVH${E}MeV -roi PTV_plan "${ROIs_plan[@]}" > trash.out
  python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV/PTV_plan${E}MeV.txt > DVH_${E}MeV/DVH_PTV.out
#  nohup  python3 ../starter_kit/mhd_viewer_RayS.py DOSE_${E}MeV.mhd -CT ../imgs/CT_plan.mhd -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -png > trash.out &
  pids+=("$!")
  cd ..
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

echo "Generating final plot"

  fail=0
  if ! wait "$pid"; then
    ((fail++))
    echo "Process PID $pid terminated with an error"
  fi

  if (( fail == 0 )); then
    echo "Process terminated successfully"
  else
    echo "$fail process terminated with an error"
  fi

mv compare_all_DVHs.py starter_kit/
echo "Enjoy!"











