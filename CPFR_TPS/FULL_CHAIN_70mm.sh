 #!/bin/bash
#++++++++++++++++++++++++++++++++++
set -e
set -o pipefail
trap 'echo "Errore: comando \"$BASH_COMMAND\" fallito alla riga $LINENO. Premi ctrl+c per uscire"; exit 1' ERR
#set -x

file="$1"
SETUP_VERIFICA_PTV="$file"  # conserva il setup originale anche dopo shift e cambi di variabile
CT="$2"
PTV="$3"
MARKER="$4"

if [ "$#" -lt 4 ]; then
  echo "Use: file CT PTV MARKER [ROI1 ROI2 ROI3 ROI4 ...]"
  exit 1
fi
echo "$1 $2 $3 $4" > input_log.out
shift 4
# Array che contiene i nomi puliti
ROIs=()

# Ciclo sugli argomenti rimanenti
for path in "$@"; do
  # Rimuove il percorso e ".mhd"
  name="${path##*/}"
  name="${name%.mhd}"
  ROIs+=("$name")
done
#++++++++++++++++++++++++++++++++++
echo "ROI: ${ROIs[@]}" >> input_log.out

# Inizializza variabili
energies=()
preD=""
preV=""
Pdose=""
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
        "PRESCRIPTION DOSE:"*)
            preD=${line#PRESCRIPTION DOSE: }
            ;;
        "PERCENTAGE PRESCRIPTION VOLUME:"*)
            preV=${line#PERCENTAGE PRESCRIPTION VOLUME: }
            ;;
        "PERCENTAGE PRESCRIPTION DOSE:"*)
            Pdose=${line#PERCENTAGE PRESCRIPTION DOSE: }
            ;;
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

# Array associativi per energia (pulses e indice della cartella AUTO_sim)
declare -A pulsesFLASH_by_E
declare -A pulsesCONV_by_E
declare -A AUTO_IDX

#CALIBRATION FACTORS (primaries/pulse)
kFLASH_7MeV=9.382238805934787e+11
kCONV_7MeV=1.42155133423254e+10
kFLASH_9MeV=7.106140035174675e+11
kCONV_9MeV=1.07668788411738e+10

#++++++++++++++++++++++++++++++++++
# HELPER
#++++++++++++++++++++++++++++++++++
W_PIDS=(); W_NAMES=()
wait_all() {  # $1 = descrizione; usa W_PIDS / W_NAMES e li svuota
  local fail=0 k=0 pid
  for pid in "${W_PIDS[@]}"; do
    if ! wait "$pid"; then
      ((++fail))
      echo "Process PID $pid (${W_NAMES[k]}) terminated with an error"
    fi
    ((++k))
  done
  if ((fail==0)); then echo "$1: OK"; else echo "$1: $fail process(es) terminated with an error"; fi
  W_PIDS=(); W_NAMES=()
}
set_rpp() {   # file  body  k  valore   (modifica in place)
  bash starter_kit/object_card_modifier.sh "$1" "$1.tmp" RPP "$2" "$3" "$4"
  mv "$1.tmp" "$1"
}
set_card() {  # file  card  k  valore  flag   (modifica in place)
  bash starter_kit/card_modifier.sh "$1" "$1.tmp" "$2" "$3" "$4" "$5"
  mv "$1.tmp" "$1"
}

#1. Get BEAM DIRECTION

mkdir -p OUTS

nohup python3 starter_kit/GetDirection.py "$CT" -PTV "$PTV" -marker "$MARKER"  > out.out &
pid=$!
echo "Computing beam direction"
if wait "$pid"; then
  echo "Process terminated successfully"
else
  echo "Process PID $pid terminated with an error"
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
  print "";
}')

echo "Beam direction: Vx:${Vx} Vy:${Vy} Vz:${Vz}"
echo "Vx:${Vx}" > OUTS/beam_direction.out
echo "Vy:${Vy}" >> OUTS/beam_direction.out
echo "Vz:${Vz}" >> OUTS/beam_direction.out


echo "Rotating PTV, ROIs and CT"
python3 starter_kit/RotROIsQuat.py "$PTV" "$Vx" "$Vy" "$Vz" imgs/PTV_ROT.mhd > ptv.out
echo "PTV rotated and saved as imgs/PTV_ROT.mhd"
for roi in "${ROIs[@]}"; do
  python3 starter_kit/RotROIsQuat.py "imgs/${roi}.mhd" "$Vx" "$Vy" "$Vz" imgs/${roi}_ROT.mhd > trash.out
  echo "${roi} rotated and saved as imgs/${roi}_ROT.mhd"
done
python3 starter_kit/RotMapsQuat.py "$CT" "$Vx" "$Vy" "$Vz" imgs/CT_ROT.mhd >    imgs/CT_ROT_matrix.txt
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
  print "";
}')
echo "X_ptv[idx]: $idxPTVx | Y_ptv[idx]: $idxPTVy | Z_ptv[idx]: $idxPTVz"

echo "X_ptv[idx]: $idxPTVx" > OUTS/PTV_ROT.out
echo "Y_ptv[idx]: $idxPTVy" >> OUTS/PTV_ROT.out
echo "Z_ptv[idx]: $idxPTVz" >> OUTS/PTV_ROT.out

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
  print "";
}')
echo "X_ptv: $PTVx | Y_ptv: $PTVy | Z_ptv: $PTVz"

echo "X_ptv: $PTVx" >> OUTS/PTV_ROT.out
echo "Y_ptv: $PTVy" >> OUTS/PTV_ROT.out
echo "Z_ptv: $PTVz" >> OUTS/PTV_ROT.out

SHx=$(echo "$PTVx * -1" | bc -l)
SHy=$(echo "$PTVy * -1" | bc -l)

python3 starter_kit/map_shift.py "imgs/CT_ROT.mhd" "${SHx}" "${SHy}" 0 imgs/CT_SH.mhd > trash.out
echo "CT shifted  and saved as imgs/CT_SH.mhd"
python3 starter_kit/map_shift.py "imgs/PTV_ROT.mhd" "${SHx}" "${SHy}" 0 imgs/PTV_SH.mhd > trash.out
echo "PTV_shifted  and saved as imgs/PTV_SH.mhd"

for roi in "${ROIs[@]}"; do
  python3 starter_kit/map_shift.py "imgs/${roi}_ROT.mhd" "${SHx}" "${SHy}" 0  imgs/${roi}_SH.mhd > trash.out
  echo "${roi} shifted and saved as imgs/${roi}_SH.mhd"
done
echo

# 3. FIELD SIZE dimensions (modify slit size with applicator)
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
if wait "$pid"; then
  echo "Process terminated successfully"
else
  echo "Process PID $pid terminated with an error"
fi

echo "Geometry optimization details reported in rectangles.out"
echo " "
echo "Cropping PTV, CT and ROIs"
cp rectangles.out OUTS/

python3 starter_kit/map_info.py imgs/CT_SH.mhd > info.out
cp info.out OUTS/
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
    print "";
}')

cropZ=$(( idxPTVz - 1 ))
if ((cropZ<0)); then
    cropZ=0
fi
echo "$cropZ" > imgs/cropZ.txt
python3 starter_kit/map_crop.py -i imgs/CT_SH.mhd -o imgs/CT_CROP.mhd -ixi 0 -ixf "$dimX" -iyi 0 -iyf "$dimY" -izi "$cropZ" -izf "$dimZ" > crop.out
cp crop.out OUTS/
python3 starter_kit/map_crop.py -i imgs/PTV_SH.mhd -o imgs/PTV_CROP.mhd -ixi 0 -ixf "$dimX" -iyi 0 -iyf "$dimY" -izi "$cropZ" -izf "$dimZ" > trash.out
for roi in "${ROIs[@]}"; do
  python3 starter_kit/map_crop.py -i imgs/${roi}_SH.mhd -o imgs/${roi}_CROP.mhd -ixi 0 -ixf "$dimX" -iyi 0 -iyf "$dimY" -izi "$cropZ" -izf "$dimZ" > trash.out
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
  print "";
}')
echo "Images cropped at [cm]: $X0crop $Y0crop $Z0crop"

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
  print "";
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
  print "";
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
  print "";
}')

Xgrid=8 #cm     -Xgrid/2 |------0------| +Xgrid/2
Ygrid=8 #cm     -Ygrid/2 |------0------| +Ygrid/2
Zgrid=6 #cm            0 |-------------| +Zgrid

python3 starter_kit/GetGridSize.py -ct imgs/CT_CROP.mhd -size ${Xgrid} ${Ygrid} ${Zgrid} > grid_size.out
cp grid_size.out OUTS/
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
  print "";
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
  print "";
}')

echo "Cropped images dimensions: $dimXcrop $dimYcrop $dimZcrop"
echo "New images length: $newLx $newLy $newLz"
echo "DOSE GRID SIZE: ${gridXidx} ${gridYidx} ${gridZidx}"
echo "X: ${Xgrid} (${gridXin} ${gridXf}) [cm]"
echo "Y: ${Ygrid} (${gridYin} ${gridYf}) [cm]"
echo "Z: ${Zgrid} (${gridZin} ${gridZf}) [cm]"

bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 4 "${gridXf}" 1
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 5 "${gridYf}" 1
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 6 "${gridZf}" 1
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 1 "${gridXin}" 2
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 2 "${gridYin}" 2
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 4 "${gridXidx}" 2
bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "USRBIN" 5 "${gridYidx}" 2
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "USRBIN" 6 "${gridZidx}" 2

bash starter_kit/card_modifier.sh starter_kit/EF70mm.inp starter_kit/tmpEF70mm.inp "VOXELS" 1 "${X0crop}"
bash starter_kit/card_modifier.sh starter_kit/tmpEF70mm.inp starter_kit/EF70mm.inp "VOXELS" 2 "${Y0crop}"

echo "FLUKA inp modified"
SHz=$(echo "${Z0crop} * -1" | bc -l)
python3 starter_kit/map_shift.py "imgs/CT_CROP.mhd" 0 0 ${SHz} imgs/CT_plan.mhd > trash.out
python3 starter_kit/map_shift.py "imgs/PTV_CROP.mhd" 0 0 ${SHz} imgs/PTV_plan.mhd > trash.out
for roi in "${ROIs[@]}"; do
  python3 starter_kit/map_shift.py "imgs/${roi}_CROP.mhd" 0 0 ${SHz} imgs/${roi}_plan.mhd > trash.out
done

echo "z-offset removed from CT; PTV and ROIs"

cp starter_kit/fredCTHU2flukaVoxels.py imgs/
cp starter_kit/hu2materials_fredEMGPU.txt imgs/
cp starter_kit/map_clamp.py imgs/
cp starter_kit/FRED_mhd_io.py imgs/
cd imgs/

file=CT_plan.mhd
filename="${file%%.*}"
filenameclamped="${filename}_clamped.mhd"

python3 map_clamp.py ${file} -o ${filenameclamped} -vmin -1000 -vmax 1376 >  conversion1.log  2>&1 \
&& python3 fredCTHU2flukaVoxels.py ${filenameclamped} hu2materials_fredEMGPU.txt >  conversion.log  2>&1
mv CTHU.vxl ${filename}.vxl
rm ${filenameclamped}

mv ${filename}.vxl CT.vxl
echo "CT .mhd file converted to .vxl as required by FLUKA"
cd ..

mv imgs/fredCTHU2flukaVoxels.py starter_kit/
mv imgs/hu2materials_fredEMGPU.txt starter_kit/
mv imgs/map_clamp.py starter_kit/

python3 starter_kit/map_astype.py imgs/CT_plan.mhd float32 > trash.out
python3 starter_kit/map_astype.py imgs/PTV_plan.mhd float32 > trash.out
for roi in "${ROIs[@]}"; do
  python3 starter_kit/map_astype.py imgs/${roi}_plan.mhd float32 > trash.out
done

echo "CT, PTV and ROIs values converted to float32"

# Angolo migliore dello shaper (usato anche per le varianti plus/minus)
read -r deg < <(grep "BEST_ANGLE:" rectangles.out | awk '{
  for(i=1; i <= NF; i++) {
    if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
      print $i; exit;
    }
  }
  print "";
}')

# ROI nel formato richiesto da ComputeDVH (percorsi relativi dalla cartella sim)
all_rois_path=()
ROIs_plan=()
for roi in "${ROIs[@]}"; do
  all_rois_path+=("../imgs/${roi}_plan.mhd")
  ROIs_plan+=("${roi}_plan")
done

#++++++++++++++++++++++++++++++++++
# SIMULAZIONI NOMINALI (tutte le energie in parallelo)
#++++++++++++++++++++++++++++++++++
W_PIDS=(); W_NAMES=()      # azzerati UNA volta, fuori dal loop sulle energie

for E in "${energies[@]}"; do
  mkdir -p "sim${E}MeV"
  rm -rf sim${E}MeV/*
  python3 starter_kit/GetBestSize_less.py -ptv imgs/PTV_plan.mhd -lutpath starter_kit/LUTs -E ${E} > "sim${E}MeV/field_size.out"

  read -r bestXsize < <(grep "Xsize:" sim${E}MeV/field_size.out | awk '{
    for(i=1; i <= NF; i++) {
      if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
        print $i; exit;
      }
    }
    print "";
  }')

  read -r bestYsize < <(grep "Ysize:" sim${E}MeV/field_size.out | awk '{
    for(i=1; i <= NF; i++) {
      if ($i ~ /^[+-]?[0-9]+(\.[0-9]+)?$/) {
        print $i; exit;
      }
    }
    print "";
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

  # lam1 / lam2: Ymin (k=3) e Ymax (k=4); lam3 / lam4: Xmin (k=1) e Xmax (k=2)
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam1 3 "${Yi_up}"
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam1 4 "${Yf_up}"
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam2 3 "${Yi_down}"
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam2 4 "${Yf_down}"
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam3 1 "${Xi_left}"
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam3 2 "${Xf_left}"
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam4 1 "${Xi_right}"
  set_rpp "sim${E}MeV/EF70mm${E}MeV.inp" lam4 2 "${Xf_right}"

  python3 starter_kit/GetBEV.py -ptv imgs/PTV_SH.mhd -slitsize "${width}" "${height}" -W "${bestXsize}" -H "${bestYsize}" -angle "${deg}" -showsave 1 -output "BEV${E}MeV.jpg"

  echo "---------------------------------------------------------------"
  echo "                      FIELD SIZE ${E}MeV                       "
  echo "---------------------------------------------------------------"
  echo " "
  echo "                    WIDTH: ${bestXsize} cm                     "
  echo "                   HEIGHT: ${bestYsize} cm                     "
  echo " "
  echo "---------------------------------------------------------------"

  degRot=$(echo "${deg} * -1" | bc -l)

  set_card "sim${E}MeV/EF70mm${E}MeV.inp" "ROT-DEFI" 3 "${degRot}" 2
  set_card "sim${E}MeV/EF70mm${E}MeV.inp" "SOURCE" 7 "${E}MeV" 1
  set_card "sim${E}MeV/EF70mm${E}MeV.inp" "START" 1 "${primaries_per_CPU}" 1
  cp imgs/CT.vxl "sim${E}MeV/CT.vxl"
  cp starter_kit/simkit${E}MeV/* "sim${E}MeV/"

  ( cd "sim${E}MeV" && bash crea_input_homemade.sh EF70mm${E}MeV.inp ${INPs} > trash.out )
  echo "FLUKA input files created for ${E}MeV treatment dose evaluation"
  mkdir -p "sim${E}MeV/logs"

  for run in $(seq 1 "${INPs}"); do
    ( cd "sim${E}MeV" && nohup ${FLUPRO}/flutil/rfluka run_${run}R.inp -e fluka_EF_${E}MeV.exe -N0 -M1 > "logs/run_${run}.log" 2>&1 ) &
    W_PIDS+=($!)
    W_NAMES+=("${E}MeV_run${run}")
  done
done

echo "Waiting all the simulations to be completed"
wait_all "Simulazioni nominali"

# bnn -> mhd
for E in "${energies[@]}"; do
  for run in $(seq 1 "${INPs}"); do
    ( cd "sim${E}MeV" && nohup ../starter_kit/bnn2mhd run_${run}R001_fort.23 dose_tot_run_${run}.mhd -Gy > trash.out 2>&1 ) &
    W_PIDS+=($!)
    W_NAMES+=("bnn2mhd_${E}MeV_run${run}")
  done
done
echo "Waiting all .bnn files to be converted to .mhd maps"
wait_all "Conversione bnn->mhd nominale"

# media dei run
for E in "${energies[@]}"; do
  (
    cd "sim${E}MeV"
    check=(dose_tot_run*)
    if (( ${#check[@]} == 1 )); then
      cp "${check[0]}" copy_dose_tot_run.mhd
    fi
    nohup python3 ../starter_kit/map_combine.py -avg dose_tot_run* > trash.out 2>&1
  ) &
  W_PIDS+=($!)
  W_NAMES+=("map_combine_${E}MeV")
done
echo "Waiting maps to be combined"
wait_all "Combinazione mappe nominale"

DVHdirs=()
for E in "${energies[@]}"; do
  cd sim${E}MeV
  cp ../starter_kit/mhd_smooth.x ./
  ./mhd_smooth.x avg.mhd -o avg_smooth.mhd > trash.out
  mv avg_smooth.mhd DOSE_${E}MeV_GRID.mhd

  python3 ../starter_kit/map_refill.py -ct ../imgs/CT_plan.mhd -dose DOSE_${E}MeV_GRID.mhd -out DOSE_${E}MeV.mhd
  python3 ../starter_kit/map_astype.py DOSE_${E}MeV.mhd float32

  echo "Creating dose map and DVH for ${E}MeV"
  mkdir -p DVH${E}MeV_FLASH
  mkdir -p DVH${E}MeV_CONV
  factor_FLASH="kFLASH_${E}MeV"
  factor_CONV="kCONV_${E}MeV"

  cp DOSE_${E}MeV.mhd DOSE_${E}MeV_FLASH.mhd
  mv DOSE_${E}MeV.mhd DOSE_${E}MeV_CONV.mhd
  python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_FLASH.mhd -multiplier "${!factor_FLASH}"   # 1pulse
  cp DOSE_${E}MeV_FLASH.mhd DOSE_${E}MeV_1pulses_FLASH.mhd
  python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_CONV.mhd -multiplier "${!factor_CONV}"   # 1pulse
  cp DOSE_${E}MeV_CONV.mhd DOSE_${E}MeV_1pulses_CONV.mhd
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_FLASH.mhd -type float -fileLabel ${E}MeV_FLASH -dir DVH${E}MeV_FLASH --binningDVH 10000 > trash.out
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_CONV.mhd -type float -fileLabel ${E}MeV_CONV -dir DVH${E}MeV_CONV --binningDVH 10000 > trash.out
  pulses_FLASH=$(python3 ../starter_kit/findPulses.py DVH${E}MeV_FLASH/PTV_plan${E}MeV_FLASH.txt ${preD} ${preV} --preD_percent ${Pdose})
  pulses_CONV=$(python3 ../starter_kit/findPulses.py DVH${E}MeV_CONV/PTV_plan${E}MeV_CONV.txt ${preD} ${preV} --preD_percent ${Pdose})

  # memorizza i pulses di QUESTA energia
  pulsesFLASH_by_E[$E]=$pulses_FLASH
  pulsesCONV_by_E[$E]=$pulses_CONV

  python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_FLASH.mhd -multiplier "${pulses_FLASH}"
  python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_CONV.mhd -multiplier "${pulses_CONV}"
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_FLASH.mhd -type float -fileLabel ${E}MeV_FLASH -dir DVH${E}MeV_FLASH --binningDVH 1000  > trash.out
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_CONV.mhd -type float -fileLabel ${E}MeV_CONV -dir DVH${E}MeV_CONV --binningDVH 1000  > trash.out
  python3 ../starter_kit/plotDVH.py -label1 ${E}MeV_FLASH -dir1 DVH${E}MeV_FLASH -roi PTV_plan "${ROIs_plan[@]}" > trash.out
  python3 ../starter_kit/plotDVH.py -label1 ${E}MeV_CONV -dir1 DVH${E}MeV_CONV -roi PTV_plan "${ROIs_plan[@]}" > trash.out
  python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV_FLASH/PTV_plan${E}MeV_FLASH.txt -preD ${preD} > DVH${E}MeV_FLASH/DVH_PTV_FLASH.out
  python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV_CONV/PTV_plan${E}MeV_CONV.txt -preD ${preD} > DVH${E}MeV_CONV/DVH_PTV_CONV.out
  cp DOSE_${E}MeV_FLASH.mhd DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH.mhd
  cp DOSE_${E}MeV_CONV.mhd DOSE_${E}MeV_${pulses_CONV}pulses_CONV.mhd
  DVHdirs+=("sim${E}MeV/DVH${E}MeV_FLASH")
  DVHdirs+=("sim${E}MeV/DVH${E}MeV_CONV")
  cd ..
  echo "Ripristino le dimensioni originali della Mappa di Dose per ${E}MeV..."

  # 1. De-crop (padding)
  python3 starter_kit/mhd_decrop.py -i sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH.mhd -o sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED.mhd -cz "$cropZ"
  python3 starter_kit/mhd_decrop.py -i sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV.mhd -o sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED.mhd -cz "$cropZ"

  # 2. De-shift
  invSHx=$PTVx
  echo "$invSHx" > imgs/invSHx.txt
  invSHy=$PTVy
  echo "$invSHy" > imgs/invSHy.txt
  invSHz=$Z0crop
  echo "$invSHz" > imgs/invSHz.txt

  python3 starter_kit/map_shift.py "sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED.mhd" "${invSHx}" "${invSHy}" "${invSHz}" sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED_deSH.mhd
  python3 starter_kit/map_shift.py "sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED.mhd" "${invSHx}" "${invSHy}" "${invSHz}" sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED_deSH.mhd
  awk '/ROTATION MATRIX/{flag=1; next} /Rotated image/{flag=0} flag' imgs/CT_ROT_matrix.txt | tr -d '[]' > imgs/rotation_matrix.txt

  # 3. Contro-rotazione
  python3 starter_kit/RotDoseQuat.py "sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED_deSH.mhd" "imgs/rotation_matrix.txt" "sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_FINAL.mhd"
  echo "Dose Map FLASH rotated and saved in sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_FINAL.mhd"

  python3 starter_kit/RotDoseQuat.py "sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED_deSH.mhd" "imgs/rotation_matrix.txt" "sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_FINAL.mhd"
  echo "Dose Map CONV rotated and saved in sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_FINAL.mhd"

  mv sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH.mhd sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_BEV.mhd
  mv sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV.mhd  sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_BEV.mhd
  echo "---------------------------------------------------------------"
  echo "                   PULSES ${E}MeV FLASH: ${pulses_FLASH}       "
  echo "                   PULSES ${E}MeV CONV: ${pulses_CONV}         "
  echo "---------------------------------------------------------------"
done

echo "DVHdirs content:"
printf "%s\n" "${DVHdirs[@]}"

python3 starter_kit/GetDVHPlot.py "${DVHdirs[@]}" --xunit "cGy" --out "./DVH_ALL.png" --label-mode file --min-volume 0.01

mkdir -p files
for E in "${energies[@]}"; do

  cp sim${E}MeV/EF70mm${E}MeV.inp files/

  # indice della cartella AUTO per QUESTA energia (salvato, non piu' variabile globale "i")
  idx=1
  while [ -d "AUTO_sim${E}MeV_$idx" ]; do
    ((++idx))
  done
  mv "sim${E}MeV" "AUTO_sim${E}MeV_$idx"
  AUTO_IDX[$E]=$idx

  pf="${pulsesFLASH_by_E[$E]}"
  pc="${pulsesCONV_by_E[$E]}"
  if [ -z "$pf" ] || [ -z "$pc" ]; then
    echo "ERRORE: pulses non trovati per E=${E} (pf='${pf}' pc='${pc}')" >&2
    exit 1
  fi

  rm -f "AUTO_sim${E}MeV_$idx"/*PADDED_deSH.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/*PADDED.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/*PADDED_deSH.raw
  rm -f "AUTO_sim${E}MeV_$idx"/*PADDED.raw
  rm -f "AUTO_sim${E}MeV_$idx"/*fort.23
  rm -f "AUTO_sim${E}MeV_$idx"/*fort.19
  rm -f "AUTO_sim${E}MeV_$idx"/*.log
  rm -f "AUTO_sim${E}MeV_$idx"/*R00*.out
  rm -f "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_1pulses_CONV.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_1pulses_FLASH.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_GRID.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_orig.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_CONV.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_FLASH.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/dose_tot*.mhd
  rm -f "AUTO_sim${E}MeV_$idx"/avg.mhd

  python3 starter_kit/create_isodosemap.py -i "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pf}pulses_FLASH_BEV.mhd -o "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pf}pulses_FLASH_BEV_isodosemap.mhd
  python3 starter_kit/create_isodosemap.py -i "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pc}pulses_CONV_BEV.mhd -o "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pc}pulses_CONV_BEV_isodosemap.mhd
  python3 starter_kit/create_isodosemap.py -i "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pf}pulses_FLASH_FINAL.mhd -o "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pf}pulses_FLASH_FINAL_isodose.mhd
  python3 starter_kit/create_isodosemap.py -i "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pc}pulses_CONV_FINAL.mhd -o "AUTO_sim${E}MeV_$idx"/DOSE_${E}MeV_${pc}pulses_CONV_FINAL_isodose.mhd
done

# isocentro del PTV in cm: una riga "x y z"
python3 starter_kit/ptv_isocenter.py imgs/PTV_plan.mhd > isocenter.out
ISO=$(cat isocenter.out)
[[ -n "$ISO" ]] || { echo "Isocentro non calcolato"; exit 1; }

export FRED_MATPLOTLIB_BACKEND=Agg

shopt -s nullglob
for E in "${energies[@]}"; do
  idx="${AUTO_IDX[$E]}"
  SIMDIR="AUTO_sim${E}MeV_${idx}"

  for MODE in FLASH CONV; do
    # il numero di impulsi nel nome viene sostituito dall'asterisco
    maps=( "${SIMDIR}"/DOSE_${E}MeV_*pulses_${MODE}_BEV_isodosemap.mhd )

    if [[ ${#maps[@]} -eq 0 ]]; then
      echo "Attenzione: mappa ${MODE} non trovata in ${SIMDIR}"
      continue
    elif [[ ${#maps[@]} -gt 1 ]]; then
      echo "Attenzione: piu' mappe ${MODE} in ${SIMDIR}, uso ${maps[0]}"
    fi
    MAP="${maps[0]}"
    OUT="${SIMDIR}/ISODOSE_BEV_${E}_${MODE}_${idx}.png"

    rm -f plot.png
    python3 -W ignore starter_kit/map_viewer_isodose.py \
        -CT imgs/CT_plan.mhd -map "$MAP" \
        -Dgoal 1 -ISO $ISO -png > /dev/null
    if [[ -f plot.png ]]; then
      mv plot.png "$OUT"
      echo "Creato $OUT"
    else
      echo "Errore: il viewer non ha prodotto l'immagine per $MAP"
    fi
  done
done

#++++++++++++++++++++++++++++++++++
# VARIANTI APERTURA +2 mm / -2 mm  (in SEQUENZA: prima plus, poi minus)
#++++++++++++++++++++++++++++++++++
run_variant() {
  local TAG="$1"                 # plus | minus
  local VAR="${TAG}2mm"
  local E SIM AUTO inp run MODE p base degRot cropZ invSHx invSHy invSHz
  local Xi_right Xf_right Xi_left Xf_left Yi_down Yf_down Yi_up Yf_up
  local pulses_FLASH pulses_CONV factor_FLASH factor_CONV
  local -A P
  local DVHdirs=()
  W_PIDS=(); W_NAMES=()

  # ---- 1. input + lancio simulazioni (tutte le energie di QUESTA variante) ----
  echo "=== ${VAR}: preparazione e lancio simulazioni ==="
  for E in "${energies[@]}"; do
    SIM="sim${E}MeV_${VAR}"
    AUTO="AUTO_sim${E}MeV_${AUTO_IDX[$E]}"
    rm -rf "$SIM"; mkdir -p "$SIM/logs"
    inp="$SIM/EF70mm${E}MeV.inp"
    cp "files/EF70mm${E}MeV.inp" "$inp"
    cp "$AUTO/aperture_${VAR}.out" "$SIM/"

    # legge SOLO il file di aperture di questa variante
    read -r Xi_right Xf_right Xi_left Xf_left Yi_down Yf_down Yi_up Yf_up < <(
      awk '/^APERTURE:/{
        for(i=2;i<=NF;i++)
          if($i ~ /^[+-]?[0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?$/)
            printf "%s ", $i;
        print ""
      }' "$SIM/aperture_${VAR}.out"
    )

    set_rpp "$inp" lam1 3 "$Yi_up";    set_rpp "$inp" lam1 4 "$Yf_up"
    set_rpp "$inp" lam2 3 "$Yi_down";  set_rpp "$inp" lam2 4 "$Yf_down"
    set_rpp "$inp" lam3 1 "$Xi_left";  set_rpp "$inp" lam3 2 "$Xf_left"
    set_rpp "$inp" lam4 1 "$Xi_right"; set_rpp "$inp" lam4 2 "$Xf_right"

    degRot=$(echo "${deg} * -1" | bc -l)
    set_card "$inp" "ROT-DEFI" 3 "$degRot" 2
    set_card "$inp" "SOURCE" 7 "${E}MeV" 1
    set_card "$inp" "START" 1 "${primaries_per_CPU}" 1

    cp imgs/CT.vxl "$SIM/CT.vxl"
    cp starter_kit/simkit${E}MeV/* "$SIM/"
    ( cd "$SIM" && bash crea_input_homemade.sh "EF70mm${E}MeV.inp" "${INPs}" > trash.out )
    echo "FLUKA input files created for ${E}MeV (${VAR})"

    for run in $(seq 1 "${INPs}"); do
      ( cd "$SIM" && nohup "${FLUPRO}/flutil/rfluka" "run_${run}R.inp" -e "fluka_EF_${E}MeV.exe" -N0 -M1 \
          > "logs/run_${run}.log" 2>&1 ) &
      W_PIDS+=($!); W_NAMES+=("${E}MeV_${VAR}_run${run}")
    done
  done

  echo "Waiting ${VAR} simulations"
  wait_all "Simulazioni ${VAR}"       # <-- plus si ferma qui finche' NON ha finito

  # ---- 2. bnn -> mhd ----
  for E in "${energies[@]}"; do
    SIM="sim${E}MeV_${VAR}"
    for run in $(seq 1 "${INPs}"); do
      ( cd "$SIM" && nohup ../starter_kit/bnn2mhd "run_${run}R001_fort.23" "dose_tot_run_${run}.mhd" -Gy \
          > trash.out 2>&1 ) &
      W_PIDS+=($!); W_NAMES+=("bnn2mhd_${E}MeV_${VAR}_run${run}")
    done
  done
  wait_all "Conversione bnn->mhd ${VAR}"

  # ---- 3. media dei run ----
  for E in "${energies[@]}"; do
    SIM="sim${E}MeV_${VAR}"
    (
      cd "$SIM"
      check=(dose_tot_run*)
      if (( ${#check[@]} == 1 )); then cp "${check[0]}" dose_tot_run_copy.mhd; fi
      nohup python3 ../starter_kit/map_combine.py -avg dose_tot_run* > trash.out 2>&1
    ) &
    W_PIDS+=($!); W_NAMES+=("map_combine_${E}MeV_${VAR}")
  done
  wait_all "Combinazione mappe ${VAR}"

  # ---- 4. DVH, impulsi, de-crop / de-shift / contro-rotazione, isodose ----
  cropZ=$(<imgs/cropZ.txt)
  invSHx=$(<imgs/invSHx.txt); invSHy=$(<imgs/invSHy.txt); invSHz=$(<imgs/invSHz.txt)

  for E in "${energies[@]}"; do
    SIM="sim${E}MeV_${VAR}"
    cd "$SIM"
    cp ../starter_kit/mhd_smooth.x ./
    ./mhd_smooth.x avg.mhd -o avg_smooth.mhd > trash.out
    mv avg_smooth.mhd DOSE_${E}MeV_GRID.mhd
    python3 ../starter_kit/map_refill.py -ct ../imgs/CT_plan.mhd -dose DOSE_${E}MeV_GRID.mhd -out DOSE_${E}MeV.mhd
    python3 ../starter_kit/map_astype.py DOSE_${E}MeV.mhd float32

    echo "Creating dose map and DVH for ${E}MeV (${VAR})"
    mkdir -p DVH${E}MeV_FLASH DVH${E}MeV_CONV
    factor_FLASH="kFLASH_${E}MeV"; factor_CONV="kCONV_${E}MeV"

    cp DOSE_${E}MeV.mhd DOSE_${E}MeV_FLASH.mhd
    mv DOSE_${E}MeV.mhd DOSE_${E}MeV_CONV.mhd
    python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_FLASH.mhd -multiplier "${!factor_FLASH}"
    cp DOSE_${E}MeV_FLASH.mhd DOSE_${E}MeV_1pulses_FLASH.mhd
    python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_CONV.mhd -multiplier "${!factor_CONV}"
    cp DOSE_${E}MeV_CONV.mhd DOSE_${E}MeV_1pulses_CONV.mhd
    ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_FLASH.mhd -type float -fileLabel ${E}MeV_FLASH -dir DVH${E}MeV_FLASH --binningDVH 10000 > trash.out
    ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_CONV.mhd  -type float -fileLabel ${E}MeV_CONV  -dir DVH${E}MeV_CONV  --binningDVH 10000 > trash.out

    # stessi argomenti del caso nominale (incluso --preD_percent)
    pulses_FLASH=$(python3 ../starter_kit/findPulses.py DVH${E}MeV_FLASH/PTV_plan${E}MeV_FLASH.txt ${preD} ${preV} --preD_percent ${Pdose})
    pulses_CONV=$(python3 ../starter_kit/findPulses.py DVH${E}MeV_CONV/PTV_plan${E}MeV_CONV.txt   ${preD} ${preV} --preD_percent ${Pdose})

    python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_FLASH.mhd -multiplier "${pulses_FLASH}"
    python3 ../starter_kit/map_rescale.py DOSE_${E}MeV_CONV.mhd  -multiplier "${pulses_CONV}"
    ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_FLASH.mhd -type float -fileLabel ${E}MeV_FLASH -dir DVH${E}MeV_FLASH --binningDVH 1000 > trash.out
    ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_CONV.mhd  -type float -fileLabel ${E}MeV_CONV  -dir DVH${E}MeV_CONV  --binningDVH 1000 > trash.out
    python3 ../starter_kit/plotDVH.py -label1 ${E}MeV_FLASH -dir1 DVH${E}MeV_FLASH -roi PTV_plan "${ROIs_plan[@]}" > trash.out
    python3 ../starter_kit/plotDVH.py -label1 ${E}MeV_CONV  -dir1 DVH${E}MeV_CONV  -roi PTV_plan "${ROIs_plan[@]}" > trash.out
    python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV_FLASH/PTV_plan${E}MeV_FLASH.txt -preD ${preD} > DVH${E}MeV_FLASH/DVH_PTV_FLASH.out
    python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV_CONV/PTV_plan${E}MeV_CONV.txt   -preD ${preD} > DVH${E}MeV_CONV/DVH_PTV_CONV.out
    cp DOSE_${E}MeV_FLASH.mhd DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH.mhd
    cp DOSE_${E}MeV_CONV.mhd  DOSE_${E}MeV_${pulses_CONV}pulses_CONV.mhd
    DVHdirs+=("${SIM}/DVH${E}MeV_FLASH" "${SIM}/DVH${E}MeV_CONV")   # locale: niente mix plus/minus
    cd ..

    echo "Ripristino le dimensioni originali della mappa di dose per ${E}MeV (${VAR})..."
    P=([FLASH]="$pulses_FLASH" [CONV]="$pulses_CONV")
    for MODE in FLASH CONV; do
      p="${P[$MODE]}"
      base="${SIM}/DOSE_${E}MeV_${p}pulses_${MODE}"
      python3 starter_kit/mhd_decrop.py -i "${base}.mhd" -o "${base}_PADDED.mhd" -cz "$cropZ"
      python3 starter_kit/map_shift.py "${base}_PADDED.mhd" "$invSHx" "$invSHy" "$invSHz" "${base}_PADDED_deSH.mhd"
      python3 starter_kit/RotDoseQuat.py "${base}_PADDED_deSH.mhd" "imgs/rotation_matrix.txt" "${base}_FINAL.mhd"
      mv "${base}.mhd" "${base}_BEV.mhd"
      python3 starter_kit/create_isodosemap.py -i "${base}_BEV.mhd"   -o "${base}_BEV_isodosemap.mhd"
      python3 starter_kit/create_isodosemap.py -i "${base}_FINAL.mhd" -o "${base}_FINAL_isodose.mhd"
  done
 

    

    echo "---------------------------------------------------------------"
    echo "          PULSES ${E}MeV FLASH (${VAR}): ${pulses_FLASH}"
    echo "          PULSES ${E}MeV CONV  (${VAR}): ${pulses_CONV}"
    echo "---------------------------------------------------------------"
  done

  python3 starter_kit/GetDVHPlot.py "${DVHdirs[@]}" --xunit "cGy" --out "./DVH_ALL_${VAR}.png" --label-mode file --min-volume 0.01

  # ---- 5. pulizia e rinomina (indice AUTO corretto per ogni energia) ----
  for E in "${energies[@]}"; do
    SIM="sim${E}MeV_${VAR}"
    rm -f "$SIM"/*PADDED.mhd "$SIM"/*PADDED_deSH.mhd "$SIM"/*PADDED.raw "$SIM"/*PADDED_deSH.raw
    rm -f "$SIM"/*fort.23 "$SIM"/*fort.19 "$SIM"/*.log "$SIM"/*R00*.out
    rm -f "$SIM"/DOSE_${E}MeV_{1pulses_CONV,1pulses_FLASH,CONV,FLASH,GRID,orig}.mhd
    rm -f "$SIM"/dose_tot*.mhd "$SIM"/avg.mhd
    mv "$SIM" "AUTO_sim${E}MeV_${AUTO_IDX[$E]}_${VAR}"
  done
        # ---- 6. PNG mappe di isodose per apertura +2 mm / -2 mm ----
python3 starter_kit/ptv_isocenter.py imgs/PTV_plan.mhd > isocenter.out
ISO=$(cat isocenter.out)
[[ -n "$ISO" ]] || { echo "Isocentro non calcolato"; exit 1; }

export FRED_MATPLOTLIB_BACKEND=Agg

shopt -s nullglob  
for E in "${energies[@]}"; do
    AUTO="AUTO_sim${E}MeV_${AUTO_IDX[$E]}_${VAR}"

    for MODE in FLASH CONV; do
      # Seleziona la mappa BEV dell'energia e della modalità corrente
      maps=( "${AUTO}"/DOSE_${E}MeV_*pulses_${MODE}_BEV_isodosemap.mhd )

      if [[ ${#maps[@]} -eq 0 ]]; then
        echo "Attenzione: mappa isodose ${MODE} non trovata in ${AUTO}"
        continue
      elif [[ ${#maps[@]} -gt 1 ]]; then
        echo "Attenzione: trovate più mappe isodose ${MODE} in ${AUTO}; uso ${maps[0]}"
      fi

      MAP="${maps[0]}"
      OUT="${AUTO}/ISODOSE_BEV_${E}_${MODE}_${AUTO_IDX[$E]}_${VAR}.png"

      rm -f plot.png
      python3 -W ignore starter_kit/map_viewer_isodose.py \
        -CT imgs/CT_plan.mhd \
        -map "$MAP" \
        -Dgoal 1 -ISO $ISO -png > /dev/null

      if [[ -f plot.png ]]; then
        mv plot.png "$OUT"
        echo "Creato $OUT"
      else
        echo "Errore: il viewer non ha prodotto il PNG per $MAP"
      fi
    done
  done

  echo "=== ${VAR} completata ==="

}
# aperture +/-2 mm per TUTTE le energie
for E in "${energies[@]}"; do
  python3 starter_kit/vary_aperture.py "AUTO_sim${E}MeV_${AUTO_IDX[$E]}/aperture.out" \
      --delta-mm 2 --outdir "AUTO_sim${E}MeV_${AUTO_IDX[$E]}"
done

run_variant plus     # lancia, ATTENDE la fine, post-processing
run_variant minus    # parte solo quando plus e' completamente concluso

# ============================================================================
# CURRENT INVOCATION MANIFEST: NO SCAN OF OLD AUTO_sim DIRECTORIES
# ============================================================================
# AUTO_IDX[E] was assigned by this exact invocation, at AUTO_sim creation.
CURRENT_AUTO_SIMS=()
CURRENT_SIM_ARGS=()
for E in "${energies[@]}"; do
  base="AUTO_sim${E}MeV_${AUTO_IDX[$E]}"
  for sim in "$base" "${base}_plus2mm" "${base}_minus2mm"; do
    [[ -d "$sim" && ! -L "$sim" ]] || { echo "ERROR: Current run simulation missing: $sim" >&2; exit 1; }
    CURRENT_AUTO_SIMS+=("$sim")
    CURRENT_SIM_ARGS+=(--sim "$sim")
  done
done
[[ ${#CURRENT_AUTO_SIMS[@]} -gt 0 ]] || { echo "ERROR: No current run AUTO_sim folders" >&2; exit 1; }
# Unique, run-scoped report. Never reused and never overwritten.
PTV_RUN_DIR="$PWD/PTV_coverage_comparison_$(date +%Y%m%d_%H%M%S)_$$"
PTV_RUN_TABLE="$PTV_RUN_DIR/PTV_coverage_table.csv"
printf 'CURRENT FULL CHAIN RUN: %s\n' "${CURRENT_AUTO_SIMS[*]}"

# ============================================================================
# FINAL PTV CHECK, YES-ONLY COVERAGE SUMMARIES, CURRENT-RUN YES COLLECTION
# ============================================================================
echo ""
echo "=============== CURRENT RUN: FINAL PTV COVERAGE CHECK / RESULTS ==============="
python3 starter_kit/ptv_final_check.py "$SETUP_VERIFICA_PTV" "$PTV_RUN_DIR" "${CURRENT_AUTO_SIMS[@]}"


echo "Current-run PTV comparison, YES-only summaries, and results collection completed."

# ============================================================================
# OAR ANALYSIS: ONLY CURRENT-RUN AUTO_sim, EXACT MATCH TO THIS RUN'S PTV TABLE
# ============================================================================
# The OARs match the same ROIs_plan passed to ComputeDVH for this patient.
# No historical AUTO_sim or historical results/ directory is re-analyzed.
echo ""
echo "=============== CURRENT RUN: OAR DMAX / D(V) COMPARISON ==============="
OAR_PLAN_ARGS=()
for roi in "${ROIs_plan[@]}"; do
  OAR_PLAN_ARGS+=(--oar-plan "imgs/${roi}.mhd")
done
python3 starter_kit/oar_comparison.py  --root "$PWD" --ptv-table "$PTV_RUN_TABLE" "${CURRENT_SIM_ARGS[@]}" "${OAR_PLAN_ARGS[@]}" --prescription-gy "${preD}" --prescription-level-percent "${Pdose}"

echo "Current-run OAR comparison and YES-only per-plan table_oar.txt completed."
echo ""
echo "================ FINAL PTV COMPARISON (THIS RUN ONLY) ================"
cat "$PTV_RUN_DIR/PTV_coverage_comparison.txt"
echo "================ FINAL OAR COMPARISON (THIS RUN ONLY) ================"
echo "Dmax and D(v) minima, including different winners, are printed above."
echo "OAR text comparison is saved in oar_analysis/run_*/OAR_comparison_table.txt"
echo "PTV YES plan reports: results/<AUTO_sim>_<MODE>/summary_coverage.txt and table_oar.txt"
echo "Historical AUTO_sim and results/ folders have not been re-analyzed or overwritten."
# ============================================================================
# FIX: OPTIMIZATION OF CONV PULSE NUMBERS — CURRENT FULL CHAIN RUN ONLY
# ============================================================================
# This block starts ONLY after the Monte Carlo, PTV and OAR post-processing.
# CURRENT_SIM_ARGS was populated from AUTO_IDX immediately after the simulation
# phase, and contains only the original, plus2mm and minus2mm simulations
# created in THIS invocation (3 cases for one energy; 6 for both energies).
#
# pulse_number_optimizer.sh optimizes only CONV; FLASH is NOT optimized.
# It creates new, unique output folders under pulse_optimization_runs/run_*
# and reads dose maps, ROI and DVHs without overwriting original .mhd/.raw,
# dose, isodose, DVH or AUTO_sim content.
echo "================ CONV PULSE OPTIMIZATION — CURRENT RUN ================"

if python3 starter_kit/pulse_number_optimizer.py \
    --setup "$SETUP_VERIFICA_PTV" \
    "${CURRENT_SIM_ARGS[@]}" \
    --require-complete-run; then

    echo "CONV pulse optimization completed."
    echo "Reports: pulse_optimization_runs/run_*/CONFRONTO_FINALE_IMPULSI.txt"

else
    status=$?
    echo "WARNING: CONV pulse optimizer returned an error (exit code $status)." >&2
    echo "Review the newest pulse_optimization_runs/run_*/errori.csv and terminal log." >&2
    echo "Original AUTO_sim and dose/isodose maps have not been changed by this block." >&2
    exit "$status"
fi


echo "Enjoy! Premi ctrl+c per uscire"
