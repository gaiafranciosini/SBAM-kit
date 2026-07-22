#!/bin/bash
#++++++++++++++++++++++++++++++++++
set -e
set -o pipefail
trap 'echo "Errore: comando \"$BASH_COMMAND\" fallito alla riga $LINENO".  Premi ctrl+c per uscire;  exit 1' ERR
#set -x

file="$1"
CT="$2"
PTV="$3"
MARKER="$4"
FIELD_ID="$5"

if [ "$#" -lt 5 ]; then
  echo "Use: FILE CT PTV MARKER FIELD_ID [ROI1 ROI2 ROI3 ROI4 ...]"
  exit 1
fi

echo "$1 $2 $3 $4 $5" > input_log.out
shift 5
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
#pulses=""
preD=""
preV=""
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
#        "PULSES:"*)
#            pulses=(${line#PULSES: })
#            ;;
        "PERCENTAGE PRESCRIPTION DOSE:"*)
            preD=${line#PERCENTAGE PRESCRIPTION DOSE: }
            ;;
        "PERCENTAGE PRESCRIPTION VOLUME:"*)
            preV=${line#PERCENTAGE PRESCRIPTION VOLUME: }
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


#CALIBRATION FACTORS (primaries/pulse)
kFLASH_7MeV=9.382238805934787e+11
kCONV_7MeV=1.42155133423254e+10
kFLASH_9MeV=7.106140035174675e+11
kCONV_9MeV=1.07668788411738e+10

for E in "${energies[@]}"; do

  if [ -d "sim${E}MeV" ]; then
      rm -rf "OLD_sim${E}MeV"
      mv sim${E}MeV OLD_sim${E}MeV
  fi

  mkdir sim${E}MeV
  cp files/EF70mm${E}MeV.inp sim${E}MeV
  bash starter_kit/card_modifier.sh sim${E}MeV/EF70mm${E}MeV.inp sim${E}MeV/tmp.inp "START" 1 "${primaries_per_CPU}" 1
  mv sim${E}MeV/tmp.inp sim${E}MeV/EF70mm${E}MeV.inp

done

# Controllo esistenza file ManualFieldSize
manual_file="ManualFieldSize_${FIELD_ID}.out"
if [ ! -f "$manual_file" ]; then
    echo "Errore: il file $manual_file non esiste!"
    exit 1
fi

while IFS=':' read -r key value; do
    [[ -z "$key" || -z "$value" ]] && continue
    # Pulizia base per evitare problemi di formattazione nascosti (spazi o carriage return)
    key=$(echo "$key" | tr -d '\r' | xargs)
    value=$(echo "$value" | tr -d '\r' | xargs)
    declare "${key}=${value}"
done < "$manual_file"

# top_cm
# bottom_cm
# left_cm
# right_cm
# opening_top_bottom_cm
# opening_left_right_cm
# angle_deg

# ATTENTION: BEV x-axis is reversed in FLUKA, left and right must be inverted
FLUKA_left_cm=${right_cm}
FLUKA_right_cm=${left_cm}


width=8.0
height=4.0
angles=(0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85)



Xgrid=8 #cm     -Xgrid/2 |------0------| +Xgrid/2
Ygrid=8 #cm     -Ygrid/2 |------0------| +Ygrid/2
Zgrid=6 #cm            0 |-------------| +Zgrid


Ws=()
Hs=()


for E in "${energies[@]}"; do


# 4. MODIFY FLUKA file - field size

  python3 starter_kit/ManualGetSlitMargins.py ${width} ${height} ${FLUKA_right_cm} ${FLUKA_left_cm} ${top_cm} ${bottom_cm}  > "sim${E}MeV/aperture.out"
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



degRot=$(echo "${rotation_angle_deg} * -1" | bc -l)

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
  (( ++i ))
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
  (( ++i ))
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
  nohup python3 ../starter_kit/mhd_combine.py -avg dose_tot_run* > trash.out 2>&1 &
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
DVHdirs=()
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
    print ""
  }')

#  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -divider ${max}
#  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV.mhd -multiplier 50
#
  echo "Creating dose map and DVH for ${E}MeV"
  mkdir -p DVH${E}MeV_FLASH
  mkdir -p DVH${E}MeV_CONV
  factor_FLASH="kFLASH_${E}MeV"
  factor_CONV="kCONV_${E}MeV"

  cp DOSE_${E}MeV.mhd DOSE_${E}MeV_FLASH.mhd
  mv DOSE_${E}MeV.mhd DOSE_${E}MeV_CONV.mhd
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV_FLASH.mhd -multiplier "${!factor_FLASH}"   # 1pulse
  cp DOSE_${E}MeV_FLASH.mhd DOSE_${E}MeV_1pulses_FLASH.mhd
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV_CONV.mhd -multiplier "${!factor_CONV}"   # 1pulse
  cp DOSE_${E}MeV_CONV.mhd DOSE_${E}MeV_1pulses_CONV.mhd
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_FLASH.mhd -type float -fileLabel ${E}MeV_FLASH -dir DVH${E}MeV_FLASH --binningDVH 10000 > trash.out
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_CONV.mhd -type float -fileLabel ${E}MeV_CONV -dir DVH${E}MeV_CONV --binningDVH 10000 > trash.out
  pulses_FLASH=$(python3 ../starter_kit/findPulses.py DVH${E}MeV_FLASH/PTV_plan${E}MeV_FLASH.txt ${preD} ${preV} )
  pulses_CONV=$(python3 ../starter_kit/findPulses.py DVH${E}MeV_CONV/PTV_plan${E}MeV_CONV.txt ${preD} ${preV} )
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV_FLASH.mhd -multiplier "${pulses_FLASH}"
  python3 ../starter_kit/mhd_rescale.py DOSE_${E}MeV_CONV.mhd -multiplier "${pulses_CONV}"
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_FLASH.mhd -type float -fileLabel ${E}MeV_FLASH -dir DVH${E}MeV_FLASH --binningDVH 1000  > trash.out
  ../starter_kit/ComputeDVH/ComputeDVH.py -Dgoal 1 -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -dose DOSE_${E}MeV_CONV.mhd -type float -fileLabel ${E}MeV_CONV -dir DVH${E}MeV_CONV --binningDVH 1000  > trash.out
  python3 ../starter_kit/plotDVH.py -label1 ${E}MeV_FLASH -dir1 DVH${E}MeV_FLASH -roi PTV_plan "${ROIs_plan[@]}" > trash.out
  python3 ../starter_kit/plotDVH.py -label1 ${E}MeV_CONV -dir1 DVH${E}MeV_CONV -roi PTV_plan "${ROIs_plan[@]}" > trash.out
  python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV_FLASH/PTV_plan${E}MeV_FLASH.txt > DVH${E}MeV_FLASH/DVH_PTV_FLASH.out
  python3 ../starter_kit/readDVH.py -DVH DVH${E}MeV_CONV/PTV_plan${E}MeV_CONV.txt > DVH${E}MeV_CONV/DVH_PTV_CONV.out
  cp DOSE_${E}MeV_FLASH.mhd DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH.mhd
  cp DOSE_${E}MeV_CONV.mhd DOSE_${E}MeV_${pulses_CONV}pulses_CONV.mhd
  DVHdirs+=("sim${E}MeV/DVH${E}MeV_FLASH")
  DVHdirs+=("sim${E}MeV/DVH${E}MeV_CONV")
#  nohup  python3 ../starter_kit/mhd_viewer_RayS.py DOSE_${E}MeV.mhd -CT ../imgs/CT_plan.mhd -roi ../imgs/PTV_plan.mhd "${all_rois_path[@]}" -png > trash.out &
#  pids+=("$!")
  cd ..

   echo "Ripristino le dimensioni originali della Mappa di Dose per ${E}MeV..."
cropZ=$(<imgs/cropZ.txt)
  # 1. De-crop (padding)
  python3 starter_kit/mhd_decrop.py -i sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH.mhd -o sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED.mhd -cz "$cropZ"
  python3 starter_kit/mhd_decrop.py -i sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV.mhd -o sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED.mhd -cz "$cropZ"

  # 2. De-shift
  invSHx=$(<imgs/invSHx.txt)
  invSHy=$(<imgs/invSHy.txt)
  invSHz=$(<imgs/invSHz.txt)
  #invSHx=$PTVx
  #invSHy=$PTVy
  #invSHz=$Z0crop

  python3 starter_kit/mhd_shift.py "sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED.mhd" "$invSHx" "$invSHy" "$invSHz" sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED_deSH.mhd
  python3 starter_kit/mhd_shift.py "sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED.mhd" "$invSHx" "$invSHy" "$invSHz" sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED_deSH.mhd
  awk '/ROTATION MATRIX/{flag=1; next} /Rotated image/{flag=0} flag' imgs/CT_ROT_matrix.txt | tr -d '[]' > imgs/rotation_matrix.txt

  # 3. Contro-rotazione
  python3 starter_kit/RotDoseQuat.py "sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_PADDED_deSH.mhd" "imgs/rotation_matrix.txt" "sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_FINAL.mhd"

echo "Dose Map CONV rotated and saved in sim${E}MeV/DOSE_${E}MeV_${pulses_FLASH}pulses_FLASH_FINAL.mhd"

python3 starter_kit/RotDoseQuat.py "sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_PADDED_deSH.mhd" "imgs/rotation_matrix.txt" "sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_FINAL.mhd"

echo "Dose Map CONV rotated and saved in sim${E}MeV/DOSE_${E}MeV_${pulses_CONV}pulses_CONV_FINAL.mhd"



  echo "---------------------------------------------------------------"
  echo "                   PULSES ${E}MeV FLASH: ${pulses_FLASH}       "
  echo "                   PULSES ${E}MeV CONV: ${pulses_CONV}         "
  echo "---------------------------------------------------------------"


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
echo "DVHdirs content:"
printf "%s\n" "${DVHdirs[@]}"

python3 starter_kit/GetDVHPlot.py "${DVHdirs[@]}" --xunit "cGy" --out "./DVH_ALL_${FIELD_ID}.png"

for E in ${energies[@]}; do
# Cartella di output univoca per evitare sovrascritture in run successivi
mv "sim${E}MeV" "MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}"
echo "L R T B = left right top bottom slit aperture in cm (ID: ${FIELD_ID})"
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/*PADDED.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/*PADDED_deSH.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/*PADDED.raw
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/*PADDED_deSH.rawi=1
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/DOSE_${E}MeV_1pulses_CONV.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/DOSE_${E}MeV_1pulses_FLASH.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/DOSE_${E}MeV_CONV.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/DOSE_${E}MeV_FLASH.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/DOSE_${E}MeV_GRID.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/DOSE_${E}MeV_orig.mhd
rm -f MANUAL_sim${E}MeV_L${FLUKA_left_cm}_R${FLUKA_right_cm}_T${top_cm}_B${bottom_cm}_cm_${FIELD_ID}/dose_tot*.mhd
done

echo "Enjoy! Premere ctrl+c per uscire"











