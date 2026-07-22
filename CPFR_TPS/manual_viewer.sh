#!/bin/bash
# Script per la SOLA visualizzazione delle mappe di dose finali con percorsi dinamici

file_energie="$1"
sim="$2"
CT="$3"
PTV="$4"
MARKER="$5"
file_setup="$6"
INPUT_MANUAL="ManualFieldSize_${sim}.out"

if [ "$#" -lt 6 ]; then
  echo "Uso: $0 <config_file> <CT> <PTV> <MARKER> [ROI1 ROI2 ROI3 ROI4 ...] <file_setup>"
  exit 1
fi

# Verifica esistenza file parametri manuali
if [ ! -f "$INPUT_MANUAL" ]; then
    echo "Errore: Il file $INPUT_MANUAL non è stato trovato."
    exit 1
fi

shift 6

# ======================================================================
# --- ESTRAZIONE PARAMETRI GEOMETRICI (da ManualFieldSize.out) ---
# Leggiamo i valori numerici dalle righe 3, 4, 5 e 6
# --- ESTRAZIONE PARAMETRI GEOMETRICI DINAMICA ---
# --- ESTRAZIONE PARAMETRI GEOMETRICI DINAMICA ---
TOP=$(grep -i "TOP" "$INPUT_MANUAL" | head -n 1 | awk '{print $2}' | tr -d '\r')
BOTTOM=$(grep -i "BOTTOM" "$INPUT_MANUAL" | head -n 1 | awk '{print $2}' | tr -d '\r')
LEFT=$(grep -i "LEFT" "$INPUT_MANUAL" | head -n 1 | awk '{print $2}' | tr -d '\r')
RIGHT=$(grep -i "RIGHT" "$INPUT_MANUAL" | head -n 1 | awk '{print $2}' | tr -d '\r')

# Costruiamo il suffisso per le cartelle
SUFFIX_DIR="_L${LEFT}_R${RIGHT}_T${TOP}_B${BOTTOM}_cm_${sim}"

# ======================================================================

# --- LETTURA ENERGIE TRAMITE GREP (da setup.out) ---
riga_energia=$(grep "^ENERGY:" "setup.out" | head -n 1)
valori_energia="${riga_energia#ENERGY: }"
read -ra energies <<< "$valori_energia"

if [ ${#energies[@]} -eq 0 ]; then
    echo "Errore: Nessuna riga 'ENERGY:' trovata nel file setup.out."
    exit 1
fi

echo "Energie trovate: ${energies[@]}"
echo "Geometria rilevata: L=$LEFT, R=$RIGHT, T=$TOP, B=$BOTTOM"

# ======================================================================
# --- GESTIONE ROI ---
ROIs=()
for path in "$@"; do
  name="${path##*/}"
  name="${name%.mhd}"
  ROIs+=("$name")
done

ALL_ROI_PATHS=("${PTV}")
ALL_ROI_LABELS=("PTV")

if [[ -f "$MARKER" ]]; then
    ALL_ROI_PATHS+=("$MARKER")
    ALL_ROI_LABELS+=("Marker")
fi

for roi in "${ROIs[@]}"; do
    if [[ -f "imgs/${roi}.mhd" ]]; then
        ALL_ROI_PATHS+=("imgs/${roi}.mhd")
        ALL_ROI_LABELS+=("$roi")
    elif [[ -f "${roi}.mhd" ]]; then
        ALL_ROI_PATHS+=("${roi}.mhd")
        ALL_ROI_LABELS+=("$roi")
    fi
done

# ======================================================================
# --- INIZIO VISUALIZZAZIONE MAPPE ---
echo "Avvio visualizzatore per le mappe di dose finali..."

for E in "${energies[@]}"; do

    # Costruiamo il nome della cartella specifica per questa energia
    # Esempio: MANUAL_sim7MeV_L 3.0_R 3.0_T 2.0_B 2.0_cm
    DIR_ENERGIA="MANUAL_sim${E}MeV${SUFFIX_DIR}"

    if [ ! -d "$DIR_ENERGIA" ]; then
        echo "Attenzione: Cartella $DIR_ENERGIA non trovata. Salto energia ${E}MeV."
        continue
    fi

    # 1. Troviamo i file corretti all'interno della cartella specifica
    DOSE_CONV_FILE=$(ls "${DIR_ENERGIA}/DOSE_${E}MeV_"*pulses_CONV_FINAL.mhd 2>/dev/null | head -n 1)
    DOSE_FLASH_FILE=$(ls "${DIR_ENERGIA}/DOSE_${E}MeV_"*pulses_FLASH_FINAL.mhd 2>/dev/null | head -n 1)

    # 2. Avviamo il viewer per la dose CONV
    if [[ -f "$DOSE_CONV_FILE" ]]; then
        echo "Apertura CONV: $DOSE_CONV_FILE"
        TITOLO_CONV="Dose_CONV_${E}_MeV"
        nohup python3 starter_kit/viewer_zoom_manual.py \
            --ct "$CT" \
            --roi "${ALL_ROI_PATHS[@]}" \
            --labels "${ALL_ROI_LABELS[@]}" \
            --dose "$DOSE_CONV_FILE" \
            --dose-alpha 0.5 \
            --dose-cmap "jet"  \
            --title "$TITOLO_CONV" > "trash_viewer_conv_${E}.out" 2>&1 &
    else
        echo "Mappa CONV non trovata in $DIR_ENERGIA"
    fi

    # 3. Avviamo il viewer per la dose FLASH
    if [[ -f "$DOSE_FLASH_FILE" ]]; then
        echo "Apertura FLASH: $DOSE_FLASH_FILE"
        TITOLO_FLASH="Dose_FLASH_${E}_MeV"
        nohup python3 starter_kit/viewer_zoom_manual.py \
            --ct "$CT" \
            --roi "${ALL_ROI_PATHS[@]}" \
            --labels "${ALL_ROI_LABELS[@]}" \
            --dose "$DOSE_FLASH_FILE" \
            --dose-alpha 0.5 \
            --dose-cmap "hot" \
            --title "$TITOLO_FLASH" > "trash_viewer_flash_${E}.out" 2>&1 &
    else
        echo "Mappa FLASH non trovata in $DIR_ENERGIA"
    fi

done

echo "Enjoy! Le finestre si stanno aprendo in background."
