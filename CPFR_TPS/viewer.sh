#!/bin/bash
# Script per la SOLA visualizzazione delle mappe di dose finali

file_energie="$1"
sim="$2"
CT="$3"
PTV="$4"
MARKER="$5"

if [ "$#" -lt 5 ]; then
  echo "Uso: $0 <config_file> <sim_id> <CT> <PTV> <MARKER> [ROI1 ROI2 ...]"
  exit 1
fi

shift 5

# ======================================================================
# --- LETTURA ENERGIE TRAMITE GREP ---
riga_energia=$(grep "^ENERGY:" "setup.out" | head -n 1)
valori_energia="${riga_energia#ENERGY: }"

# Rimuove eventuali ritorni a capo Windows (\r) per sicurezza
valori_energia=$(echo "$valori_energia" | tr -d '\r')

read -ra energies <<< "$valori_energia"

if [ ${#energies[@]} -eq 0 ]; then
    echo "Errore: Nessuna riga 'ENERGY:' trovata nel file setup.out o nessun valore specificato."
    exit 1
fi

echo "Energie trovate: ${energies[@]}"
# ======================================================================

# Array che conterrà i nomi puliti delle ROI
ROIs=()
for path in "$@"; do
  name="${path##*/}"
  name="${name%.mhd}"
  ROIs+=("$name")
done

# Ricostruiamo la lista completa dei file ROI (.mhd) originali passati in input
ALL_ROI_PATHS=("$PTV")
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
# --- INIZIO VISUALIZZAZIONE MAPPE FINALI ---
echo "Avvio visualizzatore per le mappe di dose finali..."

# Apriamo i visualizzatori per ogni energia
for E in "${energies[@]}"; do

    # Troviamo i file corretti usando l'asterisco per scavalcare il numero esatto di impulsi
    DOSE_CONV_FILE=$(ls AUTO_sim${E}MeV_${sim}/DOSE_${E}MeV_*pulses_CONV_FINAL.mhd 2>/dev/null | head -n 1)
    DOSE_FLASH_FILE=$(ls AUTO_sim${E}MeV_${sim}/DOSE_${E}MeV_*pulses_FLASH_FINAL.mhd 2>/dev/null | head -n 1)

    # 1. Avviamo il viewer per la dose CONV
    if [[ -f "$DOSE_CONV_FILE" ]]; then
        echo "Apertura visualizzatore CONV per ${E}MeV (File: $DOSE_CONV_FILE)..."
        TITOLO_CONV="Dose_CONV_${E}_MeV"
        
        nohup python3 starter_kit/viewer_zoom.py \
            --ct "$CT" \
            --roi "${ALL_ROI_PATHS[@]}" \
            --labels "${ALL_ROI_LABELS[@]}" \
            --dose "$DOSE_CONV_FILE" \
            --dose-alpha 0.5 \
            --dose-cmap "jet" \
            --title "$TITOLO_CONV" > "trash_viewer_conv_${E}.out" 2>&1 &
    else
        echo "Attenzione: Mappa CONV non trouvata per ${E}MeV."
    fi

    # 2. Avviamo il viewer per la dose FLASH
    if [[ -f "$DOSE_FLASH_FILE" ]]; then
        echo "Apertura visualizzatore FLASH per ${E}MeV (File: $DOSE_FLASH_FILE)..."
        TITOLO_FLASH="Dose_FLASH_${E}_MeV"
        
        nohup python3 starter_kit/viewer_zoom.py \
            --ct "$CT" \
            --roi "${ALL_ROI_PATHS[@]}" \
            --labels "${ALL_ROI_LABELS[@]}" \
            --dose "$DOSE_FLASH_FILE" \
            --dose-alpha 0.5 \
            --dose-cmap "hot" \
            --title "$TITOLO_FLASH" > "trash_viewer_flash_${E}.out" 2>&1 &
    else
        echo "Attenzione: Mappa FLASH non trovata per ${E}MeV."
    fi

done
# --- FINE VISUALIZZAZIONE MAPPE FINALI ---
# ======================================================================

echo "Enjoy! Le finestre si stanno aprendo in background."
