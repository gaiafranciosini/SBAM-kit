#!/bin/bash

# Lancia: ./run_fluka.sh

INP_FILE="/loc_work2/mperri/Pancreas_PZ3/Patch_Att/Simulazione_celle/simulazione_field1_field1_att_mlc.inp"
BATCH_SIZE=15

# Controlli iniziali
if [ ! -f "$INP_FILE" ]; then
    echo "Errore: il file '$INP_FILE' non esiste."
    exit 1
fi

# Raccoglie tutti i .txt nella directory corrente
TXT_FILES=(*.txt)
TOTAL=${#TXT_FILES[@]}
echo "Trovati $TOTAL file .txt da elaborare, batch da $BATCH_SIZE."
echo ""

# ── FASE 1: run FLUKA a batch da 15 ──────────────────────────────────────────

i=0
while [ $i -lt $TOTAL ]; do

    echo "============================================================"
    echo "Avvio batch: file $((i+1)) - $((i+BATCH_SIZE < TOTAL ? i+BATCH_SIZE : TOTAL)) di $TOTAL"
    echo "============================================================"

    PIDS=()
    CELL_NUMS=()
    NEW_INPS=()

    # Lancia fino a BATCH_SIZE job in parallelo
    for (( j=i; j<i+BATCH_SIZE && j<TOTAL; j++ )); do

        TXT_NAME="${TXT_FILES[$j]}"         # es. 12.txt
        CELL_NUM="${TXT_NAME%.txt}"         # es. 12
        NEW_INP="simulazione_cella${CELL_NUM}.inp"

        echo "  Creo $NEW_INP con SOURCE → $TXT_NAME"

        # Copia l'originale e sostituisce il .txt nella riga SOURCE
        sed "s/[^ ]*\.txt/${TXT_NAME}/" "$INP_FILE" > "$NEW_INP"

        if [ $? -ne 0 ]; then
            echo "  [ERRORE] sed fallito per cella ${CELL_NUM}, salto."
            continue
        fi

        # Lancia FLUKA in background
        $FLUPRO/flutil/rfluka "$NEW_INP" -N0 -M1 -e SAFEST_PIN.exe &

        PIDS+=($!)
        CELL_NUMS+=("$CELL_NUM")
        NEW_INPS+=("$NEW_INP")

    done

    # Aspetta che tutti i job del batch finiscano
    echo ""
    echo "  Attendo il completamento del batch..."
    for k in "${!PIDS[@]}"; do
        PID=${PIDS[$k]}
        CELL_NUM=${CELL_NUMS[$k]}
        wait "$PID"
        EXIT_CODE=$?
        if [ $EXIT_CODE -ne 0 ]; then
            echo "  [ERRORE] FLUKA fallito per cella ${CELL_NUM} (exit code $EXIT_CODE)"
        else
            echo "  OK: cella ${CELL_NUM} completata."
        fi
    done

    # Cancella i .inp temporanei del batch
    echo ""
    echo "  Cancello i .inp temporanei..."
    for NEW_INP in "${NEW_INPS[@]}"; do
        rm -f "$NEW_INP"
        echo "  Cancellato $NEW_INP"
    done

    i=$((i + BATCH_SIZE))

done

echo ""
echo "============================================================"
echo "Tutti i run completati."
echo ""

# ── FASE 2: conversione fort.21 → .mhd ──────────────────────────────────────

echo "Avvio conversione fort.21 → .mhd ..."
echo "============================================================"

CONVERTED=0
FAILED=0

for FORT_FILE in simulazione_cella*001_fort.21; do

    if [ ! -f "$FORT_FILE" ]; then
        echo "Nessun file fort.21 trovato."
        break
    fi

    CELL_NUM=$(echo "$FORT_FILE" | sed 's/simulazione_cella\([0-9]*\)001_fort\.21/\1/')
    MHD_OUT="dose${CELL_NUM}.mhd"

    echo "  Conversione: $FORT_FILE → $MHD_OUT"
    bnn2mhd "$FORT_FILE" "$MHD_OUT"

    if [ $? -ne 0 ]; then
        echo "  [ERRORE] Conversione fallita per $FORT_FILE"
        FAILED=$((FAILED + 1))
    else
        echo "  OK: $MHD_OUT creato."
        CONVERTED=$((CONVERTED + 1))
    fi

done

echo "============================================================"
echo "Conversione completata: $CONVERTED ok, $FAILED falliti."
