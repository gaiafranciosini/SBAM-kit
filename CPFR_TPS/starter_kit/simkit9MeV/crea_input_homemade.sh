#!/bin/bash

# Uso: ./replica_randomiz.sh file_input N
# Dove:
#   file_input = file .inp di partenza (contiene la carta RANDOMIZ)
#   N          = numero di repliche da creare

if [ "$#" -ne 2 ]; then
    echo "Uso: $0 file_input N"
    exit 1
fi

file_input="$1"
N="$2"

if [ ! -f "$file_input" ]; then
    echo "Errore: file $file_input non trovato."
    exit 1
fi

for ((i=1; i<=N; i++)); do
    file_output="run_${i}R.inp"
    new_number=$((1000 + i))

    # Copia il file di input
    cp "$file_input" "$file_output"

    # Modifica SOLO WHAT(2) (colonne 21–30) della carta RANDOMIZ
    # mantenendo TUTTO il resto identico (nome carta, WHAT(1), WHAT(3–6), SDUM...)
    awk -v num="$new_number" '
    BEGIN { OFS="" }
    {
        # Controlla se nei primi 10 caratteri c’è "RANDOMIZ"
        if (substr($0,1,10) ~ /^RANDOMIZ/) {
            name = substr($0,  1,10)         # posizione 1–10
            w1   = substr($0, 11,10)         # posizione 11–20 (WHAT(1))
            w2   = sprintf("%10d", num)      # posizione 21–30 (nuovo WHAT(2))
            rest = substr($0, 31)            # dal carattere 31 in poi (WHAT(3–6)+SDUM)

            print name, w1, w2, rest
        } else {
            print $0
        }
    }' "$file_output" > tmp && mv tmp "$file_output"

    echo "Creato: $file_output con RANDOMIZ WHAT(2) = $(printf '%d' "$new_number")"
done
