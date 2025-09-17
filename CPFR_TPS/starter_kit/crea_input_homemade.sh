#!/bin/bash

# Controlla se ci sono almeno 2 argomenti
if [ "$#" -ne 2 ]; then
    echo "Uso: $0 file_input N"
    exit 1
fi

file_input="$1"
N="$2"

# Controlla se il file esiste
if [ ! -f "$file_input" ]; then
    echo "Errore: file $file_input non trovato."
    exit 1
fi

# Loop per creare N copie
for ((i=1; i<=N; i++)); do
    file_output="EF70mm_run${i}R.inp"
    new_number=$((1000 + i))

    # Copia il file
    cp "$file_input" "$file_output"

    # Modifica la riga che inizia con RANDOMIZ
    sed -i -E "s/^(RANDOMIZ.*)1000/\1${new_number}/" "$file_output"

    echo "Creato: $file_output con RANDOMIZ $new_number"
done
