#!/usr/bin/env bash
# Uso:
#   ./fix_fluka_free.sh input.inp output.inp RPP lam1 3 12.34 [occ]
#                                                 ^k (1..6) ^nuovo  ^occorrenza (default 1)
#   Sostituisce il k-esimo numero (tra i sei) della riga FREE format:
#   CARD OBJ v1 v2 v3 v4 v5 v6 [eventuale commento...]

set -euo pipefail

in="$1"
out="$2"
card="$3"       # es. RPP
obj="$4"        # es. lam1
k="$5"          # 1..6 (indice del valore da sostituire tra i sei)
newval="$6"     # nuovo valore (stringa numerica)
occ="${7:-1}"   # quale occorrenza della coppia (card,obj) modificare (default 1)

# Controlli input
if ! [[ "$k" =~ ^[1-6]$ ]]; then
  echo "Errore: l'indice k deve essere 1..6 (seleziona quale dei sei valori sostituire)" >&2
  exit 1
fi
if ! [[ "$occ" =~ ^[0-9]+$ ]] || [ "$occ" -lt 1 ]; then
  echo "Errore: occorrenza deve essere un intero >= 1" >&2
  exit 1
fi

awk -v CARD="$card" -v OBJ="$obj" -v K="$k" -v NEWVAL="$newval" -v OCC="$occ" '
BEGIN { seen = 0; OFS = " " }

{
  raw = $0  # linea originale (le altre righe verranno stampate identiche)

  # Se la riga ha almeno 2 campi, controlla card e obj
  if (NF >= 2 && $1 == CARD && $2 == OBJ) {
    seen++
    if (seen == OCC) {
      # Assicurati che ci siano almeno i 6 valori dopo i primi due campi
      # (Se non ci sono, li consideriamo vuoti e li aggiungiamo)
      n = NF
      num_needed = 8
      if (n < num_needed) n = num_needed

      # Copia in un array locale per sicurezza
      for (i=1; i<=NF; i++) f[i]=$i

      # Prepara valori mancanti come stringhe vuote
      for (i=NF+1; i<=n; i++) f[i]=""

      # Sostituisci il K-esimo valore tra i sei numeri (campi 3..8)
      target_index = 2 + K
      f[target_index] = NEWVAL

      # Ricostruisci la riga in FREE format: CARD OBJ v1..v6 (+ eventuali campi extra)
      out = f[1] OFS f[2]
      for (i=3; i<=8; i++) out = out OFS f[i]
      if (NF > 8) {
        for (i=9; i<=NF; i++) out = out OFS f[i]
      }
      print out
      next
    }
  }

  # Tutte le altre righe (e le occorrenze non target) restano INALTERATE
  print raw
}
END {
  if (seen < OCC) {
    printf("Avviso: riga FREE con card \"%s\" e oggetto \"%s\" occorrenza %d non trovata (trovate %d)\n",
           CARD, OBJ, OCC, seen) > "/dev/stderr"
  }
}
' "$in" > "$out"
