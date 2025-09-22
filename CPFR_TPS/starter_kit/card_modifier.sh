#!/usr/bin/env bash
# Uso:
#   ./fix_fluka.sh input.inp output.inp ROT-DEFI 2 30 [occ]
#                                                ^WHAT(i) ^nuovo  ^occorrenza (default 1)

in="$1"
out="$2"
card="$3"         # es. ROT-DEFI
idx="$4"          # 1..7 (WHAT(i) or SDUM)
newval="$5"       # es. 30, -45.5, 0, ecc.
occ="${6:-1}"     # quale occorrenza della card (1=prima; 2=seconda; ...)

# Controlli veloci
if ! [[ "$idx" =~ ^[1-7]$ ]]; then
  echo "Errore: WHAT index deve essere 1..7 (1..6=WHAT(i), 7=SDUM)" >&2; exit 1
fi
if ! [[ "$occ" =~ ^[0-9]+$ ]] || [ "$occ" -lt 1 ]; then
  echo "Errore: occorrenza deve essere un intero >= 1" >&2; exit 1
fi

gawk -v card="$card" -v idx="$idx" -v newval="$newval" -v occ="$occ" '
function trim(s){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
function clip10(s){ return (length(s) > 10 ? substr(s, 1, 10) : s) }

BEGIN { seen = 0 }

{
  # Leggi solo il nome card dalle colonne canoniche (1..10)
  name = trim(substr($0, 1, 10))

  if (name == card) {
    # È una riga della card richiesta → conta occorrenze
    seen++

     if (seen == occ) {
       # SOLO sulla riga target usiamo FIELDWIDTHS per i campi a 10 char
       FIELDWIDTHS = "10 10 10 10 10 10 10 10"
       # forza il ricalcolo dei campi $1..$8 secondo FIELDWIDTHS per questo record
       $0 = $0


      # Con FIELDWIDTHS attivo, $1..$8 sono i campi a larghezza fissa
      for (i=1; i<=8; i++) {
        f[i] = $i
        c[i] = trim(f[i])
      }

      # idx 1..6 → WHAT(i); idx 7 → SDUM
      if (idx >= 1 && idx <= 6) {
        c[1+idx] = newval
      } else if (idx == 7) {
        c[8] = newval
      }

      # Tronca ogni campo a max 10 caratteri
      for (i=1; i<=8; i++) c[i] = clip10(c[i])

      # Stampa rispettando l’allineamento FLUKA:
      # - nome e SDUM a sinistra (%-10s)
      # - WHAT(i) a destra (%10s)
      printf "%-10s%10s%10s%10s%10s%10s%10s%-10s\n",
             c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8]
      next
    }
  }

  # Tutte le altre righe (e le occorrenze non target) restano INALTERATE
  print $0
}
END {
  if (seen < occ) {
    printf("Avviso: card \"%s\" occorrenza %d non trovata (trovate %d)\n", card, occ, seen) > "/dev/stderr"
  }
}
' "$in" > "$out"
