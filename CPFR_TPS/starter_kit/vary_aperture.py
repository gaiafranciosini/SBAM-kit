#!/usr/bin/env python3
"""
vary_aperture.py

Da un file aperture.out del tipo
    APERTURE: 1.7 5.7 -5.7 -1.7 -5.9 -1.9 1.9 5.9
    FIELD_SIZE: 3.4 3.8
crea due varianti, SENZA toccare l'originale:
    <nome>_plus<N>mm.<ext>   apertura allargata di N mm lungo x e y
    <nome>_minus<N>mm.<ext>  apertura ristretta di N mm lungo x e y

N (--delta-mm, default 2) e' la variazione dell'apertura TOTALE lungo ciascun asse:
ogni lamella si sposta quindi di N/2 mm (1 mm per N=2). Per spostare ogni lamella
di 2 mm (apertura totale +/- 4 mm) usare --delta-mm 4.

Convenzione (valori in cm, come nel file):
    APERTURE = Xr_in Xr_out  Xl_out Xl_in  Yb_out Yb_in  Yt_in Yt_out
Entrambi i bordi di ogni lamella ricevono lo stesso spostamento, cosi' lo spessore
resta invariato. FIELD_SIZE varia di N/10 cm per asse. Tutto il resto del file
(altre righe, commenti, a-capo) resta identico.

Uso:
    python3 vary_aperture.py aperture.out
    python3 vary_aperture.py aperture.out --delta-mm 2 --outdir varianti/
"""

import argparse
import re
import sys
from decimal import Decimal, InvalidOperation
from pathlib import Path

LABEL = re.compile(r"^(\s*)(APERTURE|FIELD_SIZE)(\s*:\s*)(.*?)(\r?\n?)$", re.I)
NUM = re.compile(r"^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$")


def read_numbers(text, n, key):
    parts = text.split(maxsplit=n)
    if len(parts) < n or not all(NUM.fullmatch(p) for p in parts[:n]):
        raise ValueError(f"Riga {key}: attesi {n} numeri, trovato {text!r}")
    extra = parts[n] if len(parts) > n else ""
    if extra and not extra.startswith(("#", "!", ";")):
        raise ValueError(f"Riga {key}: campi extra non riconosciuti: {extra!r}")
    return [Decimal(p) for p in parts[:n]], extra


def fmt(x):
    s = format(x.normalize(), "f")
    return "0" if s in ("-0", "+0") else s


def check_aperture(v):
    """Controlla il formato atteso (segni e ordine); apertura interna sempre > 0."""
    ok = (v[0] < v[1] and v[2] < v[3] and v[4] < v[5] and v[6] < v[7]
          and v[0] >= 0 and v[1] > 0 and v[2] < 0 and v[3] <= 0
          and v[4] < 0 and v[5] <= 0 and v[6] >= 0 and v[7] > 0)
    if not ok:
        raise ValueError("APERTURE non ha il formato atteso "
                         "(Xr_in Xr_out Xl_out Xl_in Yb_out Yb_in Yt_in Yt_out): "
                         "non modifico il file.")
    if not (-v[3] > 0 and v[0] > 0 and -v[5] > 0 and v[6] > 0):
        raise ValueError("Una lamella e' al centro o oltre: apertura non positiva.")


def transform(text, d):
    """d = spostamento per lamella in cm (positivo = apertura piu' larga)."""
    out, aperture, field = [], None, None
    for line in text.splitlines(keepends=True):
        m = LABEL.match(line)
        if not m:
            out.append(line)
            continue
        lead, label, sep, body, end = m.groups()
        key = label.upper()
        if key == "APERTURE":
            if aperture is not None:
                raise ValueError("APERTURE compare piu' di una volta.")
            v, extra = read_numbers(body, 8, key)
            check_aperture(v)
            new = [a + s for a, s in zip(v, (d, d, -d, -d, -d, -d, d, d))]
            check_aperture(new)
            aperture = (v, new)
        else:
            if field is not None:
                raise ValueError("FIELD_SIZE compare piu' di una volta.")
            v, extra = read_numbers(body, 2, key)
            new = [x + 2 * d for x in v]
            if not all(x > 0 for x in new):
                raise ValueError("FIELD_SIZE non positivo nella variante ristretta.")
            field = (v, new)
        out.append(lead + label + sep + " ".join(fmt(x) for x in new)
                   + (" " + extra if extra else "") + end)
    if aperture is None:
        raise ValueError("Riga APERTURE assente: non modifico il file.")
    return "".join(out), aperture, field


def dec(s):
    try:
        v = Decimal(s)
    except InvalidOperation:
        raise argparse.ArgumentTypeError(f"numero non valido: {s}")
    if v <= 0:
        raise argparse.ArgumentTypeError("deve essere > 0")
    return v


def main():
    ap = argparse.ArgumentParser(description="Crea varianti +/- di aperture.out")
    ap.add_argument("files", nargs="+", type=Path, help="file aperture.out")
    ap.add_argument("--delta-mm", type=dec, default=Decimal("2"),
                    help="variazione dell'apertura TOTALE per asse, in mm (default 2: "
                         "1 mm per lamella)")
    ap.add_argument("--outdir", type=Path, help="cartella di output (default: quella del file)")
    ap.add_argument("--force", action="store_true", help="sovrascrive varianti gia' esistenti")
    args = ap.parse_args()

    d = args.delta_mm / 20          # mm totali per asse -> cm per lamella
    tag = fmt(args.delta_mm)

    for src in args.files:
        text = src.read_text(encoding="utf-8")
        outdir = args.outdir or src.parent
        outdir.mkdir(parents=True, exist_ok=True)
        for name, sign in (("plus", 1), ("minus", -1)):
            new_text, (a0, a1), field = transform(text, sign * d)
            dest = outdir / f"{src.stem}_{name}{tag}mm{src.suffix}"
            if dest.resolve() == src.resolve():
                raise RuntimeError("Il file di output coincide con l'originale")
            if dest.exists() and not args.force:
                raise FileExistsError(f"{dest} esiste gia' (usa --force)")
            dest.write_text(new_text, encoding="utf-8")
            print(f"{dest}")
            print("  APERTURE   : " + " ".join(fmt(x) for x in a0))
            print("            -> " + " ".join(fmt(x) for x in a1))
            if field:
                print("  FIELD_SIZE : " + " ".join(fmt(x) for x in field[0])
                      + " -> " + " ".join(fmt(x) for x in field[1]))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("ERRORE:", exc, file=sys.stderr)
        sys.exit(1)
