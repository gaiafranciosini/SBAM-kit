#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
normalize_dose_to_prescription.py

Normalizza le mappe DOSE_*_BEV.mhd rispetto alla DOSE PRESCRITTA (setup.out),
NON rispetto a Dmax:  I = D_voxel[Gy] / D_prescritta[Gy]   (1.0 = 100% Rx)

Le mappe originali non vengono mai modificate: ne viene scritta una copia
in una cartella di output nuova e i risultati vanno li'.
La dose in input deve essere in Gy; se e' in Gy/primario serve un fattore
calibrato (--gy-per-map-unit), NON il numero di primari MC.
Il valore 22.5 Gy NON e' assunto: se setup.out non e' riconosciuto usare
--rx-key 'NOME_CAMPO' oppure --dose-prescription-gy (verificato a mano).

Esempi:
    python3 normalize_dose_to_prescription.py --sim-dir AUTO_sim7MeV_1 --preview
    python3 normalize_dose_to_prescription.py --root . --all-auto
"""

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import SimpleITK as sitk

LEVELS = [100, 95, 90, 80, 50, 10]
NUM_RE = re.compile(r"(?<![\w.])[-+]?(?:\d+(?:[.,]\d*)?|[.,]\d+)(?:[eEdD][-+]?\d+)?(?=$|[^\w.]|(?i:gy)\b)")
# Etichette restrittive: non prendere un numero qualsiasi (impulsi, dose/impulso...).
RX_RE = re.compile(
    r"(?i)prescribed[\s_-]*dose|dose[\s_-]*prescribed|prescription[\s_-]*dose|"
    r"dose[\s_-]*prescription|dose[\s_-]*prescritta|dose[\s_-]*di[\s_-]*prescrizione|"
    r"dosepresc(?:ritta)?|rx[\s_-]*dose|dose[\s_-]*rx|target[\s_-]*dose|"
    r"(?<![A-Za-z0-9_])prescription\b|d[\s_-]*presc(?:ritta)?"
)


def positive(name, v):
    if v is None or not np.isfinite(v) or v <= 0:
        raise ValueError(f"{name} deve essere un valore finito > 0.")
    return float(v)


def resolve(value, base):
    p = Path(value).expanduser()
    return (p if p.is_absolute() or p.exists() else base / p).resolve()


def read_prescription(setup, key=None):
    """Ritorna (riga, dose_Gy, testo). Su dubbio o conflitto: errore."""
    found = []
    for i, raw in enumerate(setup.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        line = re.split(r"\s+(?://|#|;)", raw.strip(), maxsplit=1)[0].strip()
        if not line or line[0] in "#;%" or line.startswith("//"):
            continue
        if key:
            pos = line.casefold().find(key.casefold())
            tail = line[pos + len(key):] if pos >= 0 else None
        else:
            m = RX_RE.search(line)
            tail = line[m.end():] if m else None
        if tail is None:
            continue
        nums = NUM_RE.findall(tail)  # solo testo DOPO l'etichetta
        if len(nums) != 1:
            raise ValueError(f"Campo dose ambiguo in {setup}, riga {i}: {raw!r}. "
                             "Usare --rx-key o --dose-prescription-gy.")
        if "cgy" in tail.lower():
            raise ValueError(f"Riga {i}: dose in cGy, non Gy. Usare --dose-prescription-gy.")
        val = float(nums[0].replace(",", ".").replace("D", "E").replace("d", "e"))
        found.append((i, positive(f"Dose riga {i}", val), raw.strip()))
    if not found:
        raise ValueError(f"Dose prescritta non riconosciuta in {setup}. "
                         "Usare --rx-key o --dose-prescription-gy.")
    if any(not np.isclose(f[1], found[0][1], rtol=0, atol=1e-9) for f in found):
        raise ValueError(f"Dosi di prescrizione diverse in {setup}: {[(f[0], f[1]) for f in found]}.")
    return found[0]


def same_geometry(a, b):
    return (a.GetSize() == b.GetSize()
            and np.allclose(a.GetSpacing(), b.GetSpacing(), atol=1e-6)
            and np.allclose(a.GetOrigin(), b.GetOrigin(), atol=1e-6)
            and np.allclose(a.GetDirection(), b.GetDirection(), atol=1e-6))


def preview(iso, outfile, percent):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    arr = sitk.GetArrayViewFromImage(iso)
    xy = arr[arr.shape[0] // 2]
    s = 100.0 if percent else 1.0
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(xy, origin="lower", cmap="turbo", vmin=0, vmax=max(1.1 * s, float(np.nanmax(xy))))
    lv = [t * s for t in (0.1, 0.5, 0.8, 0.9, 0.95, 1.0, 1.07) if xy.min() < t * s < xy.max()]
    if lv:
        ax.clabel(ax.contour(xy, levels=lv, colors="black", linewidths=1.0), fontsize=7)
    fig.colorbar(im, ax=ax, label="% dose prescritta" if percent else "Dose / prescrizione")
    ax.set_title("Isodosi normalizzate alla prescrizione (slice centrale XY)")
    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    plt.close(fig)


def process(sim, a):
    # --- setup.out e prescrizione
    setup = resolve(a.setup, sim) if a.setup else next(
        (p for p in (sim / "setup.out", sim.parent / "setup.out") if p.is_file()), None)
    if setup is None or not setup.is_file():
        raise FileNotFoundError("setup.out non trovato: usare --setup /percorso/setup.out")
    if a.dose_prescription_gy is not None:
        rx, rx_src = positive("--dose-prescription-gy", a.dose_prescription_gy), "override esplicito"
    else:
        ln, rx, txt = read_prescription(setup, a.rx_key)
        rx_src = f"{setup} (riga {ln}: {txt})"
    factor = positive("--gy-per-map-unit", a.gy_per_map_unit)
    tol = a.tolerance_gy
    if not np.isfinite(tol) or tol < 0 or tol >= 0.1 * rx:
        raise ValueError("--tolerance-gy deve essere >= 0 e molto minore della prescrizione (Gy?)")

    # --- PTV (opzionale)
    ptv_path = resolve(a.ptv, sim) if a.ptv else sim.parent / "imgs" / "PTV_plan.mhd"
    if a.ptv and not ptv_path.is_file():
        raise FileNotFoundError(f"PTV non trovato: {ptv_path}")
    ptv = sitk.ReadImage(str(ptv_path)) if ptv_path.is_file() else None

    # --- letture e controlli PRIMA di creare output
    doses = [resolve(a.input, sim)] if a.input else sorted(sim.glob("DOSE_*_BEV.mhd"))
    if not doses:
        raise FileNotFoundError(f"Nessuna DOSE_*_BEV.mhd in {sim}")
    maps = []
    for p in doses:
        img = sitk.ReadImage(str(p))
        gy = sitk.GetArrayFromImage(img).astype(np.float64) * factor
        if not np.isfinite(gy).all() or gy.max() <= 0:
            raise ValueError(f"Mappa non valida (non finita o dose max <= 0): {p}")
        if gy.max() < 1e-6:
            raise ValueError(f"STOP: dose max di {p.name} = {gy.max():.3g} Gy. Forse Gy/primario: "
                             "serve --gy-per-map-unit calibrato.")
        if (gy < 0).any():
            print(f"ATTENZIONE: {p.name} contiene dosi negative")
        use_ptv = ptv is not None and same_geometry(img, ptv)
        if ptv is not None and not use_ptv:
            if a.ptv:
                raise ValueError(f"Geometria non compatibile con PTV: {p}")
            print(f"ATTENZIONE: geometria PTV diversa da {p.name}; salto copertura PTV")
        maps.append((p, img, gy, use_ptv))

    # --- output
    out = sim / f"isodose_prescrizione_{datetime.now():%Y%m%d_%H%M%S}"
    (out / "dose_BEV_copie").mkdir(parents=True, exist_ok=False)
    (out / "isodose_su_prescrizione").mkdir()
    scale = 100.0 if a.percent else 1.0
    print(f"\n{'=' * 70}\nSIM: {sim}\nDose prescritta: {rx:.6f} Gy [{rx_src}]\n"
          f"Fattore Gy: {factor:g} | tolleranza: -{tol:g} Gy | output: {out}\n{'=' * 70}")
    report = {"sim_dir": str(sim), "setup": str(setup), "dose_prescribed_gy": rx, "rx_source": rx_src,
              "gy_per_map_unit": factor, "tolerance_gy": tol,
              "output_scale": "percentuale" if a.percent else "frazione (1.0 = 100% Rx)", "maps": []}

    for p, img, gy, use_ptv in maps:
        sitk.WriteImage(img, str(out / "dose_BEV_copie" / p.name), useCompression=False)  # copia
        iso = sitk.GetImageFromArray((gy / rx * scale).astype(np.float32))
        iso.CopyInformation(img)
        iso_path = out / "isodose_su_prescrizione" / f"{p.stem}_isodosemap_prescrizione.mhd"
        sitk.WriteImage(iso, str(iso_path), useCompression=False)
        mask = sitk.GetArrayViewFromImage(ptv) > 0 if use_ptv else None
        print(f"\n{p.name}: Dmax {gy.max():.6g} Gy (solo diagnostica) -> {iso_path.name}")
        levels = []
        for lv in LEVELS:
            cut = rx * lv / 100.0 - tol
            d = {"level_percent": lv, "cutoff_gy": cut, "threshold": cut / rx * scale}
            line = f"  {lv:3d}% => >= {cut:.5f} Gy"
            if mask is not None and mask.any():
                cov = int(np.count_nonzero(mask & (gy >= cut)))
                d.update(ptv_voxels=int(mask.sum()), ptv_covered=cov,
                         ptv_coverage_percent=100.0 * cov / mask.sum())
                line += f" | PTV {cov}/{int(mask.sum())} = {d['ptv_coverage_percent']:.2f}%"
            print(line)
            levels.append(d)
        if a.preview:
            preview(iso, iso_path.with_name(f"{p.stem}_anteprima_prescrizione.png"), a.percent)
        report["maps"].append({"original": str(p), "isodose": str(iso_path),
                               "dmax_gy": float(gy.max()), "levels": levels})

    (out / "report_normalizzazione.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nReport: {out / 'report_normalizzazione.json'}\nMappe originali: NON MODIFICATE.")


def main():
    ap = argparse.ArgumentParser(description="Normalizza mappe BEV alla dose prescritta in setup.out.")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--sim-dir", help="Cartella di una singola AUTO_sim/MANUAL_sim (default: .)")
    g.add_argument("--all-auto", action="store_true", help="Tutte le AUTO_sim* sotto --root")
    ap.add_argument("--root", default=".", help="Radice per --all-auto")
    ap.add_argument("--input", help="Singolo DOSE_*_BEV.mhd (default: tutti)")
    ap.add_argument("--setup", help="Percorso di setup.out (default: sim/setup.out o ../setup.out)")
    ap.add_argument("--rx-key", help="Nome del campo dose prescritta in setup.out")
    ap.add_argument("--dose-prescription-gy", type=float, help="Override esplicito della prescrizione (Gy)")
    ap.add_argument("--gy-per-map-unit", type=float, default=1.0, help="Fattore calibrato mappa -> Gy")
    ap.add_argument("--tolerance-gy", type=float, default=0.05, help="Tolleranza (Gy) sulle soglie PTV")
    ap.add_argument("--percent", action="store_true", help="Output in %% (100 = 100%% Rx)")
    ap.add_argument("--ptv", help="PTV .mhd (default: ../imgs/PTV_plan.mhd se esiste)")
    ap.add_argument("--preview", action="store_true", help="Salva PNG della slice centrale")
    a = ap.parse_args()

    if a.all_auto:
        if a.input:
            ap.error("--input e --all-auto non sono compatibili")
        root = Path(a.root).expanduser().resolve()
        sims = sorted(p for p in root.glob("AUTO_sim*") if p.is_dir() and any(p.glob("DOSE_*_BEV.mhd")))
        if not sims:
            ap.error(f"Nessuna AUTO_sim* con DOSE_*_BEV.mhd sotto {root}")
    else:
        sims = [Path(a.sim_dir or ".").expanduser().resolve()]

    failed = 0
    for sim in sims:
        try:
            process(sim, a)
        except Exception as e:
            failed += 1
            print(f"\nERRORE ({sim}): {e}", file=sys.stderr)
    if failed:
        print(f"\n{failed} simulazione/i non elaborate.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
