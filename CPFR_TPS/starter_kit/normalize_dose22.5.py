#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
normalize_dose_to_prescription.py

Normalizza mappe DOSE_*_BEV.mhd rispetto alla DOSE PRESCRITTA
(letta da setup.out), NON rispetto a Dmax.

Prima di fare i calcoli crea una COPIA INDIPENDENTE della mappa di dose
in una nuova cartella di output. Le operazioni usano ESCLUSIVAMENTE la
copia e non modificano mai le mappe originali (ne' .mhd ne' dati binari).

    I = D_voxel[Gy] / D_prescritta[Gy]

Una mappa normalizzata alla prescrizione puo' avere voxel > 1.0 (100%).
La dose in input deve essere espressa in Gy: se proviene da simulazioni
in Gy/primario serve un fattore di calibrazione/erogazione fisicamente
corretto (--gy-per-map-unit), NON necessariamente il numero di primari MC.

Esempi:
    python3 normalize_dose_to_prescription.py --sim-dir AUTO_sim7MeV_1
    python3 normalize_dose_to_prescription.py --sim-dir AUTO_sim7MeV_1 --preview
    python3 normalize_dose_to_prescription.py --root . --all-auto
    python3 normalize_dose_to_prescription.py --sim-dir AUTO_sim7MeV_1 \\
        --input DOSE_7MeV_1093pulses_CONV_BEV.mhd --tolerance-gy 0.05
    python3 normalize_dose_to_prescription.py --sim-dir AUTO_sim7MeV_1 \\
        --rx-key 'PrescribedDose' --gy-per-map-unit VALORE_CALIBRATO

Se setup.out non riporta il dato sotto una delle etichette comuni, usare
--rx-key 'NOME_DEL_CAMPO' oppure --dose-prescription-gy 22.5 (valore
controllato manualmente). Il programma NON assume 22.5 in automatico.
"""

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import SimpleITK as sitk


# Float: ammette 22.5, 22,5, 2.25e1, 2.25E+01, 2.25D+01.
NUM_RE = re.compile(r"(?<![\w.])[-+]?(?:\d+(?:[.,]\d*)?|[.,]\d+)(?:[eEdD][-+]?\d+)?(?=$|[^\w.]|(?i:gy)\b)")
# Etichette intenzionalmente restrittive: NON prendere un numero qualsiasi
# da setup.out (potrebbe essere il numero di impulsi o la dose per impulso).
RX_RE = re.compile(
    r"(?i)(?:"
    r"prescribed[\s_-]*dose|dose[\s_-]*prescribed|"
    r"prescription[\s_-]*dose|dose[\s_-]*prescription|"
    r"dose[\s_-]*prescritta|dose[\s_-]*di[\s_-]*prescrizione|"
    r"dose[\s_-]*prescrizione|dosepresc(?:ritta)?|"
    r"rx[\s_-]*dose|dose[\s_-]*rx|target[\s_-]*dose|"
    r"(?<![A-Za-z0-9_])prescription\b|"
    r"d[\s_-]*presc(?:ritta)?"
    r")"
)
LEVELS = [100, 95, 90, 80, 50, 10]


def finite_positive(name, value):
    if value is None or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} deve essere un valore finito > 0.")
    return float(value)


def parse_number(s):
    return float(s.replace(",", ".").replace("D", "E").replace("d", "e"))


def read_prescription(setup_path, rx_key=None):
    """Estrae una sola prescrizione in Gy. Su dubbio o conflitto: errore."""
    content = setup_path.read_text(encoding="utf-8", errors="replace")
    candidates = []
    for i, original_line in enumerate(content.splitlines(), start=1):
        # Commenti HTML/FRED usuali; una riga commentata non e' un parametro.
        line = original_line.strip()
        if not line or line.startswith(("#", "//", ";", "%")):
            continue
        line = re.split(r"\s+(?://|#|;)", line, maxsplit=1)[0].strip()
        if rx_key:
            pos = line.casefold().find(rx_key.casefold())
            if pos < 0:
                continue
            trailing = line[pos + len(rx_key):]
        else:
            match = RX_RE.search(line)
            if not match:
                continue
            trailing = line[match.end():]
        # Rimane solo il testo DOPO l'etichetta per evitare, ad es., numeri in 'PTV95'.
        nums = NUM_RE.findall(trailing)
        if len(nums) != 1:
            raise ValueError(
                f"Campo dose ambiguo in {setup_path}, riga {i}: {original_line!r}. "
                "Specificare --rx-key oppure --dose-prescription-gy."
            )
        if "cgy" in trailing.lower():
            raise ValueError(
                f"La dose alla riga {i} sembra espressa in cGy, non in Gy. "
                "Verificare e usare --dose-prescription-gy con il valore in Gy."
            )
        val = parse_number(nums[0])
        finite_positive(f"Dose prescritta alla riga {i}", val)
        candidates.append((i, val, original_line.strip()))

    if not candidates:
        raise ValueError(
            f"Non riconosco una dose prescritta in {setup_path}. "
            "Controllare il campo e usare --rx-key 'nome_campo' oppure "
            "--dose-prescription-gy 22.5 (solo se verificato)."
        )
    distinct = [c for c in candidates if not np.isclose(c[1], candidates[0][1], rtol=0, atol=1e-9)]
    if distinct:
        raise ValueError(
            f"Dosi di prescrizione diverse in {setup_path}: "
            f"{[(line, dose) for line, dose, _ in candidates]}. "
            "Specificare --rx-key (se univoco) o --dose-prescription-gy."
        )
    return candidates[0][1], candidates[0][2], candidates[0][0]


def resolve_cli_file(value, sim_dir):
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    if path.exists():
        return path.resolve()
    return (sim_dir / path).resolve()


def find_setup(sim_dir, supplied):
    if supplied:
        path = resolve_cli_file(supplied, sim_dir)
        if not path.is_file():
            raise FileNotFoundError(f"setup.out non trovato: {path}")
        return path
    candidates = [sim_dir / "setup.out", sim_dir.parent / "setup.out"]
    for c in candidates:
        if c.is_file():
            return c.resolve()
    raise FileNotFoundError(
        "setup.out non trovato nella cartella AUTO_sim o nella cartella superiore. "
        "Indicare il percorso con --setup /percorso/setup.out."
    )


def find_doses(sim_dir, specified):
    if specified:
        paths = [resolve_cli_file(specified, sim_dir)]
    else:
        paths = sorted(sim_dir.glob("DOSE_*_BEV.mhd"))
    if not paths:
        raise FileNotFoundError(f"Nessuna DOSE_*_BEV.mhd trovata in {sim_dir}")
    for p in paths:
        if not p.is_file() or p.suffix.lower() != ".mhd":
            raise FileNotFoundError(f"Mappa di dose .mhd non valida: {p}")
    return paths


def fresh_outdir(sim_dir):
    """Cartella univoca; non sovrascrive mai un precedente risultato."""
    base = f"isodose_prescrizione_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    for suffix in range(10000):
        candidate = sim_dir / (base if suffix == 0 else f"{base}_{suffix:03d}")
        try:
            candidate.mkdir(parents=False, exist_ok=False)
            return candidate
        except FileExistsError:
            continue
    raise RuntimeError("Impossibile creare una cartella di output univoca")


def same_geometry(a, b):
    return (
        a.GetSize() == b.GetSize()
        and np.allclose(a.GetSpacing(), b.GetSpacing(), atol=1e-6)
        and np.allclose(a.GetOrigin(), b.GetOrigin(), atol=1e-6)
        and np.allclose(a.GetDirection(), b.GetDirection(), atol=1e-6)
    )


def find_ptv(sim_dir, requested):
    if requested:
        path = resolve_cli_file(requested, sim_dir)
        if not path.is_file():
            raise FileNotFoundError(f"PTV non trovato: {path}")
        return path, True
    auto = sim_dir.parent / "imgs" / "PTV_plan.mhd"
    return (auto.resolve(), False) if auto.is_file() else (None, False)


def preview_image(iso, outfile, percent):
    """Anteprima XY con scale riferite alla prescrizione, NON a Dmax."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    arr = sitk.GetArrayViewFromImage(iso)
    k = arr.shape[0] // 2
    xy = arr[k, :, :]
    displayed = np.ma.masked_invalid(xy)
    scale = 100.0 if percent else 1.0
    fig, ax = plt.subplots(figsize=(7, 6))
    finite = xy[np.isfinite(xy)]
    max_view = float(finite.max()) if finite.size else 1.1 * scale
    vmax = max(1.1 * scale, max_view)
    im = ax.imshow(displayed, origin="lower", cmap="turbo", vmin=0, vmax=vmax)
    thresholds = np.array([0.1, 0.5, 0.8, 0.9, 0.95, 1.0, 1.07]) * scale
    present = [float(t) for t in thresholds if finite.size and np.min(finite) < t < np.max(finite)]
    if present:
        contour = ax.contour(xy, levels=present, linewidths=1.0, colors="black")
        ax.clabel(contour, inline=True, fontsize=7, fmt=lambda v: f"{v:.1f}%" if percent else f"{v:.2f}")
    fig.colorbar(im, ax=ax, label="% dose prescritta" if percent else "Dose / prescrizione")
    ax.set_title("Isodosi normalizzate alla prescrizione (slice centrale XY)")
    fig.tight_layout()
    fig.savefig(outfile, dpi=160)
    plt.close(fig)


def process_one_sim(sim_dir, args):
    setup_path = find_setup(sim_dir, args.setup)
    doses = find_doses(sim_dir, args.input)
    if args.dose_prescription_gy is not None:
        rx = finite_positive("--dose-prescription-gy", args.dose_prescription_gy)
        rx_source = f"override esplicito --dose-prescription-gy={rx}"
    else:
        rx, original_line, lineno = read_prescription(setup_path, args.rx_key)
        rx_source = f"{setup_path} (riga {lineno}: {original_line})"
    factor = finite_positive("--gy-per-map-unit", args.gy_per_map_unit)
    tolerance = args.tolerance_gy
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("--tolerance-gy deve essere >= 0")
    if tolerance >= rx * 0.10:
        raise ValueError("Tolleranza troppo ampia rispetto alla prescrizione: verificare unita' Gy")

    ptv_path, ptv_explicit = find_ptv(sim_dir, args.ptv)
    ptv_img = sitk.ReadImage(str(ptv_path)) if ptv_path else None

    # Controlli PRIMA di creare una directory o copiare file; evita output
    # fuorvianti se dose raw e' in Gy/primario e manca calibrazione.
    sources = []
    for dose_path in doses:
        src_img = sitk.ReadImage(str(dose_path))
        src_arr = sitk.GetArrayViewFromImage(src_img)
        valid = np.isfinite(src_arr)
        if not valid.any():
            raise ValueError(f"Mappa senza dosi finite: {dose_path}")
        mn = float(src_arr[valid].min())
        mx = float(src_arr[valid].max())
        if mx <= 0:
            raise ValueError(f"Dose massima non positiva in {dose_path}: {mx}")
        if mx * factor < 1e-6:
            raise ValueError(
                f"STOP: la dose massima di {dose_path.name} e' {mx:.6g} in unita' della mappa "
                f"({mx * factor:.6g} Gy con fattore {factor:g}). "
                "Potrebbe essere Gy/primario. Impossibile normalizzare in modo "
                "fisicamente corretto senza un fattore di conversione verificato. "
                "Usare --gy-per-map-unit solo se la calibrazione e' nota."
            )
        if np.any(src_arr[valid] < 0):
            print(f"ATTENZIONE: {dose_path.name} contiene dosi negative: controllare l'input")
        if ptv_img is not None and not same_geometry(src_img, ptv_img):
            if ptv_explicit:
                raise ValueError(f"Geometria non compatibile con PTV: {dose_path}")
            print(f"ATTENZIONE: geometria PTV non corrispondente a {dose_path.name}; salto analisi PTV")
        sources.append((dose_path, mn, mx))

    out_dir = fresh_outdir(sim_dir)
    copied = out_dir / "dose_BEV_copie"
    generated = out_dir / "isodose_su_prescrizione"
    copied.mkdir()
    generated.mkdir()
    # Copia del setup a fini di tracciabilita' (non modifica l'originale).
    import shutil
    shutil.copy2(setup_path, out_dir / "setup_originale_copia.out")

    report = {
        "sim_dir": str(sim_dir), "setup_source": str(setup_path),
        "dose_prescribed_gy": rx, "rx_source": rx_source,
        "gy_per_map_unit": factor, "tolerance_gy": tolerance,
        "normalization": "isodose = dose_bev_Gy / dose_prescribed_Gy",
        "output_scale": "0-100 (percentual)" if args.percent else "frazione: 1.0 = 100% prescrizione",
        "maps": [],
    }
    print(f"\n{'=' * 76}\nSIMULAZIONE: {sim_dir}\n{'=' * 76}")
    print(f"setup.out                : {setup_path}")
    print(f"Dose prescritta          : {rx:.6f} Gy  ({rx_source})")
    print(f"Conversione in Gy        : dose_mappa * {factor:g}")
    print(f"Tolleranza               : -{tolerance:.4f} Gy per soglia")
    print(f"Directory lavoro         : {out_dir}")

    for dose_path, original_min, original_max in sources:
        # Copia indipendente (incluso payload binario .raw) usando SimpleITK.
        # Anche con input .mhd che riferisce raw in cartelle esterne, l'output
        # viene riscritto in locale ed e' completamente autonomo.
        img_orig = sitk.ReadImage(str(dose_path))
        dose_copy = copied / dose_path.name
        sitk.WriteImage(img_orig, str(dose_copy), useCompression=False)
        del img_orig

        # LETTURA DALLA COPIA, mai dalla mappa di partenza.
        working_img = sitk.ReadImage(str(dose_copy))
        working_array = sitk.GetArrayViewFromImage(working_img)
        norm_array = np.asarray(working_array, dtype=np.float64) * (factor / rx)
        if not np.isfinite(norm_array).all():
            raise ValueError(f"Mappa normalizzata con valori non finiti: {dose_copy}")
        if args.percent:
            norm_array *= 100.0
        # Float32 per compatibilita' con gli script di analisi .mhd.
        iso_img = sitk.GetImageFromArray(norm_array.astype(np.float32))
        iso_img.CopyInformation(working_img)
        iso_name = f"{dose_path.stem}_isodosemap_prescrizione.mhd"
        iso_path = generated / iso_name
        sitk.WriteImage(iso_img, str(iso_path), useCompression=False)
        dmax_gy = original_max * factor
        print(f"\nInput originale   : {dose_path}")
        print(f"Copia indipendente: {dose_copy}")
        print(f"Dose massima (Gy) : {dmax_gy:.6g} (solo diagnostica; NON usata per normalizzare)")
        print(f"Isodose generata  : {iso_path}")
        print(f"Valore massimo isodose: {dmax_gy/rx*(100 if args.percent else 1):.6f}")
        print("Soglie di confronto della mappa (con -tolleranza):")
        level_summaries = []
        for level in LEVELS:
            cutoff_gy = rx * level / 100.0 - tolerance
            threshold = cutoff_gy / rx * (100.0 if args.percent else 1.0)
            dct = {"level_percent": level, "cutoff_gy_with_tolerance": cutoff_gy,
                   "normalized_threshold": threshold}
            print(f"  {level:3d}% => dose >= {cutoff_gy:.5f} Gy -> valore >= {threshold:.6f}", end="")
            if ptv_img is not None and same_geometry(working_img, ptv_img):
                mask = sitk.GetArrayViewFromImage(ptv_img) > 0
                total = int(np.count_nonzero(mask))
                if total > 0:
                    # Confronto in Gy sul valore ORIGINALE della copia (float64),
                    # per evitare errori di arrotondamento float32 vicino alla soglia.
                    covered = int(np.count_nonzero(mask & np.isfinite(working_array)
                                                      & (working_array.astype(np.float64) * factor >= cutoff_gy)))
                    cov = 100.0 * covered / total
                    dct.update({"ptv_total_voxels": total, "ptv_covered_voxels": covered,
                                "ptv_coverage_percent": cov})
                    print(f" | PTV: {covered}/{total} = {cov:.2f}%", end="")
            print()
            level_summaries.append(dct)
        if args.preview:
            preview_path = generated / f"{dose_path.stem}_anteprima_prescrizione.png"
            preview_image(iso_img, preview_path, args.percent)
            print(f"Anteprima          : {preview_path}")
        report["maps"].append({
            "original_dose_map": str(dose_path), "copied_dose_map": str(dose_copy),
            "normalized_isodose_map": str(iso_path), "dose_map_min": original_min,
            "dose_map_max": original_max, "effective_dmax_gy": dmax_gy,
            "ptv": str(ptv_path) if ptv_path is not None else None,
            "levels": level_summaries,
        })
    report_path = out_dir / "report_normalizzazione.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nReport JSON: {report_path}")
    print("Mappe originali: NON MODIFICATE.\n")
    return out_dir


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Copia mappe BEV in cartelle isolate e le normalizza alla prescrizione in setup.out."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--sim-dir", default=None, help="Cartella di una singola AUTO_sim (o MANUAL_sim)")
    mode.add_argument("--all-auto", action="store_true", help="Elabora tutte le cartelle AUTO_sim* sotto --root")
    parser.add_argument("--root", default=".", help="Radice in cui trovare AUTO_sim* con --all-auto")
    parser.add_argument("--input", help="File DOSE_*_BEV.mhd specifico (altrimenti tutti i BEV nella sim)")
    parser.add_argument("--setup", help="Percorso di setup.out (default: AUTO_sim/setup.out oppure ../setup.out)")
    parser.add_argument("--rx-key", help="Nome preciso del campo dose prescritta in setup.out se non riconosciuto")
    parser.add_argument("--dose-prescription-gy", type=float,
                        help="Override ESPLICITO dose prescritta (Gy); non esegue parsing dose da setup.out")
    parser.add_argument("--gy-per-map-unit", type=float, default=1.0,
                        help="Fattore CALIBRATO per convertire i valori della mappa in Gy; 1 se gia' in Gy")
    parser.add_argument("--tolerance-gy", type=float, default=0.05,
                        help="Tolleranza in Gy sottratta alle soglie di verifica PTV; non modifica la normalizzazione")
    parser.add_argument("--percent", action="store_true", help="Mappa in percentuale (100=100%% Rx) invece che 1=100%% Rx")
    parser.add_argument("--ptv", help="Mappa PTV .mhd per verifica copertura; cercata automaticamente in ../imgs/PTV_plan.mhd")
    parser.add_argument("--preview", action="store_true", help="Salva un'immagine PNG della slice centrale")
    args = parser.parse_args()
    if args.all_auto:
        if args.input is not None:
            parser.error("--input specifico e --all-auto non sono compatibili")
        root = Path(args.root).expanduser().resolve()
        sims = sorted(p for p in root.glob("AUTO_sim*") if p.is_dir()
                      and any(p.glob("DOSE_*_BEV.mhd")))
        if not sims:
            parser.error(f"Non trovo cartelle AUTO_sim* con mappe DOSE_*_BEV.mhd sotto {root}")
    else:
        sims = [Path(args.sim_dir or '.').expanduser().resolve()]
    errors = []
    for sim in sims:
        try:
            process_one_sim(sim, args)
        except Exception as exc:
            print(f"\nERRORE ({sim}): {exc}", file=sys.stderr)
            errors.append((sim, str(exc)))
    if errors:
        print(f"\n{len(errors)} simulazione/i non elaborate correttamente.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
