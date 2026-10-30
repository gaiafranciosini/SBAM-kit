#!/usr/bin/env python3
"""
Analizza file DVH (7 MeV, 9 MeV o altre energie).
Supporta la specifica dinamica delle ROI/strutture da linea di comando.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STATUS_EN = {
    "RAGGIUNTA": "MET",
    "NON RAGGIUNTA": "NOT MET",
    "N/A (organo a rischio)": "N/A (OAR)",
}


DEFAULT_ROIS: list[str] = []


def locate_setup(explicit: Path | None, *data_dirs: Path) -> Path:
    if explicit is not None:
        return explicit.expanduser().resolve()
 
    checked = set()
    for origin in (*data_dirs, Path.cwd()):
        base = origin.expanduser().resolve()
        for directory in (base, *base.parents):
            if directory == Path(directory.anchor):
                break
            for candidate in (directory / "prescription_inputs" / "setup.out",
                              directory / "setup.out"):
                if candidate in checked:
                    continue
                checked.add(candidate)
                if candidate.is_file():
                    return candidate
    raise FileNotFoundError("setup.out not found; refusing to use a default PTV target")


@dataclass
class Result:
    file: Path
    energy: str
    modality: str
    structure: str
    plan: str
    max_tab_gy: float
    dmax_gy: float
    mean_weighted_gy: float
    v_prescription: float
    d95_gy: float
    status: str


def read_setup_file(path: Path) -> tuple[float | None, float | None, float | None]:
    """Read prescription Gy, required PTV volume % and prescribed isodose %."""
    
    if not path.is_file():
        raise FileNotFoundError(
            f'Prescription file not found: {path}; refusing to use an implicit PTV target'
        )

    values = {}
    valid_keys = {
        'PRESCRIPTION DOSE': 'rx_gy',
        'PERCENTAGE PRESCRIPTION VOLUME': 'volume_pct',
        'PERCENTAGE PRESCRIPTION DOSE': 'isodose_pct',
    }
    for line_number, line in enumerate(path.read_text(encoding='utf-8', errors='replace').splitlines(), 1):
        if ':' not in line:
            continue
        key, raw = (part.strip() for part in line.split(':', 1))
        field = valid_keys.get(key.upper())
        if field is None:
            continue
        raw = raw.split('#', 1)[0].strip().replace(',', '.')
        try:
            number = float(raw)
        except ValueError as exc:
            raise ValueError(f'{path}:{line_number}: invalid numeric value for {key}: {raw}') from exc
        if not np.isfinite(number):
            raise ValueError(f'{path}:{line_number}: non-finite value for {key}')
        if field in values:
            raise ValueError(f'{path}:{line_number}: duplicated setup field: {key}')
        values[field] = number
    rx, vol, level = (values.get(k) for k in ('rx_gy', 'volume_pct', 'isodose_pct'))
    if rx is not None and rx <= 0:
        raise ValueError('PRESCRIPTION DOSE must be positive and expressed in Gy')
    if vol is not None and not (0 <= vol <= 100):
        raise ValueError('PERCENTAGE PRESCRIPTION VOLUME must be between 0 and 100')
    if level is not None and not (0 <= level <= 100):
        raise ValueError('PERCENTAGE PRESCRIPTION DOSE must be between 0 and 100')
    if rx is None or vol is None or level is None:
        raise ValueError(f'Missing prescription parameter(s) in {path}: '
                         'PRESCRIPTION DOSE, PERCENTAGE PRESCRIPTION VOLUME, '
                         'PERCENTAGE PRESCRIPTION DOSE are all required')
    return rx, vol, level


def parse_args() -> argparse.Namespace:
    pre_parser = argparse.ArgumentParser(add_help=False)
    
    pre_parser.add_argument("--setup-file", type=Path, default=None)
    pre_args, remaining = pre_parser.parse_known_args()
    candidates = [Path(item).parent for item in remaining
                  if item.lower().endswith(".txt") and Path(item).is_file()]
    setup_path = locate_setup(pre_args.setup_file, *candidates)

    
    setup_dose, setup_vol, setup_level = read_setup_file(setup_path)
    default_dose = setup_dose
    default_vol = setup_vol
    default_level = setup_level

    parser = argparse.ArgumentParser(description="Analizza DVH e verifica la prescrizione PTV.")

    parser.add_argument(
        "--setup-file",
        type=Path,
        default=setup_path,
        help="File di prescrizione del run (trovato automaticamente se omesso).",
    )

    parser.add_argument(
        "files",
        nargs="+",
        type=Path,
        help="File ComputeDVH .txt con due colonne: dose (cGy) e volume cumulativo (%%).",
    )
    parser.add_argument(
        "--prescription-gy",
        type=float,
        required=True,
        help=f"Dose prescritta in Gy (default: {default_dose}).",
    )
    parser.add_argument(
        "--target-volume",
        type=float,
        default=default_vol,
        help=f"Volume bersaglio per D% in %% (default: {default_vol}).",
    )
    parser.add_argument(
        "--prescription-level-percent",
        type=float,
        required=True,
        help=f"Percentuale dell'isodose di prescrizione dal setup.out (default: {default_level:g}%%).",
    )
    parser.add_argument(
        "--dmax-threshold",
        type=float,
        default=0.1,
        help="Soglia di volume %% per la stima di Dmax (default: 0.1%%).",
    )
    parser.add_argument(
        "--rois",
        nargs="+",
        default=DEFAULT_ROIS,
        help="ROI aggiuntive per file con nome non standard; di norma sono lette automaticamente dai file *_plan*.txt.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("risultati_dvh.csv"),
        help="File CSV di uscita (default: risultati_dvh.csv).",
    )
    parser.add_argument(
        "--table-file",
        type=Path,
        default=Path("dvh_table.txt"),
        help="File di testo per la tabella (default: dvh_table.txt).",
    )
    parser.add_argument("--plot", action="store_true", help="Salva i grafici PNG.")
    parser.add_argument("--show", action="store_true", help="Mostra i grafici a schermo.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Cartella di uscita per i file generati.",
    )

    args = parser.parse_args()
    # Keep legacy command-line options accepted, but never permit the PTV
    # prescription/marker to disagree with this invocation's setup.out.
    for option, entered, from_setup in (
        ('--prescription-gy', args.prescription_gy, setup_dose),
        ('--target-volume', args.target_volume, setup_vol),
        ('--prescription-level-percent', args.prescription_level_percent, setup_level),
    ):
        if not np.isclose(entered, from_setup, atol=1e-9, rtol=0):
            parser.error(f'{option}={entered:g} conflicts with setup.out '
                         f'{setup_path} ({from_setup:g}); the setup file is authoritative')
    return args


def read_dvh(path: Path) -> tuple[np.ndarray, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f"File non trovato: {path}")

    data = np.loadtxt(path, dtype=float)
    if data.ndim == 1:
        if data.size != 2:
            raise ValueError(f"{path}: expected two numeric DVH columns")
        data = data.reshape(1, 2)

    if data.ndim != 2 or data.shape[1] != 2 or data.shape[0] < 2:
        raise ValueError(f"{path}: formato non valido; servono almeno 2 righe e 2 colonne")

    if not np.all(np.isfinite(data)):
        raise ValueError(f"{path}: presenti valori non numerici o infiniti")

    order = np.argsort(data[:, 0], kind="stable")
    raw_dose = data[order, 0]
    volume = data[order, 1]

    
    dose_gy = raw_dose / 100.0

    if np.min(raw_dose) < 0 or np.min(volume) < -0.05 or np.max(volume) > 100.05:
        raise ValueError(f"{path}: invalid dose/volume values")

    # Se il volume è ordinato in senso crescente (da 0% a 100%), invertilo per DVH cumulativo
    if volume[0] < volume[-1]:
        volume = 100.0 - volume

    # Garantisce che il DVH sia monotono non crescente
    volume = np.maximum.accumulate(volume[::-1])[::-1]

    return dose_gy, volume


def detect_metadata(path: Path, known_rois: list[str]) -> tuple[str, str, str, str]:
    """Identify energy, mode and ANY patient-dependent ROI from ComputeDVH filename."""
    
    match = re.fullmatch(
        r'(?P<roi>.+?)_plan(?P<energy>\d+(?:\.\d+)?)MeV_(?P<mode>CONV|FLASH)',
        path.stem, re.IGNORECASE
    )
    if match:
        energy = match.group('energy')
        modality = match.group('mode').upper()
        structure = match.group('roi')
    else:
        name = path.name
        match_e = re.search(r'(\d+(?:\.\d+)?)MeV', name, re.IGNORECASE)
        match_m = re.search(r'(CONV|FLASH)', name, re.IGNORECASE)
        energy = match_e.group(1) if match_e else '?'
        modality = match_m.group(1).upper() if match_m else '?'
        structure = '?'
        # Compatibilita' con input storici non chiamati *_plan<E>MeV_<MODE>.
        for roi in known_rois:
            pattern = re.escape(roi).replace(r'\ ', r'[\s_\-]?')
            if re.search(pattern, name, re.IGNORECASE):
                structure = roi
                break

    # FIX 5b - La variante e la numerazione sono nella cartella AUTO_sim,
    # non necessariamente nel nome del file DVH (che contiene solo la ROI).
    ancestry = str(path)
    if re.search(r'AUTO_sim\d+(?:\.\d+)?MeV_', ancestry, re.IGNORECASE):
        plan = 'AUTO'
    elif 'MANUAL_SIM' in ancestry.upper():
        plan = 'MANUAL'
    else:
        plan = '?'
    return energy, modality, structure, plan


def volume_at_dose(dose_gy: np.ndarray, volume: np.ndarray, target_dose_gy: float) -> float:
    """Interpola il volume [%] corrispondente a una dose data."""
    return float(np.interp(target_dose_gy, dose_gy, volume, left=volume[0], right=volume[-1]))


def dose_at_volume(dose_gy: np.ndarray, volume: np.ndarray, target_volume: float) -> float:
    """Interpola la dose [Gy] corrispondente a un volume dato (es. D95 o D0.1%)."""
    vol_u, idx = np.unique(volume[::-1], return_index=True)
    dose_u = dose_gy[::-1][idx]
    return float(np.interp(target_volume, vol_u, dose_u))


def weighted_mean_dose(dose_gy: np.ndarray, volume: np.ndarray) -> float:
    """Calcola la dose media integrando la curva DVH (regola trapezoidale)."""
    total_vol = volume[0]
    if total_vol <= 0:
        return float("nan")

    # FIX 6 - NumPy >= 2.0 non ha np.trapz. Il precedente getattr(..., np.trapz)
    # falliva perche' Python valuta il valore di default PRIMA di chiamare getattr.
    # Formula trapezoidale esplicita: compatibile sia con NumPy 1.x sia con 2.x.
    integral = float(np.sum(0.5 * (volume[:-1] + volume[1:]) * np.diff(dose_gy)))
    area = dose_gy[0] * total_vol + integral
    return float(area / total_vol)


def analyze_file(
    path: Path,
    prescription_gy: float,
    target_volume: float,
    prescription_level_percent: float,
    dmax_threshold: float,
    known_rois: list[str],
) -> tuple[Result, np.ndarray, np.ndarray]:
    dose_gy, volume = read_dvh(path)
    energy, modality, structure, plan = detect_metadata(path, known_rois)

    v_cutoff_gy = prescription_gy * (prescription_level_percent / 100.0)
    v_prescription = volume_at_dose(dose_gy, volume, v_cutoff_gy)
    d95_gy = dose_at_volume(dose_gy, volume, target_volume)
    dmax_gy = dose_at_volume(dose_gy, volume, dmax_threshold)
    mean_weighted_gy = weighted_mean_dose(dose_gy, volume)

    is_ptv = structure.upper() == "PTV"
    # FIX 7 - Criterio PTV conforme a setup.out: V(Rx * livello%) >= volume%.
    # Non imporre anche D(volume%) >= 100% Rx: sarebbe un requisito diverso
    # e potrebbe contraddire la tabella PTV nominale della Full Chain.
    # La decisione PTV ufficiale della Full Chain resta quella sui voxel 3D.
    status = (
        ("RAGGIUNTA" if v_prescription >= target_volume else "NON RAGGIUNTA")
        if is_ptv
        else "N/A (organo a rischio)"
    )

    result = Result(
        file=path,
        energy=energy,
        modality=modality,
        structure=structure,
        plan=plan,
        max_tab_gy=float(np.max(dose_gy)),
        dmax_gy=dmax_gy,
        mean_weighted_gy=mean_weighted_gy,
        v_prescription=v_prescription,
        d95_gy=d95_gy,
        status=status,
    )

    return result, dose_gy, volume


def format_table(
    results: list[Result],
    prescription_gy: float,
    target_volume: float,
    prescription_level_percent: float,
) -> str:
    v_cutoff_gy = prescription_gy * (prescription_level_percent / 100.0)
    v_header = f"V{prescription_level_percent:g}% ({v_cutoff_gy:.2f}Gy) [%]"
    d_header = f"D{target_volume:g} [Gy]"

    title_lines = [
        f"PTV prescription: {prescription_gy:g} Gy | V{prescription_level_percent:g}% cutoff = {v_cutoff_gy:.3f} Gy | D{target_volume:g}% target",
        "Estimated Dmax = D0.1% by default (not the maximum voxel dose)",
        "Weighted mean dose = sum(dose x diff_vol) / total_vol",
    ]

    header = (
        f"{'Energy':<9}{'Mode':<7}{'Plan':<8}{'Structure':<15}"
        f"{'Max tab. dose [Gy]':>20}{'Est. Dmax [Gy]':>17}"
        f"{'Wtd mean dose [Gy]':>20}{v_header:>22}{d_header:>13}    Status"
    )

    rows = []
    for r in results:
        status = STATUS_EN.get(r.status, r.status)
        rows.append(
            f"{r.energy + ' MeV':<9}{r.modality:<7}{r.plan:<8}{r.structure:<15}"
            f"{r.max_tab_gy:>20.3f}{r.dmax_gy:>17.3f}{r.mean_weighted_gy:>20.3f}"
            f"{r.v_prescription:>22.3f}{r.d95_gy:>13.3f}    {status}"
        )

    width = max(len(line) for line in [header, *rows, *title_lines])

    lines = [
        "=" * width,
        *title_lines,
        "=" * width,
        header,
        "-" * width,
        *rows,
        "=" * width,
    ]

    return "\n".join(lines) + "\n"


def write_csv(
    results: list[Result],
    output: Path,
    prescription_level_percent: float,
    target_volume: float,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "file",
                "energy_MeV",
                "mode",
                "plan",
                "structure",
                "max_tabulated_dose_Gy",
                "estimated_Dmax_Gy",
                "weighted_mean_dose_Gy",
                f"V{prescription_level_percent:g}_percent",
                # FIX 8 - Se il volume richiesto cambia, cambia anche il nome Dv.
                f"D{target_volume:g}_Gy",
                "PTV_prescription",
            ]
        )
        for r in results:
            writer.writerow(
                [
                    str(r.file),
                    r.energy,
                    r.modality,
                    r.plan,
                    r.structure,
                    f"{r.max_tab_gy:.6f}",
                    f"{r.dmax_gy:.6f}",
                    f"{r.mean_weighted_gy:.6f}",
                    f"{r.v_prescription:.6f}",
                    f"{r.d95_gy:.6f}",
                    STATUS_EN.get(r.status, r.status),
                ]
            )


def make_plots(
    curves: list[tuple[Result, np.ndarray, np.ndarray]],
    output_dir: Path,
    prescription_gy: float,
    target_volume: float,
    prescription_level_percent: float,
    show: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    v_cutoff_gy = prescription_gy * (prescription_level_percent / 100.0)

    energies = sorted(
        {res.energy for res, _, _ in curves},
        key=lambda v: (float(v) if v != "?" else 999),
    )

    for energy in energies:
        plt.figure(figsize=(12, 7))

        for res, dose_gy, volume in curves:
            if res.energy != energy:
                continue
            label = f"{res.plan} - {res.structure} - {res.modality} | Mean dose={res.mean_weighted_gy:.2f} Gy"
            plt.plot(dose_gy, volume, label=label)

        plt.plot(
            v_cutoff_gy,
            target_volume,
            marker="*",
            linewidth=0,
            markersize=14,
            label=f"PTV target V({v_cutoff_gy:.3f} Gy) >= {target_volume:g}%",
        )

        plt.xlabel("Dose (Gy)")
        plt.ylabel("Volume (%)")
        # FIX 9 - Grafici e messaggi del report finale in inglese.
        plt.title(f"DVH comparison {energy} MeV")
        plt.grid(True)
        plt.legend(fontsize=8)
        plt.tight_layout()

        # Nome del PNG invariato per compatibilita' con le analisi precedenti.
        output = output_dir / f"confronto_DVH_{energy}MeV.png"
        plt.savefig(output, dpi=200, bbox_inches="tight")
        print(f"Plot saved: {output}")

        if show:
            plt.show()

        plt.close()


def main() -> int:
    args = parse_args()

    results: list[Result] = []
    curves: list[tuple[Result, np.ndarray, np.ndarray]] = []

    for path in args.files:
        try:
            res, dose_gy, volume = analyze_file(
                path,
                args.prescription_gy,
                args.target_volume,
                args.prescription_level_percent,
                args.dmax_threshold,
                args.rois,
            )
        except Exception as error:
            print(f"ERROR [{path}]: {error}", file=sys.stderr)
            return 1

        results.append(res)
        curves.append((res, dose_gy, volume))

    results.sort(
        key=lambda r: (
            float(r.energy) if r.energy != "?" else 999,
            r.modality,
            r.structure,
            r.plan,
        )
    )

    table_text = format_table(
        results,
        args.prescription_gy,
        args.target_volume,
        args.prescription_level_percent,
    )

    print("\n" + table_text)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    table_path = args.output_dir / args.table_file
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table_path.write_text(table_text, encoding="utf-8")
    print(f"Table saved: {table_path}")

    csv_path = args.output_dir / args.csv
    write_csv(results, csv_path, args.prescription_level_percent, args.target_volume)
    print(f"CSV saved: {csv_path}")

    if args.plot or args.show:
        make_plots(
            curves,
            args.output_dir,
            args.prescription_gy,
            args.target_volume,
            args.prescription_level_percent,
            args.show,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
