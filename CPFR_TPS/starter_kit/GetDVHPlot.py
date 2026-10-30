#!/usr/bin/env python3
"""Overlay ComputeDVH curves and mark the PTV prescription from setup.out.

The numeric ComputeDVH TXT dose column is in cGy. --xunit controls the
displayed unit and transforms BOTH curves and the prescription marker.

Tail trimming (display only, files on disk are never modified):
  --xmax X         do not draw doses above X (in --xunit)
  --xmax-rx F      do not draw doses above F * prescription dose (e.g. 1.2)
  --min-volume V   do not draw points with volume below V % (default 0.001)
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np


def locate_setup(explicit: Path | None, *folders: Path) -> Path:
    if explicit is not None:
        return explicit.expanduser().resolve()
    # MODIFICA: la Full Chain ORIGINALE chiama GetDVHPlot.py dalla cartella
    # del paziente, senza --setup-file. Verifichiamo che esista setup.out,
    # cercandolo in cwd o nelle directory dei dati. Nei results di un piano
    # preferiamo prescription_inputs/setup.out (copia della prescrizione).
    # PRIMA: GetDVHPlot.py non leggeva affatto setup.out e non disegnava stelle.
    checked = set()
    # Le cartelle DVH hanno precedenza sul cwd per le analisi storiche in results/.
    for origin in (*folders, Path.cwd()):
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
    raise FileNotFoundError("setup.out not found for DVH overlay; no default star is allowed")


def read_prescription(path: Path) -> tuple[float, float, float]:
    keys = {
        "PRESCRIPTION DOSE": "dose",
        "PERCENTAGE PRESCRIPTION DOSE": "isodose",
        "PERCENTAGE PRESCRIPTION VOLUME": "volume",
    }
    if not path.is_file():
        raise FileNotFoundError(f"Prescription file not found: {path}")
    values = {}
    for line_no, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if ":" not in line:
            continue
        raw_key, raw_value = line.split(":", 1)
        name = keys.get(raw_key.strip().upper())
        if name is None:
            continue
        if name in values:
            raise ValueError(f"{path}:{line_no}: duplicate {raw_key.strip()}")
        number = float(raw_value.split("#", 1)[0].strip().replace(",", "."))
        if not math.isfinite(number):
            raise ValueError(f"{path}:{line_no}: non-finite prescription parameter")
        values[name] = number
    missing = sorted(set(keys.values()) - set(values))
    if missing:
        raise ValueError(f"Missing prescription parameters in {path}: {', '.join(missing)}")
    dose, iso, volume = (values[k] for k in ("dose", "isodose", "volume"))
    if dose <= 0 or not 0 <= iso <= 100 or not 0 <= volume <= 100:
        raise ValueError(f"Invalid prescription parameters in {path}: {values}")
    return dose, iso, volume


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Overlay DVH curves from one or more folders")
    p.add_argument("folders", nargs="+", type=Path, help="Folders with ComputeDVH TXT files")
    p.add_argument("--pattern", default="*.txt")
    p.add_argument("--recursive", action="store_true")
    p.add_argument("--title", default="DVH")
    p.add_argument("--xlabel", default="Dose")
    p.add_argument("--ylabel", default="Volume")
    p.add_argument("--xunit", choices=("Gy", "cGy"), default="Gy",
                   help="Displayed dose unit; input TXT files are always in cGy")
    p.add_argument("--yunit", default="%")
    p.add_argument("--bg", choices=("white", "transparent"), default="white")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--out", type=Path, default=None)
    p.add_argument("--max-curves", type=int, default=0)
    p.add_argument("--label-mode", choices=("file", "folder_file", "folder"),
                   default="folder_file")
    # PRIMA (versione precedente): --setup-file era required=True.
    # MODIFICA: opzionale per mantenere identica la chiamata della Full Chain.
    p.add_argument("--setup-file", type=Path, default=None,
                   help="Optional exact setup.out; otherwise find it automatically")
    # MODIFICA: soglie per tagliare la coda delle curve (solo visualizzazione).
    p.add_argument("--xmax", type=float, default=None,
                   help="Max displayed dose, in --xunit; curve points beyond it are cut "
                        "(takes precedence over --xmax-rx)")
    p.add_argument("--xmax-rx", type=float, default=None,
                   help="Max displayed dose as a multiple of the prescription dose (e.g. 1.2)")
    p.add_argument("--min-volume", type=float, default=0.001,
                   help="Drop curve points with volume (%%) below this value (default 0.001)")
    return p.parse_args()


def find_txt_files(folders: List[Path], pattern: str, recursive: bool) -> List[Path]:
    files: List[Path] = []
    for folder in folders:
        if not folder.exists() or not folder.is_dir():
            raise FileNotFoundError(f"Invalid DVH directory: {folder}")
        files.extend(sorted(folder.rglob(pattern) if recursive else folder.glob(pattern)))
    return [f for f in files if f.is_file()]


def read_dvh_txt(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    xs: List[float] = []
    ys: List[float] = []
    with path.open("r", encoding="utf-8", errors="ignore") as stream:
        for line in stream:
            s = line.strip()
            if not s:
                continue
            for comment in ("#", "//"):
                if comment in s:
                    s = s.split(comment, 1)[0].strip()
            if not s:
                continue
            s = s.replace(",", " ").replace(";", " ")
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                x, y = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            if not math.isfinite(x) or not math.isfinite(y):
                raise ValueError(f"Non-finite numeric value in {path}")
            xs.append(x)
            ys.append(y)
    if len(xs) < 2:
        raise ValueError(f"Unreadable or too-short DVH file: {path}")
    xarr, yarr = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    order = np.argsort(xarr)
    return xarr[order], yarr[order]


def make_label(path: Path, mode: str) -> str:
    if mode == "file":
        return path.stem
    if mode == "folder":
        return path.parent.name
    return f"{path.parent.name}/{path.stem}"


def main() -> None:
    args = parse_args()
    setup_file = locate_setup(args.setup_file, *args.folders)
    rx_gy, iso_pct, volume_pct = read_prescription(setup_file)
    target_gy = rx_gy * iso_pct / 100.0
    # ComputeDVH numerical files contain cGy, regardless of desired plot units.
    conversion = 0.01 if args.xunit == "Gy" else 1.0
    target_x = target_gy if args.xunit == "Gy" else target_gy * 100.0

    # MODIFICA: soglia massima in dose (assoluta o relativa alla prescrizione),
    # espressa nell'unita' mostrata (--xunit).
    rx_x = rx_gy if args.xunit == "Gy" else rx_gy * 100.0
    dose_cap = None
    if args.xmax is not None:
        dose_cap = args.xmax
    elif args.xmax_rx is not None:
        dose_cap = args.xmax_rx * rx_x
    if dose_cap is not None and dose_cap <= 0:
        raise SystemExit("--xmax / --xmax-rx must be positive")
    if args.min_volume < 0:
        raise SystemExit("--min-volume must be >= 0")

    files = find_txt_files(args.folders, args.pattern, args.recursive)
    if not files:
        raise SystemExit("No matching DVH TXT files found")
    if args.max_curves > 0:
        files = files[:args.max_curves]

    xlabel, ylabel = args.xlabel, args.ylabel
    if "(" not in xlabel and "[" not in xlabel:
        xlabel = f"{xlabel} ({args.xunit})"
    if "(" not in ylabel and "[" not in ylabel:
        ylabel = f"{ylabel} ({args.yunit})"

    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    if args.bg == "transparent":
        fig.patch.set_alpha(0.0)
        ax.set_facecolor((0, 0, 0, 0))
    else:
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

    colors = plt.cm.tab20(np.linspace(0, 1, len(files)))
    plotted, errors, xmax = 0, [], target_x
    for i, path in enumerate(files):
        try:
            x, y = read_dvh_txt(path)
            # RIGA ORIGINALE: mask = y > 0.001
            # MODIFICA: soglia in volume configurabile (default identico) e
            # soglia opzionale in dose. Il DVH e' cumulativo, quindi il taglio
            # non deforma il resto della curva.
            display_all = x * conversion
            mask = y > args.min_volume
            if dose_cap is not None:
                mask &= display_all <= dose_cap
            if not np.any(mask):
                raise ValueError("No points left after dose/volume thresholds")
            display_x = display_all[mask]
            xmax = max(xmax, float(np.max(display_x)))
            ax.plot(display_x, y[mask], linewidth=2.0,
                    label=make_label(path, args.label_mode), color=colors[i])
            plotted += 1
        except Exception as exc:
            errors.append(f"{path}: {exc}")

    if plotted == 0:
        raise SystemExit("No DVH curves plotted: " + "; ".join(errors[:30]))

    # PRIMA: nessuna stella veniva disegnata in GetDVHPlot.py.
    # MODIFICA: la prescrizione verificata sopra governa la posizione della
    # stella; su asse Gy e cGy si usa la stessa conversione delle curve.
    ax.scatter([target_x], [volume_pct], marker="*", s=330,
               color="red", edgecolors="gold", linewidths=1,
               zorder=10, label=f"PTV target V({target_gy:.3f} Gy) >= {volume_pct:g}%")
    ax.set_title(args.title, fontsize=18, pad=12)
    ax.set_xlabel(xlabel, fontsize=15, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=15, labelpad=10)
    # MODIFICA: con soglia in dose l'asse X si ferma alla soglia (ma la stella
    # resta sempre visibile).
    if dose_cap is not None:
        xmax = max(dose_cap, target_x)
    ax.set_xlim(0, max(1.01 * xmax, 1e-9))
    ax.tick_params(axis="both", which="major", labelsize=13)
    ax.grid(True, which="both", alpha=0.25)
    fig.subplots_adjust(right=0.72)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
              fontsize=11)
    plt.tight_layout()

    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight",
                    transparent=(args.bg == "transparent"))
        print(f"[OK] Saved: {args.out}")
    print(f"Prescription marker: ({target_x:.3f} {args.xunit}, {volume_pct:g}%) "
          f"from {setup_file}")
    if errors:
        print("[WARN] Some DVH files were skipped:")
        for error in errors[:30]:
            print(" -", error)
    plt.close(fig)


if __name__ == "__main__":
    main()
