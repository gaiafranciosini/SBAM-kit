#!/usr/bin/env python3
"""
Overlay DVH plotter

Esempi:
  python dvh_overlay.py /path/run1 /path/run2 --title "DVH comparison" --recursive --bg white
  python dvh_overlay.py /path/a /path/b --pattern "*.txt" --bg transparent --out dvh.png --xunit Gy --yunit "%"
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot DVH curves from .txt files in one or more folders (overlaid)."
    )
    p.add_argument(
        "folders",
        nargs="+",
        type=Path,
        help="Path(s) alle cartelle contenenti i DVH .txt",
    )
    p.add_argument(
        "--pattern",
        default="*.txt",
        help='Pattern dei file (default: "*.txt")',
    )
    p.add_argument(
        "--recursive",
        action="store_true",
        help="Cerca i file anche nelle sottocartelle",
    )
    p.add_argument(
        "--title",
        default="DVH",
        help="Titolo del plot (default: DVH)",
    )
    p.add_argument(
        "--xlabel",
        default="Dose",
        help="Etichetta asse X (default: Dose)",
    )
    p.add_argument(
        "--ylabel",
        default="Volume",
        help="Etichetta asse Y (default: Volume)",
    )
    p.add_argument(
        "--xunit",
        default="Gy",
        help="Unità asse X (default: Gy). Usata solo per comporre la label, se non già presente.",
    )
    p.add_argument(
        "--yunit",
        default="%",
        help="Unità asse Y (default: %%). Usata solo per comporre la label, se non già presente.",
    )
    p.add_argument(
        "--bg",
        choices=["white", "transparent"],
        default="white",
        help="Sfondo figura: white o transparent (default: white)",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="DPI in salvataggio (default: 200)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Se specificato, salva un PNG in questo path",
    )
    p.add_argument(
        "--max-curves",
        type=int,
        default=0,
        help="Limita il numero massimo di curve (0 = nessun limite). Utile se ci sono tantissimi file.",
    )
    p.add_argument(
        "--label-mode",
        choices=["file", "folder_file", "folder"],
        default="folder_file",
        help="Come costruire la label in legenda (default: folder_file).",
    )
    return p.parse_args()


def find_txt_files(folders: List[Path], pattern: str, recursive: bool) -> List[Path]:
    files: List[Path] = []
    for folder in folders:
        if not folder.exists() or not folder.is_dir():
            raise FileNotFoundError(f"Cartella non valida: {folder}")
        if recursive:
            # Path.rglob non supporta pattern con brace/glob complessi come glob.glob,
            # ma per "*.txt" va benissimo.
            files.extend(sorted(folder.rglob(pattern)))
        else:
            files.extend(sorted(folder.glob(pattern)))
    # Tieni solo file reali
    files = [f for f in files if f.is_file()]
    return files


def read_dvh_txt(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Prova a leggere un DVH da txt con due colonne numeriche.
    - Ignora righe vuote e righe non numeriche (header/commenti).
    - Supporta separatori: spazio/tab/virgola/;.
    """
    xs: List[float] = []
    ys: List[float] = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            # rimuovi commenti inline comuni
            for c in ("#", "//"):
                if c in s:
                    s = s.split(c, 1)[0].strip()
            if not s:
                continue

            # normalizza separatori
            s = s.replace(",", " ").replace(";", " ")
            parts = s.split()
            if len(parts) < 2:
                continue

            try:
                x = float(parts[0])
                y = float(parts[1])
            except ValueError:
                continue

            xs.append(x)
            ys.append(y)

    if len(xs) < 2:
        raise ValueError(f"File DVH non leggibile o con poche righe numeriche: {path}")

    xarr = np.asarray(xs, dtype=float)
    yarr = np.asarray(ys, dtype=float)

    # Ordina per dose crescente (capita che alcuni file non siano ordinati)
    order = np.argsort(xarr)
    xarr = xarr[order]
    yarr = yarr[order]
    return xarr, yarr


def make_label(path: Path, mode: str) -> str:
    if mode == "file":
        return path.stem
    if mode == "folder":
        return path.parent.name
    # folder_file
    return f"{path.parent.name}/{path.stem}"


def main() -> None:
    args = parse_args()

    files = find_txt_files(args.folders, args.pattern, args.recursive)
    if not files:
        raise SystemExit("Nessun file trovato con il pattern richiesto nelle cartelle fornite.")

    if args.max_curves and args.max_curves > 0:
        files = files[: args.max_curves]

    # Labels asse con unità (se l'utente non le ha già messe)
    xlabel = args.xlabel
    ylabel = args.ylabel
    if args.xunit and "(" not in xlabel and "[" not in xlabel:
        xlabel = f"{xlabel} ({args.xunit})"
    if args.yunit and "(" not in ylabel and "[" not in ylabel:
        ylabel = f"{ylabel} ({args.yunit})"

    fig, ax = plt.subplots(figsize=(9.5, 6.0))

    # Background
    if args.bg == "transparent":
        fig.patch.set_alpha(0.0)
        ax.set_facecolor((0, 0, 0, 0))
    else:
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

    plotted = 0
    errors: List[str] = []
    xmax=0
    for fpath in files:
        try:
            x, y = read_dvh_txt(fpath)
            mask=y > 0.001
            if (max(x[mask])>xmax):
              xmax = max(x[mask])
              print(xmax)
            label = make_label(fpath, args.label_mode)
            ax.plot(x[mask], y[mask], linewidth=2.0, label=label)
            plotted += 1
        except Exception as e:
            errors.append(f"{fpath}: {e}")

    if plotted == 0:
        msg = "Nessun DVH plottato. Errori:\n" + "\n".join(errors[:30])
        raise SystemExit(msg)

    # Styling: assi ben visibili
    ax.set_title(args.title, fontsize=18, pad=12)
    ax.set_xlabel(xlabel, fontsize=15, labelpad=10)
    ax.set_ylabel(ylabel, fontsize=15, labelpad=10)
    ax.set_xlim(0, 1.01*xmax)
    ax.tick_params(axis="both", which="major", labelsize=13)
    ax.grid(True, which="both", alpha=0.25)

    # Legenda fuori a destra
    # Riduciamo l'area del grafico per far spazio
    fig.subplots_adjust(right=0.72)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=11,
    )

    # Tight layout (dopo adjust, per evitare tagli)
    plt.tight_layout()

    if args.out is not None:
        out_path = args.out
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            out_path,
            dpi=args.dpi,
            bbox_inches="tight",
            transparent=(args.bg == "transparent"),
        )
        print(f"[OK] Salvato: {out_path}")

    if errors:
        print("\n[WARN] Alcuni file sono stati saltati:")
        for e in errors[:30]:
            print(" -", e)
        if len(errors) > 30:
            print(f" ... (+{len(errors)-30} altri)")

#    plt.show()


if __name__ == "__main__":
    main()
