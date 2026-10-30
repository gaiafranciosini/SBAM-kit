#!/usr/bin/env python3
"""Plot a single ComputeDVH set with the prescription marker from setup.out.

ComputeDVH TXT files contain dose in cGy and cumulative volume in percent.
"""

from pathlib import Path
import argparse
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def locate_setup(explicit: Path | None, *hints: Path) -> Path:
    """Find the current plan's setup.out without a Full Chain CLI change."""
    if explicit is not None:
        return explicit.expanduser().resolve()
    # MODIFICA: la Full Chain ORIGINALE lancia plotDVH.py da sim7MeV/sim9MeV
    # e NON passa --setup-file. Cerchiamo prima nella cartella corrente e poi
    # nelle cartelle superiori. Dentro results/<piano> privilegiamo invece
    # prescription_inputs/setup.out, che conserva la prescrizione del piano.
    # PRIMA: non esisteva questa verifica; la stella usava il default 100 cGy.
    checked = set()
    # Il percorso della cartella DVH ha precedenza: utile quando si apre un
    # vecchio piano da results/ dopo che il setup.out principale e' cambiato.
    for origin in (*hints, Path.cwd()):
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
    raise FileNotFoundError("setup.out not found near the DVH plot; refusing a fixed PTV star")


def read_prescription(path: Path) -> tuple[float, float, float]:
    """Return (prescription Gy, prescribed isodose %, required PTV volume %)."""
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
            raise ValueError(f"{path}:{line_no}: duplicated {raw_key.strip()}")
        number = float(raw_value.split("#", 1)[0].strip().replace(",", "."))
        if not math.isfinite(number):
            raise ValueError(f"{path}:{line_no}: non-finite prescription parameter")
        values[name] = number
    missing = sorted(set(keys.values()) - set(values))
    if missing:
        raise ValueError(f"Missing prescription parameters in {path}: {', '.join(missing)}")
    dose, isodose, volume = (values[k] for k in ("dose", "isodose", "volume"))
    if dose <= 0 or not 0 <= isodose <= 100 or not 0 <= volume <= 100:
        raise ValueError(f"Invalid prescription values in {path}: {values}")
    return dose, isodose, volume


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot ComputeDVH and setup.out PTV target")
    parser.add_argument("-label1", type=str, default="", help="Label of DVH TXT files")
    parser.add_argument("-dir1", type=Path, default=Path("."), help="DVH directory")
    parser.add_argument("-roi", type=str, nargs="+", required=True, help="ROI names")
    # PRIMA (versione corretta precedente): --setup-file era required=True;
    # la Full Chain originale non lo passa e l'esecuzione si sarebbe fermata.
    # MODIFICA: argomento opzionale, ricerca automatica verificata sotto.
    parser.add_argument("--setup-file", type=Path, default=None,
                        help="Optional explicit setup.out; otherwise find it automatically")
    # Accept the historic parameter so existing direct calls do not error.
    # Its numeric value is deliberately ignored: prescription comes only from setup.out.
    parser.add_argument("-Dgoal", type=float, default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    setup_file = locate_setup(args.setup_file, args.dir1)
    rx_gy, isodose_pct, volume_pct = read_prescription(setup_file)
    target_gy = rx_gy * isodose_pct / 100.0
    target_cgy = target_gy * 100.0

    palette = ["r", "g", "b", "k", "y", "c", "m", "tab:orange", "tab:purple",
               "tab:pink", "tab:brown", "yellowgreen", "slategray", "hotpink",
               "tan", "lightcoral", "tomato"]
    plt.rcParams.update({"font.size": 11})
    fig, ax = plt.subplots(figsize=(6, 6))
    for i, roi in enumerate(args.roi):
        path = args.dir1 / f"{roi}{args.label1}.txt"
        x, y = np.loadtxt(path, comments=["#", "idx"], usecols=(0, 1), unpack=True)
        ax.plot(x, y, color=palette[i % len(palette)], label=roi, linewidth=1.8)

    # RIGHE ORIGINALI (non usarle: numeri fissi non legati al setup.out):
    # plt.plot([args.Dgoal*1],[95],'r*',lw=5,markerfacecolor='r', markeredgewidth=1, markeredgecolor='yellow',markersize=17)
    # star = Line2D([], [],marker='*',linewidth='0',markerfacecolor='r',markeredgewidth=1,markeredgecolor='yellow' ,markersize=17, label='V95% 95%')
    # MODIFICA: X = dose prescritta * isodose% (convertita in cGy), Y = volume%.
    target_label = f"PTV target V({target_gy:.3f} Gy) >= {volume_pct:g}%"
    ax.plot([target_cgy], [volume_pct], "r*", lw=0, markerfacecolor="r",
            markeredgecolor="yellow", markeredgewidth=1, markersize=17,
            zorder=10)
    star = Line2D([], [], marker="*", linestyle="None", markerfacecolor="r",
                  markeredgecolor="yellow", markeredgewidth=1, markersize=17,
                  label=target_label)

    ax.set_xlabel("Dose [cGy]", fontsize=11)
    ax.set_ylabel("Volume [%]", fontsize=11)
    ax.grid(True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.19,
                     box.width, box.height * 0.73])
    legend_curves = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20),
                              fancybox=True, shadow=True, ncol=3)
    ax.add_artist(legend_curves)
    ax.legend(handles=[star], bbox_to_anchor=(0, 1.02, 1, 0.2),
              loc="lower left", mode="expand", borderaxespad=0, ncol=1)

    fig.savefig(args.label1 + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Prescription marker: ({target_cgy:.3f} cGy, {volume_pct:g}%) from {setup_file}")


if __name__ == "__main__":
    main()




