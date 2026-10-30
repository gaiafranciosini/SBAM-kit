 #!/usr/bin/env python3
"""
viewer_zoom_isodose.py

Viewer CT + ROI + DOSE con:
- CT in scala di grigi
- ROI come soli contorni
- dose sovrapposta con trasparenza regolabile
- modalità TPS-like con colori discreti
- curve di isodose sovrapposte
- tre viste ortogonali: assiale, coronale, sagittale
- slider X/Y/Z e trasparenza dose
- zoom con rotella, pan con tasto destro, reset con tasto 0
- resampling automatico di ROI e dose sulla griglia CT

La modalità TPS-like usa livelli di dose RELATIVA rispetto a Dref.
Se --dose-reference non viene specificato:
  * se la mappa sembra già normalizzata (max <= 2), Dref = 1
  * altrimenti Dref = dose massima della mappa

Esempi:
  # mappa già normalizzata 0..1
  python3 viewer_zoom_isodose.py --ct CT_plan.mhd \
      --roi PTV_plan.mhd Osso_plan.mhd Muscolo_plan.mhd \
      --labels PTV Osso Muscolo \
      --dose isodosemap.mhd --dose-reference 1

  # dose assoluta in Gy, isodosi riferite alla prescrizione 22.5 Gy
  python3 viewer_zoom_isodose.py --ct CT_plan.mhd \
      --roi PTV_plan.mhd --labels PTV \
      --dose dose.mhd --dose-reference 22.5
"""

import argparse
from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.colors import Normalize, ListedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from skimage import measure


# =============================================================================
# DATI / IO
# =============================================================================

@dataclass
class Volume:
    img: sitk.Image
    arr: np.ndarray                     # shape (Z, Y, X)
    spacing: Tuple[float, float, float] # SITK: (X,Y,Z), mm
    size: Tuple[int, int, int]          # (X,Y,Z)


def read_mhd(path: str) -> Volume:
    img = sitk.ReadImage(path)
    arr = sitk.GetArrayFromImage(img)
    sx, sy, sz = img.GetSpacing()
    return Volume(
        img=img,
        arr=arr,
        spacing=(sx, sy, sz),
        size=img.GetSize(),
    )


def resample_like(
    moving_img: sitk.Image,
    reference_img: sitk.Image,
    is_label: bool = False,
) -> sitk.Image:
    """Resample moving_img sulla griglia della CT/reference."""
    interp = sitk.sitkNearestNeighbor if is_label else sitk.sitkLinear

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(reference_img)
    resampler.SetInterpolator(interp)
    resampler.SetDefaultPixelValue(0)
    resampler.SetTransform(sitk.Transform())
    return resampler.Execute(moving_img)


def same_grid(a: sitk.Image, b: sitk.Image) -> bool:
    return (
        a.GetSize() == b.GetSize()
        and a.GetSpacing() == b.GetSpacing()
        and a.GetOrigin() == b.GetOrigin()
        and a.GetDirection() == b.GetDirection()
    )


# =============================================================================
# GEOMETRIA / SLICE
# =============================================================================

def window_auto(
    arr: np.ndarray,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    p_low: float = 1,
    p_high: float = 99,
) -> Tuple[float, float]:
    if vmin is not None and vmax is not None:
        return float(vmin), float(vmax)

    arr_f = arr[np.isfinite(arr)]
    if arr_f.size == 0:
        return 0.0, 1.0

    lo, hi = np.percentile(arr_f, [p_low, p_high])
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def build_extent(vol: Volume, plane: str) -> Tuple[float, float, float, float]:
    """
    Mantiene la convenzione del viewer originale:
    coordinate visualizzate in mm a partire da 0.
    """
    sx, sy, sz = vol.spacing
    nx, ny, nz = vol.size

    if plane == "axial":
        return (0, nx * sx, 0, ny * sy)       # X,Y
    if plane == "coronal":
        return (0, nx * sx, 0, nz * sz)       # X,Z
    if plane == "sagittal":
        return (0, ny * sy, 0, nz * sz)       # Y,Z

    raise ValueError(f"Piano non riconosciuto: {plane}")


def extract_slices(arr: np.ndarray, i: int, j: int, k: int):
    """arr ha ordine (Z,Y,X)."""
    axial = arr[k, :, :]      # XY a Z=k
    coronal = arr[:, j, :]    # XZ a Y=j
    sagittal = arr[:, :, i]   # YZ a X=i
    return axial, coronal, sagittal


# =============================================================================
# ROI
# =============================================================================

def generate_colors(n: int) -> List[Tuple[float, float, float]]:
    """Colori ben separati per i contorni delle ROI."""
    if n <= 0:
        return []

    cmap = plt.get_cmap("tab10" if n <= 10 else "tab20")
    return [cmap(i % cmap.N)[:3] for i in range(n)]


def contours_phys(
    mask2d: np.ndarray,
    plane: str,
    spacing: Tuple[float, float, float],
):
    """Estrae i contorni della ROI e li converte da pixel a mm."""
    sx, sy, sz = spacing
    cs = measure.find_contours(mask2d.astype(float), level=0.5)

    out = []
    for c in cs:
        rows = c[:, 0]
        cols = c[:, 1]

        if plane == "axial":
            x_mm = cols * sx
            y_mm = rows * sy
        elif plane == "coronal":
            x_mm = cols * sx
            y_mm = rows * sz
        elif plane == "sagittal":
            x_mm = cols * sy
            y_mm = rows * sz
        else:
            raise ValueError(f"Piano non riconosciuto: {plane}")

        out.append(np.vstack([x_mm, y_mm]).T)

    return out


def draw_roi_contours(ax, contours_xy, color, lw=1.8):
    lines = []
    for xy in contours_xy:
        line, = ax.plot(
            xy[:, 0],
            xy[:, 1],
            color=color,
            linewidth=lw,
            zorder=4,
        )
        lines.append(line)
    return lines


# =============================================================================
# COLORMAP TPS-LIKE E ISODOSE
# =============================================================================

# Colori presi dalla logica TPS-like del secondo viewer.
TPS_COLORS = [
    [102/255, 255/255, 255/255],
    [0/255,   255/255,   0/255],
    [197/255, 248/255, 190/255],
    [0/255,     0/255, 255/255],
    [122/255,   0/255, 255/255],
    [255/255, 125/255,   0/255],
    [255/255, 255/255,   0/255],
    [255/255,   0/255,   0/255],
    [255/255, 184/255, 205/255],
    [197/255,   0/255, 190/255],
    [100/255, 144/255, 255/255],
]

DEFAULT_ISODOSE_LEVELS = [
    0.25, 0.50, 0.60, 0.70, 0.80,
    0.90, 0.95, 1.00, 1.05, 1.07, 1.10,
]


def build_tps_colormap(levels: np.ndarray):
    """
    Crea una colormap discreta TPS-like.
    I livelli sono frazioni di Dref: 0.95 = 95% di Dref.
    """
    if np.any(np.diff(levels) <= 0):
        raise ValueError("--isodose-levels deve contenere valori strettamente crescenti")

    if len(levels) > len(TPS_COLORS):
        # Se l'utente chiede più livelli dei colori originali,
        # campioniamo 'turbo' mantenendo una buona separazione cromatica.
        base = plt.get_cmap("turbo")
        cols = [base(x) for x in np.linspace(0.05, 0.95, len(levels))]
    else:
        cols = TPS_COLORS[:len(levels)]

    cmap = ListedColormap(cols, name="tps_like")
    cmap.set_bad((0, 0, 0, 0))  # voxel mascherati completamente trasparenti

    # Stessa idea del secondo codice: intervalli di dose associati a colori fissi.
    norm = BoundaryNorm(levels, cmap.N, clip=True)
    return cmap, norm


def relative_dose(arr: np.ndarray, dref: float) -> np.ndarray:
    if dref <= 0:
        raise ValueError("Dref deve essere > 0")
    return arr.astype(float) / float(dref)


def mask_relative_dose(rel: np.ndarray, cutoff: float):
    out = np.ma.masked_invalid(rel)
    if cutoff is not None:
        out = np.ma.masked_less(out, cutoff)
    return out


def draw_isodose_curves(
    ax,
    rel_slice: np.ndarray,
    extent,
    levels: np.ndarray,
    cmap,
    norm,
):
    """
    Disegna curve di isodose ai livelli richiesti.
    95% e 100% vengono leggermente ispessite.
    """
    finite = rel_slice[np.isfinite(rel_slice)]
    if finite.size == 0:
        return None

    # contour() richiede almeno un livello dentro il range della slice.
    lo, hi = float(np.min(finite)), float(np.max(finite))
    active_levels = np.array([x for x in levels if lo <= x <= hi], dtype=float)
    if active_levels.size == 0:
        return None

    curve_colors = [cmap(norm(x)) for x in active_levels]
    linewidths = [
        2.5 if (np.isclose(x, 0.95) or np.isclose(x, 1.00)) else 1.35
        for x in active_levels
    ]

    return ax.contour(
        rel_slice,
        levels=active_levels,
        colors=curve_colors,
        linewidths=linewidths,
        alpha=1.0,
        origin="lower",
        extent=extent,
        zorder=3,
    )


def remove_contour_set(cs):
    if cs is None:
        return

    # Matplotlib recente
    try:
        cs.remove()
        return
    except Exception:
        pass

    # Compatibilità con versioni più vecchie
    try:
        for coll in cs.collections:
            coll.remove()
    except Exception:
        pass


# =============================================================================
# CROSSHAIR
# =============================================================================

def draw_crosshair(ax, x, y):
    for art in list(ax.lines):
        if getattr(art, "_is_cross", False):
            art.remove()

    h = ax.axhline(y, linewidth=0.8, color="white", alpha=0.55, zorder=5)
    v = ax.axvline(x, linewidth=0.8, color="white", alpha=0.55, zorder=5)
    setattr(h, "_is_cross", True)
    setattr(v, "_is_cross", True)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Viewer CT + ROI + dose con colormap TPS-like e curve di isodose"
    )

    parser.add_argument("--ct", required=True, help="CT .mhd")
    parser.add_argument(
        "--roi",
        nargs="*",
        default=[],
        help="ROI .mhd; vengono mostrate come contorni",
    )
    parser.add_argument(
        "--labels",
        nargs="*",
        default=None,
        help="Etichette ROI, nello stesso ordine di --roi",
    )

    parser.add_argument(
        "--window",
        type=float,
        nargs=2,
        metavar=("VMIN", "VMAX"),
        default=None,
        help="Finestra CT",
    )

    parser.add_argument("--dose", default=None, help="Mappa dose .mhd")
    parser.add_argument(
        "--dose-reference",
        type=float,
        default=None,
        help=(
            "Dref per le isodosi. Esempio: 1 per mappa normalizzata; "
            "22.5 per dose assoluta riferita a 22.5 Gy"
        ),
    )
    parser.add_argument(
        "--dose-alpha",
        type=float,
        default=0.45,
        help="Trasparenza dose iniziale [0-1], default 0.45",
    )
    parser.add_argument(
        "--dose-cutoff",
        type=float,
        default=0.25,
        help=(
            "Nasconde la dose sotto questa frazione di Dref. "
            "Default 0.25, come la modalità TPS-like del secondo viewer"
        ),
    )

    parser.add_argument(
        "--dose-mode",
        choices=["tps", "continuous"],
        default="tps",
        help="tps = colori discreti; continuous = colormap continua",
    )
    parser.add_argument(
        "--dose-cmap",
        default="turbo",
        help="Colormap usata soltanto con --dose-mode continuous",
    )
    parser.add_argument(
        "--dose-window",
        type=float,
        nargs=2,
        metavar=("DMIN", "DMAX"),
        default=None,
        help="Range nativo dose usato soltanto con --dose-mode continuous",
    )

    parser.add_argument(
        "--isodose-levels",
        type=float,
        nargs="+",
        default=DEFAULT_ISODOSE_LEVELS,
        help=(
            "Livelli relativi a Dref. Default: "
            "0.25 0.50 0.60 0.70 0.80 0.90 0.95 1.00 1.05 1.07 1.10"
        ),
    )
    parser.add_argument(
        "--no-isodose-lines",
        action="store_true",
        help="Non disegnare le curve di isodose",
    )

    parser.add_argument("--title", default="Viewer CT + ROI + Isodose")

    args = parser.parse_args()

    if not (0.0 <= args.dose_alpha <= 1.0):
        raise SystemExit("--dose-alpha deve essere tra 0 e 1")

    if args.labels is not None and len(args.labels) != len(args.roi):
        raise SystemExit(
            f"Numero labels ({len(args.labels)}) diverso dal numero ROI ({len(args.roi)})."
        )

    labels = (
        args.labels
        if args.labels is not None
        else [f"ROI {i+1}" for i in range(len(args.roi))]
    )

    # -------------------------------------------------------------------------
    # CT
    # -------------------------------------------------------------------------
    ct = read_mhd(args.ct)
    vmin, vmax = window_auto(
        ct.arr,
        *(args.window if args.window else (None, None)),
    )

    # -------------------------------------------------------------------------
    # ROI -> griglia CT
    # -------------------------------------------------------------------------
    roi_arrays = []
    for path in args.roi:
        rimg = sitk.ReadImage(path)
        if not same_grid(rimg, ct.img):
            print(f"[ROI] resampling NN -> CT: {path}")
            rimg = resample_like(rimg, ct.img, is_label=True)
        roi_arrays.append(sitk.GetArrayFromImage(rimg))

    # -------------------------------------------------------------------------
    # DOSE -> griglia CT
    # -------------------------------------------------------------------------
    dose_vol: Optional[Volume] = None
    dose_rel = None
    dref = None

    tps_levels = np.array(args.isodose_levels, dtype=float)
    tps_cmap = tps_norm = None

    continuous_cmap = continuous_norm = None

    if args.dose is not None:
        dimg = sitk.ReadImage(args.dose)
        if not same_grid(dimg, ct.img):
            print(f"[DOSE] resampling lineare -> CT: {args.dose}")
            dimg = resample_like(dimg, ct.img, is_label=False)

        darr = sitk.GetArrayFromImage(dimg).astype(float)
        dose_vol = Volume(dimg, darr, ct.spacing, ct.size)

        finite = darr[np.isfinite(darr)]
        if finite.size == 0:
            raise SystemExit("La mappa di dose non contiene valori finiti.")

        dmax_native = float(np.max(finite))
        dmin_native = float(np.min(finite))

        if args.dose_reference is not None:
            dref = float(args.dose_reference)
        else:
            # Auto: una isodose map normalizzata tipicamente vive circa in [0,1].
            if dmin_native >= 0 and dmax_native <= 2.0:
                dref = 1.0
            else:
                dref = dmax_native

        if dref <= 0:
            raise SystemExit("Dref non valido: deve essere > 0.")

        dose_rel = relative_dose(darr, dref)

        print(f"[DOSE] min/max nativi = {dmin_native:.6g} / {dmax_native:.6g}")
        print(f"[DOSE] Dref usato = {dref:.6g}")
        print(
            "[ISO] livelli = "
            + ", ".join(f"{100*x:g}%" for x in tps_levels)
        )

        tps_cmap, tps_norm = build_tps_colormap(tps_levels)

        if args.dose_mode == "continuous":
            continuous_cmap = plt.get_cmap(args.dose_cmap).copy()
            continuous_cmap.set_bad((0, 0, 0, 0))

            if args.dose_window is not None:
                dmin, dmax = map(float, args.dose_window)
            else:
                dmin, dmax = dmin_native, dmax_native
                if dmax <= dmin:
                    dmax = dmin + 1.0

            continuous_norm = Normalize(vmin=dmin, vmax=dmax)

    # -------------------------------------------------------------------------
    # FIGURA
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(14, 5),
        gridspec_kw={"width_ratios": [1, 1, 1.05]},
    )
    ax_axial, ax_cor, ax_sag = axes

    ax_axial.set_title("Trasversale (assiale) — Z")
    ax_cor.set_title("Coronale — Y")
    ax_sag.set_title("Sagittale — X")

    ext_ax = build_extent(ct, "axial")
    ext_co = build_extent(ct, "coronal")
    ext_sa = build_extent(ct, "sagittal")

    nx, ny, nz = ct.size
    i, j, k = nx // 2, ny // 2, nz // 2

    # CT
    ct_ax, ct_co, ct_sa = extract_slices(ct.arr, i, j, k)

    im_ct_ax = ax_axial.imshow(
        ct_ax, cmap="gray", vmin=vmin, vmax=vmax,
        origin="lower", extent=ext_ax, zorder=0
    )
    im_ct_co = ax_cor.imshow(
        ct_co, cmap="gray", vmin=vmin, vmax=vmax,
        origin="lower", extent=ext_co, zorder=0
    )
    im_ct_sa = ax_sag.imshow(
        ct_sa, cmap="gray", vmin=vmin, vmax=vmax,
        origin="lower", extent=ext_sa, zorder=0
    )

    # -------------------------------------------------------------------------
    # HEATMAP DOSE
    # -------------------------------------------------------------------------
    im_dose_ax = im_dose_co = im_dose_sa = None

    if dose_vol is not None:
        if args.dose_mode == "tps":
            r_ax, r_co, r_sa = extract_slices(dose_rel, i, j, k)

            r_ax_m = mask_relative_dose(r_ax, args.dose_cutoff)
            r_co_m = mask_relative_dose(r_co, args.dose_cutoff)
            r_sa_m = mask_relative_dose(r_sa, args.dose_cutoff)

            im_dose_ax = ax_axial.imshow(
                r_ax_m, cmap=tps_cmap, norm=tps_norm,
                alpha=args.dose_alpha, origin="lower", extent=ext_ax, zorder=1
            )
            im_dose_co = ax_cor.imshow(
                r_co_m, cmap=tps_cmap, norm=tps_norm,
                alpha=args.dose_alpha, origin="lower", extent=ext_co, zorder=1
            )
            im_dose_sa = ax_sag.imshow(
                r_sa_m, cmap=tps_cmap, norm=tps_norm,
                alpha=args.dose_alpha, origin="lower", extent=ext_sa, zorder=1
            )

            # Colorbar relativa
            cbar_ax = fig.add_axes([0.92, 0.36, 0.015, 0.49])
            sm = ScalarMappable(norm=tps_norm, cmap=tps_cmap)
            sm.set_array([])
            cbar = fig.colorbar(sm, cax=cbar_ax, ticks=tps_levels)
            cbar.ax.set_yticklabels([f"{100*x:g}%" for x in tps_levels])
            cbar.set_label("Dose relativa a Dref")

        else:
            d_ax, d_co, d_sa = extract_slices(dose_vol.arr, i, j, k)

            # Nel modo continuo, cutoff espresso come frazione di Dref.
            threshold_native = args.dose_cutoff * dref
            d_ax_m = np.ma.masked_less(np.ma.masked_invalid(d_ax), threshold_native)
            d_co_m = np.ma.masked_less(np.ma.masked_invalid(d_co), threshold_native)
            d_sa_m = np.ma.masked_less(np.ma.masked_invalid(d_sa), threshold_native)

            im_dose_ax = ax_axial.imshow(
                d_ax_m, cmap=continuous_cmap, norm=continuous_norm,
                alpha=args.dose_alpha, origin="lower", extent=ext_ax, zorder=1
            )
            im_dose_co = ax_cor.imshow(
                d_co_m, cmap=continuous_cmap, norm=continuous_norm,
                alpha=args.dose_alpha, origin="lower", extent=ext_co, zorder=1
            )
            im_dose_sa = ax_sag.imshow(
                d_sa_m, cmap=continuous_cmap, norm=continuous_norm,
                alpha=args.dose_alpha, origin="lower", extent=ext_sa, zorder=1
            )

            cbar_ax = fig.add_axes([0.92, 0.36, 0.015, 0.49])
            cbar = fig.colorbar(im_dose_ax, cax=cbar_ax)
            cbar.set_label("Dose (unità native)")

    # -------------------------------------------------------------------------
    # CURVE DI ISODOSE
    # -------------------------------------------------------------------------
    iso_ax = iso_co = iso_sa = None

    def redraw_isodose(ii, jj, kk):
        nonlocal iso_ax, iso_co, iso_sa

        remove_contour_set(iso_ax)
        remove_contour_set(iso_co)
        remove_contour_set(iso_sa)
        iso_ax = iso_co = iso_sa = None

        if dose_rel is None or args.no_isodose_lines:
            return

        r_ax, r_co, r_sa = extract_slices(dose_rel, ii, jj, kk)

        iso_ax = draw_isodose_curves(
            ax_axial, r_ax, ext_ax, tps_levels, tps_cmap, tps_norm
        )
        iso_co = draw_isodose_curves(
            ax_cor, r_co, ext_co, tps_levels, tps_cmap, tps_norm
        )
        iso_sa = draw_isodose_curves(
            ax_sag, r_sa, ext_sa, tps_levels, tps_cmap, tps_norm
        )

    redraw_isodose(i, j, k)

    # -------------------------------------------------------------------------
    # ROI
    # -------------------------------------------------------------------------
    roi_colors = generate_colors(len(roi_arrays))

    roi_ax_lines: List[List[plt.Line2D]] = []
    roi_co_lines: List[List[plt.Line2D]] = []
    roi_sa_lines: List[List[plt.Line2D]] = []

    def clear_roi_lines(lines_nested):
        for lines in lines_nested:
            for ln in lines:
                try:
                    ln.remove()
                except Exception:
                    pass
        lines_nested.clear()

    def redraw_rois(ii, jj, kk):
        clear_roi_lines(roi_ax_lines)
        clear_roi_lines(roi_co_lines)
        clear_roi_lines(roi_sa_lines)

        for ridx, roi_arr in enumerate(roi_arrays):
            r_ax, r_co, r_sa = extract_slices(roi_arr, ii, jj, kk)

            c_ax = contours_phys(r_ax, "axial", ct.spacing)
            c_co = contours_phys(r_co, "coronal", ct.spacing)
            c_sa = contours_phys(r_sa, "sagittal", ct.spacing)

            roi_ax_lines.append(
                draw_roi_contours(ax_axial, c_ax, roi_colors[ridx])
            )
            roi_co_lines.append(
                draw_roi_contours(ax_cor, c_co, roi_colors[ridx])
            )
            roi_sa_lines.append(
                draw_roi_contours(ax_sag, c_sa, roi_colors[ridx])
            )

    redraw_rois(i, j, k)

    # -------------------------------------------------------------------------
    # CROSSHAIR
    # -------------------------------------------------------------------------
    sx, sy, sz = ct.spacing

    def redraw_crosshair(ii, jj, kk):
        x_phys = ii * sx
        y_phys = jj * sy
        z_phys = kk * sz

        draw_crosshair(ax_axial, x_phys, y_phys)
        draw_crosshair(ax_cor, x_phys, z_phys)
        draw_crosshair(ax_sag, y_phys, z_phys)

    redraw_crosshair(i, j, k)

    for ax in axes:
        ax.set_xlabel("mm")
        ax.set_ylabel("mm")
        ax.set_aspect("equal")

    # -------------------------------------------------------------------------
    # LAYOUT / SLIDER
    # -------------------------------------------------------------------------
    plt.subplots_adjust(
        left=0.07, right=0.90, bottom=0.36, top=0.85, wspace=0.15
    )

    axcolor = "lightgoldenrodyellow"
    ax_i = plt.axes([0.07, 0.27, 0.83, 0.03], facecolor=axcolor)
    ax_j = plt.axes([0.07, 0.23, 0.83, 0.03], facecolor=axcolor)
    ax_k = plt.axes([0.07, 0.19, 0.83, 0.03], facecolor=axcolor)
    ax_da = plt.axes([0.07, 0.15, 0.83, 0.03], facecolor=axcolor)

    s_i = Slider(ax_i, "X (i)", 0, nx - 1, valinit=i, valfmt="%0.0f")
    s_j = Slider(ax_j, "Y (j)", 0, ny - 1, valinit=j, valfmt="%0.0f")
    s_k = Slider(ax_k, "Z (k)", 0, nz - 1, valinit=k, valfmt="%0.0f")
    s_da = Slider(
        ax_da,
        "Dose α",
        0.0,
        1.0,
        valinit=(args.dose_alpha if dose_vol is not None else 0.0),
    )

    def update(_):
        nonlocal i, j, k

        i = int(round(s_i.val))
        j = int(round(s_j.val))
        k = int(round(s_k.val))

        # CT
        c_ax, c_co, c_sa = extract_slices(ct.arr, i, j, k)
        im_ct_ax.set_data(c_ax)
        im_ct_co.set_data(c_co)
        im_ct_sa.set_data(c_sa)

        # Dose
        if dose_vol is not None:
            alpha = float(s_da.val)

            if args.dose_mode == "tps":
                r_ax, r_co, r_sa = extract_slices(dose_rel, i, j, k)
                im_dose_ax.set_data(mask_relative_dose(r_ax, args.dose_cutoff))
                im_dose_co.set_data(mask_relative_dose(r_co, args.dose_cutoff))
                im_dose_sa.set_data(mask_relative_dose(r_sa, args.dose_cutoff))
            else:
                d_ax, d_co, d_sa = extract_slices(dose_vol.arr, i, j, k)
                threshold_native = args.dose_cutoff * dref
                im_dose_ax.set_data(
                    np.ma.masked_less(np.ma.masked_invalid(d_ax), threshold_native)
                )
                im_dose_co.set_data(
                    np.ma.masked_less(np.ma.masked_invalid(d_co), threshold_native)
                )
                im_dose_sa.set_data(
                    np.ma.masked_less(np.ma.masked_invalid(d_sa), threshold_native)
                )

            im_dose_ax.set_alpha(alpha)
            im_dose_co.set_alpha(alpha)
            im_dose_sa.set_alpha(alpha)

        # Prima isodose, poi ROI: le ROI restano graficamente sopra.
        redraw_isodose(i, j, k)
        redraw_rois(i, j, k)
        redraw_crosshair(i, j, k)

        fig.canvas.draw_idle()

    s_i.on_changed(update)
    s_j.on_changed(update)
    s_k.on_changed(update)
    s_da.on_changed(update)

    # -------------------------------------------------------------------------
    # ZOOM / PAN
    # -------------------------------------------------------------------------
    base_scale = 1.2

    def on_scroll(event):
        if event.inaxes not in {ax_axial, ax_cor, ax_sag}:
            return
        if event.xdata is None or event.ydata is None:
            return

        ax = event.inaxes
        xdata, ydata = event.xdata, event.ydata
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()

        scale = (1 / base_scale) if event.button == "up" else base_scale

        new_xlim = (
            xdata - (xdata - cur_xlim[0]) * scale,
            xdata + (cur_xlim[1] - xdata) * scale,
        )
        new_ylim = (
            ydata - (ydata - cur_ylim[0]) * scale,
            ydata + (cur_ylim[1] - ydata) * scale,
        )

        ax.set_xlim(new_xlim)
        ax.set_ylim(new_ylim)
        fig.canvas.draw_idle()

    pan_state = {
        "active": False,
        "ax": None,
        "x0": None,
        "y0": None,
        "xlim0": None,
        "ylim0": None,
    }

    def on_button_press(event):
        if event.button != 3:
            return
        if event.inaxes not in {ax_axial, ax_cor, ax_sag}:
            return
        if event.xdata is None or event.ydata is None:
            return

        pan_state["active"] = True
        pan_state["ax"] = event.inaxes
        pan_state["x0"] = event.xdata
        pan_state["y0"] = event.ydata
        pan_state["xlim0"] = event.inaxes.get_xlim()
        pan_state["ylim0"] = event.inaxes.get_ylim()

    def on_button_release(event):
        if event.button == 3:
            pan_state["active"] = False
            pan_state["ax"] = None

    def on_motion(event):
        if not pan_state["active"] or event.inaxes != pan_state["ax"]:
            return
        if event.xdata is None or event.ydata is None:
            return

        ax = pan_state["ax"]
        dx = event.xdata - pan_state["x0"]
        dy = event.ydata - pan_state["y0"]

        x0, x1 = pan_state["xlim0"]
        y0, y1 = pan_state["ylim0"]

        ax.set_xlim(x0 - dx, x1 - dx)
        ax.set_ylim(y0 - dy, y1 - dy)
        fig.canvas.draw_idle()

    def reset_zoom():
        ax_axial.set_xlim(ext_ax[0], ext_ax[1])
        ax_axial.set_ylim(ext_ax[2], ext_ax[3])

        ax_cor.set_xlim(ext_co[0], ext_co[1])
        ax_cor.set_ylim(ext_co[2], ext_co[3])

        ax_sag.set_xlim(ext_sa[0], ext_sa[1])
        ax_sag.set_ylim(ext_sa[2], ext_sa[3])

        fig.canvas.draw_idle()

    # -------------------------------------------------------------------------
    # TASTIERA
    # -------------------------------------------------------------------------
    step = 1

    def on_key(event):
        if event.key == "0":
            reset_zoom()
        elif event.key == "left":
            s_i.set_val(np.clip(s_i.val - step, s_i.valmin, s_i.valmax))
        elif event.key == "right":
            s_i.set_val(np.clip(s_i.val + step, s_i.valmin, s_i.valmax))
        elif event.key == "down":
            s_j.set_val(np.clip(s_j.val - step, s_j.valmin, s_j.valmax))
        elif event.key == "up":
            s_j.set_val(np.clip(s_j.val + step, s_j.valmin, s_j.valmax))
        elif event.key == "pageup":
            s_k.set_val(np.clip(s_k.val + step, s_k.valmin, s_k.valmax))
        elif event.key == "pagedown":
            s_k.set_val(np.clip(s_k.val - step, s_k.valmin, s_k.valmax))

    fig.canvas.mpl_connect("scroll_event", on_scroll)
    fig.canvas.mpl_connect("button_press_event", on_button_press)
    fig.canvas.mpl_connect("button_release_event", on_button_release)
    fig.canvas.mpl_connect("motion_notify_event", on_motion)
    fig.canvas.mpl_connect("key_press_event", on_key)

    # -------------------------------------------------------------------------
    # LEGENDA / TITOLI
    # -------------------------------------------------------------------------
    if roi_arrays:
        legend_handles = [
            plt.Line2D(
                [0], [0],
                color=roi_colors[idx],
                lw=2.5,
                label=labels[idx],
            )
            for idx in range(len(roi_arrays))
        ]
        ncol = min(len(legend_handles), 6)
        fig.legend(
            handles=legend_handles,
            loc="lower center",
            ncol=ncol,
            bbox_to_anchor=(0.5, 0.06),
        )

    titolo_pulito = args.title.replace("_", " ")
    fig.suptitle(titolo_pulito, fontsize=14, fontweight="bold", y=0.96)

    if fig.canvas.manager is not None:
        try:
            fig.canvas.manager.set_window_title(titolo_pulito)
        except Exception:
            pass

    plt.show()


if __name__ == "__main__":
    main()
