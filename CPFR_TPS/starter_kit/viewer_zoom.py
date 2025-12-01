#!/usr/bin/env python3
"""
Visualizzatore TAC + N ROI (file .mhd) + DOSE (mhd) con piani ortogonali, cursori, leggenda colori, controllo trasparenza dose

- CT in scala di grigi
- ROI visualizzate **solo come contorni** (niente riempimento)
- Dose (mhd) come heatmap sovrapposta con slider per la **trasparenza**
- Vista a 3 pannelli: coronale, trasversale (assiale), sagittale
- Cursori (slider) per X (sagittale), Y (coronale), Z (assiale)
- Slider aggiuntivo per **Dose Alpha**; finestratura dose opzionale
- Frecce: ←/→ (X), ↓/↑ (Y), PgDn/PgUp (Z)
- Resampling automatico (ROI: NN, Dose: lineare) verso la CT
- Leggenda dinamica per tutte le ROI + **colorbar** dose posizionata di fianco (a destra)
- **Zoom/Pan**: rotella per zoom, tasto destro per pan, tasto `0` per reset

Nota: la dose è visualizzata con normalizzazione **lineare** su [dmin, dmax].
"""

import argparse
from dataclasses import dataclass
from typing import List, Tuple, Optional

import math
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.colors import Normalize
from skimage import measure  # per i contorni delle ROI


@dataclass
class Volume:
    img: sitk.Image
    arr: np.ndarray  # shape (Z, Y, X)
    spacing: Tuple[float, float, float]  # (sx, sy, sz) in SITK order (X,Y,Z)
    size: Tuple[int, int, int]  # (X, Y, Z) in voxel


def read_mhd(path: str) -> Volume:
    img = sitk.ReadImage(path)
    arr = sitk.GetArrayFromImage(img)
    sx, sy, sz = img.GetSpacing()
    size = img.GetSize()
    return Volume(img=img, arr=arr, spacing=(sx, sy, sz), size=size)


def resample_like(moving_img: sitk.Image, reference_img: sitk.Image, is_label: bool = False) -> sitk.Image:
    interp = sitk.sitkNearestNeighbor if is_label else sitk.sitkLinear
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(reference_img)
    resampler.SetInterpolator(interp)
    resampler.SetDefaultPixelValue(0)
    resampler.SetTransform(sitk.Transform())
    return resampler.Execute(moving_img)


def window_auto(arr: np.ndarray, vmin: Optional[float] = None, vmax: Optional[float] = None,
                p_low: float = 1, p_high: float = 99) -> Tuple[float, float]:
    if vmin is None or vmax is None:
        arr_f = arr[np.isfinite(arr)]
        if arr_f.size == 0:
            return 0.0, 1.0
        lo, hi = np.percentile(arr_f, [p_low, p_high])
        return float(lo), float(hi)
    return float(vmin), float(vmax)


def build_extent(vol: Volume, plane: str) -> Tuple[float, float, float, float]:
    sx, sy, sz = vol.spacing
    nx, ny, nz = vol.size
    if plane == 'axial':
        return (0, nx * sx, 0, ny * sy)
    elif plane == 'coronal':
        return (0, nx * sx, 0, nz * sz)
    elif plane == 'sagittal':
        return (0, ny * sy, 0, nz * sz)
    else:
        raise ValueError('Unknown plane')


def extract_slices(arr: np.ndarray, i: int, j: int, k: int):
    axial = arr[k, :, :]
    coronal = arr[:, j, :]
    sagittal = arr[:, :, i]
    return axial, coronal, sagittal


def draw_crosshair(ax, x, y):
    for art in list(ax.lines):
        if getattr(art, "_is_cross", False):
            art.remove()
    h = ax.axhline(y, linewidth=0.8)
    v = ax.axvline(x, linewidth=0.8)
    setattr(h, "_is_cross", True)
    setattr(v, "_is_cross", True)


def generate_colors(n: int) -> List[Tuple[float, float, float]]:
    if n <= 0:
        return []
    hues = np.linspace(0, 1, n, endpoint=False)
    colors = []
    for h in hues:
        i = int(h * 6)
        f = h * 6 - i
        p = 0.1
        q = 0.9 * (1 - f * 0.9)
        t = 0.9 * (1 - (1 - f) * 0.9)
        i = i % 6
        if i == 0:
            r, g, b = 0.9, t, p
        elif i == 1:
            r, g, b = q, 0.9, p
        elif i == 2:
            r, g, b = p, 0.9, t
        elif i == 3:
            r, g, b = p, q, 0.9
        elif i == 4:
            r, g, b = t, p, 0.9
        else:
            r, g, b = 0.9, p, q
        colors.append((r, g, b))
    return colors


def contours_phys(mask2d: np.ndarray, plane: str, spacing: Tuple[float, float, float]):
    sx, sy, sz = spacing
    cs = measure.find_contours(mask2d.astype(float), level=0.5)
    out = []
    for c in cs:
        rows = c[:, 0]
        cols = c[:, 1]
        if plane == 'axial':
            x_mm = cols * sx
            y_mm = rows * sy
        elif plane == 'coronal':
            x_mm = cols * sx
            y_mm = rows * sz
        elif plane == 'sagittal':
            x_mm = cols * sy
            y_mm = rows * sz
        else:
            raise ValueError('Unknown plane')
        out.append(np.vstack([x_mm, y_mm]).T)
    return out


def draw_contours(ax, contours_xy, color, lw=1.8):
    lines = []
    for xy in contours_xy:
        line, = ax.plot(xy[:, 0], xy[:, 1], color=color, linewidth=lw)
        lines.append(line)
    return lines


# === Tacche "belle" per colorbar ===
def nice_ticks(vmin: float, vmax: float, max_ticks: int = 12):
    if not np.isfinite([vmin, vmax]).all():
        return []
    if vmax <= vmin:
        vmax = vmin + 1.0
    span = vmax - vmin
    raw_step = span / max(1, max_ticks)
    pow10 = 10 ** math.floor(math.log10(raw_step))
    for m in [1, 2, 5, 10]:
        step = m * pow10
        if span / step <= max_ticks:
            break
    start = math.floor(vmin / step) * step
    stop = math.ceil(vmax / step) * step
    ticks = np.arange(start, stop + 0.5 * step, step)
    ticks = ticks[(ticks >= vmin - 1e-9) & (ticks <= vmax + 1e-9)]
    # preferisci multipli da 50 se l'intervallo è ampio
    if span >= 500:
        step50 = max(50, int(round(step / 50.0)) * 50)
        start50 = math.ceil(vmin / step50) * step50
        ticks = np.arange(start50, vmax + 1e-9, step50)
    return ticks.tolist()


def main():
    parser = argparse.ArgumentParser(description='Viewer TAC+ROI+DOSE per file .mhd con piani ortogonali e cursori')
    parser.add_argument('--ct', required=True, help='Percorso file CT .mhd')
    parser.add_argument('--roi', required=True, nargs='+', help='Percorsi per N ROI .mhd (etichettate)')
    parser.add_argument('--labels', nargs='*', default=None, help='Etichette per le ROI (opzionale, stessa lunghezza di --roi)')
    parser.add_argument('--window', type=float, nargs=2, metavar=('VMIN', 'VMAX'), default=None, help='Finestra intensità CT per visualizzazione')
    parser.add_argument('--dose', type=str, default=None, help='Percorso file dose .mhd (opzionale)')
    parser.add_argument('--dose-window', type=float, nargs=2, metavar=('DMIN', 'DMAX'), default=None, help='Finestra dose per visualizzazione (stessa unità della dose, es. Gy)')
    parser.add_argument('--dose-cmap', type=str, default='inferno', help='Colormap Matplotlib per la dose (default: inferno)')
    parser.add_argument('--dose-alpha', type=float, default=0.5, help='Trasparenza iniziale dose [0-1] (default: 0.5)')

    args = parser.parse_args()

    if args.labels is not None and len(args.labels) != len(args.roi):
        raise SystemExit(f"Il numero di etichette ({len(args.labels)}) deve coincidere con il numero di ROI ({len(args.roi)}).")

    labels = args.labels if args.labels is not None else [f'ROI {i+1}' for i in range(len(args.roi))]

    # --- CT ---
    ct = read_mhd(args.ct)

    # --- ROI (resample NN verso CT) ---
    roi_imgs = []
    for p in args.roi:
        roi_img = sitk.ReadImage(p)
        need_resample = (
            roi_img.GetSize() != ct.img.GetSize() or
            roi_img.GetSpacing() != ct.img.GetSpacing() or
            roi_img.GetOrigin() != ct.img.GetOrigin() or
            roi_img.GetDirection() != ct.img.GetDirection()
        )
        if need_resample:
            roi_img = resample_like(roi_img, ct.img, is_label=True)
        roi_imgs.append(roi_img)
    roi_arrays = [sitk.GetArrayFromImage(r) for r in roi_imgs]

    # --- Dose (opzionale; resample lineare verso CT) ---
    dose_vol: Optional[Volume] = None
    if args.dose is not None:
        d_img = sitk.ReadImage(args.dose)
        need_resample_d = (
            d_img.GetSize() != ct.img.GetSize() or
            d_img.GetSpacing() != ct.img.GetSpacing() or
            d_img.GetOrigin() != ct.img.GetOrigin() or
            d_img.GetDirection() != ct.img.GetDirection()
        )
        if need_resample_d:
            d_img = resample_like(d_img, ct.img, is_label=False)
        d_arr = sitk.GetArrayFromImage(d_img)
        dose_vol = Volume(img=d_img, arr=d_arr, spacing=ct.spacing, size=ct.size)

    # --- Windowing (CT: percentili; Dose: LINEARE su [dmin,dmax]) ---
    vmin, vmax = window_auto(ct.arr, *(args.window if args.window else (None, None)))

    dose_norm = None
    dmin = dmax = None
    if dose_vol is not None:
        if args.dose_window:
            dmin, dmax = float(args.dose_window[0]), float(args.dose_window[1])
        else:
            arrf = dose_vol.arr[np.isfinite(dose_vol.arr)]
            if arrf.size == 0:
                dmin, dmax = 0.0, 1.0
            else:
                # Fondo a 0 se la dose è non-negativa; altrimenti usa il minimo reale
                dmin = 0.0 if arrf.min() >= 0 else float(arrf.min())
                dmax = float(arrf.max())
                if dmax <= dmin:
                    dmax = dmin + 1.0
        dose_norm = Normalize(vmin=dmin, vmax=dmax)

    # --- Figure ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), gridspec_kw={'width_ratios': [1, 1, 1.05]})
    ax_axial, ax_cor, ax_sag = axes
    ax_axial.set_title('Trasversale (assiale) — Z')
    ax_cor.set_title('Coronale — Y')
    ax_sag.set_title('Sagittale — X')

    ext_ax = build_extent(ct, 'axial')
    ext_co = build_extent(ct, 'coronal')
    ext_sa = build_extent(ct, 'sagittal')

    nx, ny, nz = ct.size
    i = nx // 2
    j = ny // 2
    k = nz // 2

    colors = generate_colors(len(roi_arrays))

    # --- Slice iniziali ---
    ct_ax, ct_co, ct_sa = extract_slices(ct.arr, i, j, k)
    im_ct_ax = ax_axial.imshow(ct_ax, cmap='gray', vmin=vmin, vmax=vmax, origin='lower', extent=ext_ax)
    im_ct_co = ax_cor.imshow(ct_co, cmap='gray', vmin=vmin, vmax=vmax, origin='lower', extent=ext_co)
    im_ct_sa = ax_sag.imshow(ct_sa, cmap='gray', vmin=vmin, vmax=vmax, origin='lower', extent=ext_sa)

    # Dose overlays (se presente) + colorbar con tacche “belle”
    im_dose_ax = im_dose_co = im_dose_sa = None
    cbar = None
    if dose_vol is not None:
        d_ax, d_co, d_sa = extract_slices(dose_vol.arr, i, j, k)
        im_dose_ax = ax_axial.imshow(d_ax, cmap=args.dose_cmap, norm=dose_norm,
                                     alpha=args.dose_alpha, origin='lower', extent=ext_ax)
        im_dose_co = ax_cor.imshow(d_co, cmap=args.dose_cmap, norm=dose_norm,
                                   alpha=args.dose_alpha, origin='lower', extent=ext_co)
        im_dose_sa = ax_sag.imshow(d_sa, cmap=args.dose_cmap, norm=dose_norm,
                                   alpha=args.dose_alpha, origin='lower', extent=ext_sa)
        # Colorbar laterale destra
        cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.6])
        cbar = fig.colorbar(im_dose_ax, cax=cbar_ax)
        cbar.set_label('Dose (unità native)')
        ticks = nice_ticks(dose_norm.vmin, dose_norm.vmax, max_ticks=12)
        if ticks:
            cbar.set_ticks(ticks)

    # ROI contours (liste di Line2D per poterle rimuovere/aggiornare)
    roi_ax_lines: List[List[plt.Line2D]] = []
    roi_co_lines: List[List[plt.Line2D]] = []
    roi_sa_lines: List[List[plt.Line2D]] = []

    for ridx, roi_arr in enumerate(roi_arrays):
        r_ax, r_co, r_sa = extract_slices(roi_arr, i, j, k)
        c_ax = contours_phys(r_ax, 'axial', ct.spacing)
        c_co = contours_phys(r_co, 'coronal', ct.spacing)
        c_sa = contours_phys(r_sa, 'sagittal', ct.spacing)
        roi_ax_lines.append(draw_contours(ax_axial, c_ax, colors[ridx]))
        roi_co_lines.append(draw_contours(ax_cor, c_co, colors[ridx]))
        roi_sa_lines.append(draw_contours(ax_sag, c_sa, colors[ridx]))

    # Crosshair
    sx, sy, sz = ct.spacing
    x_phys = i * sx
    y_phys = j * sy
    z_phys = k * sz
    draw_crosshair(ax_axial, x_phys, y_phys)
    draw_crosshair(ax_cor, x_phys, z_phys)
    draw_crosshair(ax_sag, y_phys, z_phys)

    for ax in axes:
        ax.set_xlabel('mm')
        ax.set_ylabel('mm')
        ax.set_aspect('equal')

    # Layout per slider e legenda
    plt.subplots_adjust(left=0.07, right=0.9, bottom=0.36, top=0.92, wspace=0.15)

    axcolor = 'lightgoldenrodyellow'
    ax_i = plt.axes([0.07, 0.27, 0.86, 0.03], facecolor=axcolor)
    ax_j = plt.axes([0.07, 0.23, 0.86, 0.03], facecolor=axcolor)
    ax_k = plt.axes([0.07, 0.19, 0.86, 0.03], facecolor=axcolor)
    ax_da = plt.axes([0.07, 0.15, 0.86, 0.03], facecolor=axcolor)

    s_i = Slider(ax_i, 'X (i)', 0, nx - 1, valinit=i, valfmt='%0.0f')
    s_j = Slider(ax_j, 'Y (j)', 0, ny - 1, valinit=j, valfmt='%0.0f')
    s_k = Slider(ax_k, 'Z (k)', 0, nz - 1, valinit=k, valfmt='%0.0f')
    s_da = Slider(ax_da, 'Dose α', 0.0, 1.0, valinit=(args.dose_alpha if dose_vol is not None else 0.0))

    def clear_lines(lines_nested):
        for lines in lines_nested:
            for ln in lines:
                ln.remove()
        lines_nested.clear()

    def update(_):
        nonlocal i, j, k
        i = int(round(s_i.val))
        j = int(round(s_j.val))
        k = int(round(s_k.val))

        ct_ax, ct_co, ct_sa = extract_slices(ct.arr, i, j, k)
        im_ct_ax.set_data(ct_ax)
        im_ct_co.set_data(ct_co)
        im_ct_sa.set_data(ct_sa)

        if dose_vol is not None:
            d_ax, d_co, d_sa = extract_slices(dose_vol.arr, i, j, k)
            im_dose_ax.set_data(d_ax)
            im_dose_co.set_data(d_co)
            im_dose_sa.set_data(d_sa)
            a = float(s_da.val)
            im_dose_ax.set_alpha(a)
            im_dose_co.set_alpha(a)
            im_dose_sa.set_alpha(a)
            # Norm e ticks restano invariati: mappatura lineare fissa su [dmin, dmax]

        # Aggiorna contorni ROI
        clear_lines(roi_ax_lines)
        clear_lines(roi_co_lines)
        clear_lines(roi_sa_lines)
        for ridx, roi_arr in enumerate(roi_arrays):
            r_ax, r_co, r_sa = extract_slices(roi_arr, i, j, k)
            c_ax = contours_phys(r_ax, 'axial', ct.spacing)
            c_co = contours_phys(r_co, 'coronal', ct.spacing)
            c_sa = contours_phys(r_sa, 'sagittal', ct.spacing)
            roi_ax_lines.append(draw_contours(ax_axial, c_ax, colors[ridx]))
            roi_co_lines.append(draw_contours(ax_cor, c_co, colors[ridx]))
            roi_sa_lines.append(draw_contours(ax_sag, c_sa, colors[ridx]))

        # Crosshair aggiornate
        x_phys = i * sx
        y_phys = j * sy
        z_phys = k * sz
        draw_crosshair(ax_axial, x_phys, y_phys)
        draw_crosshair(ax_cor, x_phys, z_phys)
        draw_crosshair(ax_sag, y_phys, z_phys)

        fig.canvas.draw_idle()

    # Collega gli slider all'update
    s_i.on_changed(update)
    s_j.on_changed(update)
    s_k.on_changed(update)
    s_da.on_changed(update)

    # ===== Zoom & Pan interattivi =====
    base_scale = 1.2  # fattore di zoom per step rotella

    def on_scroll(event):
        if event.inaxes not in {ax_axial, ax_cor, ax_sag}:
            return
        ax = event.inaxes
        if event.xdata is None or event.ydata is None:
            return
        xdata, ydata = event.xdata, event.ydata
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        scale = (1 / base_scale) if event.button == 'up' else base_scale
        new_xlim = (xdata - (xdata - cur_xlim[0]) * scale,
                    xdata + (cur_xlim[1] - xdata) * scale)
        new_ylim = (ydata - (ydata - cur_ylim[0]) * scale,
                    ydata + (cur_ylim[1] - ydata) * scale)
        ax.set_xlim(new_xlim)
        ax.set_ylim(new_ylim)
        fig.canvas.draw_idle()

    pan_state = {'active': False, 'ax': None, 'x0': None, 'y0': None, 'xlim0': None, 'ylim0': None}

    def on_button_press(event):
        if event.button != 3:  # tasto destro = pan
            return
        if event.inaxes not in {ax_axial, ax_cor, ax_sag}:
            return
        pan_state['active'] = True
        pan_state['ax'] = event.inaxes
        pan_state['x0'] = event.xdata
        pan_state['y0'] = event.ydata
        pan_state['xlim0'] = event.inaxes.get_xlim()
        pan_state['ylim0'] = event.inaxes.get_ylim()

    def on_button_release(event):
        if event.button != 3:
            return
        pan_state['active'] = False
        pan_state['ax'] = None

    def on_motion(event):
        if not pan_state['active'] or event.inaxes != pan_state['ax']:
            return
        if event.xdata is None or event.ydata is None:
            return
        ax = pan_state['ax']
        dx = event.xdata - pan_state['x0']
        dy = event.ydata - pan_state['y0']
        x0, x1 = pan_state['xlim0']
        y0, y1 = pan_state['ylim0']
        ax.set_xlim(x0 - dx, x1 - dx)
        ax.set_ylim(y0 - dy, y1 - dy)
        fig.canvas.draw_idle()

    def on_key_extra(event):
        if event.key == '0':
            # reset limiti ai extent fisici originali
            ax_axial.set_xlim(ext_ax[0], ext_ax[1])
            ax_axial.set_ylim(ext_ax[2], ext_ax[3])
            ax_cor.set_xlim(ext_co[0], ext_co[1])
            ax_cor.set_ylim(ext_co[2], ext_co[3])
            ax_sag.set_xlim(ext_sa[0], ext_sa[1])
            ax_sag.set_ylim(ext_sa[2], ext_sa[3])
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('scroll_event', on_scroll)
    fig.canvas.mpl_connect('button_press_event', on_button_press)
    fig.canvas.mpl_connect('button_release_event', on_button_release)
    fig.canvas.mpl_connect('motion_notify_event', on_motion)
    fig.canvas.mpl_connect('key_press_event', on_key_extra)

    # Scorciatoie tastiera per slicing
    step = 1
    def on_key(event):
        if event.key == 'left':
            s_i.set_val(np.clip(s_i.val - step, s_i.valmin, s_i.valmax))
        elif event.key == 'right':
            s_i.set_val(np.clip(s_i.val + step, s_i.valmin, s_i.valmax))
        elif event.key == 'down':
            s_j.set_val(np.clip(s_j.val - step, s_j.valmin, s_j.valmax))
        elif event.key == 'up':
            s_j.set_val(np.clip(s_j.val + step, s_j.valmin, s_j.valmax))
        elif event.key == 'pageup':
            s_k.set_val(np.clip(s_k.val + step, s_k.valmin, s_k.valmax))
        elif event.key == 'pagedown':
            s_k.set_val(np.clip(s_k.val - step, s_k.valmin, s_k.valmax))

    fig.canvas.mpl_connect('key_press_event', on_key)

    # Leggenda ROI
    legend_handles = [plt.Line2D([0], [0], color=colors[i], lw=2.5, label=labels[i]) for i in range(len(roi_arrays))]
    ncol = min(len(legend_handles), 6)
    fig.legend(handles=legend_handles, loc='lower center', ncol=ncol, bbox_to_anchor=(0.5, 0.06))

    plt.show()


if __name__ == '__main__':
    main()
