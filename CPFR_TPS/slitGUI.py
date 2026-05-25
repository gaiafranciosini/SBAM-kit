#!/usr/bin/env python3
"""
BEV Collimator GUI
==================
Beam Eye View visualizer for MHD volumes with a 4-leaf/slit rectangular collimator.
Supports multiple overlaid MHD volumes with individual solid-color tints and labels.

Features
--------
- Independent slit motion: top, bottom, left, right are NOT coupled.
- Direct numeric input for each slit distance from isocentre.
- Distances are expressed in cm from the centre.
- Overlay volumes are shown as solid-color masks where slice != 0.
- Layout optimized for smaller windows / SSH X-forwarding.

Usage examples
--------------
# Single volume
python slitGUI.py --mhd ct.mhd

# Multiple volumes
python slitGUI.py \
    --mhd imgs/CT_plan.mhd imgs/PTV_plan.mhd imgs/Muscolo_plan.mhd \
    --label CT PTV MUSCOLO \
    --color gray red cyan \
    --alpha 0.7 0.5 0.5 \
    --angle 60

# Set asymmetric slit coordinates
python slitGUI.py \
    --mhd ct.mhd ptv.mhd \
    --top 1.5 --bottom 2.0 --left 3.0 --right 4.0
"""

import argparse
import sys
import math
import json
from pathlib import Path

import numpy as np
import SimpleITK as sitk
import matplotlib
matplotlib.use("TkAgg")   # change to "Qt5Agg" if Tk not available
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, Slider, TextBox
from matplotlib.colors import LinearSegmentedColormap, to_rgba, is_color_like
import matplotlib.transforms as transforms
from matplotlib.lines import Line2D


# ──────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────

_AUTO_COLORS = [
    None,        # 0: primary volume -> grayscale
    '#FF4040',   # red
    '#40FF80',   # green
    '#4080FF',   # blue
    '#FFD700',   # gold
    '#FF40FF',   # magenta
    '#00FFFF',   # cyan
    '#FF8C00',   # orange
    '#BF40FF',   # violet
    '#80FF40',   # lime
]

LEAF_COLORS = {
    'top':    '#00CFFF',
    'bottom': '#00CFFF',
    'left':   '#FF9F00',
    'right':  '#FF9F00',
}


# ──────────────────────────────────────────────────────────────────────────────
# GEOMETRY HELPERS
# ──────────────────────────────────────────────────────────────────────────────

def find_iso_pixel(volume_shape, origin, spacing):
    """
    Compute isocentre position in pixel/slice coordinates.

    SimpleITK layout
    ----------------
    array shape : (NZ, NY, NX)
    spacing     : (sx, sy, sz) [mm]
    origin      : (ox, oy, oz) [mm]

    The isocentre is approximated as the voxel minimizing |physical_coordinate|
    independently along each axis.
    """
    nz, ny, nx = volume_shape
    sx, sy, sz = spacing
    ox, oy, oz = origin

    xs = np.arange(nx) * sx + ox
    ys = np.arange(ny) * sy + oy
    zs = np.arange(nz) * sz + oz

    return (
        int(np.argmin(np.abs(xs))),
        int(np.argmin(np.abs(ys))),
        int(np.argmin(np.abs(zs))),
    )


def voxel_z_to_cm(iz, origin_z, spacing_z):
    return (iz * spacing_z + origin_z) / 10.0


def rotation_matrix_2d(angle_deg):
    a = math.radians(angle_deg)
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s], [s, c]])


def _make_cmap(color):
    if color is None:
        return 'gray'
    return LinearSegmentedColormap.from_list(f'vol_{color}', ['#000000', color], N=256)


def _solid_rgba_overlay(slice2d, color, alpha, mask_nonzero=True):
    rgba = np.zeros(slice2d.shape + (4,), dtype=float)

    if color is None:
        color = '#ffffff'

    r, g, b, _ = to_rgba(color)
    mask = (slice2d != 0) if mask_nonzero else np.ones_like(slice2d, dtype=bool)

    rgba[..., 0] = r
    rgba[..., 1] = g
    rgba[..., 2] = b
    rgba[..., 3] = alpha * mask.astype(float)

    return rgba


def _legend_color(vol):
    return '#aaaaaa' if vol.color is None else vol.color


def _resolve_color(raw):
    if raw is None:
        return None

    raw = raw.strip()
    if raw.lower() in ('none', 'gray', 'grey', ''):
        return None

    if not is_color_like(raw):
        raise ValueError(f"Invalid color: {raw}")

    return raw


def _validate_list_lengths(n_volumes, labels, colors, alphas):
    if len(labels) > n_volumes:
        raise ValueError(f"Received {len(labels)} labels but only {n_volumes} MHD volumes.")
    if len(colors) > n_volumes:
        raise ValueError(f"Received {len(colors)} colors but only {n_volumes} MHD volumes.")
    if len(alphas) > n_volumes:
        raise ValueError(f"Received {len(alphas)} alpha values but only {n_volumes} MHD volumes.")


def _validate_alphas(alphas):
    for a in alphas:
        if a is not None and not (0.0 <= a <= 1.0):
            raise ValueError(f"Alpha must be in [0, 1], got {a}")


# ──────────────────────────────────────────────────────────────────────────────
# VOLUME DESCRIPTOR
# ──────────────────────────────────────────────────────────────────────────────

class VolumeDesc:
    def __init__(self, path, label=None, color=None, alpha=1.0):
        self.path = Path(path)
        self.label = label if label else self.path.stem
        self.color = color
        self.alpha = alpha
        self.cmap = _make_cmap(color)

        img = sitk.ReadImage(str(self.path))
        self.array = sitk.GetArrayFromImage(img)   # (NZ, NY, NX)
        self.spacing = img.GetSpacing()            # (sx, sy, sz) mm
        self.origin = img.GetOrigin()              # (ox, oy, oz) mm

    def get_slice(self, iz):
        iz = max(0, min(iz, self.array.shape[0] - 1))
        return self.array[iz, :, :]

    @property
    def nz(self):
        return self.array.shape[0]


# ──────────────────────────────────────────────────────────────────────────────
# LEAF / SLIT
# ──────────────────────────────────────────────────────────────────────────────

class Leaf:
    """
    One rectangular jaw/slit of the collimator.

    Distances are independent and measured from isocentre to the INNER edge.

    top    : inner edge at  y = +d
    bottom : inner edge at  y = -d
    left   : inner edge at  x = -d
    right  : inner edge at  x = +d
    """

    def __init__(self, side, inner_dist_cm, width_cm, height_cm, angle_deg, sp_x, sp_y):
        self.side = side
        self.inner_dist_cm = inner_dist_cm
        self.width_cm = width_cm
        self.height_cm = height_cm
        self.angle_deg = angle_deg
        self.sp_x = sp_x
        self.sp_y = sp_y
        self.patch = None

    def _rect_cm(self):
        """
        Return (x0, y0, w, h) in cm before rotation, with isocentre at (0, 0).
        """
        w = self.height_cm
        h = self.width_cm
        d = self.inner_dist_cm

        if self.side == 'top':
            return -w / 2.0, d, w, h
        if self.side == 'bottom':
            return -w / 2.0, -(d + h), w, h
        if self.side == 'left':
            return -(d + h), -w / 2.0, h, w
        return d, -w / 2.0, h, w  # right

    def make_patch(self, iso_px, ax, color='cyan', alpha=0.40):
        x0c, y0c, wc, hc = self._rect_cm()

        x0 = iso_px[0] + x0c / self.sp_x
        y0 = iso_px[1] + y0c / self.sp_y
        w = wc / self.sp_x
        h = hc / self.sp_y

        if self.patch is None:
            self.patch = Rectangle(
                (x0, y0), w, h,
                linewidth=2,
                edgecolor=color,
                facecolor=color,
                alpha=alpha,
                picker=True,
                zorder=10
            )
            ax.add_patch(self.patch)
        else:
            self.patch.set_xy((x0, y0))
            self.patch.set_width(w)
            self.patch.set_height(h)
            self.patch.set_edgecolor(color)
            self.patch.set_facecolor(color)
            self.patch.set_alpha(alpha)

        t = (
            transforms.Affine2D()
            .rotate_deg_around(iso_px[0], iso_px[1], self.angle_deg)
            + ax.transData
        )
        self.patch.set_transform(t)

    def contains_point(self, event):
        return self.patch is not None and self.patch.contains(event)[0]


# ──────────────────────────────────────────────────────────────────────────────
# APPLICATION
# ──────────────────────────────────────────────────────────────────────────────

class BEVApp:

    def __init__(self, volumes, args):
        self.volumes = volumes
        self.primary = volumes[0]

        self.iso_x, self.iso_y, self.iso_z = find_iso_pixel(
            self.primary.array.shape,
            self.primary.origin,
            self.primary.spacing
        )
        self.current_z = self.iso_z

        self.sp_x = self.primary.spacing[0] / 10.0
        self.sp_y = self.primary.spacing[1] / 10.0

        self.rect_w = args.rect_w
        self.rect_h = args.rect_h

        self.top_cm = args.top
        self.bottom_cm = args.bottom
        self.left_cm = args.left
        self.right_cm = args.right

        self.angle = args.angle
        self.output_file = args.output

        self._drag_leaf = None
        self._drag_start = None
        self._drag_coord_start = None

        self._textbox_guard = False
        self.coord_boxes = {}

        self._build_leaves()
        self._build_gui()

    # ── geometry ──────────────────────────────────────────────────────

    def _iso_px(self):
        return (self.iso_x, self.iso_y)

    def _build_leaves(self):
        kw = dict(
            width_cm=self.rect_w,
            height_cm=self.rect_h,
            angle_deg=self.angle,
            sp_x=self.sp_x,
            sp_y=self.sp_y,
        )
        self.leaves = {
            'top':    Leaf('top',    self.top_cm, **kw),
            'bottom': Leaf('bottom', self.bottom_cm, **kw),
            'left':   Leaf('left',   self.left_cm, **kw),
            'right':  Leaf('right',  self.right_cm, **kw),
        }

    def _sync_leaf_values_from_state(self):
        self.leaves['top'].inner_dist_cm = self.top_cm
        self.leaves['bottom'].inner_dist_cm = self.bottom_cm
        self.leaves['left'].inner_dist_cm = self.left_cm
        self.leaves['right'].inner_dist_cm = self.right_cm

    def _field_opening(self):
        return self.top_cm + self.bottom_cm, self.left_cm + self.right_cm

    # ── GUI construction ───────────────────────────────────────────────

    def _build_gui(self):
        self.fig = plt.figure(figsize=(14, 8), facecolor='#0d1117')
        self.fig.canvas.manager.set_window_title("BEV Collimator Viewer")

        # Main image axes
        self.ax = self.fig.add_axes([0.06, 0.14, 0.58, 0.80])
        self.ax.set_facecolor('black')
        self.ax.tick_params(colors='#aaaaaa')
        for sp in self.ax.spines.values():
            sp.set_edgecolor('#333333')

        # Right-top info panel
        self.ax_info = self.fig.add_axes([0.67, 0.50, 0.30, 0.42])
        self.ax_info.axis('off')
        self.ax_info.set_facecolor('#0d1117')

        # Right-bottom controls panel
        self.ax_controls = self.fig.add_axes([0.67, 0.18, 0.30, 0.24])
        self.ax_controls.axis('off')
        self.ax_controls.set_facecolor('#0d1117')

        # Render volumes
        self.im_handles = []
        for i, vol in enumerate(self.volumes):
            sl = vol.get_slice(self.current_z)

            if i == 0:
                im = self.ax.imshow(
                    sl,
                    cmap=vol.cmap,
                    origin='lower',
                    aspect='equal',
                    alpha=vol.alpha,
                    zorder=i + 1
                )
            else:
                overlay_rgba = _solid_rgba_overlay(sl, vol.color, vol.alpha)
                im = self.ax.imshow(
                    overlay_rgba,
                    origin='lower',
                    aspect='equal',
                    zorder=i + 1
                )

            self.im_handles.append(im)

        self.ax.set_title("Beam Eye View", color='#00CFFF', fontsize=13, fontweight='bold', pad=8)
        self.ax.set_xlabel("X (pixel)", color='#aaaaaa')
        self.ax.set_ylabel("Y (pixel)", color='#aaaaaa')

        legend_handles = [
            Line2D([0], [0], color=_legend_color(vol), linewidth=5, label=vol.label)
            for vol in self.volumes
        ]
        self.ax.legend(
            handles=legend_handles,
            loc='upper left',
            facecolor='#1c2333',
            edgecolor='#444',
            labelcolor='#dddddd',
            fontsize=9
        )

        # Crosshair + isocentre
        iso = self._iso_px()
        self.ax.axhline(iso[1], color='#ff4444', lw=0.8, ls='--', alpha=0.9, zorder=20)
        self.ax.axvline(iso[0], color='#ff4444', lw=0.8, ls='--', alpha=0.9, zorder=20)
        self.ax.plot(iso[0], iso[1], '+', color='#ff4444', ms=14, mew=2, zorder=21)

        # Leaves
        for side, leaf in self.leaves.items():
            leaf.make_patch(iso, self.ax, color=LEAF_COLORS[side], alpha=0.40)

        # Navigation buttons
        ax_prev = self.fig.add_axes([0.12, 0.05, 0.08, 0.05])
        ax_next = self.fig.add_axes([0.21, 0.05, 0.08, 0.05])
        self.btn_prev = Button(ax_prev, '◀ Prev', color='#1c2333', hovercolor='#2a3a55')
        self.btn_next = Button(ax_next, 'Next ▶', color='#1c2333', hovercolor='#2a3a55')
        self.btn_prev.label.set_color('#aaaaaa')
        self.btn_next.label.set_color('#aaaaaa')
        self.btn_prev.on_clicked(self._on_prev)
        self.btn_next.on_clicked(self._on_next)

        # Slice label
        ax_sl = self.fig.add_axes([0.31, 0.045, 0.18, 0.06])
        ax_sl.axis('off')
        ax_sl.set_facecolor('#0d1117')
        self._slice_txt = ax_sl.text(
            0.5, 0.5, self._slice_label(),
            ha='center', va='center',
            color='#00CFFF', fontsize=9,
            transform=ax_sl.transAxes,
            fontfamily='monospace'
        )

        # Horizontal angle slider
        ax_ang = self.fig.add_axes([0.62, 0.08, 0.18, 0.03])
        ax_ang.set_facecolor('#0d1117')
        self.angle_slider = Slider(
            ax_ang, 'Angle', -180, 180,
            valinit=self.angle,
            color='#00CFFF',
            initcolor='#ff4444'
        )
        self.angle_slider.label.set_color('#aaaaaa')
        self.angle_slider.valtext.set_color('#00CFFF')
        self.angle_slider.valtext.set_fontsize(9)
        self.angle_slider.on_changed(self._on_angle)

        # Save button
        ax_save = self.fig.add_axes([0.82, 0.05, 0.12, 0.05])
        self.btn_save = Button(ax_save, '💾 Save', color='#1c4433', hovercolor='#2a6a4a')
        self.btn_save.label.set_color('#aaffaa')
        self.btn_save.on_clicked(self._on_save)

        # Textboxes
        self._build_textboxes()

        # Events
        self.fig.canvas.mpl_connect('button_press_event', self._on_press)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.fig.canvas.mpl_connect('button_release_event', self._on_release)
        self.fig.canvas.mpl_connect('key_press_event', self._on_key)

        self._update_textboxes()
        self._update_info()
        plt.show()

    def _build_textboxes(self):
        self.coord_boxes = {}

        self.ax_controls.text(
            0.02, 0.95, "SLIT POSITIONS [cm]",
            color='#dddddd', fontsize=10, fontweight='bold',
            va='top', ha='left', fontfamily='monospace',
            transform=self.ax_controls.transAxes
        )

        specs = [
            ('top',    0.02, 0.55, '#00CFFF'),
            ('bottom', 0.52, 0.55, '#00CFFF'),
            ('left',   0.02, 0.12, '#FF9F00'),
            ('right',  0.52, 0.12, '#FF9F00'),
        ]

        panel_left = 0.67
        panel_bottom = 0.18
        panel_w = 0.30
        panel_h = 0.24

        for name, x, y, color in specs:
            self.ax_controls.text(
                x, y + 0.20, f"{name.upper()}",
                color=color, fontsize=9, fontweight='bold',
                va='bottom', ha='left', fontfamily='monospace',
                transform=self.ax_controls.transAxes
            )

            ax_box = self.fig.add_axes([
                panel_left + panel_w * x,
                panel_bottom + panel_h * y,
                panel_w * 0.40,
                panel_h * 0.18
            ])

            tb = TextBox(ax_box, '', initial="0.000", color='#111827', hovercolor='#1f2937')
            ax_box.set_facecolor('#111827')
            for sp in ax_box.spines.values():
                sp.set_edgecolor('#444444')
            tb.text_disp.set_color('#dddddd')
            tb.text_disp.set_fontfamily('monospace')

            tb.on_submit(lambda text, side=name: self._on_submit_coord(side, text))
            self.coord_boxes[name] = tb

        self.ax_controls.text(
            0.02, 0.02, "Edit value and press Enter",
            color='#666666', fontsize=8,
            va='bottom', ha='left', fontfamily='monospace',
            transform=self.ax_controls.transAxes
        )

    # ── slice helpers ──────────────────────────────────────────────────

    def _slice_label(self):
        z_cm = voxel_z_to_cm(self.current_z, self.primary.origin[2], self.primary.spacing[2])
        return f"Slice  {self.current_z} / {self.primary.nz - 1}\nz = {z_cm:.2f} cm"

    def _refresh_slice(self):
        for i, (vol, im) in enumerate(zip(self.volumes, self.im_handles)):
            sl = vol.get_slice(self.current_z)

            if i == 0:
                im.set_data(sl)
                im.set_cmap(vol.cmap)
                im.set_alpha(vol.alpha)
                im.autoscale()
            else:
                im.set_data(_solid_rgba_overlay(sl, vol.color, vol.alpha))

        self._slice_txt.set_text(self._slice_label())
        self.fig.canvas.draw_idle()

    def _redraw_leaves(self):
        self._sync_leaf_values_from_state()
        iso = self._iso_px()
        for side, leaf in self.leaves.items():
            leaf.make_patch(iso, self.ax, color=LEAF_COLORS[side], alpha=0.40)
        self.fig.canvas.draw_idle()

    # ── textbox handling ───────────────────────────────────────────────

    def _format_cm(self, value):
        return f"{value:.3f}"

    def _update_textboxes(self):
        if self._textbox_guard:
            return

        self._textbox_guard = True
        try:
            mapping = {
                'top': self.top_cm,
                'bottom': self.bottom_cm,
                'left': self.left_cm,
                'right': self.right_cm,
            }
            for side, value in mapping.items():
                self.coord_boxes[side].set_val(self._format_cm(value))
        finally:
            self._textbox_guard = False

    def _on_submit_coord(self, side, text):
        if self._textbox_guard:
            return

        try:
            value = float(text.strip())
        except ValueError:
            self._update_textboxes()
            return

        value = max(0.0, value)
        self._set_side_distance(side, value, redraw=True, save=False)

    def _set_side_distance(self, side, value, redraw=True, save=False):
        value = max(0.0, float(value))

        if side == 'top':
            self.top_cm = value
        elif side == 'bottom':
            self.bottom_cm = value
        elif side == 'left':
            self.left_cm = value
        elif side == 'right':
            self.right_cm = value
        else:
            raise ValueError(f"Unknown side: {side}")

        if redraw:
            self._update_textboxes()
            self._redraw_leaves()
            self._update_info()

        if save:
            self._write_output()

    # ── callbacks ─────────────────────────────────────────────────────

    def _on_prev(self, _):
        if self.current_z > 0:
            self.current_z -= 1
            self._refresh_slice()

    def _on_next(self, _):
        if self.current_z < self.primary.nz - 1:
            self.current_z += 1
            self._refresh_slice()

    def _on_key(self, event):
        if event.key in ('left', 'down'):
            self._on_prev(None)
        elif event.key in ('right', 'up'):
            self._on_next(None)

    def _on_angle(self, val):
        self.angle = float(val)
        for leaf in self.leaves.values():
            leaf.angle_deg = self.angle
        self._redraw_leaves()
        self._update_info()

    def _on_save(self, _):
        self._write_output()
        self.btn_save.label.set_text('✓ Saved!')
        self.fig.canvas.draw_idle()
        self.fig.canvas.start_event_loop(0.9)
        self.btn_save.label.set_text('💾 Save')
        self.fig.canvas.draw_idle()

    # ── drag ──────────────────────────────────────────────────────────

    def _on_press(self, event):
        if event.inaxes != self.ax or event.button != 1:
            return

        for side in ('top', 'bottom', 'left', 'right'):
            leaf = self.leaves[side]
            if leaf.contains_point(event):
                self._drag_leaf = leaf
                self._drag_start = (event.xdata, event.ydata)
                self._drag_coord_start = leaf.inner_dist_cm
                return

    def _on_motion(self, event):
        if self._drag_leaf is None or event.inaxes != self.ax:
            return

        dx = event.xdata - self._drag_start[0]
        dy = event.ydata - self._drag_start[1]
        d_r = rotation_matrix_2d(-self.angle) @ np.array([dx, dy])

        side = self._drag_leaf.side
        if side == 'top':
            new_d = self._drag_coord_start + d_r[1] * self.sp_y
        elif side == 'bottom':
            new_d = self._drag_coord_start - d_r[1] * self.sp_y
        elif side == 'left':
            new_d = self._drag_coord_start - d_r[0] * self.sp_x
        else:  # right
            new_d = self._drag_coord_start + d_r[0] * self.sp_x

        self._set_side_distance(side, max(0.0, new_d), redraw=True, save=False)

    def _on_release(self, event):
        if self._drag_leaf is not None:
            self._write_output()

        self._drag_leaf = None
        self._drag_start = None
        self._drag_coord_start = None

    # ── info panel ────────────────────────────────────────────────────

    def _update_info(self):
        ax = self.ax_info
        ax.cla()
        ax.axis('off')
        ax.set_facecolor('#0d1117')

        tb_open, lr_open = self._field_opening()
        p = self.primary
        iso_x_cm = (self.iso_x * p.spacing[0] + p.origin[0]) / 10.0
        iso_y_cm = (self.iso_y * p.spacing[1] + p.origin[1]) / 10.0
        iso_z_cm = (self.iso_z * p.spacing[2] + p.origin[2]) / 10.0

        rows = [
            ("COLLIMATOR", '#00CFFF', 11, 'bold'),
            ("", '#aaa', 5, 'normal'),
            ("Leaf size", '#888', 8, 'normal'),
            (f"W = {self.rect_w:.2f} cm   H = {self.rect_h:.2f} cm", '#ddd', 10, 'normal'),
            ("", '#aaa', 5, 'normal'),
            ("Slit distances from centre", '#888', 8, 'normal'),
            (f"Top    : {self.top_cm:.3f} cm", '#00CFFF', 10, 'normal'),
            (f"Bottom : {self.bottom_cm:.3f} cm", '#00CFFF', 10, 'normal'),
            (f"Left   : {self.left_cm:.3f} cm", '#FF9F00', 10, 'normal'),
            (f"Right  : {self.right_cm:.3f} cm", '#FF9F00', 10, 'normal'),
            ("", '#aaa', 5, 'normal'),
            ("Openings", '#888', 8, 'normal'),
            (f"Top + Bottom : {tb_open:.3f} cm", '#00CFFF', 10, 'normal'),
            (f"Left + Right : {lr_open:.3f} cm", '#FF9F00', 10, 'normal'),
            ("", '#aaa', 5, 'normal'),
            ("Angle", '#888', 8, 'normal'),
            (f"{self.angle:.1f}°", '#ff4444', 14, 'bold'),
            ("", '#aaa', 5, 'normal'),
            ("Isocentre", '#888', 8, 'normal'),
            (f"x = {iso_x_cm:.3f} cm", '#ddd', 9, 'normal'),
            (f"y = {iso_y_cm:.3f} cm", '#ddd', 9, 'normal'),
            (f"z = {iso_z_cm:.3f} cm", '#ddd', 9, 'normal'),
            (f"vox ({self.iso_x}, {self.iso_y}, {self.iso_z})", '#666', 8, 'normal'),
        ]

        y = 0.98
        dy = 0.050
        for text, color, size, weight in rows:
            ax.text(
                0.03, y, text,
                color=color,
                fontsize=size,
                fontweight=weight,
                va='top',
                fontfamily='monospace',
                transform=ax.transAxes
            )
            y -= dy

        self.fig.canvas.draw_idle()

    # ── output ────────────────────────────────────────────────────────
    def _write_output(self):
      tb_open, lr_open = self._field_opening()

      out_path = Path(self.output_file).with_suffix(".out")

      with open(out_path, 'w') as f:
          f.write(f"leaf_width_cm: {round(self.rect_w, 4)}\n")
          f.write(f"leaf_height_cm: {round(self.rect_h, 4)}\n")
          f.write(f"top_cm: {round(self.top_cm, 4)}\n")
          f.write(f"bottom_cm: {round(self.bottom_cm, 4)}\n")
          f.write(f"left_cm: {round(self.left_cm, 4)}\n")
          f.write(f"right_cm: {round(self.right_cm, 4)}\n")
          f.write(f"field_opening_top_bottom_cm: {round(tb_open, 4)}\n")
          f.write(f"field_opening_left_right_cm: {round(lr_open, 4)}\n")
          f.write(f"rotation_angle_deg: {round(self.angle, 2)}\n")
          f.write(f"current_slice_z: {int(self.current_z)}\n")
          f.write(f"isocentre_voxel_x: {self.iso_x}\n")
          f.write(f"isocentre_voxel_y: {self.iso_y}\n")
          f.write(f"isocentre_voxel_z: {self.iso_z}\n")

      print(f"[saved] {out_path}")
  
      # file separato come prima
      grid_path = out_path.parent / "grid_size.out"
      with open(grid_path, 'w') as f:
          f.write(f"{round(self.top_cm, 4)}\n")
          f.write(f"{round(self.bottom_cm, 4)}\n")
          f.write(f"{round(self.left_cm, 4)}\n")
          f.write(f"{round(self.right_cm, 4)}\n")
          f.write(f"{round(self.angle, 2)}\n")

      print(f"[saved] {grid_path}")    


# CLI

def parse_args():
    p = argparse.ArgumentParser(
        description="BEV Collimator GUI — multi-volume MHD viewer with 4 independent slits.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python slitGUI.py --mhd ct.mhd

  python slitGUI.py \
      --mhd imgs/CT_plan.mhd imgs/PTV_plan.mhd imgs/Muscolo_plan.mhd \
      --label CT PTV MUSCOLO \
      --color gray red cyan \
      --alpha 0.7 0.5 0.5 \
      --angle 60

  python slitGUI.py \
      --mhd ct.mhd ptv.mhd \
      --top 1.5 --bottom 2.0 --left 3.0 --right 4.0
"""
    )

    p.add_argument(
        '--mhd',
        action='extend',
        nargs='+',
        required=True,
        metavar='PATH',
        help="One or more paths to MHD files."
    )

    p.add_argument(
        '--label',
        action='extend',
        nargs='+',
        default=None,
        metavar='LABEL',
        help="One or more display labels in volume order."
    )

    p.add_argument(
        '--color',
        action='extend',
        nargs='+',
        default=None,
        metavar='COLOR',
        help="One or more colors in volume order. Use 'gray' or 'none' for grayscale."
    )

    p.add_argument(
        '--alpha',
        action='extend',
        nargs='+',
        type=float,
        default=None,
        metavar='ALPHA',
        help="One or more opacities in [0, 1] in volume order."
    )

    p.add_argument('--rect-w', type=float, default=4.0,
                   help="Leaf width in cm (beam direction). Default: 4.0")
    p.add_argument('--rect-h', type=float, default=8.0,
                   help="Leaf height in cm (lateral). Default: 8.0")

    p.add_argument('--top', type=float, default=2.0,
                   help="Top slit inner-edge distance from isocentre [cm]. Default: 2.0")
    p.add_argument('--bottom', type=float, default=2.0,
                   help="Bottom slit inner-edge distance from isocentre [cm]. Default: 2.0")
    p.add_argument('--left', type=float, default=3.0,
                   help="Left slit inner-edge distance from isocentre [cm]. Default: 3.0")
    p.add_argument('--right', type=float, default=3.0,
                   help="Right slit inner-edge distance from isocentre [cm]. Default: 3.0")

    p.add_argument('--angle', type=float, default=0.0,
                   help="Initial collimator rotation [degrees]. Default: 0")
    p.add_argument('--output', default='collimator_output.json',
                   help="Output JSON file. Default: collimator_output.json")

    return p.parse_args()


def main():
    args = parse_args()
    n = len(args.mhd)

    labels = list(args.label or [])
    colors = list(args.color or [])
    alphas = list(args.alpha or [])

    try:
        _validate_list_lengths(n, labels, colors, alphas)
        _validate_alphas(alphas)

        for name in ('top', 'bottom', 'left', 'right'):
            value = getattr(args, name)
            if value < 0:
                raise ValueError(f"{name} must be >= 0, got {value}")
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    labels += [None] * (n - len(labels))
    colors += [None] * (n - len(colors))
    alphas += [None] * (n - len(alphas))

    volumes = []
    for i, path_str in enumerate(args.mhd):
        path = Path(path_str)
        if not path.exists():
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

        label = labels[i]

        raw_color = colors[i]
        try:
            if raw_color is not None:
                color = _resolve_color(raw_color)
            else:
                color = None if i == 0 else (_AUTO_COLORS[i] if i < len(_AUTO_COLORS) else '#ffffff')
        except ValueError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            sys.exit(1)

        alpha = alphas[i] if alphas[i] is not None else (1.0 if i == 0 else 0.45)

        print(f"Loading [{label or path.stem}]  {path}  color={color}  alpha={alpha} …")
        vol = VolumeDesc(path=str(path), label=label, color=color, alpha=alpha)
        print(f"  shape={vol.array.shape}  spacing={vol.spacing}  origin={vol.origin}")
        volumes.append(vol)

    BEVApp(volumes, args)


if __name__ == '__main__':
    main()
