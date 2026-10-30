#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GetBestSize_LUTZ_nearest_PTV.py (versione compatta)

Sceglie l'apertura X/Y migliore tra le LUT (inviluppo rettangolare dell'isodose
95%) confrontandole con il PTV_plan in BEV, slice per slice lungo Z.

Per ogni slice PTV si usa la riga LUT con Z piu' vicina. Una LUT e' ammessa se la
copertura geometrica 3D del PTV e' >= --min-coverage; tra le ammesse vince quella
con minima differenza volumetrica (volume fuori PTV + volume PTV scoperto).
Se nessuna LUT e' ammessa: fallback 5x5 cm2 (NON validato, avviso su stderr).

OUTPUT: su stdout SOLO due righe (tutto il resto va su stderr):
    Xsize: <cm>
    Ysize: <cm>

Esempio:
    python3 GetBestSize_LUTZ_nearest_PTV.py -ptv imgs/PTV_plan.mhd \\
        -lutpath LUTs -E 7 > aperture.out

Unita': il PTV e' letto con FRED_mhd_io.readMap; spacing e origine sono gia' in cm.

NB: confronto GEOMETRICO e morfologico, non una DVH. Con --z-reference relative
Z=0 e' la prima slice occupata del PTV, non la profondita' fisica dalla superficie.
"""
import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from FRED_mhd_io import readMap  # spacing e origine gia in cm, Map indicizzata [x, y, z]


TOL = 1e-6  # cm: assorbe l'errore float32 di spacing/origine di readMap


def log(*args):
    print(*args, file=sys.stderr)


@dataclass
class Slice:
    z_cm: float
    x: np.ndarray
    y: np.ndarray
    area_cm2: float
    bbox_x_cm: float
    bbox_y_cm: float
    n_voxel: int


def parse_args():
    p = argparse.ArgumentParser(description='Apertura ottimale LUT vs PTV (output: Xsize/Ysize)')
    p.add_argument('-ptv', required=True, help='PTV_plan.mhd in BEV')
    p.add_argument('-lutpath', required=True, help='Directory LUTs')
    p.add_argument('-E', required=True, choices=('7', '9'), help='Energia [MeV]')
    p.add_argument('-min_size', nargs=2, type=float, metavar=('X', 'Y'), default=None,
                   help='Minima apertura X/Y nominale in cm')
    p.add_argument('--z-reference', choices=('relative', 'absolute'), default='relative')
    p.add_argument('--min-coverage', type=float, default=95.0,
                   help='Copertura geometrica 3D minima [%%]')
    p.add_argument('-v', '--verbose', action='store_true', help='Tabella per LUT su stderr')
    a = p.parse_args()
    if not math.isfinite(a.min_coverage) or not 0 <= a.min_coverage <= 100:
        p.error('--min-coverage richiede un valore fra 0 e 100')
    if a.min_size and not all(math.isfinite(v) and v > 0 for v in a.min_size):
        p.error('-min_size richiede due numeri positivi')
    return a


def candidate_files(folder, E):
    """LUT_<E>MeV_<X>x<Y>cm2 con X in 2..5 passo 0.2 e Y in 2,3,4,5."""
    pat = re.compile(r'^LUT_' + re.escape(E) + r'MeV_(\d+(?:\.\d+)?)x(\d+(?:\.\d+)?)cm2(?:\.[^.]+)?$')
    files = []
    for path in Path(folder).iterdir():
        m = pat.fullmatch(path.name) if path.is_file() else None
        if m is None:
            continue
        x, y = float(m[1]), float(m[2])
        ix = round((x - 2.0) / 0.2)
        iy = round((y - 2.0) / 0.2)
        if (0 <= ix <= 15 and abs(x - (2 + ix * 0.2)) < 1e-7
                and 0 <= iy <= 15 and abs(y - (2 + iy * 0.2)) < 1e-7):
            files.append((x, y, path))
    return sorted(files, key=lambda v: (v[0], v[1]))


def load_ptv(path, z_reference):
    # readMap (FRED_mhd_io): nn = n. voxel, hs = spacing [cm], x0 = origine [cm],
    # Map indicizzata [x, y, z], direction = matrice di orientamento.
    nn, hs, x0, Map, direction = readMap(path)
    arr = np.asarray(Map)
    mask = np.isfinite(arr) & (arr > 0)
    if mask.ndim != 3 or not mask.any():
        raise ValueError('PTV non valido: mappa non 3D o priva di voxel PTV')
    if direction is not None and np.size(direction) == 9:
        if not np.allclose(np.asarray(direction, dtype=float).reshape(3, 3), np.eye(3), atol=1e-5):
            raise ValueError('PTV non in assi XYZ BEV canonici (Direction non identita)')
    sp = np.round(np.asarray(hs, dtype=np.float64), 6)  # float32 -> 0.1f diventa 0.1
    origin = np.asarray(x0, dtype=np.float64)
    if sp.shape != (3,) or np.any(sp <= 0) or not np.all(np.isfinite(sp)):
        raise ValueError('Spacing non valido')

    coords = np.argwhere(mask)  # colonne: (x, y, z)
    ks = np.unique(coords[:, 2])
    center_x = origin[0] + 0.5 * (coords[:, 0].min() + coords[:, 0].max()) * sp[0]
    center_y = origin[1] + 0.5 * (coords[:, 1].min() + coords[:, 1].max()) * sp[1]
    z0 = origin[2] + int(ks[0]) * sp[2]
    slices = []
    for k in ks:
        xs, ys = np.nonzero(mask[:, :, int(k)])
        z = origin[2] + int(k) * sp[2]
        if z_reference == 'relative':
            z -= z0
        slices.append(Slice(float(z), origin[0] + xs * sp[0], origin[1] + ys * sp[1],
                            float(len(xs) * sp[0] * sp[1]),
                            float((xs.max() - xs.min() + 1) * sp[0]),
                            float((ys.max() - ys.min() + 1) * sp[1]), int(len(xs))))
    return slices, (float(center_x), float(center_y)), float(sp[0] * sp[1]), float(sp[2])


def load_lut(file):
    data = np.atleast_2d(np.loadtxt(file, comments='#', dtype=float))
    if data.shape[1] != 5:
        raise ValueError('LUT richiede 5 colonne: Z Xmin Xmax Ymin Ymax')
    if not np.isfinite(data[:, 0]).all():
        raise ValueError('Profondita Z non finita')
    data = data[np.argsort(data[:, 0], kind='stable')]
    if np.any(np.diff(data[:, 0]) <= 0):
        raise ValueError('Profondita Z duplicate nelle righe LUT')
    return data


def rectangle(row, sl, dx, dy, voxel_area, min_size):
    """Confronta una riga LUT con una slice PTV.
    Ritorna None se non utilizzabile, altrimenti (voxel coperti, area fuori, area scoperta)."""
    _, xmin, xmax, ymin, ymax = map(float, row)
    if not np.all(np.isfinite([xmin, xmax, ymin, ymax])) or xmin >= xmax or ymin >= ymax:
        return None
    w, h = xmax - xmin, ymax - ymin
    if min_size and (w + TOL < min_size[0] or h + TOL < min_size[1]):
        return None
    if w * h + TOL < sl.area_cm2 or w + TOL < sl.bbox_x_cm or h + TOL < sl.bbox_y_cm:
        return None
    inside = ((sl.x >= xmin + dx - TOL) & (sl.x <= xmax + dx + TOL)
              & (sl.y >= ymin + dy - TOL) & (sl.y <= ymax + dy + TOL))
    n = int(np.count_nonzero(inside))
    inter = n * voxel_area
    return n, max(0., w * h - inter), max(0., sl.area_cm2 - inter)


def evaluate(lut, slices, center, voxel_area, dz, min_size):
    """Ritorna (copertura % 3D, differenza volumetrica cm3) oppure None."""
    valid = np.isfinite(lut).all(axis=1) & (lut[:, 1] < lut[:, 2]) & (lut[:, 3] < lut[:, 4])
    if not valid.any():
        return None
    # Centro LUT unico per tutta la tabella (nessun recentramento per slice)
    ctr = np.median(np.column_stack(((lut[valid, 1] + lut[valid, 2]) * 0.5,
                                     (lut[valid, 3] + lut[valid, 4]) * 0.5)), axis=0)
    dx, dy = center[0] - float(ctr[0]), center[1] - float(ctr[1])
    tol = dz * 0.5 + TOL  # mezza slice di tolleranza agli estremi
    covered, extra, missed = 0, 0., 0.
    for sl in slices:
        res = None
        if lut[0, 0] - tol <= sl.z_cm <= lut[-1, 0] + tol:
            row = lut[np.argmin(np.abs(lut[:, 0] - sl.z_cm))]
            res = rectangle(row, sl, dx, dy, voxel_area, min_size)
        if res is None:  # slice non coperta: conta tutta come scoperta (conservativo)
            missed += sl.area_cm2 * dz
        else:
            n, outside, uncovered = res
            covered += n
            extra += outside * dz
            missed += uncovered * dz
    total = sum(s.n_voxel for s in slices)
    return 100. * covered / total, extra + missed


def main():
    a = parse_args()
    slices, center, voxel_area, dz = load_ptv(a.ptv, a.z_reference)
    vptv = sum(s.n_voxel for s in slices) * voxel_area * dz
    if vptv > 100. or vptv < 0.001:
        log(f'ATTENZIONE: volume PTV insolito ({vptv:.4f} cm3). Verifica unita e Map di FRED_mhd_io!')

    admitted = []  # (dV, area nominale, X, Y, copertura)
    for x, y, file in candidate_files(a.lutpath, a.E):
        if a.min_size and (x + TOL < a.min_size[0] or y + TOL < a.min_size[1]):
            continue
        try:
            res = evaluate(load_lut(file), slices, center, voxel_area, dz, a.min_size)
        except (OSError, ValueError) as exc:
            log(f'LUT {file.name} non analizzabile: {exc}')
            continue
        if res is None:
            continue
        cov, dv = res
        ok = cov + TOL >= a.min_coverage
        if a.verbose:
            log(f'{x:g}x{y:g}  cop={cov:6.2f}%  dV={dv:9.3f} cm3  {"AMMESSA" if ok else "scartata"}')
        if ok:
            admitted.append((dv, x * y, x, y, cov))

    if admitted:
        dv, _, x, y, cov = min(admitted)
        log(f'Migliore: {x:g}x{y:g} cm2 (copertura {cov:.2f}%, dV {dv:.3f} cm3)')
    else:
        x = y = 5.0
        log(f'ATTENZIONE: nessuna LUT con copertura >= {a.min_coverage:g}%. '
            'Fallback 5x5 cm2 NON validato: verificare con FRED/DVH.')
        if a.min_size and (a.min_size[0] > 5 or a.min_size[1] > 5):
            log('ATTENZIONE: il 5x5 e inferiore a -min_size richiesto.')
    print(f'Xsize: {x:g}')
    print(f'Ysize: {y:g}')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as ex:
        log('ERRORE:', ex)
        sys.exit(1)
