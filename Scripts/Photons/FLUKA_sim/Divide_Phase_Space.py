#!/usr/bin/env python
import argparse
import numpy as np
import struct
import os
import matplotlib.path as mpath
import matplotlib.patches as patches


# ─────────────────────────────────────────────────────────────────────────────
# Geometria: intersezione cella ↔ contorno
# ─────────────────────────────────────────────────────────────────────────────

def segments_intersect(p1, p2, p3, p4):
    def cross2d(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def on_segment(p, q, r):
        return (min(p[0], r[0]) <= q[0] <= max(p[0], r[0]) and
                min(p[1], r[1]) <= q[1] <= max(p[1], r[1]))

    d1 = cross2d(p3, p4, p1)
    d2 = cross2d(p3, p4, p2)
    d3 = cross2d(p1, p2, p3)
    d4 = cross2d(p1, p2, p4)

    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and \
       ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True

    if d1 == 0 and on_segment(p3, p1, p4): return True
    if d2 == 0 and on_segment(p3, p2, p4): return True
    if d3 == 0 and on_segment(p1, p3, p2): return True
    if d4 == 0 and on_segment(p1, p4, p2): return True

    return False


def cell_intersects_contour_exact(ix, iy, grid_xmin, grid_ymin, step,
                                   path_ptv, contour_pts):
    x0 = grid_xmin + ix * step
    y0 = grid_ymin + iy * step
    x1 = x0 + step
    y1 = y0 + step

    corners = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
    if np.any(path_ptv.contains_points(corners)):
        return True

    cell_sides = [
        ((x0, y0), (x1, y0)),
        ((x1, y0), (x1, y1)),
        ((x1, y1), (x0, y1)),
        ((x0, y1), (x0, y0)),
    ]
    n = len(contour_pts)
    for i in range(n):
        p3 = contour_pts[i]
        p4 = contour_pts[(i + 1) % n]
        for (p1, p2) in cell_sides:
            if segments_intersect(p1, p2, p3, p4):
                return True

    return False


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Grid partitioning + rototraslazione → file .txt per cella')
    parser.add_argument("-i",        required=True,
                        help="File phase space binario")
    parser.add_argument("-p",        nargs=3, type=float, required=True,
                        metavar=('px','py','pz'),
                        help="Posizione target (traslazione)")
    parser.add_argument("-l",        nargs=3, type=float, required=True,
                        metavar=('lx','ly','lz'),
                        help="Vettore l (colonna 0 matrice R)")
    parser.add_argument("-u",        nargs=3, type=float, required=True,
                        metavar=('ux','uy','uz'),
                        help="Vettore u (colonna 1 matrice R)")
    parser.add_argument("-f",        nargs=3, type=float, required=True,
                        metavar=('fx','fy','fz'),
                        help="Vettore f (colonna 2 matrice R)")
    parser.add_argument("-contour",  required=True,
                        help="File contorno PTV (coppie X Y per riga)")
    parser.add_argument("-grid",     type=float, default=0.3,
                        help="Passo griglia in cm (default: 0.3)")
    parser.add_argument("-outdir",   default="output_cells",
                        help="Cartella di output (default: output_cells)")
    parser.add_argument("-V",        type=int, default=0,
                        help="Livello di verbosità [0-5]")
    args = parser.parse_args()

    # ── 1. Matrice di rotazione e vettore di traslazione ──────────────────────
    l_vec = np.array(args.l, dtype=np.float64)
    u_vec = np.array(args.u, dtype=np.float64)
    f_vec = np.array(args.f, dtype=np.float64)

    l_vec /= np.linalg.norm(l_vec)
    u_vec /= np.linalg.norm(u_vec)
    f_vec /= np.linalg.norm(f_vec)

    R        = np.column_stack([l_vec, u_vec, f_vec])   # 3×3 ortonormale

    # Punto di sparo + 39 cm lungo la direzione f
    p_shoot  = np.array(args.p, dtype=np.float64)
    z0_phsp = -8.0
    T_target = p_shoot + (44) * f_vec
    print(f"[ROT] Punto di sparo:   {p_shoot}")
    print(f"[ROT] Punto traslazione: {T_target}  (p + 39*f)")

    # ── 2. Lettura phase space binario ────────────────────────────────────────
    print("[PS] Lettura phase space...")
    with open(args.i, 'rb') as fin:
        header = fin.read(12)
        ivers, nprim_per_pulse, numpulses = struct.unpack('iii', header)
        pulse_offsets = list(struct.unpack(f'{numpulses}Q',
                                           fin.read(8 * numpulses)))
        fin.seek(0, 2)
        pulse_offsets.append(fin.tell())

        rec_size = struct.calcsize('f' * 8)
        total_npart = sum(
            (pulse_offsets[i+1] - pulse_offsets[i]) // rec_size
            for i in range(numpulses)
        )

        print(f"[PS] Pulse: {numpulses}  |  Particelle totali: {total_npart}")

        phsp = np.zeros((total_npart, 8), dtype=np.float32)
        curr = 0
        for i in range(numpulses):
            fin.seek(pulse_offsets[i])
            n_in_pulse = (pulse_offsets[i+1] - pulse_offsets[i]) // rec_size
            chunk = struct.unpack('f' * 8 * n_in_pulse,
                                  fin.read(pulse_offsets[i+1] - pulse_offsets[i]))
            phsp[curr:curr + n_in_pulse, :] = np.array(chunk).reshape(-1, 8)
            curr += n_in_pulse

#   phsp[:, 3] += 47.8
    # ── 3. Scambio ID: 1 → 3 e 3 → 1 ────────────────────────────────────────
    mask_1 = phsp[:, 0] == 1.0
    mask_3 = phsp[:, 0] == 3.0
    phsp[mask_1, 0] = 3.0
    phsp[mask_3, 0] = 7.0
    print(f"[ID] Scambio ID: {mask_1.sum()} particelle 1→3, "
          f"{mask_3.sum()} particelle 3→1")

    # ── 4. Griglia 2D con maschera contorno ───────────────────────────────────
    print("[GRID] Caricamento contorno PTV...")
    contour_pts = np.loadtxt(args.contour)
    path_ptv    = mpath.Path(contour_pts)

    grid_xmin = np.min(contour_pts[:, 0])
    grid_xmax = np.max(contour_pts[:, 0])
    grid_ymin = np.min(contour_pts[:, 1])
    grid_ymax = np.max(contour_pts[:, 1])

    step    = args.grid
    nx_side = int(np.ceil((grid_xmax - grid_xmin) / step))
    ny_side = int(np.ceil((grid_ymax - grid_ymin) / step))

    print(f"[GRID] Area PTV: X[{grid_xmin:.3f}, {grid_xmax:.3f}]  "
          f"Y[{grid_ymin:.3f}, {grid_ymax:.3f}]")
    print(f"[GRID] Griglia: {nx_side}×{ny_side} = {nx_side*ny_side} celle  "
          f"(passo={step} cm)")
    print("[GRID] Calcolo celle attive...")

    active_cell_ids = []
    for cell_id in range(nx_side * ny_side):
        ix = cell_id % nx_side
        iy = cell_id // nx_side
        if cell_intersects_contour_exact(ix, iy, grid_xmin, grid_ymin,
                                          step, path_ptv, contour_pts):
            active_cell_ids.append(cell_id)

    active_cell_ids = np.array(active_cell_ids)
    print(f"[GRID] Celle attive: {len(active_cell_ids)} / {nx_side*ny_side}")

    # ── 5. Assegnazione particelle → celle ────────────────────────────────────
    ix_p = np.floor((phsp[:, 1] - grid_xmin) / step).astype(int)
    iy_p = np.floor((phsp[:, 2] - grid_ymin) / step).astype(int)

    mask_bounds = ((ix_p >= 0) & (ix_p < nx_side) &
                   (iy_p >= 0) & (iy_p < ny_side))
    phsp  = phsp[mask_bounds]
    ix_p  = ix_p[mask_bounds]
    iy_p  = iy_p[mask_bounds]

    particle_cell_idx = ix_p + iy_p * nx_side

    mask_active = np.isin(particle_cell_idx, active_cell_ids)
    phsp              = phsp[mask_active]
    particle_cell_idx = particle_cell_idx[mask_active]

    unique_cells = np.unique(particle_cell_idx)
    print(f"[GRID] Celle con almeno una particella: {len(unique_cells)}")

    # ── 6. Cartella di output ─────────────────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)

    # ── 7. Loop celle: rototraslazione + scrittura file ───────────────────────
    print(f"[OUT] Scrittura file .txt in: {args.outdir}")

    for cell_id in unique_cells:
        mask        = (particle_cell_idx == cell_id)
        phsp_cell   = phsp[mask]
        n_part      = len(phsp_cell)

        pos_loc = phsp_cell[:, 1:4].astype(np.float64)
        dir_loc = phsp_cell[:, 4:7].astype(np.float64)
        pos_rot = pos_loc @ R.T + T_target
        dir_rot = dir_loc @ R.T

        part_ids = phsp_cell[:, 0].astype(np.float64)
        ekin     = phsp_cell[:, 7].astype(np.float64)

        out_path = os.path.join(args.outdir, f'{int(cell_id)}.txt')
        with open(out_path, 'w') as fout:
            for j in range(n_part):
                fout.write(f"{part_ids[j]:.0f}  "
                           f"{pos_rot[j,0]:.6f}  {pos_rot[j,1]:.6f}  "
                           f"{pos_rot[j,2]:.6f}  "
                           f"{dir_rot[j,0]:.6f}  {dir_rot[j,1]:.6f}  "
                           f"{dir_rot[j,2]:.6f}  {ekin[j]:.6f}\n")

        if args.V >= 1:
            print(f"  cell {int(cell_id):5d}: {n_part:8d} particelle  →  {out_path}")

    print(f"\n[DONE] File scritti: {len(unique_cells)}  |  "
          f"Cartella: {args.outdir}")


if __name__ == "__main__":
    main()
