#!/usr/bin/env python
import argparse
import numpy as np
import struct
import random
import os
from fredemlib import *
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def segments_intersect(p1, p2, p3, p4):
    """
    Verifica se il segmento p1-p2 interseca il segmento p3-p4.
    Usa il metodo delle orientazioni (cross product).
    """
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
    """
    Logica B: una cella e' attiva se:
      (1) almeno un angolo della cella e' dentro il contorno, OPPURE
      (2) almeno un lato della cella interseca almeno un segmento del contorno.
    """
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


def plot_grid(contour_pts, grid_xmin, grid_xmax, grid_ymin, grid_ymax,
              step, active_cell_ids, nx_side, outdir):
    """Produce un grafico della griglia sovrapposta al contorno PTV."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    for cell_id in active_cell_ids:
        ix = cell_id % nx_side
        iy = cell_id // nx_side
        x0 = grid_xmin + ix * step
        y0 = grid_ymin + iy * step
        rect = patches.Rectangle(
            (x0, y0), step, step,
            linewidth=0.5,
            edgecolor='#0F6E56',
            facecolor='#1D9E75',
            alpha=0.25
        )
        ax.add_patch(rect)

    for x in np.arange(grid_xmin, grid_xmax + step, step):
        ax.axvline(x, color='gray', linewidth=0.2, alpha=0.4)
    for y in np.arange(grid_ymin, grid_ymax + step, step):
        ax.axhline(y, color='gray', linewidth=0.2, alpha=0.4)

    ptv_patch = patches.Polygon(
        contour_pts, closed=True,
        linewidth=2, edgecolor='#085041',
        facecolor='none', zorder=5
    )
    ax.add_patch(ptv_patch)

    ax.set_title(
        f'Griglia SPB sul PTV  —  passo={step} cm  |  celle attive={len(active_cell_ids)}',
        fontsize=11
    )
    ax.set_xlabel('X [cm]')
    ax.set_ylabel('Y [cm]')
    ax.set_xlim(grid_xmin - step * 2, grid_xmax + step * 2)
    ax.set_ylim(grid_ymin - step * 2, grid_ymax + step * 2)
    ax.set_aspect('equal')

    legend_elements = [
        patches.Patch(facecolor='#1D9E75', edgecolor='#0F6E56',
                      alpha=0.4, label='Cella attiva (interseca PTV)'),
        patches.Patch(facecolor='none', edgecolor='#085041',
                      linewidth=2, label='Contorno PTV')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)

    plot_path = os.path.join(outdir, 'grid_plot.png')
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"[PLOT] Griglia salvata in: {plot_path}")


def main():
    parser = argparse.ArgumentParser(description='FredEM: Grid Partitioning + Manual Rototranslation')
    parser.add_argument("-i", help="Input phasespace binary file", required=True)
    parser.add_argument("-npulse", help="Number of pulses to process", type=int, default=100000)
    parser.add_argument("-p", nargs=3, type=float, required=True, help="Target Position P (x y z)")
    parser.add_argument("-l", nargs=3, type=float, required=True, help="Vector l (left/column 0)")
    parser.add_argument("-u", nargs=3, type=float, required=True, help="Vector u (up/column 1)")
    parser.add_argument("-f", nargs=3, type=float, required=True, help="Vector f (forward/column 2)")
    parser.add_argument("-random_selection", action="store_true", help="Select pulses randomly")
    parser.add_argument("-contour", help="File del contorno field_X_contour.txt", required=True)
    parser.add_argument("-grid", help="Grid step size in cm", type=float, default=0.3)
    parser.add_argument("-V", help="Verbose level [0,5]", type=int, default=0)
    parser.add_argument("-outdir", help="Percorso della cartella di output", default="output_sim")
    args = parser.parse_args()

    # --- 1. Definizione Matrice di Rotazione e Traslazione ---
    l_vec = np.array(args.l, dtype=np.float64) / np.linalg.norm(args.l)
    u_vec = np.array(args.u, dtype=np.float64) / np.linalg.norm(args.u)
    f_vec = np.array(args.f, dtype=np.float64) / np.linalg.norm(args.f)

    R = np.column_stack([l_vec, u_vec, f_vec])
    T_target = np.array(args.p, dtype=np.float64)

    # --- 2. Caricamento Dati Phase Space ---
    fred = libFredEM()
    fred.verbose = args.V

    with open(args.i, 'rb') as fin:
        header = fin.read(12)
        ivers, nprim_per_pulse, numpulses = struct.unpack('iii', header)
        pulse_offsets = list(struct.unpack(f'{numpulses}Q', fin.read(8 * numpulses)))
        fin.seek(0, 2)
        pulse_offsets.append(fin.tell())

        idx_available = list(range(numpulses))
        if args.npulse > numpulses:
            print(f"[WARN] Richieste {args.npulse} pulse ma il file ne contiene solo {numpulses}.")

        if args.random_selection:
            selected_pulses = random.sample(idx_available, min(args.npulse, numpulses))
        else:
            selected_pulses = idx_available[:min(args.npulse, numpulses)]

        total_npart = 0
        for i in selected_pulses:
            total_npart += (pulse_offsets[i+1] - pulse_offsets[i]) // struct.calcsize('f' * 8)

        phsp = np.zeros((total_npart, 8), dtype=np.float32)
        curr = 0
        for i in selected_pulses:
            fin.seek(pulse_offsets[i])
            n_in_pulse = (pulse_offsets[i+1] - pulse_offsets[i]) // struct.calcsize('f' * 8)
            data_chunk = struct.unpack('f' * 8 * n_in_pulse, fin.read(pulse_offsets[i+1] - pulse_offsets[i]))
            phsp[curr:curr + n_in_pulse, :] = np.array(data_chunk).reshape(-1, 8)
            curr += n_in_pulse

    # --- 3. Definizione Griglia con Maschera Contorno (Logica B) ---
    contour_pts = np.loadtxt(args.contour)
    path_ptv = mpath.Path(contour_pts)

    grid_xmin = np.min(contour_pts[:, 0])
    grid_xmax = np.max(contour_pts[:, 0])
    grid_ymin = np.min(contour_pts[:, 1])
    grid_ymax = np.max(contour_pts[:, 1])

    step = args.grid
    nx_side = int(np.ceil((grid_xmax - grid_xmin) / step))
    ny_side = int(np.ceil((grid_ymax - grid_ymin) / step))

    print(f"Area PTV: X[{grid_xmin:.3f}, {grid_xmax:.3f}] Y[{grid_ymin:.3f}, {grid_ymax:.3f}]")
    print(f"Griglia totale: {nx_side}x{ny_side} = {nx_side*ny_side} celle  (passo={step} cm)")
    print(f"Calcolo celle attive con intersezione geometrica esatta...")

    active_cell_ids = []
    for cell_id in range(nx_side * ny_side):
        ix = cell_id % nx_side
        iy = cell_id // nx_side
        if cell_intersects_contour_exact(ix, iy, grid_xmin, grid_ymin,
                                          step, path_ptv, contour_pts):
            active_cell_ids.append(cell_id)

    active_cell_ids = np.array(active_cell_ids)
    print(f"Celle attive: {len(active_cell_ids)} su {nx_side*ny_side} totali")

    # Assegnazione particelle alle celle
    ix_p = np.floor((phsp[:, 1] - grid_xmin) / step).astype(int)
    iy_p = np.floor((phsp[:, 2] - grid_ymin) / step).astype(int)

    mask_bounds = (ix_p >= 0) & (ix_p < nx_side) & (iy_p >= 0) & (iy_p < ny_side)
    phsp = phsp[mask_bounds]
    ix_p = ix_p[mask_bounds]
    iy_p = iy_p[mask_bounds]

    particle_grid_idx = ix_p + iy_p * nx_side

    mask_in_active_cells = np.isin(particle_grid_idx, active_cell_ids)
    phsp = phsp[mask_in_active_cells]
    grid_indices = particle_grid_idx[mask_in_active_cells]

    unique_regions = np.unique(grid_indices)
    print(f"Regioni attive con particelle: {len(unique_regions)}")

    # --- Output dir ---
    output_dir = args.outdir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # --- Plot diagnostico griglia ---
    plot_grid(contour_pts, grid_xmin, grid_xmax, grid_ymin, grid_ymax,
              step, active_cell_ids, nx_side, output_dir)

    # --- EXIT: rimuovere dopo aver verificato la griglia ---
    #print("\n[EXIT] Uscita prima della simulazione per verifica griglia.")
    #exit(0)

    # --- 4. Setup Phantom ---
    fred.loadHUCalibration('/NAS_arpg/gaia/GITFRED/fredem/src/data/mat/HU2Materials.txt')
    HUCalMin, HUCalMax = fred.getHUCalRange()
    nn, hs, x0, CT, _ = fred.readMap('PZ1/Pancreas_PZ3/imrt_mhd_int16/Phantom.mhd')
    nn_m, hs_m, x0_m, MASK, _ = fred.readMap('PZ1/Pancreas_PZ3/imrt_mhd_int16/mask.mhd')
    clampedCT = np.clip(CT, HUCalMin, HUCalMax).astype(np.int32)
    fred.setPhantom(nn, hs, x0, material=clampedCT)
    fred.setMask(MASK)

    # --- 5. Preparazione Output e Metadata ---
    metadata_path = os.path.join(output_dir, "regions_metadata.txt")
    f_meta = open(metadata_path, "w")
    f_meta.write("SPB_ID\tN_Part\t"
                 "Xmin_cell\tXmax_cell\tYmin_cell\tYmax_cell\n")

    print(f"Particelle: {len(phsp)} | Regioni attive con particelle: {len(unique_regions)}")

    # --- 6. Loop sulle Regioni e Simulazione ---
    distanza = 0.0  # cm a monte del target lungo la direzione di ogni particella

    for r_id in unique_regions:
        mask = (grid_indices == r_id)
        phsp_region = phsp[mask].copy()

        if len(phsp_region) < 5:
            print(f"  [SKIP] Region {int(r_id)}: solo {len(phsp_region)} particelle, saltata.")
            continue

        # Bordi geometrici esatti della cella
        ix_cell = int(r_id) % nx_side
        iy_cell = int(r_id) // nx_side
        cell_xmin = grid_xmin + ix_cell * step
        cell_xmax = grid_xmin + (ix_cell + 1) * step
        cell_ymin = grid_ymin + iy_cell * step
        cell_ymax = grid_ymin + (iy_cell + 1) * step

        f_meta.write(f"{int(r_id)}\t{len(phsp_region)}\t"
                     f"{cell_xmin:.3f}\t{cell_xmax:.3f}\t{cell_ymin:.3f}\t{cell_ymax:.3f}\n")
        f_meta.flush()

        # --- TRASFORMAZIONE ---
        pos_loc = phsp_region[:, 1:4].astype(np.float64)
        dir_loc = phsp_region[:, 4:7].astype(np.float64)

        # Ruota posizioni e direzioni nel sistema del phantom
        pos_rot = pos_loc @ R.T + T_target
        dir_rot = dir_loc @ R.T

        # Sposta ogni particella a monte di `distanza` cm lungo la SUA direzione.
        # Le particelle di bordo (direzione inclinata) vengono posizionate
        # correttamente a monte del phantom, non con un offset fisso uguale per tutte.
        pos_entry = pos_rot - distanza * dir_rot

        phsp_region[:, 1:4] = pos_entry.astype(np.float32)
        phsp_region[:, 4:7] = dir_rot.astype(np.float32)

        # --- Esecuzione FredEM ---
        print(f"Region {int(r_id)}: Processing {len(phsp_region)} particles...", flush=True)
        fred.loadPhaseSpace(phsp_region)
        fred.clearScorers()
        fred.closeSetup()
        fred.trackPhaseSpace()

        nnD, hsD, x0D, Dose = fred.getMaskedDose(MASK)
        Dose /= len(phsp_region)#(nprim_per_pulse * len(selected_pulses))

        out_path = os.path.join(output_dir, f'Dose_SPB_{int(r_id)}.mhd')
        fred.writeMap(out_path, hsD, x0D, Dose)

    f_meta.close()
    print(f"\nSimulazione completata. Metadata salvati in: {metadata_path}")


if __name__ == "__main__":
    main()
