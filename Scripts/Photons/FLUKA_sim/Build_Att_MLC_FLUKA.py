#!/usr/bin/env python3

"""

Uso:

  python Build_Att_MLC.py -txt thicknesses.txt -inp geometry.inp -field N
 gemetry.inp sono gli inp di FLUKA in FLUKAsim con dentro già rototraslazione, bisogna aggiungere la ct e spostare end trasfrom. Thickness.txt lo fa MakeAttenuator.py
"""
 
import argparse

import os
 
# ===========================================================================

# COSTANTI GEOMETRIA (cm)

# ===========================================================================

Z_ALDISC_END  = 2.006

Z_ATT_OFFSET  = 34.0

MLC_THICKNESS = 3.0

MLC_MARGIN    = 10.0
 
Z_ATT_START = Z_ALDISC_END + Z_ATT_OFFSET

Z_MLC_END   = Z_ATT_START

Z_MLC_START = Z_MLC_END - MLC_THICKNESS
 
CHUNK = 8
 
# ===========================================================================

# CARICAMENTO FILE SPESSORI

# ===========================================================================
 
def load_thicknesses(path):

    cells = {}

    with open(path) as f:

        for line in f:

            s = line.strip()

            if not s or s.startswith('PIN_ID'):

                continue

            p = s.split()

            if len(p) < 6:

                continue

            try:

                pin_id    = int(p[0])

                xmin      = float(p[1])

                xmax      = float(p[2])

                ymin      = float(p[3])

                ymax      = float(p[4])

                thickness = float(p[5])

                cells[pin_id] = (xmin, xmax, ymin, ymax, thickness)

            except ValueError:

                continue

    return cells
 
# ===========================================================================

# GEOMETRIA MLC

# ===========================================================================
 
def make_mlc_frame(cells):

    """

    Restituisce il dizionario per MLCOUT (il rettangolo esterno in tungsteno)

    e le info del bounding box per il riepilogo a schermo.

    I buchi non sono piu' un singolo MLCIN ma uno per cella (make_hole_bodies).

    """

    xmin_glob = min(v[0] for v in cells.values())

    xmax_glob = max(v[1] for v in cells.values())

    ymin_glob = min(v[2] for v in cells.values())

    ymax_glob = max(v[3] for v in cells.values())
 
    mlcout = {

        'xmin': xmin_glob - MLC_MARGIN,

        'xmax': xmax_glob + MLC_MARGIN,

        'ymin': ymin_glob - MLC_MARGIN,

        'ymax': ymax_glob + MLC_MARGIN,

        'zmin': Z_MLC_START,

        'zmax': Z_MLC_END,

        'bname': 'MLCOUT',

    }

    bbox = {

        'xmin': xmin_glob, 'xmax': xmax_glob,

        'ymin': ymin_glob, 'ymax': ymax_glob,

    }

    return mlcout, bbox
 
 
def make_hole_bodies(cells):

    """

    Per ogni cella crea un RPP alla z dell'MLC con la stessa xy della cella.

    Questi RPP sono i buchi nel telaio MLC -> il buco segue esattamente

    il perimetro dell'attenuatore, senza aperture dove non ci sono celle.
 
    Body name : H{pin_id}   (es. H2, H250)

    Region name: AH{pin_id} (es. AH2, AH250) -> aria nel buco

    """

    holes = []

    for pin_id in sorted(cells):

        xmin, xmax, ymin, ymax, _ = cells[pin_id]

        holes.append({

            'xmin': xmin, 'xmax': xmax,

            'ymin': ymin, 'ymax': ymax,

            'zmin': Z_MLC_START,

            'zmax': Z_MLC_END,

            'bname': f'H{pin_id}',

            'rname': f'AH{pin_id}',

        })

    return holes
 
# ===========================================================================

# UTILITIES PER MODIFICA .INP

# ===========================================================================
 
def find_geometry_bounds(lines):

    geobegin = next(i for i, l in enumerate(lines) if 'GEOBEGIN' in l)

    geoend   = next(i for i, l in enumerate(lines) if 'GEOEND'   in l)

    ends = [i for i, l in enumerate(lines)

            if l.strip() == 'END' and geobegin < i < geoend]

    return geobegin, geoend, ends[0], ends[1]
 
 
def build_vac_block(lines, att_body_names):

    """

    VAC sottrae MLCOUT (il telaio esterno) e tutte le celle attenuatore.

    I buchi H{pid} sono dentro MLCOUT quindi sono gia' esclusi dalla VAC

    quando si sottrae MLCOUT -> non servono sottrazioni aggiuntive per i buchi.

    """

    vac_start = None

    for i, line in enumerate(lines):

        if line.strip().startswith('VAC'):

            vac_start = i

            j = i + 1

            while j < len(lines) and (lines[j].startswith(' ') or

                                       lines[j].startswith('\t')):

                j += 1

            vac_end = j

            break
 
    if vac_start is None:

        raise RuntimeError('Regione VAC non trovata nel file .inp')
 
    base = ('VAC          5 +cnero  -TaDisc -AlDisc '

            '-LinacTi -VacLinac -AirBet -C1 -C2 -MLCOUT\n')
 
    conts = []

    for k in range(0, len(att_body_names), CHUNK):

        chunk = att_body_names[k:k + CHUNK]

        conts.append('   ' + ' '.join(f'-{n}' for n in chunk) + '\n')
 
    return lines[:vac_start] + [base] + conts + lines[vac_end:]
 
# ===========================================================================

# MODIFICA PRINCIPALE DEL .INP

# ===========================================================================
 
def modify_inp(inp_path, cells, mlcout, holes):

    with open(inp_path) as f:

        lines = f.readlines()
 
    _, _, bodies_end, regions_end = find_geometry_bounds(lines)
 
    # ------------------------------------------------------------------

    # 1. Bodies: attenuatore (V{pid})

    # ------------------------------------------------------------------

    att_bodies = []

    for pin_id in sorted(cells):

        xmin, xmax, ymin, ymax, sp = cells[pin_id]

        if sp <= 0.0:

            continue

        bname = f'V{pin_id}'

        att_bodies.append(

            f'RPP {bname:<8} '

            f'{xmin:.4f} {xmax:.4f} '

            f'{ymin:.4f} {ymax:.4f} '

            f'{Z_ATT_START:.4f} {Z_ATT_START + sp:.4f}\n'

        )
 
    # ------------------------------------------------------------------

    # 2. Bodies: MLCOUT + buchi per cella H{pid}

    # ------------------------------------------------------------------

    def rpp_line(b):

        return (f'RPP {b["bname"]:<8} '

                f'{b["xmin"]:.4f} {b["xmax"]:.4f} '

                f'{b["ymin"]:.4f} {b["ymax"]:.4f} '

                f'{b["zmin"]:.4f} {b["zmax"]:.4f}\n')
 
    mlc_bodies = [rpp_line(mlcout)] + [rpp_line(h) for h in holes]
 
    lines = lines[:bodies_end] + att_bodies + mlc_bodies + lines[bodies_end:]
 
    # ------------------------------------------------------------------

    # 3. Aggiorna la regione VAC

    # ------------------------------------------------------------------

    att_body_names = [f'V{pid}' for pid in sorted(cells) if cells[pid][4] > 0.0]

    lines = build_vac_block(lines, att_body_names)
 
    # ------------------------------------------------------------------

    # 4. Regions (indici ricalcolati)

    # ------------------------------------------------------------------

    _, _, _, regions_end = find_geometry_bounds(lines)
 
    # Regioni celle attenuatore

    att_regions = [f'R{pid:<11} 5 +V{pid}\n'

                   for pid in sorted(cells) if cells[pid][4] > 0.0]
 
    # Regione MLC: MLCOUT meno tutti i buchi per cella

    hole_names = [h['bname'] for h in holes]

    rmlc_base = 'RMLC         5 +MLCOUT'

    # Prima riga: fino a CHUNK buchi

    first_chunk = hole_names[:CHUNK]

    rmlc_line = rmlc_base + ' ' + ' '.join(f'-{n}' for n in first_chunk) + '\n'

    rmlc_conts = []

    for k in range(CHUNK, len(hole_names), CHUNK):

        chunk = hole_names[k:k + CHUNK]

        rmlc_conts.append('   ' + ' '.join(f'-{n}' for n in chunk) + '\n')

    rmlc_lines = [rmlc_line] + rmlc_conts
 
    # Regioni aria nei buchi AH{pid}

    air_regions = [f'{h["rname"]:<12} 5 +{h["bname"]}\n' for h in holes]
 
    lines = (lines[:regions_end]

             + att_regions

             + rmlc_lines

             + air_regions

             + lines[regions_end:])
 
    # ------------------------------------------------------------------

    # 5. ASSIGNMA

    # ------------------------------------------------------------------

    last_am = max(i for i, l in enumerate(lines) if 'ASSIGNMA' in l)
 
    att_am  = [f'ASSIGNMA    TUNGSTEN  R{pid}\n'

               for pid in sorted(cells) if cells[pid][4] > 0.0]

    mlc_am  = ['ASSIGNMA    BLCKHOLE  RMLC\n']

    air_am  = [f'ASSIGNMA         AIR  {h["rname"]}\n' for h in holes]
 
    lines = lines[:last_am + 1] + att_am + mlc_am + air_am + lines[last_am + 1:]
 
    return lines
 
# ===========================================================================

# MAIN

# ===========================================================================
 
def main():

    parser = argparse.ArgumentParser(

        description='Aggiunge attenuatore passivo e MLC a telaio a un file FLUKA .inp'

    )

    parser.add_argument('-txt', required=True,

                        help='File .txt con spessori celle')

    parser.add_argument('-inp', required=True,

                        help='File geometria FLUKA .inp da modificare')

    parser.add_argument('-field', required=True, type=int,

                        help='Numero del field (usato nel nome del file di output)')

    args = parser.parse_args()
 
    base     = os.path.splitext(os.path.basename(args.inp))[0]

    out_name = f'{base}_field{args.field}_att_mlc.inp'
 
    cells = load_thicknesses(args.txt)

    n_att = sum(1 for v in cells.values() if v[4] > 0.0)

    print(f'Caricate {len(cells)} celle, di cui {n_att} con spessore > 0')
 
    mlcout, bbox = make_mlc_frame(cells)

    holes = make_hole_bodies(cells)

    print(f'Generati {len(holes)} buchi nell\'MLC (uno per cella)')
 
    print(f'\nGeometria attenuatore:')

    print(f'  Z_ATT_START  = {Z_ATT_START:.3f} cm')

    print(f'  Bounding box x: [{bbox["xmin"]:.3f}, {bbox["xmax"]:.3f}] cm')

    print(f'  Bounding box y: [{bbox["ymin"]:.3f}, {bbox["ymax"]:.3f}] cm')
 
    print(f'\nGeometria MLC:')

    print(f'  Z_MLC_START  = {Z_MLC_START:.3f} cm')

    print(f'  Z_MLC_END    = {Z_MLC_END:.3f} cm')

    print(f'  MLCOUT x: [{mlcout["xmin"]:.3f}, {mlcout["xmax"]:.3f}] cm')

    print(f'  MLCOUT y: [{mlcout["ymin"]:.3f}, {mlcout["ymax"]:.3f}] cm')

    print(f'  Margine:  {MLC_MARGIN} cm per lato')
 
    lines_out = modify_inp(args.inp, cells, mlcout, holes)

    with open(out_name, 'w') as f:

        f.writelines(lines_out)
 
    print(f'\nFile di output: {out_name}')
 
if __name__ == '__main__':

    main()
 
