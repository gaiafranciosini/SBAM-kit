from __future__ import annotations

import argparse
import importlib.util
import csv
from contextlib import redirect_stdout, redirect_stderr
import datetime as dt
import json
import math
import os
import re
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

import numpy as np
try:
    import SimpleITK as sitk
except ImportError as exc:
    raise SystemExit('ERRORE: SimpleITK non disponibile: attiva CNR_tps_venv.') from exc

ROOT = Path.cwd().resolve()
SCRIPT_DIR = Path(__file__).resolve().parent
SIM_NAME = re.compile(r'AUTO_sim(7|9)MeV_(\d+)(?:_(plus2mm|minus2mm))?\Z')
MODES = ('CONV',)  # FIX: ottimizzazione SOLO CONV; FLASH resta escluso
CHECK_TOLERANCE_GY = 0.05  # Solo per report diagnostico, mai per lo status nominale
UPPER_PERCENT = 107.0  # Dmax 107%: limite di sicurezza, non soglia di copertura


def fail(msg):
    raise RuntimeError(msg)


def parse_setup(path):
    keys = {'PRESCRIPTION DOSE': 'rx_gy',
            'PERCENTAGE PRESCRIPTION DOSE': 'iso_pct',
            'PERCENTAGE PRESCRIPTION VOLUME': 'volume_pct'}
    result = {}
    for line_no, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        if ':' not in line:
            continue
        label, raw = line.split(':', 1)
        key = keys.get(label.strip().upper())
        if key is None:
            continue
        if key in result:
            fail(f'{path}:{line_no}: parametro duplicato: {label}')
        value = raw.split('#', 1)[0].strip().replace(',', '.')
        try:
            result[key] = float(value)
        except ValueError:
            fail(f'{path}:{line_no}: valore non numerico: {label} = {value}')
    missing = set(keys.values()) - result.keys()
    if missing:
        fail('setup.out incompleto: ' + ', '.join(sorted(missing)))
    if (not all(math.isfinite(v) for v in result.values())
            or result['rx_gy'] <= 0
            or not 0 < result['iso_pct'] <= 100
            or not 0 < result['volume_pct'] <= 100):
        fail(f'Prescrizione non valida: {result}')
    return result


def requested_energies(path):
    for line in path.read_text(encoding='utf-8').splitlines():
        if re.match(r'^\s*ENERGY\s*:', line, re.I):
            tail = line.split(':', 1)[1].split('#', 1)[0].strip()
            arr = tail.split()
            if not arr or len(arr) != len(set(arr)) or any(e not in ('7', '9') for e in arr):
                fail(f'ENERGY non valido in {path}: {tail!r}')
            return set(arr)
    fail(f'Manca ENERGY nel setup: {path}')


def geometry_equal(a, b):
    return (a.GetSize() == b.GetSize()
            and np.allclose(a.GetSpacing(), b.GetSpacing(), rtol=0, atol=1e-6)
            and np.allclose(a.GetOrigin(), b.GetOrigin(), rtol=0, atol=1e-6)
            and np.allclose(a.GetDirection(), b.GetDirection(), rtol=0, atol=1e-7))


def write_image_exclusive(image, destination):
    # Importante: mkdir della directory padre e verifica prima di generare .mhd/.raw.
    # Ogni cartella e' nuova, quindi nemmeno .raw puo' sovrascrivere file preesistenti.
    if destination.exists() or destination.with_suffix('.raw').exists():
        fail(f'PROTEZIONE SOVRASCRITTURA: file gia presente: {destination}')
    sitk.WriteImage(image, str(destination), useCompression=False)


def logged(command, report):
    if report.exists():
        fail(f'PROTEZIONE SOVRASCRITTURA: log gia presente: {report}')
    with report.open('x', encoding='utf-8') as log:
        log.write('COMMAND: ' + ' '.join(map(str, command)) + '\n\n')
        log.flush()
        rc = subprocess.run([str(v) for v in command], cwd=ROOT,
                            stdout=log, stderr=subprocess.STDOUT, check=False,
                            env={**os.environ, 'MPLBACKEND': 'Agg'}).returncode
    if rc:
        fail(f'Comando fallito (exit={rc}), vedi: {report}')


def header_table(headers, rows):
    data = [[str(v) for v in headers]] + [[str(v) for v in r] for r in rows]
    widths = [max(len(row[i]) for row in data) for i in range(len(headers))]
    lines = []
    for idx, row in enumerate(data):
        lines.append(' | '.join(item.ljust(widths[i]) for i, item in enumerate(row)))
        if idx == 0:
            lines.append('-+-'.join('-' * x for x in widths))
    return '\n'.join(lines)


def baseline_files(sim, energy, mode):
    if mode != 'CONV':
        fail('Ottimizzazione FLASH disabilitata: e consentita soltanto CONV.')
    pattern = f'DOSE_{energy}MeV_*pulses_{mode}_BEV.mhd'
    maps = sorted(p for p in sim.glob(pattern) if p.is_file() and not p.is_symlink())
    if len(maps) != 1:
        fail(f'{sim.name} {mode}: attesa una sola mappa BEV, trovate {len(maps)}: {maps}')
    match = re.fullmatch(fr'DOSE_{energy}MeV_(\d+)pulses_{mode}_BEV\.mhd', maps[0].name)
    if not match:
        fail(f'Formato del nome della dose non valido: {maps[0]}')
    p0 = int(match.group(1))
    if p0 < 1:
        fail('Numero impulsi iniziali non valido')
    final = sim / maps[0].name.replace('_BEV.mhd', '_FINAL.mhd')
    if not final.is_file() or final.is_symlink():
        fail(f'Mappa finale della Full Chain assente: {final}')
    dvh = sim / f'DVH{energy}MeV_{mode}'
    files = sorted(dvh.glob(f'*_plan{energy}MeV_{mode}.txt'))
    if not files or not any(p.name.lower() == f'ptv_plan{energy}mev_{mode}.txt'.lower() for p in files):
        fail(f'Curve DVH iniziali PTV/ROI non trovate: {dvh}')
    return p0, maps[0], final, files


# FIX: carica una volta sola il dvh_analyze.py aggiornato e usa le sue funzioni
# reali. Avviare un nuovo interprete con matplotlib ad OGNI impulso (6/12 piani)
# era molto lento; ora ogni DVH mantiene gli stessi calcoli, output e colonne.
_DVH_CACHE = {}


def get_dvh_module(script):
    key = str(Path(script).resolve())
    if key not in _DVH_CACHE:
        spec = importlib.util.spec_from_file_location('_dvh_analyze_pulse_run', key)
        if spec is None or spec.loader is None:
            fail(f'Impossibile importare dvh_analyze.py: {key}')
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        _DVH_CACHE[key] = module
    return _DVH_CACHE[key]


def save_scaled_dvh(copies, case_dir, energy, mode, ratio, scripts, setup_copy,
                    params, *, plot=False):
    output = case_dir / 'DVH'
    output.mkdir(exist_ok=False)
    generated = []
    for source in copies:
        arr = np.loadtxt(source, dtype=float, ndmin=2)
        if arr.ndim != 2 or arr.shape[1] != 2 or arr.shape[0] < 2 or not np.isfinite(arr).all():
            fail(f'DVH sorgente non valido: {source}')
        arr[:, 0] *= ratio  # Dose in cGy nel formato ComputeDVH
        dest = output / source.name
        with dest.open('x', encoding='utf-8') as fp:
            np.savetxt(fp, arr, fmt='%.12g')
        generated.append(dest)
    analysis = case_dir / 'dvh_analyze'
    analysis.mkdir(exist_ok=False)
    mod = get_dvh_module(scripts['dvh'])
    results, curves = [], []
    for src in generated:
        res, dose_gy, cumulative_volume = mod.analyze_file(
            src, params['rx_gy'], params['volume_pct'], params['iso_pct'], 0.1, [])
        results.append(res)
        curves.append((res, dose_gy, cumulative_volume))
    results.sort(key=lambda r: (float(r.energy) if r.energy != '?' else 999,
                                r.modality, r.structure, r.plan))
    tab = mod.format_table(results, params['rx_gy'],
                           params['volume_pct'], params['iso_pct'])
    with (analysis / 'dvh_table.txt').open('x', encoding='utf-8') as f:
        f.write(tab)
    csv_file = analysis / 'risultati_dvh.csv'
    if csv_file.exists():
        fail(f'PROTEZIONE SOVRASCRITTURA DVH: {csv_file}')
    mod.write_csv(results, csv_file, params['iso_pct'], params['volume_pct'])
    if plot:
        # FIX SOLO PRESENTAZIONE: conserva il log completo del plot nel caso,
        # ma non affolla il terminale con "Plot saved" e warning Matplotlib.
        with (analysis / 'plot_output.log').open('x', encoding='utf-8') as plot_log:
            with redirect_stdout(plot_log), redirect_stderr(plot_log):
                mod.make_plots(curves, analysis, params['rx_gy'],
                               params['volume_pct'], params['iso_pct'], show=False)
    with (case_dir / 'dvh_analyze.log').open('x', encoding='utf-8') as f:
        f.write('Funzioni del dvh_analyze.py copiato, eseguite nello stesso '
                'processo dell ottimizzatore per evitare riavvii ad ogni impulso.\n')
        f.write('TABELLA DVH (stesso calcolo dello script standalone):\n')
        f.write(tab)



def export_case(p, p0, dose_img, final_img, ptv_copy, source_dvhs, copy_scripts,
                setup_copy, params, mask, dose_voxels, plan_dir, energy, mode):
    case = plan_dir / 'casi' / f'{p}pulses'
    case.mkdir(exist_ok=False)
    ratio = p / p0
    base = f'DOSE_{energy}MeV_{p}pulses_{mode}'
    bev = case / f'{base}_BEV.mhd'
    final = case / f'{base}_FINAL.mhd'
    write_image_exclusive(sitk.Multiply(dose_img, ratio), bev)
    write_image_exclusive(sitk.Multiply(final_img, ratio), final)
    # Normalizzazione Rx su una directory di sola esportazione, mai sull'AUTO_sim.
    norm_log = case / 'normalize_dose22.5.log'
    logged([sys.executable, copy_scripts['normalize'], '--sim-dir', case,
            '--input', bev.name, '--setup', setup_copy, '--ptv', ptv_copy,
            '--dose-prescription-gy', params['rx_gy'],
            '--tolerance-gy', CHECK_TOLERANCE_GY], norm_log)
    normalized = sorted(case.glob('isodose_prescrizione_*/isodose_su_prescrizione/'
                                  + bev.stem + '_isodosemap_prescrizione.mhd'))
    if len(normalized) != 1:
        fail(f'Isodose di prescrizione assente o ambigua: {case}')
    # Il check di copertura con tolleranza e' solo DIAGNOSTICO.
    check_log = case / 'check_isodose_ptv22.5.txt'
    logged([sys.executable, copy_scripts['check'], '-i', normalized[0], '-p', ptv_copy,
            '--normalization', 'prescription', '--isodose-scale', 'fraction',
            '--dose-prescription-gy', params['rx_gy'],
            '--tolerance-gy', CHECK_TOLERANCE_GY,
            '--slice-threshold', params['volume_pct']], check_log)
    iso_image = sitk.ReadImage(str(normalized[0]))
    ptv_image = sitk.ReadImage(str(ptv_copy))
    if not geometry_equal(iso_image, ptv_image):
        fail(f'Geometria PTV/isodose non coincidente nel caso {case}')
    frac = sitk.GetArrayFromImage(iso_image)
    iso_threshold = params['iso_pct'] / 100.0
    coverage_real = 100.0 * np.count_nonzero(frac[mask] >= iso_threshold) / int(mask.sum())
    check_yes = coverage_real + 1e-9 >= params['volume_pct']
    # La ricerca e la normalizzazione devono concordare sulla copertura nominale.
    predicted = 100.0 * np.count_nonzero(dose_voxels * ratio >= params['rx_gy'] * iso_threshold) / dose_voxels.size
    if abs(predicted - coverage_real) > 100.0 / dose_voxels.size + 1e-8:
        fail(f'Copertura non coerente con la mappa normalizzata: {case}; '
             f'previsione {predicted:.8f}% vs reale {coverage_real:.8f}%')
    save_scaled_dvh(source_dvhs, case, energy, mode, ratio, copy_scripts, setup_copy, params, plot=True)
    return {'p': p, 'coverage_real': coverage_real, 'yes': check_yes,
            'Dmax_ptv': float(np.max(dose_voxels) * ratio),
            'dose_bev': str(bev), 'dose_final': str(final),
            'isodose': str(normalized[0]), 'check_report': str(check_log),
            'dvh_report': str(case / 'dvh_analyze' / 'dvh_table.txt')}


def optimize_plan(sim, mode, output, ptv_source, setup_copy, params, scripts,
                  max_extra):
    if mode != 'CONV':
        fail('Ottimizzazione FLASH disabilitata: e consentita soltanto CONV.')
    energy = SIM_NAME.fullmatch(sim.name).group(1)
    plan_name = f'{sim.name}_{mode}'
    plan_dir = output / 'piani' / plan_name
    plan_dir.mkdir(exist_ok=False)
    inp = plan_dir / 'input_copiati'
    inp.mkdir(exist_ok=False)
    source_dvhs_dir = inp / 'DVH_originali'
    source_dvhs_dir.mkdir(exist_ok=False)
    (plan_dir / 'casi').mkdir(exist_ok=False)
    p0, dose_orig, final_orig, dvhs = baseline_files(sim, energy, mode)
    # Snapshot isolato: nessun calcolo di output sulla cartella della Full Chain.
    src_dose = inp / dose_orig.name
    src_final = inp / final_orig.name
    write_image_exclusive(sitk.ReadImage(str(dose_orig)), src_dose)
    write_image_exclusive(sitk.ReadImage(str(final_orig)), src_final)
    dose = sitk.ReadImage(str(src_dose))
    final = sitk.ReadImage(str(src_final))
    local_dvhs = []
    for d in dvhs:
        target = source_dvhs_dir / d.name
        if target.exists():
            fail(f'Copia DVH gia presente: {target}')
        shutil.copy2(d, target)
        local_dvhs.append(target)
    ptv = sitk.ReadImage(str(ptv_source))
    if not geometry_equal(dose, ptv):
        fail(f'Geometrie BEV/PTV_plan diverse per {plan_name}. Non applico resampling impliciti.')
    ptv_mask = sitk.GetArrayFromImage(ptv) > 0
    if not np.any(ptv_mask):
        fail('PTV senza voxel')
    base_voxels = np.asarray(sitk.GetArrayFromImage(dose)[ptv_mask], dtype=np.float64)
    if not np.isfinite(base_voxels).all() or (base_voxels < 0).any():
        fail(f'Dosi nel PTV non valide: {plan_name}')
    base_max = float(base_voxels.max())
    if base_max <= 0:
        fail(f'Dmax PTV zero: {plan_name}')
    rx = params['rx_gy']
    iso_gy = rx * params['iso_pct'] / 100.0
    target_upper = rx * UPPER_PERCENT / 100.0
    max_p = p0 + max_extra
    memo = {}
    tested_pulses = []
    test_root = plan_dir / 'test_impulsi'
    test_root.mkdir(exist_ok=False)

    def metrics(p):
        # FIX: ricerca esclusivamente per CONV; FLASH non viene analizzata.
        # Le misure cliniche usate per DECIDERE sono sulla mappa voxel 3D,
        # con soglia nominale da setup.out, senza applicare 0.05 Gy.
        if p not in memo:
            ratio = p / p0
            cov = 100.0 * np.count_nonzero(base_voxels * ratio >= iso_gy) / base_voxels.size
            dmax = base_max * ratio
            memo[p] = {'p': p, 'cov': cov, 'dmax': dmax,
                       'yes': cov + 1e-9 >= params['volume_pct'],
                       'upper_ok': dmax <= target_upper + 1e-9,
                       'distance': abs(dmax - target_upper)}
        return memo[p]

    def test(p, phase):
        # FIX: ogni nuovo numero di impulsi produce la SUA tabella originale
        # dvh_analyze.py, inclusi PTV/OAR, in una directory ESCLUSIVA.
        # Non e' mai richiamato changePulse.sh e nessuna AUTO_sim e' modificata.
        if p in tested_pulses:
            return metrics(p)
        m = metrics(p)
        case = test_root / f'{p}pulses'
        case.mkdir(exist_ok=False)
        save_scaled_dvh(local_dvhs, case, energy, mode, p / p0,
                        scripts, setup_copy, params, plot=False)
        csv_file = case / 'dvh_analyze' / 'risultati_dvh.csv'
        with csv_file.open(newline='', encoding='utf-8') as f:
            ptv_rows = [r for r in csv.DictReader(f) if r['structure'].strip().upper() == 'PTV']
        if len(ptv_rows) != 1:
            fail(f'{plan_name}: manca o e duplicato PTV nel CSV per {p} impulsi')
        ptv_row = ptv_rows[0]
        level = f'V{params["iso_pct"]:g}_percent'
        vol_col = f'D{params["volume_pct"]:g}_Gy'
        if level not in ptv_row or vol_col not in ptv_row:
            fail(f'{plan_name}: colonne DVH attese non presenti: {level}, {vol_col}')
        m.update({'max_tab': float(ptv_row['max_tabulated_dose_Gy']),
                  'est_dmax': float(ptv_row['estimated_Dmax_Gy']),
                  'mean': float(ptv_row['weighted_mean_dose_Gy']),
                  'v_dvh': float(ptv_row[level]),
                  'dvol_dvh': float(ptv_row[vol_col]),
                  'phase': phase,
                  'ptv_dvh_status': ptv_row['PTV_prescription']})
        tested_pulses.append(p)
        # Righe simili alla tabella originale dvh_analyze.py, ma con numero
        # impulsi, Dmax REALE del PTV e la copertura nominale 3D ben visibili.
        row = [str(p), energy + ' MeV', mode, 'AUTO', 'PTV',
               f'{m["max_tab"]:.3f}', f'{m["est_dmax"]:.3f}',
               f'{m["mean"]:.3f}', f'{m["v_dvh"]:.3f}',
               f'{m["dvol_dvh"]:.3f}', f'{m["cov"]:.3f}',
               f'{m["dmax"]:.3f}', m['ptv_dvh_status'],
               'YES' if m['yes'] else 'NO', 'OK' if m['upper_ok'] else '>107%']
        headers = ['P', 'Energy', 'Mode', 'Plan', 'Structure', 'Max tab [Gy]',
                   'Est Dmax [Gy]', 'Wtd mean [Gy]', f'V{params["iso_pct"]:g}% DVH [%]',
                   f'D{params["volume_pct"]:g} [Gy]', 'PTV 3D [%]',
                   'Dmax PTV [Gy]', 'PTV DVH', 'PTV 3D', '107%']
        # FIX SOLO PRESENTAZIONE: lascia INVARIATI intestazioni, colonne e valori.
        label = 'INITIAL FULL CHAIN' if p == p0 else f'PULSE TEST ({p-p0:+d} vs initial)'
        print(f'\n{label} — {p} pulses', flush=True)
        print(header_table(headers, [row]), flush=True)
        # Copia della stessa tabella anche nella directory del test specifico.
        with (case / 'tabella_pulses.txt').open('x', encoding='utf-8') as f:
            f.write(header_table(headers, [row]) + '\n')
        return m

    initial = test(p0, 'INIZIO FULL CHAIN')
    # FASE 1: cerco MINIMO numero di impulsi con copertura PTV sufficiente.
    # Se p0 copre: -1 fino al primo NO. Se p0 non copre: +1 fino al primo YES.
    # Tutti i passi intermedi producono una tabella indipendente.
    if initial['yes']:
        p = p0
        while p > 1:
            if p0 - p >= max_extra:
                fail(f'{plan_name}: ricerca verso il basso incompleta dopo {max_extra} impulsi; '
                     'aumentare --max-extra-pulses senza modificare gli originali')
            # FIX SOLO PRESENTAZIONE: motivo del passaggio in UNA riga, in inglese.
            print(f'Coverage {metrics(p)["cov"]:.3f}% meets the {params["volume_pct"]:g}% target; '
                  f'test {p-1} pulses (-1).', flush=True)
            candidate = test(p - 1, 'FASE 1: SCENDO -1')
            if not candidate['yes']:
                print(f'Coverage {candidate["cov"]:.3f}% is below target; '
                      f'keep {p} pulses.', flush=True)
                break
            p -= 1
        min_p = p
    else:
        p = p0
        while not metrics(p)['yes']:
            if p - p0 >= max_extra:
                fail(f'{plan_name}: copertura PTV non raggiunta dopo +{max_extra} impulsi; '
                     'aumentare --max-extra-pulses, oppure controllare la dose')
            print(f'Coverage {metrics(p)["cov"]:.3f}% is below the {params["volume_pct"]:g}% target; '
                  f'test {p+1} pulses (+1).', flush=True)
            p += 1
            new = test(p, 'FASE 1: SALGO +1')
            if new['yes']:
                print(f'Coverage {new["cov"]:.3f}% meets target; '
                      f'stop at {p} pulses.', flush=True)
        min_p = p

    min_case = metrics(min_p)
    # FASE 2: la dose massima del PTV non deve superare il 107% di Rx.
    # Dmax aumenta monotonicamente con p: se il MINIMO che copre supera 107%,
    # non esiste alcun numero di impulsi che soddisfi entrambi i vincoli.
    # Se invece il minimo e gia <=107%, non aumentare inutilmente p:
    # l'utente ha richiesto il numero MINIMO che copre con Dmax accettabile.
    final_p = min_p
    final_case = min_case
    feasible = min_case['upper_ok']
    status = ('COPERTURA_OK_DMAX_OK' if feasible
              else 'INCOMPATIBILE: COPERTURA PTV OK MA DMAX > 107%')
    # FIX SOLO PRESENTAZIONE: la fase 2 resta nel riepilogo finale e nei file,
    # senza ripetere una spiegazione fra le tabelle dei diversi piani.

    phases = {'inizio_full_chain': p0,
              'fase1_copertura_minima': min_p,
              'fase2_verifica_107_finale': final_p}
    result_cases = {}
    for p in sorted(set(phases.values())):
        # Copie salvate come prima, ma senza messaggi intermedi a terminale.
        result_cases[p] = export_case(p, p0, dose, final, ptv_source, local_dvhs,
                                     scripts, setup_copy, params, ptv_mask, base_voxels,
                                     plan_dir, energy, mode)
        if result_cases[p]['yes'] != metrics(p)['yes']:
            fail(f'{plan_name}: esito nominale dopo normalizzazione diverso da ricerca a {p} impulsi')
    if not result_cases[min_p]['yes']:
        fail(f'{plan_name}: il minimo selezionato NON rispetta la copertura PTV')

    # Tutti i numeri provati, anche quelli NO, in una tabella ordinata.
    head = ['IMPULSI', 'FASE', 'Max tab [Gy]', 'Est Dmax [Gy]', 'Media [Gy]',
            f'V{params["iso_pct"]:g}% DVH [%]', f'D{params["volume_pct"]:g} [Gy]',
            'Copertura 3D [%]', 'PTV DVH', 'PTV 3D', 'Dmax PTV [Gy]', '107%', 'Dmax-107 [Gy]']
    rows = []
    for p in tested_pulses:
        m = metrics(p)
        rows.append([p, m['phase'], f'{m["max_tab"]:.3f}',
                     f'{m["est_dmax"]:.3f}', f'{m["mean"]:.3f}',
                     f'{m["v_dvh"]:.3f}', f'{m["dvol_dvh"]:.3f}',
                     f'{m["cov"]:.3f}', m['ptv_dvh_status'],
                     'YES' if m['yes'] else 'NO', f'{m["dmax"]:.3f}', 'YES' if m['upper_ok'] else 'NO',
                     f'{m["dmax"]-target_upper:+.3f}'])
    detailed = header_table(head, rows)
    lines = [f'PIANO: {plan_name}',
             f'SETUP: Rx={rx:g} Gy; isodose={params["iso_pct"]:g}%; '
             f'PTV >= {params["volume_pct"]:g}%',
             f'Vincolo PTV: V({iso_gy:.5f} Gy) >= {params["volume_pct"]:g}% (3D nominale)',
             f'Vincolo Dmax PTV: <= {target_upper:.5f} Gy (107% della Rx)',
             'Dmax PTV = massimo dei voxel appartenenti al PTV nella mappa BEV.',
             'Est Dmax [Gy] della tabella dvh_analyze.py = D0.1%, NON Dmax del voxel.',
             f'Full Chain: {p0} impulsi; copertura={initial["cov"]:.5f}%; '
             f'Dmax={initial["dmax"]:.5f} Gy.',
             f'Fase 1: {min_p} impulsi (minimo che soddisfa la copertura).',
             f'Fase 2: {status}',
             f'Scelta finale: {final_p} impulsi; PTV={"YES" if final_case["yes"] else "NO"}; '
             f'Dmax={final_case["dmax"]:.5f} Gy.',
             f'Differenza rispetto alla Full Chain: {final_p-p0:+d} impulsi (negativo = meno impulsi).',
             '', 'TUTTI I NUMERI DI IMPULSI PROVATI', detailed, '',
             'Le tabelle COMPLETE PTV/OAR per ogni test sono in test_impulsi/<P>pulses/dvh_analyze/dvh_table.txt.',
             'Le mappe BEV/FINAL, isodosi normalizzate e check ufficiale sono salvati SOLO come copie.',
             'PTV DVH = stato ricavato dai DVH; PTV 3D = esito nominale delle mappe normalizzate.',
             'Se DVH e PTV 3D discordano, controllare geometria e binning prima di valutare il piano.',
             'La copertura che decide YES/NO e nominale (senza 0.05 Gy); il check esterno produce anche la diagnostica con tolleranza.',
             'I DVH sono riscalati linearmente dai DVH iniziali (cGy), non rigenerati da FLUKA/ComputeDVH.',
             'Nessuna mappa o ROI originale della Full Chain viene modificata.']
    if not feasible:
        lines += ['', 'ATTENZIONE: NESSUN PIANO CONFORME AD ENTRAMBI I VINCOLI NEL MODELLO LINEARE.',
                  'Il piano indicato copre il PTV ma la Dmax supera il 107%.',
                  'Non usare il risultato come piano conforme alla dose massima.']
    with (plan_dir / 'tabella_piano.txt').open('x', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    with (plan_dir / 'tutti_impulsi.csv').open('x', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(head); w.writerows(rows)
    # FIX SOLO PRESENTAZIONE: la tabella dettagliata per piano rimane in
    # tabella_piano.txt e tutti_impulsi.csv, senza ripeterla a schermo.
    return dict(piano=plan_name, energia_MeV=energy, modalita=mode,
                variante=(SIM_NAME.fullmatch(sim.name).group(3) or 'originale'),
                impulsi_full_chain=p0, copertura_full_chain_percent=round(initial['cov'], 6),
                PTV_full_chain='YES' if initial['yes'] else 'NO',
                Dmax_full_chain_Gy=round(initial['dmax'], 6),
                impulsi_minimo_copertura=min_p,
                copertura_minima_percent=round(min_case['cov'], 6),
                Dmax_minima_Gy=round(min_case['dmax'], 6),
                delta_impulsi_minimo=min_p-p0,
                impulsi_vicino_107=final_p,
                copertura_vicino_107_percent=round(final_case['cov'], 6),
                Dmax_vicino_107_Gy=round(final_case['dmax'], 6),
                delta_impulsi_107=final_p-p0,
                scostamento_107_Gy=round(final_case['dmax']-target_upper, 6),
                PTV_vicino_107='YES' if final_case['yes'] else 'NO',
                PTV_e_Dmax_conformi='YES' if feasible else 'NO',
                fase2_status=status, impulsi_testati=len(tested_pulses),
                mappa_FullChain_sorgente=str(dose_orig),
                cartella_risultati=str(plan_dir))


def main():
    parser = argparse.ArgumentParser(description='Ottimizzatore impulsi SOLO CONV: 3 piani per energia (7/9 MeV), copie senza sovrascrittura')
    parser.add_argument('--sim', action='append', required=True,
                        help='Cartella AUTO_sim esatta della Full Chain; ripetere per ciascuna variante')
    parser.add_argument('--setup', type=Path, default=ROOT / 'setup.out', help='File setup.out di QUESTA Full Chain')
    parser.add_argument('--max-extra-pulses', type=int, default=500,
                        help='Numero massimo di passi (+1 o -1) da testare rispetto agli impulsi iniziali')
    parser.add_argument('--require-complete-run', action='store_true',
                        help='Richiede 3 varianti CONV per ogni energia richiesta dal setup (3 o 6 piani)')
    args = parser.parse_args()
    if args.max_extra_pulses < 0:
        parser.error('Valori negativi per la ricerca non consentiti')
    setup = args.setup.expanduser().resolve()
    if not setup.is_file():
        fail(f'File setup.out non trovato: {setup}')
    params = parse_setup(setup)
    energies = requested_energies(setup)
    names = args.sim
    if len(set(names)) != len(names):
        fail('Simulazione passata due volte: elenco non valido')
    sims = []
    for name in names:
        if '/' in name or name in ('.', '..') or not SIM_NAME.fullmatch(name):
            fail(f'Nome AUTO_sim non valido: {name}')
        sim = ROOT / name
        if not sim.is_dir() or sim.is_symlink() or sim.resolve().parent != ROOT:
            fail(f'Simulazione inesistente o collegamento simbolico: {sim}')
        if SIM_NAME.fullmatch(name).group(1) not in energies:
            fail(f'La simulazione {name} ha energia non richiesta da setup.out')
        sims.append(sim)
    if args.require_complete_run:
        expected_sims = 3 * len(energies)
        if len(sims) != expected_sims:
            fail(f'Attese {expected_sims} simulazioni (3 varianti per energia), trovate {len(sims)}')
        for energy in energies:
            groups = {}
            for sim in sims:
                m = SIM_NAME.fullmatch(sim.name)
                if m.group(1) == energy:
                    groups.setdefault(m.group(2), set()).add(m.group(3) or 'originale')
            if len(groups) != 1 or next(iter(groups.values())) != {'originale','plus2mm','minus2mm'}:
                fail(f'Per {energy} MeV sono necessari originale, plus2mm e minus2mm dello stesso indice: {groups}')

    required_scripts = {'normalize': SCRIPT_DIR / 'dose_to_prescription.py',
                    'check': SCRIPT_DIR / 'check_isodose_PTV.py',
                    'dvh': SCRIPT_DIR / 'dvh_analyze.py'}
    ptv_src = ROOT / 'imgs' / 'PTV_plan.mhd'
    for key, path in [('PTV', ptv_src), *required_scripts.items()]:
        if not path.is_file():
            fail(f'Manca {key}: {path}')
    # Run id con precisione microsecondo + UUID; mkdir esclusivo.
    run_id = dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '_' + uuid.uuid4().hex[:8]
    output = ROOT / 'pulse_optimization_runs' / ('run_' + run_id)
    output.parent.mkdir(exist_ok=True)
    output.mkdir(exist_ok=False)
    (output / 'piani').mkdir(exist_ok=False)
    inputs = output / 'input_copiati'
    inputs.mkdir(exist_ok=False)
    scripts_dir = inputs / 'scripts'
    scripts_dir.mkdir(exist_ok=False)
    setup_copy = inputs / 'setup.out'
    shutil.copy2(setup, setup_copy)
    copy_scripts = {}
    for key, path in required_scripts.items():
        dst = scripts_dir / path.name
        shutil.copy2(path, dst)
        copy_scripts[key] = dst
    ptv_copy = inputs / 'PTV_plan.mhd'
    write_image_exclusive(sitk.ReadImage(str(ptv_src)), ptv_copy)
    print('\n' + '=' * 85, flush=True)
    print(f'PULSE OPTIMIZATION CONV — SOLO COPIE — {run_id}', flush=True)
    print(f'Setup: Rx={params["rx_gy"]:g} Gy, isodose={params["iso_pct"]:g}%, '
          f'PTV minimo={params["volume_pct"]:g}%; Dmax target={UPPER_PERCENT:g}% Rx', flush=True)
    print(f'Simulazioni correnti: {len(sims)}, piani: {len(sims)} (SOLO CONV; FLASH esclusa)', flush=True)
    print(f'Risultati nuovi: {output}', flush=True)
    # I file originali restano protetti; omettiamo annunci intermedi nel log.
    results, errors = [], []
    for sim in sims:
        for mode in MODES:
            name = f'{sim.name}_{mode}'
            print(f'\n----- {name} -----', flush=True)
            try:
                result = optimize_plan(sim, mode, output, ptv_copy, setup_copy, params,
                                       copy_scripts, args.max_extra_pulses)
                results.append(result)
            except Exception as exc:
                errors.append({'piano': name, 'errore': str(exc)})
                print(f'  ERRORE {name}: {exc}', flush=True)
    columns = ['piano','energia_MeV','modalita','variante','impulsi_full_chain',
               'copertura_full_chain_percent','PTV_full_chain','Dmax_full_chain_Gy',
               'impulsi_minimo_copertura','copertura_minima_percent','Dmax_minima_Gy',
               'delta_impulsi_minimo','impulsi_vicino_107','copertura_vicino_107_percent',
               'Dmax_vicino_107_Gy','delta_impulsi_107','scostamento_107_Gy',
               'PTV_vicino_107','PTV_e_Dmax_conformi','fase2_status','impulsi_testati','mappa_FullChain_sorgente','cartella_risultati']
    with (output / 'confronto_impulsi_CONV_3_6_piani.csv').open('x', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        w.writerows(results)
    with (output / 'errori.csv').open('x', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['piano','errore'])
        w.writeheader()
        w.writerows(errors)
    small_rows = [[r['piano'], r['impulsi_full_chain'],
                   f'{r["copertura_full_chain_percent"]:.3f}', r['PTV_full_chain'],
                   f'{r["Dmax_full_chain_Gy"]:.3f}',
                   r['impulsi_minimo_copertura'],
                   f'{r["copertura_minima_percent"]:.3f}',
                   f'{r["Dmax_minima_Gy"]:.3f}',
                   r['PTV_e_Dmax_conformi'], r['impulsi_testati']]
                  for r in results]
    shown = header_table(['PIANO','FC p0','FC cop%','PTV FC','FC Dmax Gy',
                          'p MIN','Min cop%','Min Dmax Gy','2 VINCOLI','Test'],
                         small_rows) if results else 'NESSUN PIANO ANALIZZATO'
    summary = ['PULSE NUMBER OPTIMIZER — SOLO CONV (3 O 6 PIANI)',
               f'Result directory: {output}',
               f'Prescrizione: {params["rx_gy"]:g} Gy; isodose {params["iso_pct"]:g}%; '
               f'volume PTV >= {params["volume_pct"]:g}%; Dmax PTV <=107% Rx.',
               'FASE 1: +1 o -1 impulso dal piano Full Chain per cercare il minimo con copertura nominale.',
               'FASE 2: verifica della Dmax PTV sul minimo: se supera 107%, segnala incompatibilita.',
               'Una tabella dvh_analyze.py per OGNI impulso CONV provato (PTV ed OAR).',
               'Le tabelle per piano sono in piani/<PIANO>/tabella_piano.txt e tutti_impulsi.csv.',
               'Mappe e isodosi solo in copie isolate dentro la cartella nuova di run.',
               'Originali e .raw originali: MAI sovrascritti da questo programma.',
               '', shown, '',
               f'Piani analizzati: {len(results)}/{len(sims)}',
               f'Piani conformi a PTV + Dmax: {sum(r["PTV_e_Dmax_conformi"] == "YES" for r in results)}',
               f'Errori tecnici: {len(errors)}',
               'Modalita ottimizzata: CONV. FLASH mai analizzata o modificata.',
               'Modello di dose: riscalamento lineare, senza nuove simulazioni Monte Carlo.']
    for error in errors:
        summary.append('ERRORE: ' + error['piano'] + ': ' + error['errore'])
    with (output / 'CONFRONTO_FINALE_IMPULSI.txt').open('x', encoding='utf-8') as f:
        f.write('\n'.join(summary) + '\n')
    with (output / 'parametri_e_sorgenti.json').open('x', encoding='utf-8') as f:
        json.dump({'params': params, 'energies': sorted(energies),
                   'simulations': [s.name for s in sims], 'errors': errors,
                   'source_fullchain_setup': str(setup),
                   'source_scripts': {k: str(v) for k, v in required_scripts.items()},
                   'target_Dmax_percent': UPPER_PERCENT,
                   'selection_rule': 'minimum_CONV_pulses_with_nominal_PTV_coverage_and_Dmax_upper_check',
                   'optimized_modes': ['CONV'],
                   'dose_model': 'linear_pulse_rescale'}, f, ensure_ascii=False, indent=2)
    # FIX SOLO PRESENTAZIONE: tabella conclusiva IDENTICA, senza il testo
    # descrittivo intermedio. Il report TXT continua a contenere ogni dettaglio.
    print('\n' + '=' * 85, flush=True)
    print('FINAL COMPARISON — CONV ONLY', flush=True)
    print(shown, flush=True)
    print(f'Analyzed: {len(results)}/{len(sims)} | PTV+Dmax compliant: '
          f'{sum(r["PTV_e_Dmax_conformi"] == "YES" for r in results)} | '
          f'Errors: {len(errors)}', flush=True)
    print(f'Results: {output}', flush=True)
    # Legenda aggiunta SOLO al log (non cambia nomi o contenuto delle tabelle).
    legend = [
        '', 'TABLE LEGEND',
        'P / IMPULSI = numero di impulsi CONV.',
        'Energy / Mode / Plan / Structure = energia, modalita, tipo di piano, ROI.',
        'Max tab [Gy] = dose massima riportata nella tabella DVH.',
        'Est Dmax [Gy] = D0.1% del DVH; NON e la Dmax effettiva del PTV.',
        'Wtd mean [Gy] = dose media pesata sul volume della ROI.',
        f'V{params["iso_pct"]:g}% DVH [%] = percentuale PTV che riceve almeno '
        f'{params["iso_pct"]:g}% della dose prescritta.',
        f'D{params["volume_pct"]:g} [Gy] = dose ricevuta da almeno '
        f'{params["volume_pct"]:g}% del PTV.',
        'PTV 3D [%] = copertura nominale calcolata sui voxel del PTV.',
        'Dmax PTV [Gy] = MASSIMO nei voxel del PTV; controllato rispetto al 107%.',
        'PTV DVH = MET/NOT MET secondo il DVH; PTV 3D = YES/NO da mappa 3D.',
        '107% = OK o >107% nelle tabelle dei test; 2 VINCOLI = YES/NO nel confronto finale.',
        'FC p0 = impulsi iniziali della Full Chain; FC cop% = copertura iniziale.',
        'p MIN = minimo che copre il PTV; Min cop% e Min Dmax Gy = relativi valori.',
        'Test = numero di valori di impulsi provati.',
        f'Vincoli: copertura >= {params["volume_pct"]:g}% della dose-isodose '
        f'{params["iso_pct"]:g}%; Dmax PTV <= {params["rx_gy"]*UPPER_PERCENT/100:.3f} Gy.',
        'Dose riscalata linearmente; nessuna nuova simulazione Monte Carlo.'
    ]
    print('\n'.join(legend), flush=True)
    # La legenda viene salvata anche accanto al riepilogo del run.
    with (output / 'LEGENDA_TABELLE.txt').open('x', encoding='utf-8') as f:
        f.write('\n'.join(legend).lstrip('\n') + '\n')
    with (output / 'CONFRONTO_FINALE_IMPULSI.txt').open('a', encoding='utf-8') as f:
        f.write('\n'.join(legend) + '\n')
    return 2 if errors else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'ERRORE GENERALE pulse_number_optimizer: {exc}', file=sys.stderr)
        sys.exit(1)

