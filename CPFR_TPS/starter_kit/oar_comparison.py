#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""OAR Dmax and D(v) for only the exact AUTO_sim folders created in one run.

Required: --sim SIM (repeated), --ptv-table THIS_RUN_CSV.
Dynamic OAR masks from --oar-plan; call the real dvh_analyze.py for every plan.
Save individual tables only in results/<SIM>_<MODE> for PTV YES.
Originals and previous runs are preserved; reports are in English.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

LEVELS = (0.1, 1.0, 5.0, 10.0, 20.0, 50.0, 80.0, 95.0)
EPS_GY = 1e-6
SIM_RX = re.compile(r'^AUTO_sim(7|9)MeV_(\d+)(?:_(plus2mm|minus2mm))?$', re.I)
MODES = ('CONV', 'FLASH')


def csv_out(path, columns, rows):
    with path.open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def geometry_equal(a, b):
    return (a.GetSize() == b.GetSize()
            and np.allclose(a.GetSpacing(), b.GetSpacing(), rtol=0, atol=1e-5)
            and np.allclose(a.GetOrigin(), b.GetOrigin(), rtol=0, atol=1e-5)
            and np.allclose(a.GetDirection(), b.GetDirection(), rtol=0, atol=1e-6))


def cumulative_curve(path):
    """Read native ComputeDVH text: numeric dose in cGy, cumulative volume %."""
    data = np.atleast_2d(np.loadtxt(path, dtype=np.float64, comments='#'))
    if data.shape[1] != 2 or len(data) < 2 or not np.isfinite(data).all():
        raise ValueError(f'Invalid cumulative DVH (two numeric columns expected): {path}')
    doses, volumes = data[:, 0] / 100.0, data[:, 1]
    if (np.min(doses) < 0 or np.max(doses) > 10000 or
            np.min(volumes) < -0.05 or np.max(volumes) > 100.05):
        raise ValueError(f'Unreasonable DVH dose/volume units: {path}')
    order = np.argsort(doses)
    doses, volumes = doses[order], volumes[order]
    if np.max(np.diff(volumes)) > 1:
        raise ValueError(f'DVH volume is not cumulative/monotone: {path}')
    volumes = np.minimum.accumulate(np.clip(volumes, 0, 100))
    return doses, volumes


def dose_at_volume(sorted_doses, volume_percent):
    """D(v) in Gy: minimum dose received by AT LEAST v% of OAR voxels.

    In a discrete voxel population, select the lowest of the highest
    ceil(v*n/100) voxels. This avoids interpolation/binning artifacts.
    """
    n = len(sorted_doses)
    count = max(1, min(n, math.ceil(volume_percent * n / 100.0)))
    return float(sorted_doses[n - count])


def oar_paths(root, explicitly_given):
    if explicitly_given:
        selected = [Path(p) for p in explicitly_given]
    else:
        # input_log.out is written by the Full Chain from the actual ROI
        # arguments, not from a fixed patient-specific set of organs.
        input_log = root / 'input_log.out'
        if not input_log.is_file():
            raise RuntimeError('No OAR inputs: supply --oar-plan repeatedly or provide input_log.out')
        roi_lines = [line.strip()[4:].strip() for line in input_log.read_text(
            encoding='utf-8', errors='replace').splitlines()
                     if line.lstrip().startswith('ROI:')]
        if len(roi_lines) != 1:
            raise RuntimeError('input_log.out must contain exactly one ROI: line; '
                               'or supply --oar-plan repeatedly')
        selected = [Path('imgs') / (name + '_plan.mhd') for name in roi_lines[0].split()]
    if not selected:
        raise RuntimeError('No OARs specified; cannot compare organ sparing')
    roi_files = []
    for path in selected:
        path = path if path.is_absolute() else root / path
        path = path.resolve()
        if not path.is_file() or not path.stem.lower().endswith('_plan'):
            raise RuntimeError(f'Expected OAR_plan.mhd mask not found: {path}')
        if path.stem.lower() in ('ptv_plan', 'ct_plan', 'marker_plan'):
            raise RuntimeError(f'Non-OAR mask was included: {path}')
        if path in roi_files:
            raise RuntimeError(f'Duplicate OAR mask: {path}')
        roi_files.append(path)
    return roi_files


def find_dvh_analyzer(root):
    for base in (root, root / 'starter_kit'):
        for name in ('dvh_analyze.py', 'DVH_analyze.py'):
            candidate = base / name
            if candidate.is_file():
                return candidate
    raise FileNotFoundError('dvh_analyze.py not found in patient root or starter_kit')


def run_native_analyzer(script, files, out, root, rx, rx_level_percent):
    native_dir = out / 'dvh_analyze'
    native_dir.mkdir(exist_ok=False)
    command = [
        sys.executable,
        str(script),
        *map(str, files),
        '--prescription-gy', str(rx),
        '--prescription-level-percent', str(rx_level_percent),
        '--plot',
        '--output-dir', str(native_dir),
    ]
    log = out / 'dvh_analyze.log'
    with log.open('x', encoding='utf-8') as stream:
        stream.write('COMMAND: ' + ' '.join(command) + '\n\n')
        stream.flush()
        done = subprocess.run(command, cwd=root, stdout=stream,
                              stderr=subprocess.STDOUT,
                              env={**os.environ, 'MPLBACKEND': 'Agg'}, check=False)
    if done.returncode:
        raise RuntimeError(f'dvh_analyze.py failed (exit {done.returncode}): {log}')
    return log


def ptv_statuses(current_ptv_table, expected_names):
    """Only the PTV report of THIS Full Chain invocation is authoritative."""
    if not current_ptv_table.is_file():
        raise RuntimeError(f'Current-run PTV coverage table not found: {current_ptv_table}')
    with current_ptv_table.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        required = {'simulation', 'mode', 'nominal_result'}
        if not required <= set(reader.fieldnames or []):
            raise RuntimeError(f'Unexpected current PTV table columns: {current_ptv_table}')
        rows = list(reader)
    statuses = {}
    for row in rows:
        key = f'{row["simulation"]}_{row["mode"]}'
        if key in statuses:
            raise RuntimeError(f'Duplicate PTV entry: {key}')
        if row['nominal_result'] not in ('YES', 'NO'):
            raise RuntimeError(f'Invalid or incomplete PTV result for {key}: {row["nominal_result"]}')
        statuses[key] = row['nominal_result']
    if set(statuses) != set(expected_names):
        raise RuntimeError('Current-run PTV table differs from current-run AUTO_sim list: '
                           f'missing={sorted(set(expected_names)-set(statuses))}, '
                           f'unexpected={sorted(set(statuses)-set(expected_names))}')
    return statuses, str(current_ptv_table)


def selected_simulations(root, supplied):
    """Never scan pre-existing AUTO_sim folders; require an explicit run manifest."""
    if not supplied:
        raise RuntimeError('Supply at least one --sim from the current Full Chain invocation')
    sims, seen = [], set()
    for name in supplied:
        if not SIM_RX.fullmatch(name) or name in seen:
            raise RuntimeError(f'Invalid or duplicate current-run AUTO_sim: {name}')
        seen.add(name)
        candidate = root / name
        if not candidate.is_dir() or candidate.is_symlink() or candidate.resolve().parent != root:
            raise RuntimeError(f'Current-run AUTO_sim directory missing or unsafe: {candidate}')
        sims.append(candidate)
    return sims


def inspect_plan(sim, energy, mode, root, output, rois, roi_images, sitk, analyzer, statuses, rx, rx_level_percent):
    name = f'{sim.name}_{mode}'
    target = output / 'per_plan' / name
    target.mkdir(parents=True, exist_ok=False)
    dvh_dir = sim / f'DVH{energy}MeV_{mode}'
    doses = sorted(sim.glob(f'DOSE_{energy}MeV_*pulses_{mode}_BEV.mhd'))
    if len(doses) != 1:
        raise RuntimeError(f'{name}: expected exactly one absolute Gy BEV map, found {len(doses)}')
    if not dvh_dir.is_dir():
        raise RuntimeError(f'{name}: missing numeric DVH folder: {dvh_dir}')
    dose_img = sitk.ReadImage(str(doses[0]))
    dose_data = sitk.GetArrayFromImage(dose_img)
    selected_dvhs = []
    for roi in rois:
        numeric = dvh_dir / f'{roi.stem}{energy}MeV_{mode}.txt'
        if not numeric.is_file():
            raise RuntimeError(f'{name}: missing OAR DVH file: {numeric}')
        selected_dvhs.append(numeric)
    # Include the native PTV curve in the dvh_analyze output for context.
    ptv_dvh = dvh_dir / f'PTV_plan{energy}MeV_{mode}.txt'
    if ptv_dvh.is_file():
        selected_dvhs.insert(0, ptv_dvh)
    analyzer_log = run_native_analyzer(analyzer, selected_dvhs, target, root, rx, rx_level_percent)
    status = statuses.get(name, 'UNKNOWN')
    basic = dict(plan=name, simulation=sim.name, energy_MeV=int(energy), mode=mode,
                 aperture=('plus2mm' if sim.name.endswith('_plus2mm') else
                           'minus2mm' if sim.name.endswith('_minus2mm') else 'nominal'),
                 PTV_result=status, absolute_dose_map=str(doses[0]),
                 dvh_analyze_log=str(analyzer_log))
    rows, curves = [], {}
    for roi, numeric in zip(rois, selected_dvhs[-len(rois):]):
        organ = roi.stem[:-5]
        if roi not in roi_images:
            roi_images[roi] = sitk.ReadImage(str(roi))
        mask_img = roi_images[roi]
        if not geometry_equal(dose_img, mask_img):
            raise RuntimeError(f'{name}: ROI {organ} and BEV dose grids differ; refusing implicit resampling')
        mask = sitk.GetArrayFromImage(mask_img) > 0
        values = np.asarray(dose_data[mask], dtype=np.float64)
        if len(values) == 0 or not np.isfinite(values).all():
            raise RuntimeError(f'{name}: empty OAR / nonfinite OAR dose: {organ}')
        sorted_values = np.sort(values)
        # Read native numerical curves: plotted as originally computed, not digitised.
        x, y = cumulative_curve(numeric)
        curves[organ] = (x, y)
        for fraction in LEVELS:
            rows.append({**basic, 'OAR': organ,
                         'Dmax_Gy': float(sorted_values[-1]),
                         'volume_percent': fraction,
                         'D_v_Gy': dose_at_volume(sorted_values, fraction),
                         'OAR_voxel_count': len(values),
                         'numeric_DVH': str(numeric)})
    csv_out(target / 'OAR_metrics.csv', list(rows[0]), rows)
    (target / 'metrics_definition.txt').write_text(
        f'PLAN: {name}\nPTV result: {status}\n'
        'Dmax_Gy: maximum OAR voxel dose from ABSOLUTE-Gy BEV dose map.\n'
        'D_v_Gy: dose received by at least v% of OAR voxels; derived from BEV map and OAR_plan.\n'
        'Plots: original numeric ComputeDVH curves; dvh_analyze.py reports are also saved.\n'
        'The raw-map Dmax may differ from Dmax estimates reported by dvh_analyze.py.\n',
        encoding='utf-8')
    return basic, rows, curves


def minima(rows, metric, group_fields):
    buckets = defaultdict(list)
    for row in rows:
        key = tuple(row[col] for col in group_fields)
        buckets[key].append(row)
    result = []
    for group in sorted(buckets):
        rs = buckets[group]
        low = min(r[metric] for r in rs)
        winners = sorted(set(r['plan'] for r in rs if abs(r[metric] - low) <= EPS_GY))
        result.append({**dict(zip(group_fields, group)),
                       'minimum_dose_Gy': low,
                       'best_plan_s': '; '.join(winners),
                       'compared_plans': len(set(r['plan'] for r in rs))})
    return result


def is_nondominated(vectors):
    """True Pareto front: different organs and D(v) curves can trade off."""
    out = {}
    for name, vec in vectors.items():
        dominated = False
        for other, vals in vectors.items():
            if name == other:
                continue
            if all(v <= w + EPS_GY for v, w in zip(vals, vec)) and any(
                    v < w - EPS_GY for v, w in zip(vals, vec)):
                dominated = True
                break
        out[name] = not dominated
    return out


def csv_rankings(output, prefix, rows):
    """All plans and PTV-YES only; each split into global and by-energy tables."""
    dmax_rows = [r for r in rows if r['volume_percent'] == LEVELS[0]]
    for table, selected, fields, metric in (
        ('Best_Dmax_per_OAR', dmax_rows, ['OAR'], 'Dmax_Gy'),
        ('Best_Dv_per_OAR_and_volume', rows, ['OAR', 'volume_percent'], 'D_v_Gy'),
        ('Best_Dmax_per_energy_and_OAR', dmax_rows, ['energy_MeV', 'OAR'], 'Dmax_Gy'),
        ('Best_Dv_per_energy_OAR_and_volume', rows,
         ['energy_MeV', 'OAR', 'volume_percent'], 'D_v_Gy'),
    ):
        winner_rows = minima(selected, metric, fields)
        csv_out(output / (prefix + table + '.csv'),
                fields + ['minimum_dose_Gy', 'best_plan_s', 'compared_plans'], winner_rows)
    return (minima(dmax_rows, 'Dmax_Gy', ['OAR']),
            minima(rows, 'D_v_Gy', ['OAR', 'volume_percent']))


def main(argv=None):
    parser = argparse.ArgumentParser(description='OAR Dmax / equal-volume dose comparison: only the CURRENT Full Chain run')
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--sim', action='append', required=True, help='Exact current-run AUTO_sim directory; repeat as needed')
    parser.add_argument('--ptv-table', type=Path, required=True, help='Exact PTV_coverage_table.csv created by this run')
    parser.add_argument('--oar-plan', action='append', default=[], metavar='imgs/Organ_plan.mhd',
                        help='May be repeated; Full Chain passes its live ROIs_plan array')
    parser.add_argument(
        '--prescription-gy',
        type=float,
        required=True,
        help='Prescription dose in Gy',
    )
    parser.add_argument(
        '--prescription-level-percent',
        type=float,
        required=True,
        help='Prescription isodose level in percent',
    )
    args = parser.parse_args(argv)
    root = args.root.expanduser().resolve()
    try:
        import SimpleITK as sitk
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError('Missing SimpleITK/numpy/matplotlib in planning environment') from exc
    rois = oar_paths(root, args.oar_plan)
    analyzer = find_dvh_analyzer(root)
    sims = selected_simulations(root, args.sim)
    rx = args.prescription_gy
    rx_level_percent = args.prescription_level_percent
    current_names = {f'{sim.name}_{mode}' for sim in sims for mode in MODES}
    ptv_table = args.ptv_table if args.ptv_table.is_absolute() else root / args.ptv_table
    ptv_status, ptv_source = ptv_statuses(ptv_table.resolve(), current_names)
    stamp = dt.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    output = root / 'oar_analysis' / ('run_' + stamp)
    output.mkdir(parents=True, exist_ok=False)
    print('=== OAR ANALYSIS: CURRENT FULL CHAIN RUN ONLY (7/9 MeV, CONV/FLASH) ===', flush=True)
    print(f'Root: {root}\nOutput: {output}\nOARs: {", ".join(p.stem[:-5] for p in rois)}', flush=True)
    print(f'OARs sourced from: {"Full Chain ROIs_plan arguments" if args.oar_plan else "input_log.out"}', flush=True)
    print(f'PTV status source: {ptv_source}', flush=True)
    print(f'dvh_analyze.py: {analyzer}', flush=True)
    details, observations, curves_by_oar, errors = {}, [], defaultdict(dict), []
    roi_images = {}
    for sim in sims:
        energy = SIM_RX.fullmatch(sim.name).group(1)
        for mode in MODES:
            plan = f'{sim.name}_{mode}'
            try:
                info, rows, curves = inspect_plan(sim, energy, mode, root, output,
                                                   rois, roi_images, sitk, analyzer, ptv_status, rx, rx_level_percent)
                details[plan] = info
                observations.extend(rows)
                for organ, curve in curves.items():
                    curves_by_oar[organ][plan] = curve
                print(f'  OK: {plan} (PTV {info["PTV_result"]})', flush=True)
            except Exception as exc:
                errors.append({'plan': plan, 'error': str(exc)})
                print(f'  ERROR: {plan}: {exc}', flush=True)
    csv_out(output / 'analysis_errors.csv', ['plan', 'error'], errors)
    if not observations:
        raise RuntimeError(f'No successful OAR analysis; see {output / "analysis_errors.csv"}')
    csv_out(output / 'All_OAR_measurements.csv', list(observations[0]), observations)
    maxima = [r for r in observations if r['volume_percent'] == LEVELS[0]]
    csv_out(output / 'Dmax_all_plans.csv', list(maxima[0]), maxima)
    csv_out(output / 'Dose_at_equal_OAR_volume.csv', list(observations[0]), observations)
    # Comparison matrices include ONLY the explicitly listed current-run AUTO_sim.
    names = sorted(details)
    organ_names = sorted({r['OAR'] for r in observations})
    dmax_index = {(r['OAR'], r['plan']): r['Dmax_Gy'] for r in maxima}
    csv_out(output / 'Dmax_comparison_matrix.csv', ['OAR', *names], [
        {'OAR': organ, **{name: dmax_index.get((organ, name), '') for name in names}}
        for organ in organ_names])
    dv_index = {(r['OAR'], r['volume_percent'], r['plan']): r['D_v_Gy']
                for r in observations}
    csv_out(output / 'Dv_comparison_matrix.csv', ['OAR', 'volume_percent', *names], [
        {'OAR': organ, 'volume_percent': v,
         **{name: dv_index.get((organ, v, name), '') for name in names}}
        for organ in organ_names for v in LEVELS])

    all_dmax, all_dv = csv_rankings(output, 'ALL_', observations)
    qualified = [r for r in observations if r['PTV_result'] == 'YES']
    if qualified:
        yes_dmax, yes_dv = csv_rankings(output, 'PTV_YES_', qualified)
    else:
        yes_dmax, yes_dv = [], []
    # Non-dominated plans are not an asserted single clinical optimum.
    pareto_rows = []
    for label, subset in [('ALL', observations), ('PTV_YES', qualified)]:
        pool = sorted({r['plan'] for r in subset})
        dmax_vec = {}
        dv_vec = {}
        for plan in pool:
            p = [r for r in subset if r['plan'] == plan]
            dmax_by_oar = {r['OAR']: r['Dmax_Gy'] for r in p if r['volume_percent'] == LEVELS[0]}
            dv_by_oar = {(r['OAR'], r['volume_percent']): r['D_v_Gy'] for r in p}
            if len(dmax_by_oar) == len(organ_names):
                dmax_vec[plan] = [dmax_by_oar[o] for o in organ_names]
            if len(dv_by_oar) == len(organ_names) * len(LEVELS):
                dv_vec[plan] = [dv_by_oar[o, v] for o in organ_names for v in LEVELS]
        max_p = is_nondominated(dmax_vec)
        dv_p = is_nondominated(dv_vec)
        for plan in pool:
            pareto_rows.append({'comparison_group': label, 'plan': plan,
                                'PTV_status': details[plan]['PTV_result'],
                                'Pareto_for_OAR_Dmax': max_p.get(plan, ''),
                                'Pareto_for_OAR_Dv': dv_p.get(plan, '')})
    csv_out(output / 'OAR_Pareto_tradeoffs.csv',
            ['comparison_group', 'plan', 'PTV_status', 'Pareto_for_OAR_Dmax', 'Pareto_for_OAR_Dv'], pareto_rows)

    plots = output / 'DVH_comparison_plots'
    plots.mkdir(exist_ok=False)
    for organ, source in sorted(curves_by_oar.items()):
        # Show 7 and 9 MeV together for a given organ, only variants generated in this invocation.
        fig, ax = plt.subplots(figsize=(12, 7))
        for plan, (x, y) in sorted(source.items()):
            ax.plot(x, y, linewidth=1.5,
                    label=f'{plan} (PTV {details[plan]["PTV_result"]})')
        ax.set(xlabel='Dose (Gy)', ylabel='Volume receiving at least this dose (%)',
               title=f'Cumulative DVH comparison — {organ} (7 and 9 MeV)',
               ylim=(-1, 101))
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, bbox_to_anchor=(1.01, 1), loc='upper left')
        fig.tight_layout()
        safe = re.sub(r'[^a-zA-Z0-9_-]', '_', organ)
        fig.savefig(plots / f'DVH_comparison_{safe}_7and9MeV.png', dpi=150, bbox_inches='tight')
        plt.close(fig)

    summary = [
        'CURRENT FULL CHAIN RUN — 7/9 MeV OAR DOSE COMPARISON',
        f'Successfully analyzed plans: {len(details)}; failed: {len(errors)}',
        f'OAR list (from Full Chain): {", ".join(organ_names)}',
        f'PTV source: {ptv_source}',
        'Dmax = actual maximum Gy in OAR_plan voxel mask on absolute-dose BEV map.',
        'D(v) = minimum Gy received by at least v% of the OAR; identical volume fractions in every plan.',
        'Native dvh_analyze.py outputs and the numeric DVH curves are saved per plan.',
        'Current-run PTV NO plans appear in ALL comparisons, but are excluded from PTV_YES recommendations.',
        'Neither OAR prescription limits nor a single clinical winner are assumed.',
        '', 'LOWEST DMAX PER OAR — ALL PLANS (7 + 9 MeV)',
    ]
    for r in all_dmax:
        summary.append(f'  {r["OAR"]}: {r["minimum_dose_Gy"]:.5f} Gy | {r["best_plan_s"]}')
    summary.extend(['', 'LOWEST DOSE AT EQUAL OAR VOLUME — ALL PLANS (7 + 9 MeV)'])
    for r in all_dv:
        summary.append(f'  {r["OAR"]} D{r["volume_percent"]:g}%: '
                       f'{r["minimum_dose_Gy"]:.5f} Gy | {r["best_plan_s"]}')
    summary.extend(['', 'LOWEST DMAX PER OAR — PTV YES PLANS ONLY'])
    for r in yes_dmax:
        summary.append(f'  {r["OAR"]}: {r["minimum_dose_Gy"]:.5f} Gy | {r["best_plan_s"]}')
    if not yes_dmax:
        summary.append('  No PTV-YES eligible plans with complete OAR metrics')
    summary.extend(['', 'LOWEST DOSE AT EQUAL OAR VOLUME — PTV YES PLANS ONLY'])
    for r in yes_dv:
        summary.append(f'  {r["OAR"]} D{r["volume_percent"]:g}%: '
                       f'{r["minimum_dose_Gy"]:.5f} Gy | {r["best_plan_s"]}')
    if not yes_dv:
        summary.append('  No PTV-YES eligible plans with complete OAR metrics')
    summary.extend(['', 'For energy-specific winners see Best_*_per_energy_*.csv.',
                    'Pareto summaries allow multiple optimal trade-off plans.',
                    'A plan with the lowest Dmax need not also have the lowest D(v).',
                    'All source AUTO_sim folders are preserved. Only current-run YES results receive OAR reports.',
                    f'Output directory: {output}'])
    # Human-readable Dmax/D(v) tables: current-run plans only.
    # The report is placed in each current-run PTV=YES results/<AUTO_sim>_<MODE>/.
    # Existing runs are neither opened for writing nor re-analyzed.
    fractions = LEVELS
    fmt = lambda x: f'{float(x):.5f}'
    by_organ = defaultdict(dict)
    for r in observations:
        item = by_organ[r['OAR']].setdefault(r['plan'],
              {'Dmax': r['Dmax_Gy'], 'PTV': r['PTV_result'], 'Dv': {}})
        item['Dv'][r['volume_percent']] = r['D_v_Gy']
    table = [
        'OAR DOSE COMPARISON — CURRENT FULL CHAIN RUN ONLY',
        f'Current-run AUTO_sim folders: {", ".join(p.name for p in sims)}',
        f'Compared OARs (patient-dependent): {", ".join(organ_names)}',
        'Dmax: highest absolute Gy within each OAR_plan mask.',
        'D(v): Gy delivered to at least v% of the SAME OAR volume in each plan.',
        'PTV=NO plans are shown for reference; winners eligible for selection must have PTV=YES.',
    ]
    for organ in organ_names:
        table.extend(['', '=' * 85, f'OAR: {organ}', '=' * 85])
        head = ['PLAN', 'PTV', 'Dmax Gy', *[f'D{v:g}% Gy' for v in fractions]]
        widths = [max(35, max(len(n) for n in names)), 5, 11] + [12] * len(fractions)
        def format_line(items):
            return ' | '.join(f'{str(v):<{width}}' for v, width in zip(items, widths))
        table.append(format_line(head))
        table.append('-' * len(table[-1]))
        for plan in names:
            data = by_organ[organ].get(plan)
            if data is None:
                continue
            table.append(format_line([plan, data['PTV'], fmt(data['Dmax']),
                                     *[fmt(data['Dv'][v]) for v in fractions]]))
    table.extend(['', 'LOWEST DMAX PER OAR — CURRENT RUN (ALL PLANS)'])
    for r in all_dmax:
        table.append(f'  {r["OAR"]}: {r["minimum_dose_Gy"]:.5f} Gy — {r["best_plan_s"]}')
    table.extend(['', 'LOWEST D(V) AT EQUAL OAR VOLUME — CURRENT RUN (ALL PLANS)'])
    for r in all_dv:
        table.append(f'  {r["OAR"]} D{r["volume_percent"]:g}%: '
                     f'{r["minimum_dose_Gy"]:.5f} Gy — {r["best_plan_s"]}')
    table.extend(['', 'LOWEST DMAX PER OAR — PTV-YES PLANS ONLY'])
    for r in yes_dmax:
        table.append(f'  {r["OAR"]}: {r["minimum_dose_Gy"]:.5f} Gy — {r["best_plan_s"]}')
    table.extend(['', 'LOWEST D(V) AT EQUAL OAR VOLUME — PTV-YES PLANS ONLY'])
    for r in yes_dv:
        table.append(f'  {r["OAR"]} D{r["volume_percent"]:g}%: '
                     f'{r["minimum_dose_Gy"]:.5f} Gy — {r["best_plan_s"]}')
    if not yes_dmax:
        table.append('  No PTV=YES plans available.')
    table.extend(['', 'NOTE: Dmax and D(v) minima may occur in different plans.',
                  'The comparison does not establish a single clinical winner.',
                  'Unanalyzable plans are excluded from minima; see analysis_errors.csv.'])
    common_report = '\n'.join(table) + '\n'
    with (output / 'OAR_comparison_table.txt').open('x', encoding='utf-8') as stream:
        stream.write(common_report)

    # Only current-run PTV-YES plans get an individual table_oar.txt.
    # Never overwrite already existing results, including prior Full Chain runs.
    import shutil
    for plan, info in sorted(details.items()):
        local = output / 'per_plan' / plan
        individual = [f'CURRENT PLAN: {plan}', f'PTV result: {info["PTV_result"]}',
                      '', 'PLAN-SPECIFIC OAR METRICS']
        for organ in organ_names:
            measures = by_organ[organ].get(plan)
            if measures is None:
                continue
            individual.append(f'{organ}: Dmax={fmt(measures["Dmax"])} Gy; ' +
                              '; '.join(f'D{v:g}%={fmt(measures["Dv"][v])} Gy' for v in fractions))
        individual.extend(['', common_report])
        # No per-plan table is created for PTV NO, and original AUTO_sim
        # directories are not modified. Only compliant plans get table_oar.txt.
        if info['PTV_result'] == 'YES':
            destination = root / 'results' / plan
            if not destination.is_dir() or destination.is_symlink():
                raise RuntimeError(f'Current-run PTV-YES destination missing: {destination}')
            with (destination / 'table_oar.txt').open('x', encoding='utf-8') as stream:
                stream.write('\n'.join(individual))
            # Keep native dvh_analyze.py outputs alongside the winning plan.
            # Per-plan native and numerical results remain available even for NO
            # in the run-specific oar_analysis/per_plan folder.
            (destination / 'OAR_metrics.csv').write_bytes((local / 'OAR_metrics.csv').read_bytes())
            shutil.copytree(local / 'dvh_analyze', destination / 'dvh_analyze', dirs_exist_ok=False)
            shutil.copy2(local / 'dvh_analyze.log', destination / 'dvh_analyze.log')
            print(f'  OAR TABLE SAVED: {destination / "table_oar.txt"}', flush=True)

    (output / 'OAR_comparison_summary.txt').write_text('\n'.join(summary) + '\n', encoding='utf-8')
    # Include the FULL current-run OAR comparison table in full_log.txt.
    # This table also reports all-plans and PTV-YES-only Dmax / D(v) minima.
    print('\n' + common_report, flush=True)
    if errors:
        print(f'WARNING: {len(errors)} plan(s) could not be assessed. '
              'See analysis_errors.csv; those plans were excluded from minima.', flush=True)
        return 2
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'ERROR IN OAR COMPARISON: {exc}', file=sys.stderr)
        sys.exit(1)
