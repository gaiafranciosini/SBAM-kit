from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import SimpleITK as sitk

ROOT = Path.cwd().resolve()
SETUP = Path(sys.argv[1]).expanduser()
REPORT_DIR = Path(sys.argv[2])
CURRENT_SIM_NAMES = sys.argv[3:]
if not SETUP.is_absolute():
    SETUP = ROOT / SETUP
SETUP = SETUP.resolve()
NORMALIZE = ROOT / 'starter_kit' / 'dose_to_prescription.py'
CHECK = ROOT / 'starter_kit' / 'check_isodose_PTV.py'
PTV = ROOT / 'imgs' / 'PTV_plan.mhd'


def read_setup():
    parameters = {
        'PRESCRIPTION DOSE': 'rx',
        'PERCENTAGE PRESCRIPTION VOLUME': 'volume',
        'PERCENTAGE PRESCRIPTION DOSE': 'isodose',
    }
    values = {}
    if not SETUP.is_file():
        raise RuntimeError(f'Setup file not found: {SETUP}')
    for line_num, row in enumerate(SETUP.read_text(encoding='utf-8').splitlines(), 1):
        if ':' not in row:
            continue
        key, raw = row.split(':', 1)
        field = parameters.get(key.strip().upper())
        if field is None:
            continue
        if field in values:
            raise ValueError(f'{SETUP}:{line_num}: duplicate setup field {key}')
        raw = raw.split('#', 1)[0].strip().replace(',', '.')
        if not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', raw):
            raise ValueError(f'{SETUP}:{line_num}: invalid numeric value: {key}={raw!r}')
        value = float(raw)
        if not math.isfinite(value):
            raise ValueError(f'{SETUP}:{line_num}: non-finite value')
        values[field] = value
    if set(parameters.values()) - values.keys():
        raise ValueError('Missing prescription parameters in setup.out')
    if values['rx'] <= 0 or not 0 <= values['volume'] <= 100 or not 0 <= values['isodose'] <= 100:
        raise ValueError(f'Invalid prescription parameters: {values}')
    # Technical parameter needed by the normalization/check tools.
    # It does NOT change nominal coverage or the YES/NO decision.
    values['tolerance_tool'] = 0.05
    return values


def same_geometry(a, b):
    return (a.GetSize() == b.GetSize()
            and np.allclose(a.GetSpacing(), b.GetSpacing(), atol=1e-6)
            and np.allclose(a.GetOrigin(), b.GetOrigin(), atol=1e-6)
            and np.allclose(a.GetDirection(), b.GetDirection(), atol=1e-6))


def run_and_log(command, log):
    with log.open('x', encoding='utf-8') as stream:
        stream.write('COMMAND: ' + ' '.join(map(str, command)) + '\n\n')
        stream.flush()
        done = subprocess.run([str(part) for part in command], cwd=ROOT,
                              stdout=stream, stderr=subprocess.STDOUT, check=False)
    if done.returncode:
        raise RuntimeError(f'Command failed ({done.returncode}): {log}')


def normalize_rx(sim, dose, rx, tolerance_tool, outdir, tag):
    old = {p.name for p in sim.glob('isodose_prescrizione_*') if p.is_dir()}
    run_and_log([sys.executable, NORMALIZE, '--sim-dir', sim, '--input', dose.name,
                 '--setup', SETUP, '--ptv', PTV, '--dose-prescription-gy', f'{rx:.12g}',
                 '--tolerance-gy', f'{tolerance_tool:.12g}'],
                outdir / f'{tag}_normalize.log')
    new = sorted(p for p in sim.glob('isodose_prescrizione_*')
                 if p.is_dir() and p.name not in old)
    if len(new) != 1:
        raise RuntimeError(f'Ambiguous prescription normalization output for {dose}: {new}')
    folder = new[0]
    report = folder / 'report_normalizzazione.json'
    if not report.is_file():
        raise RuntimeError(f'Normalization report missing: {report}')
    data = json.loads(report.read_text(encoding='utf-8'))
    if abs(float(data['dose_prescribed_gy']) - rx) > 1e-8:
        raise RuntimeError('Prescription dose in normalization report does not match setup.out')
    maps = data.get('maps', [])
    if len(maps) != 1:
        raise RuntimeError('Expected exactly one normalized prescription map')
    normalized = folder / 'isodose_su_prescrizione' / (dose.stem + '_isodosemap_prescrizione.mhd')
    if not normalized.is_file() or Path(maps[0]['isodose']).resolve() != normalized.resolve():
        raise RuntimeError(f'Prescription-normalized map missing or inconsistent with report: {normalized}')
    return normalized


def nominal_coverage(normalized_path, ptv_image, isodose_percentage):
    dose_img = sitk.ReadImage(str(normalized_path))
    if not same_geometry(dose_img, ptv_image):
        raise RuntimeError('Prescription-normalized map and PTV geometry do not match')
    dose = sitk.GetArrayFromImage(dose_img).astype(np.float64, copy=False)
    ptv = sitk.GetArrayFromImage(ptv_image) > 0
    total = int(np.count_nonzero(ptv))
    if not total:
        raise RuntimeError('PTV contains no voxels')
    covered = int(np.count_nonzero(ptv & np.isfinite(dose) &
                                   (dose >= isodose_percentage / 100.0)))
    return total, covered, 100.0 * covered / total



# ================================================================
# An English individual summary for each PTV-YES plan ONLY.
# The global current-run comparison still includes every PTV status.
# ================================================================
def write_coverage_summary(row, rx, requested_iso, requested_volume, ptv_image, path):
    sim, mode = row['simulation'], row['mode']
    result = row['nominal_result']
    text = [
        'PTV COVERAGE SUMMARY — PRESCRIPTION-NORMALIZED ISODOSES',
        '=' * 84,
        f'Simulation: {sim}',
        f'Mode: {mode}',
        f'Prescription dose: {rx:g} Gy',
        f'Requested prescription isodose: {requested_iso:g}% ({rx*requested_iso/100:g} Gy)',
        f'Required PTV volume coverage: {requested_volume:g}%',
        'Normalization: D(x) / prescribed dose (not Dmax)',
        'Coverage rule: PTV voxels with dose >= requested isodose; no dose tolerance.',
        f'Nominal PTV status: {result}',
        '',
    ]
    if result == 'ERROR':
        text.extend([f'Analysis error: {row["diagnostics"]}',
                     'No nominal coverage values are available for this plan.'])
        with path.open('x', encoding='utf-8') as stream:
            stream.write('\n'.join(text) + '\n')
        return

    normalized_path = Path(row['normalized_map'])
    image = sitk.ReadImage(str(normalized_path))
    if not same_geometry(image, ptv_image):
        raise RuntimeError(f'PTV and prescription-normalized dose geometry differ: {sim} {mode}')
    iso = sitk.GetArrayFromImage(image).astype(np.float64, copy=False)
    ptv = sitk.GetArrayFromImage(ptv_image) > 0
    n_ptv = int(ptv.sum())
    if not n_ptv:
        raise RuntimeError('Empty PTV mask')

    # Report both the requested threshold and the common check isodose levels;
    # all thresholds refer to prescribed dose, never Dmax.
    levels = sorted({float(requested_iso), 100., 95., 90., 80., 50., 10.}, reverse=True)
    text.extend([
        f'PTV voxels (3D): {n_ptv}',
        f'Covered at requested level: {row["nominal_covered_voxels"]} voxels',
        f'Nominal coverage: {row["nominal_coverage_percent"]:.5f}%',
        f'Nominal decision: {result}',
        f'Prescription-normalized map: {normalized_path}',
        f'Full check_isodose_PTV.py report: {row["diagnostics"]}',
        '',
        '3D PTV COVERAGE BY PRESCRIPTION ISODOSE (NOMINAL, NO TOLERANCE)',
        f'{"Isodose %Rx":>12} | {"Threshold Gy":>13} | {"Covered voxels":>14} | {"Coverage %":>11} | Setup level',
        '-' * 82,
    ])
    finite = np.isfinite(iso)
    for level in levels:
        n = int(np.count_nonzero(ptv & finite & (iso >= level/100.0)))
        pct = 100.0 * n / n_ptv
        is_requested = abs(level - requested_iso) <= 1e-10
        text.append(f'{level:12.3f} | {rx*level/100:13.5f} | {n:14d} | {pct:11.5f} | {"*" if is_requested else ""}')
        if is_requested and abs(n - int(row['nominal_covered_voxels'])) > 0:
            raise RuntimeError(f'Coverage count mismatch in summary for {sim} {mode}')

    text.append('* marks the isodose level requested in setup.out.')

    # Slice-by-slice reporting at the exact dynamic threshold from setup.out.
    # z indices refer to the PTV_plan volume array, not the CT original slices.
    slices = []
    for z in range(ptv.shape[0]):
        mask = ptv[z]
        total = int(np.count_nonzero(mask))
        if total == 0:
            continue
        covered = int(np.count_nonzero(mask & finite[z] & (iso[z] >= requested_iso/100.0)))
        cov = covered * 100.0 / total
        slices.append((z, total, covered, cov, cov >= requested_volume - 1e-10))
    text.extend(['',
        f'SLICE-BY-SLICE PTV COVERAGE AT {requested_iso:g}% OF PRESCRIPTION',
        'These per-slice values are descriptive; the YES/NO decision above uses the 3D volume.',
        f'Slices containing PTV: {len(slices)}',
        f'Slices meeting {requested_volume:g}% coverage: {sum(x[4] for x in slices)}/{len(slices)}',
        f'{"Slice z":>8} | {"PTV voxels":>12} | {"Covered":>10} | {"Coverage %":>11} | {"Meets volume":>12}',
        '-' * 70,
    ])
    for z, total, covered, cov, good in slices:
        text.append(f'{z:8d} | {total:12d} | {covered:10d} | {cov:11.5f} | {"YES" if good else "NO":>12}')
    with path.open('x', encoding='utf-8') as stream:
        stream.write('\n'.join(text) + '\n')


def write_run_coverage_summary(records, outdir, rx, requested_iso, requested_volume):
    lines = [
        'CURRENT FULL CHAIN RUN — PTV COVERAGE SUMMARY (ALL NEW SIMULATIONS)',
        'Only plans generated in this invocation. Historical AUTO_sim are excluded.',
        f'Prescription: {rx:g} Gy; requested isodose: {requested_iso:g}%',
        f'Required covered PTV volume: {requested_volume:g}%',
        'Coverage is nominal; no tolerance is used for YES/NO.',
        '',
        f'{"SIMULATION":<42} {"MODE":<6} {"COV NOM %":>12} {"RESULT":>8}',
        '-' * 76,
    ]
    for row in records:
        cov = row['nominal_coverage_percent']
        val = f'{cov:.5f}' if isinstance(cov, float) else '-'
        lines.append(f'{row["simulation"]:<42} {row["mode"]:<6} {val:>12} {row["nominal_result"]:>8}')
    lines.extend(['', 'Per-plan summaries:',
                  '  Only PTV YES plans receive an individual summary_coverage.txt in results/.',
                  '  PTV NO plans do not receive per-plan summary files.',
                  '', 'Detailed native check output is saved separately for every plan.'])
    with (outdir / 'PTV_coverage_comparison.txt').open('x', encoding='utf-8') as stream:
        stream.write('\n'.join(lines) + '\n')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def copy_verified(src, dst):
    """Make an independent copy and verify SHA-256; never hardlink originals."""
    if not src.is_file() or src.is_symlink():
        raise RuntimeError(f'Source file missing or symbolic link not allowed: {src}')
    if dst.exists() or dst.is_symlink():
        raise FileExistsError(f'Destination file already exists: {dst}')
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    if src.stat().st_size != dst.stat().st_size or digest(src) != digest(dst):
        raise RuntimeError(f'Copy verification failed: {src} -> {dst}')


def mhd_payloads(mhd):
    """Read ElementDataFile to retain all RAW payloads required by MHD."""
    if mhd.suffix.lower() != '.mhd':
        return []
    names = []
    for line in mhd.read_text(encoding='utf-8', errors='replace').splitlines():
        if re.match(r'^\s*ElementDataFile\s*=', line, re.IGNORECASE):
            names.append(line.split('=', 1)[1].strip().strip('"'))
    if len(names) != 1:
        raise RuntimeError(f'MHD does not contain exactly one ElementDataFile: {mhd}')
    name = names[0]
    if name.upper() == 'LOCAL':
        return []
    if not name or name.upper().startswith('LIST') or '%' in name:
        raise RuntimeError(f'Unsupported ElementDataFile: {mhd}: {name}')
    if Path(name).is_absolute():
        raise RuntimeError(f'Absolute ElementDataFile path is not supported: {mhd}: {name}')
    payload = (mhd.parent / name).resolve(strict=True)
    if not payload.is_file():
        raise RuntimeError(f'RAW payload not found: {payload}')
    return [payload]


def is_other_mode(path, mode):
    opposite = 'CONV' if mode == 'FLASH' else 'FLASH'
    # Match exact mode tokens (including DVH7MeV_CONV and 7MeV_CONV.png).
    return any(re.search(r'(?:^|[_.-])' + opposite + r'(?:[_.-]|$)', part, re.I)
               for part in path.parts)


def choose_aperture(sim):
    """Use ONLY the aperture coordinates corresponding to this variant."""
    if sim.name.endswith('_plus2mm'):
        expected = sim / 'aperture_plus2mm.out'
    elif sim.name.endswith('_minus2mm'):
        expected = sim / 'aperture_minus2mm.out'
    else:
        expected = sim / 'aperture.out'
    if not expected.is_file():
        raise RuntimeError(f'Variant-specific aperture file not found: {expected}')
    return expected


def validate_mhd_tree(root):
    """Ensure all copied MHD maps have self-contained RAW payloads."""
    for mhd in root.rglob('*.mhd'):
        for payload in mhd_payloads(mhd):
            if not payload.is_relative_to(root.resolve()):
                raise RuntimeError(f'MHD depends on RAW outside results/: {mhd} -> {payload}')
            if not payload.is_file():
                raise RuntimeError(f'MHD payload missing: {mhd} -> {payload}')


def copy_winner(row, destination_root, outdir, rx, iso_level, volume, ptv_image):
    sim = ROOT / row['simulation']
    mode = row['mode']
    energy = row['energy_MeV']
    name = f'{sim.name}_{mode}'
    dest = destination_root / name
    if dest.exists() or dest.is_symlink():
        raise FileExistsError(f'Result already exists; refusing to overwrite: {dest}')
    # Copy into a temporary folder; publish only after validation.
    temp = Path(tempfile.mkdtemp(prefix=f'.{name}_INCOMPLETE_', dir=destination_root))
    try:
        selected = {}
        # Copy files for the selected mode. Shared files support reproducibility
        # (input, .dat, CT.vxl, scripts, logs, etc.).
        for source in sorted(sim.rglob('*')):
            if source.is_symlink():
                raise RuntimeError(f'Symbolic link found in simulation: {source}')
            if not source.is_file():
                continue
            relative = source.relative_to(sim)
            if any(part.startswith('isodose_prescrizione_') for part in relative.parts):
                # Never mix previous CONV/FLASH normalization runs.
                continue
            if is_other_mode(relative, mode):
                continue
            if source.name.lower().startswith('aperture') and source.suffix.lower() == '.out':
                # Aperture files are handled separately and unambiguously.
                continue
            # Exclude the other modality's DVH and plots.
            selected[relative] = source

        required_dvh = sim / f'DVH{energy}MeV_{mode}'
        if not required_dvh.is_dir() or not any(required_dvh.glob('*.txt')):
            raise RuntimeError(f'Missing or incomplete DVH for {name}: {required_dvh}')
        for expected_glob in (f'DOSE_{energy}MeV_*pulses_{mode}_BEV.mhd',
                              f'DOSE_{energy}MeV_*pulses_{mode}_FINAL.mhd'):
            if len(list(sim.glob(expected_glob))) != 1:
                raise RuntimeError(f'Expected exactly one dose map for {name}: {expected_glob}')
        for rel, source in selected.items():
            copy_verified(source, temp / rel)

        # Also provide aperture.out for a uniform folder structure, while
        # retaining the variant-specific filename.
        aperture = choose_aperture(sim)
        copy_verified(aperture, temp / 'aperture.out')
        if aperture.name != 'aperture.out':
            copy_verified(aperture, temp / aperture.name)

        # Include prescription and PTV inputs (and the RAW payload) for traceability.
        copy_verified(SETUP, temp / 'prescription_inputs' / 'setup.out')
        copy_verified(PTV, temp / 'prescription_inputs' / 'PTV_plan.mhd')
        for source in mhd_payloads(PTV):
            rel = source.relative_to(PTV.parent.resolve())
            copy_verified(source, temp / 'prescription_inputs' / rel)

        # Include the NEW prescription-normalized map, not a Dmax-normalized map.
        normalized = Path(row['normalized_map'])
        normalized_dest = temp / 'prescription_normalized_isodose' / normalized.name
        copy_verified(normalized, normalized_dest)
        for source in mhd_payloads(normalized):
            rel = source.relative_to(normalized.parent.resolve())
            copy_verified(source, normalized_dest.parent / rel)

        # Generate individual PTV summary ONLY for a plan with nominal PTV YES.
        # It is written directly into this YES plan's temporary results/ folder.
        # NO plans do not receive a summary in their AUTO_sim or elsewhere.
        write_coverage_summary(row, rx, iso_level, volume, ptv_image,
                               temp / 'summary_coverage.txt')

        report = Path(row['diagnostics'])
        copy_verified(report, temp / 'PTV_coverage_check.txt')
        normlog = outdir / f'{name}_normalize.log'
        if normlog.is_file():
            copy_verified(normlog, temp / 'prescription_normalization.log')
        metadata = {
            'source_simulation': str(sim),
            'mode': mode,
            'energy_MeV': energy,
            'aperture_file': aperture.name,
            'isodose_percent': iso_level,
            'prescription_dose_Gy': rx,
            'required_PTV_volume_percent': volume,
            'nominal_coverage_percent': row['nominal_coverage_percent'],
            'nominal_result': 'YES',
            'note': 'Only selected-mode files and shared files were copied; originals were preserved.'
        }
        (temp / 'plan_summary.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        validate_mhd_tree(temp)
        temp.rename(dest)
        return dest
    except Exception:
        (temp / 'INCOMPLETE_COPY.txt').write_text('Copy/validation not completed; original files preserved.\n', encoding='utf-8')
        raise


def main():
    params = read_setup()
    rx, required_volume, required_iso, tolerance_tool = (
        params[k] for k in ('rx', 'volume', 'isodose', 'tolerance_tool'))
    for necessary in (NORMALIZE, CHECK, PTV):
        if not necessary.is_file():
            raise RuntimeError(f'Required file not found: {necessary}')
    ptv_image = sitk.ReadImage(str(PTV))
    # Exact, ordered manifest from AUTO_IDX of THIS execution only.
    if not CURRENT_SIM_NAMES or len(CURRENT_SIM_NAMES) != len(set(CURRENT_SIM_NAMES)):
        raise RuntimeError('Empty or duplicate current-run AUTO_sim manifest')
    sims = []
    for name in CURRENT_SIM_NAMES:
        if not re.fullmatch(r'AUTO_sim(?:7|9)MeV_\d+(?:_(?:plus|minus)2mm)?', name):
            raise RuntimeError(f'Invalid current-run AUTO_sim name: {name}')
        sim = ROOT / name
        if not sim.is_dir() or sim.is_symlink() or sim.resolve().parent != ROOT:
            raise RuntimeError(f'Current-run AUTO_sim missing or unsafe: {sim}')
        sims.append(sim)
    outdir = REPORT_DIR
    if outdir.parent.resolve() != ROOT:
        raise RuntimeError(f'Invalid report directory: {outdir}')
    outdir.mkdir(exist_ok=False)
    print(f'Patient directory: {ROOT}', flush=True)
    print(f'Setup file: {SETUP}', flush=True)
    print(f'Prescription: {rx:g} Gy; isodose level: {required_iso:g}%; '
          f'minimum PTV volume coverage: {required_volume:g}%', flush=True)
    print(f'Current-run simulations only: {len(sims)}; reports: {outdir}', flush=True)

    records = []
    for sim in sims:
        energy = re.search(r'AUTO_sim(7|9)MeV_', sim.name).group(1)
        for mode in ('CONV', 'FLASH'):
            tag = f'{sim.name}_{mode}'
            row = {
                'simulation': sim.name, 'energy_MeV': energy, 'mode': mode,
                'prescription_dose_Gy': rx, 'isodose_level_percent': required_iso,
                'dose_threshold_Gy': rx * required_iso / 100.,
                'required_PTV_volume_percent': required_volume,
                'PTV_voxels': '', 'nominal_covered_voxels': '',
                'nominal_coverage_percent': '', 'nominal_result': 'ERROR',
                'normalized_map': '', 'diagnostics': '',
            }
            try:
                doses = sorted(p for p in sim.glob(f'DOSE_{energy}MeV_*pulses_{mode}_BEV.mhd')
                               if p.is_file())
                if len(doses) != 1:
                    raise RuntimeError(f'Expected one BEV dose map, found {len(doses)}: {[p.name for p in doses]}')
                normalized = normalize_rx(sim, doses[0], rx, tolerance_tool, outdir, tag)
                row['normalized_map'] = str(normalized)
                report = outdir / f'{tag}_check_isodose_ptv.txt'
                run_and_log([sys.executable, CHECK, '-i', normalized, '-p', PTV,
                             '--normalization', 'prescription', '--isodose-scale', 'fraction',
                             '--dose-prescription-gy', f'{rx:.12g}',
                             '--tolerance-gy', f'{tolerance_tool:.12g}',
                             '--slice-threshold', f'{required_volume:.12g}'], report)
                total, covered, percent = nominal_coverage(normalized, ptv_image, required_iso)
                row.update(PTV_voxels=total, nominal_covered_voxels=covered,
                           nominal_coverage_percent=round(percent, 8),
                           nominal_result='YES' if percent >= required_volume - 1e-10 else 'NO',
                           diagnostics=str(report))
            except Exception as exc:
                row['diagnostics'] = str(exc)
            records.append(row)

    # Report ALL current-run PTV checks in one run-scoped comparison file.
    # Individual summary_coverage.txt is created ONLY in results/ for PTV YES.
    write_run_coverage_summary(records, outdir, rx, required_iso, required_volume)
    print(f'Current-run coverage comparison: {outdir / "PTV_coverage_comparison.txt"}', flush=True)

    # Build the full table before copying plans; never select plans from
    # incomplete checks.
    with (outdir / 'PTV_coverage_table.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)
    print('\n' + '-' * 100, flush=True)
    print(f'{"SIMULATION":<42} {"MODE":<6} {"NOMINAL COV %":>14} {"YES/NO":>7}', flush=True)
    print('-' * 100, flush=True)
    for row in records:
        value = row['nominal_coverage_percent']
        display = f'{value:.3f}' if isinstance(value, float) else '-'
        print(f'{row["simulation"]:<42} {row["mode"]:<6} {display:>10} '
              f'{row["nominal_result"]:>7}', flush=True)
        if row['nominal_result'] == 'ERROR':
            print('  ERROR:', row['diagnostics'], flush=True)
    print('-' * 100, flush=True)
    print(f'Nominal criterion: at least {required_volume:g}% of the PTV must receive '
          f'at least {required_iso:g}% of {rx:g} Gy.', flush=True)
    print(f'Full table: {outdir / "PTV_coverage_table.csv"}', flush=True)

    errors = [r for r in records if r['nominal_result'] == 'ERROR']
    if errors:
        raise RuntimeError(f'{len(errors)} incomplete checks: results/ has NOT been updated; see CSV report')

    winners = [r for r in records if r['nominal_result'] == 'YES']
    results_dir = ROOT / 'results'
    # Preflight: never overwrite existing result directories.
    destinations = [results_dir / f'{r["simulation"]}_{r["mode"]}' for r in winners]
    already = [p for p in destinations if p.exists() or p.is_symlink()]
    if already:
        raise FileExistsError('Results already exist; refusing to overwrite: ' + ', '.join(map(str, already)))
    results_dir.mkdir(exist_ok=True)
    print(f'\nNominally compliant plans: {len(winners)}/{len(records)}', flush=True)
    for row in winners:
        dest = copy_winner(row, results_dir, outdir, rx, required_iso, required_volume, ptv_image)
        print('  COPIED AND VERIFIED:', dest, flush=True)
    no_count = sum(r['nominal_result'] == 'NO' for r in records)
    print(f'NO plans (not copied): {no_count}', flush=True)
    print('Original AUTO_sim directories preserved unchanged.', flush=True)
    print('RESULTS DIRECTORY:', results_dir, flush=True)


try:
    main()
except Exception as exc:
    print('ERROR IN FINAL PTV COVERAGE CHECK / RESULTS:', exc, file=sys.stderr)
    sys.exit(1)
