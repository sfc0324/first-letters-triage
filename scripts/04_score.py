#!/usr/bin/env python3
"""Combine inventory + zarrsweep integrity + low-res metrics into a transparent triage score.

Every input number is read from results/*.json (produced by scripts 01-03) or from the zarrsweep
full-sweep report; nothing is typed in by hand. Each factor is converted to a percentile rank among
the 21 eligible volumes (0 = worst, 1 = best) so units do not matter, then weighted.

This is a *data-readiness* heuristic for choosing where to spend ink-hunting effort.
It does NOT measure ink and cannot predict whether letters exist.
"""
import json, csv, pathlib, collections, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = pathlib.Path(__file__).resolve().parents[1]
SWEEP = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else R/'data/zarrsweep_full_sweep.json'
inv = {o['scroll']: o for o in json.load(open(R/'results/inventory.json'))}
m4 = json.load(open(R/'results/lowres_metrics.json'))
m3 = json.load(open(R/'results/l3_slabs.json'))
sweep = {r['store'].rstrip('/'): r for r in json.load(open(SWEEP))}

WEIGHTS = {  # factor: (weight, description)
    'sheet_contrast_l3':     (0.30, 'L3 raw-intensity contrast surface vs non-surface papyrus (5 slabs)'),
    'sheet_contrast_l3_min': (0.10, 'worst of the 5 L3 slabs (penalises locally poor regions)'),
    'surface_clean':         (0.20, '1 - fraction of predicted surface voxels outside the scroll mask (L4)'),
    'z_coverage':            (0.10, 'fraction of papyrus z-slices with normal surface density (L4)'),
    'tooling':               (0.20, 'mean of: published segments on this volume, lasagna, umbilicus'),
    'papyrus_volume_cm3':    (0.10, 'papyrus volume (more surface to search)'),
}
ALT = {'equal': {k: 1 / len(WEIGHTS) for k in WEIGHTS},
       'data_only': {'sheet_contrast_l3': .4, 'sheet_contrast_l3_min': .2, 'surface_clean': .25, 'z_coverage': .15,
                     'tooling': 0, 'papyrus_volume_cm3': 0}}

rows = []
for sid, o in inv.items():
    t = next(v for v in o['volumes'] if v['is_eligible_volume']); a = t['assets']
    stores = [p.rstrip('/') for k in ('ome-zarr', 'surface-prediction-zarr') for e in a.get(k, []) for p in e['paths']]
    sev = collections.Counter(f['severity'] for s in stores for f in sweep.get(s, {'findings': []})['findings'])
    swept = all(s in sweep for s in stores)
    objs = sum(l['objects'] for s in stores if s in sweep for l in sweep[s]['levels'])
    seg = o['segments_by_volume'].get(o['eligible_volume'], 0)
    hires = [f"{v['volume_id']} ({v['pixel_size_um']} um)" for v in o['volumes']
             if not v['is_eligible_volume'] and v['pixel_size_um'] and v['pixel_size_um'] < t['pixel_size_um']]
    r4, r3 = m4[sid], m3[sid]
    rows.append({
        'scroll': sid, 'reference': o['reference'], 'volume': o['eligible_volume'],
        'pixel_size_um': t['pixel_size_um'], 'energy_keV': t['energy_keV'],
        'surface_model': (r4['surface_model'] or {}).get('model_id'),
        'zarrsweep_swept': swept, 'zarrsweep_errors': sev.get('error', 0), 'zarrsweep_warnings': sev.get('warning', 0),
        'zarrsweep_objects': objs,
        'papyrus_volume_cm3': r4['papyrus_volume_cm3'], 'papyrus_z_extent_mm': r4['papyrus_z_extent_mm'],
        'surface_frac_l4': r4['surface_frac'], 'surface_outside_frac_l4': r4['surface_outside_frac'],
        'surface_clean': 1 - r4['surface_outside_frac'], 'z_coverage': r4['z_coverage'],
        'sheet_contrast_l4': r4['sheet_contrast'], 'sheet_contrast_l3': r3['sheet_contrast'],
        'sheet_contrast_l3_min': r3['sheet_contrast_l3_min'] if 'sheet_contrast_l3_min' in r3 else r3['sheet_contrast_min_slab'],
        'segments_on_volume': seg, 'lasagna': int(bool(a.get('lasagna'))), 'umbilicus': int(bool(a.get('umbilicus'))),
        'tooling': (int(seg > 0) + int(bool(a.get('lasagna'))) + int(bool(a.get('umbilicus')))) / 3,
        'higher_res_volumes': '; '.join(hires),
    })

elig = [r for r in rows if not r['reference']]
def pct(vals, x):  # percentile rank among eligible (ties averaged), 0..1
    v = np.array(vals); return float(((v < x).sum() + 0.5 * ((v == x).sum() - 1)) / (len(v) - 1))
for scheme_name, W in [('score', {k: w for k, (w, _) in WEIGHTS.items()})] + list(ALT.items()):
    for r in rows:
        if r['zarrsweep_errors'] > 0 or not r['zarrsweep_swept']:
            r[scheme_name] = 0.0; continue
        r[scheme_name] = sum(W[k] * pct([e[k] for e in elig], r[k]) for k in W)
for name in ['score', 'equal', 'data_only']:
    for i, r in enumerate(sorted(elig, key=lambda r: -r[name]), 1): r['rank_' + name] = i
elig.sort(key=lambda r: r['rank_score'])
refs = [r for r in rows if r['reference']]

cols = ['rank_score', 'scroll', 'volume', 'score', 'rank_equal', 'rank_data_only', 'pixel_size_um', 'energy_keV',
        'surface_model', 'zarrsweep_errors', 'zarrsweep_warnings', 'zarrsweep_objects', 'sheet_contrast_l3',
        'sheet_contrast_l3_min', 'sheet_contrast_l4', 'surface_frac_l4', 'surface_outside_frac_l4', 'z_coverage',
        'papyrus_volume_cm3', 'papyrus_z_extent_mm', 'segments_on_volume', 'lasagna', 'umbilicus', 'higher_res_volumes']
with open(R/'results/ranking.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols + ['reference'], extrasaction='ignore'); w.writeheader()
    for r in elig + refs: w.writerow(r)
json.dump({'weights': {k: {'weight': w, 'meaning': d} for k, (w, d) in WEIGHTS.items()},
           'alternative_weightings': ALT, 'eligible': elig, 'references': refs},
          open(R/'results/ranking.json', 'w'), indent=1)

def f(x, n=3): return f'{x:.{n}f}'
md = ['| # | scroll | volume | score | rank (equal / data-only) | sheet contrast L3 (min slab) | surface outside mask | z-coverage | papyrus cm³ | segments | lasagna | umbilicus | zarrsweep errors |',
      '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for r in elig + refs:
    md.append(f"| {r['rank_score'] if not r['reference'] else 'ref'} | {r['scroll']} | {r['volume']} | {f(r['score'])} | "
              f"{r.get('rank_equal', '-')} / {r.get('rank_data_only', '-')} | {f(r['sheet_contrast_l3'])} ({f(r['sheet_contrast_l3_min'])}) | "
              f"{f(r['surface_outside_frac_l4'])} | {f(r['z_coverage'])} | {r['papyrus_volume_cm3']:.0f} | {r['segments_on_volume']} | "
              f"{r['lasagna']} | {r['umbilicus']} | {r['zarrsweep_errors']} |")
open(R/'results/ranking.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))

# figures
fig, ax = plt.subplots(figsize=(9, 6))
names = [r['scroll'] for r in elig][::-1]
contrib = {k: [w * pct([e[k] for e in elig], r[k]) for r in elig][::-1] for k, (w, _) in WEIGHTS.items()}
left = np.zeros(len(names))
for k, v in contrib.items():
    ax.barh(names, v, left=left, label=k); left += np.array(v)
ax.set_xlabel('triage score (weighted percentile ranks)'); ax.legend(fontsize=7, loc='lower right')
ax.set_title('First Letters eligible volumes: data-readiness triage'); fig.tight_layout()
fig.savefig(R/'figures/ranking.png', dpi=110); plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 5))
for r in elig:
    ax.scatter(r['surface_outside_frac_l4'], r['sheet_contrast_l3'], s=r['papyrus_volume_cm3'] / 2, c='tab:blue', alpha=.6)
    ax.annotate(r['scroll'].replace('PHerc', ''), (r['surface_outside_frac_l4'], r['sheet_contrast_l3']), fontsize=7)
for r in refs:
    ax.scatter(r['surface_outside_frac_l4'], r['sheet_contrast_l3'], s=r['papyrus_volume_cm3'] / 2, c='tab:red', marker='*')
    ax.annotate(r['scroll'].replace('PHerc', '') + ' (ink found)', (r['surface_outside_frac_l4'], r['sheet_contrast_l3']), fontsize=7, color='tab:red')
ax.set_xlabel('surface predictions outside scroll mask (L4 fraction, lower = cleaner)')
ax.set_ylabel('sheet contrast at L3 (higher = sheets stand out)')
ax.set_title('Eligible volumes (blue, size = papyrus volume) vs ink-found references (red)')
fig.tight_layout(); fig.savefig(R/'figures/contrast_vs_clean.png', dpi=110); plt.close(fig)
