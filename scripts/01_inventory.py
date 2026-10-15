#!/usr/bin/env python3
"""Inventory every public asset of the 21 First Letters-eligible scrolls.

Inputs (all fetched by scripts/00_fetch.sh, nothing hand-typed):
  data/prizeEligibility.json  - official eligible list (villa repo, source of scrollprize.org/prizes)
  data/metadata.json          - s3://vesuvius-challenge-open-data/metadata.json
Outputs: results/inventory.json, results/inventory.csv
"""
import json, csv, collections, pathlib
R = pathlib.Path(__file__).resolve().parents[1]
elig = json.load(open(R/'data/prizeEligibility.json'))['first-letters-2027']
# Cross-check against the rendered scrollprize.org/prizes page (2nd list on the page = First Letters).
import re
for snap in sorted(R.glob('data/prizes_snapshot_*.html')):
    page = re.findall(r'>(PHerc\w+)</a> — <code>(\d+)</code>', snap.read_text())
    gp_n = len(json.load(open(R/'data/prizeEligibility.json'))['grand-prize-2027'])
    assert page[gp_n:] == [(e['scroll'], e['volume']) for e in elig], f'{snap.name} disagrees with prizeEligibility.json'
    print(f'OK: {snap.name} lists the same {len(elig)} First Letters volumes as prizeEligibility.json')
# Calibration references: scrolls where ink/letters HAVE been found at a comparable ~8.6 um resolution
# (PHerc0343: First Letters awarded 2026-10-08; PHerc1447: ink found, per scrollprize.substack.com).
REFERENCES = [{'scroll': 'PHerc0343', 'volume': '20250521140437', 'reference': True},
              {'scroll': 'PHerc1447', 'volume': '20250521151220', 'reference': True}]
elig = [dict(e, reference=False) for e in elig] + REFERENCES
md = json.load(open(R/'data/metadata.json'))['samples']

def paths(item):
    return [o['path'] for o in item.get('origins', [])]

out = []
for e in elig:
    sid, vid = e['scroll'], e['volume']
    s = md[sid]
    vols = []
    for v in s['volumes'].values():
        p = v.get('properties') or {}
        dtypes = collections.defaultdict(list)
        for dd in v['data']:
            entry = {'paths': paths(dd)}
            if dd.get('parameters'): entry['parameters'] = dd['parameters']
            dtypes[dd['type']].append(entry)
        vols.append({'volume_id': v['id'], 'long_id': v['long_id'], 'scan_id': v.get('scan_id'),
                     'pixel_size_um': p.get('pixel_size_um'), 'energy_keV': p.get('energy_keV'),
                     'detector_distance_mm': p.get('detector_distance_mm'), 'shape_zyx': p.get('shape'),
                     'superseded_by': p.get('superseded_by'), 'is_eligible_volume': v['id'] == vid,
                     'assets': dtypes})
    segs = collections.Counter()
    for g in s.get('segments', {}).values():
        segs[g.get('original_volume_id')] += 1
    scans = [{'scan_id': k, 'long_id': sc['long_id'],
              'beamline': sc['creation']['metadata'].get('beamline'),
              'location': sc['creation']['metadata'].get('location')} for k, sc in s['scans'].items()]
    target = next((v for v in vols if v['is_eligible_volume']), None)
    out.append({'reference': e['reference'], 'scroll': sid, 'eligible_volume': vid, 'eligible_volume_found_in_metadata': target is not None,
                'scans': scans, 'volumes': vols, 'segments_by_volume': dict(segs),
                'segments_total': sum(segs.values())})

json.dump(out, open(R/'results/inventory.json', 'w'), indent=1)
cols = ['reference', 'scroll', 'eligible_volume', 'pixel_size_um', 'energy_keV', 'detector_distance_mm', 'shape_zyx',
        'n_scans', 'n_volumes', 'other_volumes', 'surface_predictions', 'ink_predictions', 'lasagna',
        'normal_grids', 'umbilicus', 'segments_on_eligible_volume', 'segments_total']
with open(R/'results/inventory.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(cols)
    for o in out:
        t = next(v for v in o['volumes'] if v['is_eligible_volume']); a = t['assets']
        others = ['%s@%sum/%skeV' % (v['volume_id'], v['pixel_size_um'], v['energy_keV']) for v in o['volumes'] if not v['is_eligible_volume']]
        w.writerow([o['reference'], o['scroll'], o['eligible_volume'], t['pixel_size_um'], t['energy_keV'], t['detector_distance_mm'],
                    'x'.join(map(str, t['shape_zyx'] or [])), len(o['scans']), len(o['volumes']), ';'.join(others),
                    len(a.get('surface-prediction-zarr', [])), len(a.get('ink-detection-3d-zarr', [])),
                    len(a.get('lasagna', [])), len(a.get('normal-grids', [])), len(a.get('umbilicus', [])),
                    o['segments_by_volume'].get(o['eligible_volume'], 0), o['segments_total']])
print(open(R/'results/inventory.csv').read())
