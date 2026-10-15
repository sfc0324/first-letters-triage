#!/usr/bin/env python3
"""Finer check at OME-Zarr level 3 (8x downsampled, ~69-75 um/voxel) on 5 evenly spaced z-slabs
(64 slices each, at 10/30/50/70/90 % of the papyrus z-extent, measured on level 5).
Level 4 voxels (~140-150 um) are as thick as a papyrus sheet, so sheet_contrast is re-measured here.
Writes results/l3_slabs.json. Same definitions as scripts/02_lowres_metrics.py.
"""
import json, pathlib, numpy as np, zarr
R = pathlib.Path(__file__).resolve().parents[1]
B = 'https://vesuvius-challenge-open-data.s3.amazonaws.com/'
m4 = json.load(open(R/'results/lowres_metrics.json'))
outp = R/'results/l3_slabs.json'
res = json.load(open(outp)) if outp.exists() else {}
for sid, r in m4.items():
    if sid in res: continue
    raw = zarr.open_array(B + r['raw_store'] + '/3', mode='r')
    srf = zarr.open_array(B + r['surface_store'] + '/3', mode='r')
    Z = min(raw.shape[0], srf.shape[0]); Y = min(raw.shape[1], srf.shape[1]); X = min(raw.shape[2], srf.shape[2])
    # papyrus z-range from level 5 (x4 for level 3)
    a4 = zarr.open_array(B + r['raw_store'] + '/5', mode='r')[:]
    zz = np.flatnonzero((a4 > 0).any(axis=(1, 2))) * 4
    z0, z1 = int(zz.min()), int(min(zz.max() + 4, Z))
    slabs = []
    S = dict(npap=0, nsin=0, nsurf=0, sum_s=0.0, sum_ns=0.0, n_ns=0, sum_p=0.0, sum_p2=0.0)
    for q in (0.1, 0.3, 0.5, 0.7, 0.9):
        zc = int(z0 + q * (z1 - z0)); a, b = max(0, zc - 32), min(Z, zc + 32)
        rw = raw[a:b, :Y, :X]; sf = srf[a:b, :Y, :X] > 0
        p = rw > 0; sin = sf & p; rp = rw[p].astype(np.float64)
        d = {'z_range_level3': [a, b], 'papyrus_voxels': int(p.sum()), 'surface_frac': float(sin.sum() / max(p.sum(), 1)),
             'surface_outside_frac': float((sf & ~p).sum() / max(sf.sum(), 1)),
             'sheet_contrast': float((rw[sin].mean() - rw[p & ~sf].mean()) / rp.std()) if sin.any() else None}
        slabs.append(d)
        S['npap'] += int(p.sum()); S['nsin'] += int(sin.sum()); S['nsurf'] += int(sf.sum())
        S['sum_s'] += float(rw[sin].sum()); S['sum_ns'] += float(rw[p & ~sf].sum()); S['n_ns'] += int((p & ~sf).sum())
        S['sum_p'] += float(rp.sum()); S['sum_p2'] += float((rp ** 2).sum())
    mp = S['sum_p'] / S['npap']; sd = (S['sum_p2'] / S['npap'] - mp ** 2) ** 0.5
    res[sid] = {'level': 3, 'slabs': slabs,
                'surface_frac': S['nsin'] / S['npap'],
                'surface_outside_frac': (S['nsurf'] - S['nsin']) / S['nsurf'],
                'sheet_contrast': (S['sum_s'] / S['nsin'] - S['sum_ns'] / S['n_ns']) / sd,
                'sheet_contrast_min_slab': min(s['sheet_contrast'] for s in slabs if s['sheet_contrast'] is not None)}
    json.dump(res, open(outp, 'w'), indent=1)
    print(sid, {k: round(v, 4) for k, v in res[sid].items() if isinstance(v, float)}, flush=True)
