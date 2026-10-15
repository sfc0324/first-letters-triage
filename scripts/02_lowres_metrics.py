#!/usr/bin/env python3
"""CPU-only low-resolution metrics for each eligible volume.

Reads OME-Zarr level L (default 4 = 16x downsampled) of the eligible CT volume and of its
published surface prediction directly over HTTPS from the public bucket. No GPU, no auth.
Writes results/lowres_metrics.json and per-scroll preview PNGs in figures/previews/.

Definitions (all on level L voxels):
  papyrus        raw > 0   (eligible volumes are the *-masked.zarr exports: outside-scroll voxels are 0)
  surface        surface prediction > 0 (stores are binary 0/255, threshold 0.2 baked in)
  surface_frac   |surface & papyrus| / |papyrus|
  surface_outside_frac  |surface & ~papyrus| / |surface|   (predictions outside the mask)
  sheet_contrast (mean raw on surface&papyrus - mean raw on papyrus&~surface) / std raw on papyrus
  z_coverage     fraction of z-slices containing papyrus whose surface_frac >= 0.5 * median slice surface_frac
"""
import json, sys, pathlib, time
import numpy as np, zarr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = pathlib.Path(__file__).resolve().parents[1]
B = 'https://vesuvius-challenge-open-data.s3.amazonaws.com/'
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
inv = json.load(open(R/'results/inventory.json'))
outp = R/'results/lowres_metrics.json'
res = json.load(open(outp)) if outp.exists() else {}
(R/'figures/previews').mkdir(parents=True, exist_ok=True)

for o in inv:
    sid = o['scroll']
    if sid in res and res[sid].get('level') == L: continue
    t = next(v for v in o['volumes'] if v['is_eligible_volume'])
    rawp = t['assets']['ome-zarr'][0]['paths'][0].rstrip('/')
    surfp = t['assets']['surface-prediction-zarr'][0]['paths'][0].rstrip('/')
    t0 = time.time()
    raw = zarr.open_array(B + rawp + f'/{L}', mode='r')[:]
    srf = zarr.open_array(B + surfp + f'/{L}', mode='r')[:]
    if raw.shape != srf.shape:
        z = min(raw.shape[0], srf.shape[0]); y = min(raw.shape[1], srf.shape[1]); x = min(raw.shape[2], srf.shape[2])
        shape_note = f'shape mismatch raw {raw.shape} vs surface {srf.shape}; cropped to common'
        raw, srf = raw[:z, :y, :x], srf[:z, :y, :x]
    else:
        shape_note = None
    pap = raw > 0; s = srf > 0
    npap = int(pap.sum()); nsurf = int(s.sum()); sin = s & pap; nsin = int(sin.sum())
    rp = raw[pap].astype(np.float64)
    mean_s = float(raw[sin].mean()) if nsin else None
    mean_ns = float(raw[pap & ~s].mean())
    std_p = float(rp.std())
    pz = pap.sum(axis=(1, 2)); sz = sin.sum(axis=(1, 2))
    zmask = pz > 0
    fz = sz[zmask] / pz[zmask]
    med = float(np.median(fz))
    vox_um = t['pixel_size_um'] * 2 ** L
    r = {'level': L, 'raw_store': rawp, 'surface_store': surfp,
         'surface_model': t['assets']['surface-prediction-zarr'][0].get('parameters'),
         'shape_zyx_level': list(raw.shape), 'voxel_um_level': vox_um, 'shape_note': shape_note,
         'papyrus_voxels': npap, 'papyrus_volume_cm3': npap * (vox_um * 1e-4) ** 3,
         'papyrus_z_extent_mm': float(zmask.sum() * vox_um / 1000),
         'surface_voxels': nsurf, 'surface_frac': nsin / npap,
         'surface_outside_frac': (nsurf - nsin) / nsurf if nsurf else None,
         'raw_mean_papyrus': float(rp.mean()), 'raw_std_papyrus': std_p,
         'raw_mean_surface': mean_s, 'raw_mean_nonsurface': mean_ns,
         'sheet_contrast': (mean_s - mean_ns) / std_p if mean_s is not None else None,
         'z_median_surface_frac': med, 'z_coverage': float((fz >= 0.5 * med).mean()),
         'seconds': round(time.time() - t0, 1)}
    # previews: middle axial slice and middle coronal slice, surface overlaid in red
    zc = int(np.flatnonzero(zmask).mean())
    yc = raw.shape[1] // 2
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [1, raw.shape[0] / raw.shape[1] * 0 + 1]})
    for a, img, sm, title in [(ax[0], raw[zc], s[zc], f'axial z={zc}'), (ax[1], raw[:, yc, :], s[:, yc, :], f'coronal y={yc}')]:
        a.imshow(img, cmap='gray', vmin=0, vmax=255)
        ov = np.zeros(img.shape + (4,)); ov[sm] = [1, 0, 0, 0.45]; a.imshow(ov)
        a.set_title(title, fontsize=9); a.axis('off')
    fig.suptitle(f'{sid} {o["eligible_volume"]}  level {L} ({vox_um:.0f} um/vox)  red = surface prediction', fontsize=10)
    fig.tight_layout(); fig.savefig(R/f'figures/previews/{sid}.png', dpi=90); plt.close(fig)
    res[sid] = r
    json.dump(res, open(outp, 'w'), indent=1)
    print(sid, {k: r[k] for k in ('papyrus_volume_cm3', 'surface_frac', 'sheet_contrast', 'z_coverage', 'seconds')}, flush=True)
