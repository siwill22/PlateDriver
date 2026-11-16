"""Export one PlateDriver run folder into the viewer's format.

    conda run -n pygmt17 python viewer/prep/export_run.py <run_dir> [--name NAME] [--label LABEL]

Writes viewer/public/data/runs/<name>/ and adds the run to viewer/public/data/runs.json.
Works on any run folder; run_SDP.py is not changed.  See docs/viewer_plan.md.

V1 exports boundaries, plates and velocities:

  run.json                 name, force model and modes, key params, saved steps, git commit
  boundaries/manifest.json petrify BoundarySeries manifest: one frame per saved step
  boundaries/NNNN.geojson  boundary polylines with type, trench polarity and plate IDs
  plates/NNNN.json         per-plate pole, rate, RMS speed, area, continental fraction and
                           label position (plate_stats.csv, plate_id_NNNN.nc)
  velocities.json          petrify VelocityField: fixed points, east/north km/Myr per step
  grids/continents/NNNN.bin  uint8 class grid (Archive layout, below): 0 ocean, 1 continent,
                           2 strong core, 3 belt
  grids/strength/NNNN.bin  uint8 rift strength (below): value x STRENGTH_SCALE, 255 no data
  diagnostics.json         every step's whole-run numbers for the plot above the time bar
                           (step_stats.csv, plate_stats.csv), and the rift and weld times

Grids use the Geode Archive byte layout: longitude fastest, -180 + i * 360 / nlon
(no repeated seam), then latitude ascending, -90 + j * 180 / (nlat - 1) (both poles).
At 0.5 deg (--grid-deg) that is 720 x 361, every second node of PlateDriver's 0.25 deg grid.

Continents: PlateDriver's continental mask is continent_NNNN.nc > 0 (craton_NNNN.nc in runs
made before November 2026; the same cells as land and as no seafloor age).  With
continent_strength: cores, each continental node takes the core/belt flag of its nearest core
point (cores_NNNN.csv, before November 2026 craton_cores_NNNN.csv), the same lookup the model
uses for continental strength (self_driving_plates._continent_factor).

Strength: the plate's strength against rifting at each node, from the same function the
rift test uses (self_driving_plates.rift_strength): strength per unit length of cut, in units
of full slab pull per unit length -- the model's vertically integrated strength.  Continent:
rift_strength_continent x the core / belt factor; ocean: oceanic_strength_ratio x
rift_strength_continent x max(oceanic_strength_floor, sqrt(min(age, 80) / 80)), age from
seafloor_age_NNNN.nc.  Suture weakening (sutures: tracers) is not included: suture points
are not saved.

Everything is on the sphere (lon/lat, WGS84 degrees):
- boundary polylines are densified along great circles to at most --densify-deg;
- trench polarity is read from the plate-ID raster 0.7 deg either side of each segment,
  as diagnose_run.py's map sequence does;
- velocities are omega x r from each plate's pole and rate, as diagnose_run.py does;
- plate labels sit at the area-weighted centroid of the plate's cells (diagnose_run.py).
"""

import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import xarray as xr
import yaml

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, REPO)
import self_driving_plates as sdp          # noqa: E402
from diagnose_run import _plate_labels     # noqa: E402

DATA = os.path.join(REPO, 'viewer', 'public', 'data')
STRENGTH_SCALE = 0.02        # strength = byte * STRENGTH_SCALE (0 to 5.08); 255 = no data

# PlateDriver's 'conservative' means "neither ridge nor subduction by the divergence
# thresholds", not a checked transform, so it is drawn as Geode's 'other'.
TYPE_MAP = {'ridge': 'ridge', 'subduction': 'subduction', 'conservative': 'other'}

POLARITY_PROBE_DEG = 0.7     # diagnose_run._overriding_side_markers offset_deg


def saved_steps(run_dir):
    """Steps with a plate-ID grid.  A step with a single plate has no boundaries file; it is
    exported with no boundaries."""
    steps = []
    for f in os.listdir(run_dir):
        m = re.fullmatch(r'plate_id_(\d{4})\.nc', f)
        if m:
            steps.append(int(m.group(1)))
    return sorted(steps)


def slerp_densify(coords, max_deg):
    """Insert great-circle points so no piece of the line is longer than max_deg."""
    xyz = sdp.lonlat_to_xyz(coords[:, 0], coords[:, 1])
    out = [xyz[0]]
    for a, b in zip(xyz[:-1], xyz[1:]):
        ang = np.degrees(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0)))
        n = int(np.ceil(ang / max_deg))
        if n > 1:
            th = np.radians(ang)
            for k in range(1, n):
                f = k / n
                out.append((np.sin((1 - f) * th) * a + np.sin(f * th) * b) / np.sin(th))
        out.append(b)
    lon, lat = sdp.xyz_to_lonlat(np.array(out))
    return np.stack([lon, lat], axis=1)


def trench_polarity(p0, p1, overriding, plate_id):
    """'Left' or 'Right' (side of travel p0 -> p1 the overriding plate is on), or None.

    Same test as diagnose_run._overriding_side_markers: step POLARITY_PROBE_DEG from the
    segment midpoint along the great-circle normal each way, and look up the plate there.
    """
    pa = sdp.lonlat_to_xyz(p0[:, 0], p0[:, 1])
    pb = sdp.lonlat_to_xyz(p1[:, 0], p1[:, 1])
    mid_lon, mid_lat = sdp.xyz_to_lonlat(pa + pb)
    normal = np.cross(pa, pb)                          # left of the direction a -> b
    lon_l, lat_l = sdp.offset_along_normal(mid_lon, mid_lat, normal, POLARITY_PROBE_DEG)
    lon_r, lat_r = sdp.offset_along_normal(mid_lon, mid_lat, -normal, POLARITY_PROBE_DEG)
    left = sdp.raster_lookup(plate_id, lon_l, lat_l) == overriding
    right = sdp.raster_lookup(plate_id, lon_r, lat_r) == overriding
    return np.where(left, 'Left', np.where(right, 'Right', ''))


def export_boundaries(src, plate_id, max_deg):
    """PlateDriver's 2-point segments, chained into polylines.

    petrify spaces trench triangles along each feature and restarts the spacing at every
    feature, so 2-point segments would each get their own triangle.  Consecutive segments
    are joined when one ends where the next starts and type, plates and polarity match.
    Returns (geojson dict, counts).
    """
    feats = src['features']
    n = len(feats)
    p0 = np.array([f['geometry']['coordinates'][0] for f in feats], float).reshape(n, 2)
    p1 = np.array([f['geometry']['coordinates'][-1] for f in feats], float).reshape(n, 2)
    props = [f['properties'] for f in feats]
    btype = [TYPE_MAP.get(p.get('boundary_type'), 'other') for p in props]
    over = np.array([int(p.get('overriding_plate') or 0) for p in props])
    pol = np.full(n, '', dtype=object)
    sub = np.array([t == 'subduction' for t in btype]) & (over > 0)
    if sub.any():
        pol[sub] = trench_polarity(p0[sub], p1[sub], over[sub], plate_id)

    def key(i):
        p = props[i]
        return (btype[i], int(p.get('plateid1_poly1') or 0), int(p.get('plateid1_poly2') or 0),
                int(p.get('downgoing_plate') or 0), int(p.get('overriding_plate') or 0), pol[i])

    out, chain, chain_key = [], [], None

    def flush():
        if not chain:
            return
        t, a, b, down, ov, side = chain_key
        coords = slerp_densify(np.array(chain), max_deg)
        prop = {'boundary_type': t, 'plates': [a, b]}
        if t == 'subduction':
            prop.update(downgoing_plate=down, overriding_plate=ov)
            if side:
                prop['polarity'] = side
        out.append({'type': 'Feature', 'properties': prop,
                    'geometry': {'type': 'LineString',
                                 'coordinates': np.round(coords, 4).tolist()}})

    for i, f in enumerate(feats):
        c = f['geometry']['coordinates']
        k = key(i)
        if chain and k == chain_key and np.allclose(chain[-1], c[0], atol=1e-9):
            chain.extend(c[1:])
        else:
            flush()
            chain, chain_key = [list(x) for x in c], k
    flush()

    counts = {t: sum(1 for f in out if f['properties']['boundary_type'] == t)
              for t in ('ridge', 'subduction', 'other')}
    counts['subduction_undecorated'] = sum(
        1 for f in out if f['properties']['boundary_type'] == 'subduction'
        and 'polarity' not in f['properties'])
    counts['segments_in'] = n
    return {'type': 'FeatureCollection', 'features': out}, counts


def arrow_points(spacing_deg):
    """Roughly equal-area points (fewer along a parallel towards the poles), as
    diagnose_run._velocity_arrows samples them."""
    lon, lat = [], []
    for la in np.arange(-78 + spacing_deg / 2, 78, spacing_deg):
        n = max(1, int(round(360 / spacing_deg * np.cos(np.radians(la)))))
        lon.append(-180 + (np.arange(n) + 0.5) * 360 / n)
        lat.append(np.full(n, la))
    return np.concatenate(lon), np.concatenate(lat)


def velocities_en(lon, lat, plate_id, rows):
    """East and north velocity (km/Myr) at each point from its plate's rotation, omega x r
    (diagnose_run._velocity_arrows).  0 where the point is on no plate."""
    pid = sdp.raster_lookup(plate_id, lon, lat)
    poles = {int(r.plate_id): sdp.lonlat_to_xyz(r.pole_lon, r.pole_lat) * np.radians(r.rate_deg_per_myr)
             for r in rows.itertuples()}
    r = sdp.lonlat_to_xyz(lon, lat)
    omega = np.array([poles.get(int(p), np.zeros(3)) if np.isfinite(p) else np.zeros(3) for p in pid])
    v = np.cross(omega, r) * sdp.EARTH_RADIUS_KM
    east = np.stack([-np.sin(np.radians(lon)), np.cos(np.radians(lon)), np.zeros_like(lon)], axis=1)
    north = np.cross(r, east)
    return (v * east).sum(axis=1), (v * north).sum(axis=1)


CONTINENT_CLASSES = [
    {'value': 0, 'name': 'ocean'},
    {'value': 1, 'name': 'continent'},
    {'value': 2, 'name': 'strong core'},
    {'value': 3, 'name': 'belt'},
]


def grid_nodes(deg):
    """Archive-layout node lon/lat (1-D) for a global grid at deg spacing."""
    nlon, nlat = int(round(360 / deg)), int(round(180 / deg)) + 1
    return -180 + np.arange(nlon) * 360 / nlon, -90 + np.arange(nlat) * 180 / (nlat - 1)


def continent_classes(craton, cores, lons, lats):
    """uint8 (nlat, nlon): 0 ocean, 1 continent, 2 craton core, 3 belt."""
    lon2, lat2 = np.meshgrid(lons, lats)
    cont = np.nan_to_num(sdp.raster_lookup(craton, lon2.ravel(), lat2.ravel())) > 0
    out = np.where(cont, 1, 0).astype(np.uint8)
    if cores is not None and len(cores) and cont.any():
        from scipy.spatial import cKDTree
        tree = cKDTree(sdp.lonlat_to_xyz(cores.lon.values.astype(float), cores.lat.values.astype(float)))
        near = tree.query(sdp.lonlat_to_xyz(lon2.ravel()[cont], lat2.ravel()[cont]))[1]
        out[cont] = np.where(cores.belt.values.astype(bool)[near], 3, 2)
    return out.reshape(lat2.shape)


def strength_bytes(craton, age, cores, params, lons, lats):
    """uint8 (nlat, nlon) rift strength (rift_strength), byte = strength / STRENGTH_SCALE."""
    lon2, lat2 = np.meshgrid(lons, lats)
    lo, la = lon2.ravel(), lat2.ravel()
    xyz = sdp.lonlat_to_xyz(lo, la)
    cont = np.nan_to_num(sdp.raster_lookup(craton, lo, la)) > 0
    a = None if age is None else sdp.raster_lookup(age, lo, la).astype(float)
    tree = None
    if cores is not None and len(cores):
        from scipy.spatial import cKDTree
        tree = cKDTree(sdp.lonlat_to_xyz(cores.lon.values.astype(float), cores.lat.values.astype(float)))
    sig = sdp.rift_strength(xyz, cont, a, params, cores if tree is not None else None, tree)
    b = np.clip(np.rint(sig / STRENGTH_SCALE), 0, 254)
    return np.where(np.isfinite(sig), b, 255).astype(np.uint8).reshape(lat2.shape)


def diagnostics(run_dir, step_stats, plate_stats):
    """Per-step series for the viewer's diagnostics plot: {time, series: [{id, label, units,
    lines: [{name, values}]}], events: {rift: [t], weld: [t]}}.  Values are rounded; NaN is
    written as null."""
    t = step_stats.time_myr.values
    by_t = plate_stats.groupby('time_myr')

    def per_t(x):
        return x.reindex(t).values

    def r(v, nd=2):
        return [None if not np.isfinite(x) else round(float(x), nd) for x in np.asarray(v, float)]

    cont = step_stats.continent_area_km2 + step_stats.continent_hidden_km2 \
        if 'continent_hidden_km2' in step_stats else None
    series = [
        dict(id='plates', label='Plates', units='count', lines=[
            dict(name='all', values=r(step_stats.n_plates, 0)),
            dict(name='< 0.5% of sphere', values=r(per_t(by_t.area_frac.apply(lambda a: (a < 0.005).sum())), 0))]),
        dict(id='rms', label='RMS speed', units='km/Myr', lines=[
            dict(name='RMS speed', values=r(step_stats.rms_speed_km_per_myr, 1))]),
        dict(id='boundaries', label='Boundary length', units='1000 km', lines=[
            dict(name='trench', values=r(step_stats.trench_km / 1e3, 1)),
            dict(name='ridge', values=r(step_stats.ridge_km / 1e3, 1))]),
        dict(id='largest', label='Largest plate', units='% of sphere', lines=[
            dict(name='largest plate', values=r(100 * per_t(by_t.area_frac.max()), 1))]),
        dict(id='netrot', label='Net rotation', units='deg/Myr', lines=[
            dict(name='net rotation', values=r(step_stats.net_rotation_deg_per_myr, 3))]),
    ]
    if cont is not None:
        series.append(dict(id='hidden', label='Continent hidden by overlap', units='% of continent', lines=[
            dict(name='hidden', values=r(100 * step_stats.continent_hidden_km2 / cont, 2))]))
    events = {}
    for key, fname in (('rift', 'rift_events.csv'), ('weld', 'welds.csv')):
        f = os.path.join(run_dir, fname)
        events[key] = [float(x) for x in pd.read_csv(f).time_myr] if os.path.exists(f) and os.path.getsize(f) else []
    return {'time': r(t, 1), 'series': series, 'events': events}


def clean(x):
    """JSON-safe float (NaN -> None)."""
    return None if x is None or not np.isfinite(x) else float(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('run_dir')
    ap.add_argument('--name', help='folder name under runs/ (default: run folder name)')
    ap.add_argument('--label', help='label shown in the viewer (default: name)')
    ap.add_argument('--densify-deg', type=float, default=1.0)
    ap.add_argument('--grid-deg', type=float, default=0.5, help='grid spacing of exported grids')
    ap.add_argument('--arrow-deg', type=float, default=10.0,
                    help='velocity arrow spacing (diagnose_run sequence uses 10)')
    args = ap.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    name = args.name or os.path.basename(run_dir.rstrip('/'))
    label = args.label or name
    out = os.path.join(DATA, 'runs', name)
    for sub in ('boundaries', 'plates', 'grids/continents', 'grids/strength'):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    info = yaml.safe_load(open(os.path.join(run_dir, 'run_info.yaml')))
    step_stats = pd.read_csv(os.path.join(run_dir, 'step_stats.csv'))
    plate_stats = pd.read_csv(os.path.join(run_dir, 'plate_stats.csv'))
    steps = saved_steps(run_dir)
    time_of = dict(zip(step_stats.step, step_stats.time_myr))

    alon, alat = arrow_points(args.arrow_deg)
    glons, glats = grid_nodes(args.grid_deg)
    has_cores = False
    has_strength = False
    params = (info.get('config') or {}).get('params') or {}
    strength_max = 0.0
    vel_frames, frames, step_list = {}, [], []
    for step in steps:
        t = float(time_of[step])
        tag = '{:04d}'.format(step)
        plate_id = xr.open_dataarray(os.path.join(run_dir, 'plate_id_{}.nc'.format(tag)))
        rows = plate_stats[plate_stats.step == step]

        bfile = os.path.join(run_dir, 'boundaries_{}.geojson'.format(tag))
        if os.path.exists(bfile):
            gj, counts = export_boundaries(json.load(open(bfile)), plate_id, args.densify_deg)
        else:
            gj = {'type': 'FeatureCollection', 'features': []}
            counts = dict(ridge=0, subduction=0, other=0, subduction_undecorated=0, segments_in=0)
        with open(os.path.join(out, 'boundaries', tag + '.geojson'), 'w') as f:
            json.dump(gj, f, separators=(',', ':'))
        frames.append({'time': t, 'file': tag + '.geojson'})

        labels = {int(p): (lo, la) for lo, la, p in _plate_labels(plate_id)}
        plates = []
        for r in rows.itertuples():
            lo, la = labels.get(int(r.plate_id), (None, None))
            plates.append({
                'id': int(r.plate_id), 'pole_lat': clean(r.pole_lat), 'pole_lon': clean(r.pole_lon),
                'rate_deg_per_myr': clean(r.rate_deg_per_myr),
                'rms_speed_km_per_myr': clean(r.rms_speed_km_per_myr),
                'area_frac': clean(r.area_frac), 'continental_frac': clean(r.continental_frac),
                'label_lon': clean(lo), 'label_lat': clean(la)})
        with open(os.path.join(out, 'plates', tag + '.json'), 'w') as f:
            json.dump({'step': step, 'time': t, 'plates': plates}, f, separators=(',', ':'))

        craton = sdp.continent_grid(run_dir, tag)
        if craton is not None:
            core_file = sdp.cores_file(run_dir, tag)
            cores = pd.read_csv(core_file) if core_file else None
            has_cores = has_cores or cores is not None
            cls = continent_classes(craton, cores, glons, glats)
            cls.tofile(os.path.join(out, 'grids', 'continents', tag + '.bin'))
            w = np.cos(np.radians(glats))[:, None]
            counts['continent_area_frac'] = round(float((w * (cls > 0)).sum() / (w * np.ones_like(cls)).sum()), 4)
            counts['belt_share_of_continent'] = round(float((w * (cls == 3)).sum() / max((w * (cls > 0)).sum(), 1e-12)), 4)
            afile = os.path.join(run_dir, 'seafloor_age_{}.nc'.format(tag))
            age = xr.open_dataset(afile).z if os.path.exists(afile) else None
            sb = strength_bytes(craton, age, cores, params, glons, glats)
            sb.tofile(os.path.join(out, 'grids', 'strength', tag + '.bin'))
            has_strength = True
            strength_max = max(strength_max, float(sb[sb < 255].max()) * STRENGTH_SCALE)

        ve, vn = velocities_en(alon, alat, plate_id, rows)
        vel_frames[str(int(round(t)))] = [np.round(ve, 2).tolist(), np.round(vn, 2).tolist()]
        step_list.append({'step': step, 'time': t, 'n_plates': len(plates), **counts})
        print('step {:4d}  t {:5.0f}  {:3d} plates  {ridge} ridge / {subduction} trench '
              '({subduction_undecorated} undecorated) / {other} other polylines from {segments_in} segments'
              .format(step, t, len(plates), **counts))

    with open(os.path.join(out, 'diagnostics.json'), 'w') as f:
        json.dump(diagnostics(run_dir, step_stats, plate_stats), f, separators=(',', ':'))
    with open(os.path.join(out, 'boundaries', 'manifest.json'), 'w') as f:
        json.dump({'frames': frames}, f, indent=1)
    with open(os.path.join(out, 'velocities.json'), 'w') as f:
        json.dump({'lon': np.round(alon, 4).tolist(), 'lat': np.round(alat, 4).tolist(),
                   'frames': vel_frames}, f, separators=(',', ':'))

    cfg = info.get('config', {})
    run = {
        'name': name, 'label': label, 'source': run_dir,
        'git_commit': info.get('git_commit'), 'started': info.get('started'),
        'force_model': info.get('force_model'),
        'force_modes': info.get('force_modes_resolved'),
        'forces': cfg.get('forces'), 'params': cfg.get('params'),
        'steps': step_list,
        'grid': {'deg': args.grid_deg, 'nlon': len(glons), 'nlat': len(glats), 'layout': 'archive'},
        'layers': [{
            'id': 'continents', 'kind': 'categorical', 'path': 'grids/continents',
            'encode_min': 0, 'encode_max': 255, 'classes': CONTINENT_CLASSES,
            'has_cores': has_cores,
        }] + ([{
            'id': 'strength', 'kind': 'scalar', 'path': 'grids/strength',
            'scale': STRENGTH_SCALE, 'nodata': 255, 'max': round(strength_max, 2),
            'label': 'Rift strength', 'units': 'slab pull per unit length',
        }] if has_strength else []),
    }
    with open(os.path.join(out, 'run.json'), 'w') as f:
        json.dump(run, f, indent=1, default=str)

    runs_file = os.path.join(DATA, 'runs.json')
    runs = json.load(open(runs_file)) if os.path.exists(runs_file) else []
    runs = [r for r in runs if r['name'] != name] + [{'name': name, 'label': label, 'path': 'runs/' + name}]
    with open(runs_file, 'w') as f:
        json.dump(runs, f, indent=1)
    print('wrote', out)


if __name__ == '__main__':
    main()
