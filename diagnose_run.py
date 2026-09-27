"""Diagnostic figures for one or more PlateDriver runs.

    python diagnose_run.py <output_dir> [<output_dir> ...] [--labels A B] [--out figs]
                           [--steps 0 50 100] [--every 5] [--min-area-frac 0.002]

Reads step_stats.csv and plate_stats.csv (written every step by run_SDP.py) and the
saved grids, and writes PyGMT figures to --out (default: the first output_dir):

  sequence.png        the main figure: the plate configuration through time as small
                      multiples (seafloor age, ridges, trenches with their polarity,
                      velocity arrows, plate IDs).  Several runs: one column per run.
  fu_plot.png         plate RMS speed against trench fraction (downgoing trench length
                      / perimeter) and against continental fraction, one point per plate
                      per sampled step.  Forsyth & Uyeda (1975) found speed rises with
                      trench fraction and falls with continental fraction.
  timeseries.png      global RMS speed (and fastest plate), net rotation, number of
                      plates and wall-clock seconds per step.
  maps_NNNN.png       for each requested step: seafloor age, plate boundaries by type
                      and plate velocity arrows, one map per run.

Several runs are overlaid (fu_plot, timeseries) or stacked (maps) for comparison.
"""

import argparse
import os

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import pygmt

import self_driving_plates as sdp

COLOURS = ['dodgerblue3', 'orange', 'seagreen', 'purple', 'firebrick', 'gray30']
R = sdp.EARTH_RADIUS_KM


def load_run(path):
    steps = pd.read_csv(os.path.join(path, 'step_stats.csv'))
    plates = pd.read_csv(os.path.join(path, 'plate_stats.csv'))
    return steps, plates


def fu_plot(runs, labels, out, every, min_area_frac):
    fig = pygmt.Figure()
    samples = []
    for (steps, plates) in runs:
        p = plates[(plates.step % every == 0) & (plates.area_frac >= min_area_frac)]
        samples.append(p)
    vmax = max(float(p.rms_speed_km_per_myr.max()) for p in samples if len(p)) * 1.1
    with fig.subplot(nrows=1, ncols=2, figsize=('18c', '8c'), margins=['1.2c', '0.5c']):
        for panel, (col, xlabel) in enumerate([('trench_frac', 'Trench fraction (downgoing / perimeter)'),
                                                ('continental_frac', 'Continental fraction of plate area')]):
            with fig.set_panel(panel=panel):
                fig.basemap(region=[0, 1, 0, vmax], projection='X?',
                            frame=['WSne', 'xaf+l' + xlabel, 'yaf+lPlate RMS speed (km/Myr)'])
                for k, (p, label) in enumerate(zip(samples, labels)):
                    if not len(p):
                        continue
                    r = np.corrcoef(p[col], p.rms_speed_km_per_myr)[0, 1] if len(p) > 2 else np.nan
                    fig.plot(x=p[col], y=p.rms_speed_km_per_myr, style='c0.12c',
                             fill=COLOURS[k % len(COLOURS)], pen='0.2p,black',
                             transparency=30, label='{} (r = {:.2f})'.format(label, r))
                fig.legend(position='jTR+o0.1c', box='+gwhite+p0.5p')
    fig.savefig(os.path.join(out, 'fu_plot.png'), dpi=200)


def timeseries(runs, labels, out):
    panels = [
        ('rms_speed_km_per_myr', 'RMS speed (km/Myr)', 'max_plate_rms_km_per_myr'),
        ('net_rotation_deg_per_myr', 'Net rotation (deg/Myr)', None),
        ('n_plates', 'Number of plates', None),
        ('seconds_step', 'Seconds per step', None),
    ]
    tmax = max(float(s.time_myr.max()) for s, _ in runs)
    fig = pygmt.Figure()
    with fig.subplot(nrows=len(panels), ncols=1, figsize=('16c', '20c'), sharex='b',
                     margins=['0.5c', '0.4c']):
        for panel, (col, ylabel, col2) in enumerate(panels):
            vals = [s[col].dropna() for s, _ in runs if col in s]
            if col2:
                vals += [s[col2].dropna() for s, _ in runs if col2 in s]
            vals = [v for v in vals if len(v)]
            ymax = max(float(v.max()) for v in vals) * 1.1 if vals else 1.0
            with fig.set_panel(panel=panel):
                bottom = panel == len(panels) - 1
                fig.basemap(region=[0, max(tmax, 1), 0, max(ymax, 1e-6)], projection='X?',
                            frame=['WSne' if bottom else 'Wsne',
                                   'xaf+lTime (Myr)' if bottom else 'xaf', 'yaf+l' + ylabel])
                for k, ((steps, _), label) in enumerate(zip(runs, labels)):
                    if col not in steps:
                        continue
                    colour = COLOURS[k % len(COLOURS)]
                    ok = steps[col].notna()
                    fig.plot(x=steps.time_myr[ok], y=steps[col][ok], pen='1.2p,' + colour,
                             label=label if panel == 0 else None)
                    if col2 and col2 in steps:
                        fig.plot(x=steps.time_myr, y=steps[col2], pen='0.8p,{},--'.format(colour),
                                 label=(label + ' fastest plate') if panel == 0 else None)
                if panel == 0:
                    fig.legend(position='jTR+o0.1c', box='+gwhite+p0.5p')
    fig.savefig(os.path.join(out, 'timeseries.png'), dpi=200)


def _velocity_arrows(plate_id, plate_rows, spacing_deg=6.0):
    """Arrow positions, azimuths and speeds (km/Myr) from plate IDs and plate poles."""
    # Roughly equal-area sampling: fewer arrows along a parallel towards the poles
    lon, lat = [], []
    for la in np.arange(-78 + spacing_deg / 2, 78, spacing_deg):
        n = max(1, int(round(360 / spacing_deg * np.cos(np.radians(la)))))
        lon.append(-180 + (np.arange(n) + 0.5) * 360 / n)
        lat.append(np.full(n, la))
    lon, lat = np.concatenate(lon), np.concatenate(lat)
    pid = sdp.raster_lookup(plate_id, lon, lat)
    poles = {int(r.plate_id): sdp.lonlat_to_xyz(r.pole_lon, r.pole_lat) * np.radians(r.rate_deg_per_myr)
             for r in plate_rows.itertuples()}
    r = sdp.lonlat_to_xyz(lon, lat)
    omega = np.array([poles.get(int(p), np.zeros(3)) if np.isfinite(p) else np.zeros(3) for p in pid])
    v = np.cross(omega, r) * R
    east = np.stack([-np.sin(np.radians(lon)), np.cos(np.radians(lon)), np.zeros_like(lon)], axis=1)
    north = np.cross(r, east)
    ve, vn = (v * east).sum(axis=1), (v * north).sum(axis=1)
    speed = np.hypot(ve, vn)
    keep = speed > 0.5
    return lon[keep], lat[keep], np.degrees(np.arctan2(ve, vn))[keep] % 360, speed[keep]


def _plate_labels(plate_id):
    lons, lats = sdp.raster_coords(plate_id)
    lon2d, lat2d = np.meshgrid(lons, lats)
    vals = plate_id.values
    w = sdp.cell_area_weights(lats, vals.shape)
    out = []
    for pid in np.unique(vals[np.isfinite(vals) & (vals > 0)]).astype(int):
        m = vals == pid
        c = (sdp.lonlat_to_xyz(lon2d[m], lat2d[m]) * w[m][:, None]).sum(axis=0)
        out.append((*sdp.xyz_to_lonlat(c), str(pid)))
    return out


def _overriding_side_markers(b, plate_id, offset_deg=0.7, every=3):
    """Points just inside the overriding plate along subduction segments (thinned), so
    subduction polarity can be read off a map: markers sit on the overriding side."""
    sub = b[(b.boundary_type == 'subduction') & (b.overriding_plate > 0)]
    if not len(sub):
        return np.array([]), np.array([])
    sub = sub.iloc[::every]
    coords = sub.geometry.get_coordinates()
    first = coords.groupby(level=0).head(1)
    last = coords.groupby(level=0).tail(1)
    pa = sdp.lonlat_to_xyz(first.x.values, first.y.values)
    pb = sdp.lonlat_to_xyz(last.x.values, last.y.values)
    mid_lon, mid_lat = sdp.xyz_to_lonlat(pa + pb)
    normal = np.cross(pa, pb)                      # left of the direction a -> b
    lon_l, lat_l = sdp.offset_along_normal(mid_lon, mid_lat, normal, offset_deg)
    lon_r, lat_r = sdp.offset_along_normal(mid_lon, mid_lat, -normal, offset_deg)
    over = sub.overriding_plate.values
    left = sdp.raster_lookup(plate_id, lon_l, lat_l) == over
    right = sdp.raster_lookup(plate_id, lon_r, lat_r) == over
    lon = np.where(left, lon_l, np.where(right, lon_r, np.nan))
    lat = np.where(left, lat_l, np.where(right, lat_r, np.nan))
    ok = np.isfinite(lon)
    return lon[ok], lat[ok]


def draw_map_panel(fig, path, steps, plates, step, title, arrow_spacing=6.0,
                   label_font='9p', annotate=True):
    """One map of a run's configuration at a saved step: seafloor age (continents grey),
    ridges red, trenches blue with triangles on the overriding side, transforms and other
    boundaries black, plate velocity arrows and plate IDs.  Draws into the current panel."""
    # grey map interior shows through where the age grid is NaN (continents)
    fig.basemap(region='d', projection='N?', frame=['f', '+ggray75+t' + title])
    age_file = os.path.join(path, 'seafloor_age_{:04d}.nc'.format(step))
    if os.path.exists(age_file):
        pygmt.makecpt(cmap='batlow', series=[0, 150], reverse=True)
        fig.grdimage(age_file, cmap=True, nan_transparent=True)
    pid_file = os.path.join(path, 'plate_id_{:04d}.nc'.format(step))
    plate_id = xr.open_dataarray(pid_file) if os.path.exists(pid_file) else None
    bfile = os.path.join(path, 'boundaries_{:04d}.geojson'.format(step))
    if os.path.exists(bfile):
        b = gpd.read_file(bfile)
        for kind, pen in [('conservative', '0.6p,black'), ('ridge', '1.2p,red'),
                          ('subduction', '1.2p,blue')]:
            sel = b[b.boundary_type == kind]
            if len(sel):
                fig.plot(data=sel, pen=pen)
        if plate_id is not None:
            lon, lat = _overriding_side_markers(b, plate_id)
            if len(lon):
                fig.plot(x=lon, y=lat, style='t0.13c', fill='blue', pen='0.2p,white')
    rows = plates[plates.step == step]
    if plate_id is not None and len(rows):
        x, y, az, sp = _velocity_arrows(plate_id, rows, spacing_deg=arrow_spacing)
        fig.plot(x=x, y=y, direction=[az, sp / 250.0],
                 style='V0.1c+e+a45+n', pen='0.5p,black', fill='black')
        for lon, lat, text in _plate_labels(plate_id):
            fig.text(x=lon, y=lat, text=text, font=label_font + ',Helvetica-Bold,white',
                     fill='black@40', clearance='0.04c')
    srow = steps[steps.step == step]
    if annotate and len(srow):
        fig.text(position='BL', offset='0.1c/0.1c', font='7p,Helvetica,black',
                 fill='white@20', justify='BL',
                 text='RMS {:.0f} km/Myr, net rot {:.2f} deg/Myr, {} plates'.format(
                     float(srow.rms_speed_km_per_myr.iloc[0]),
                     float(srow.net_rotation_deg_per_myr.iloc[0]), int(srow.n_plates.iloc[0])))


def _map_key(fig, width='9c'):
    pygmt.makecpt(cmap='batlow', series=[0, 150], reverse=True)
    fig.colorbar(frame='xaf+lSeafloor age (Myr); grey = continent', position='JBC+w{}/0.3c+o0/0.7c'.format(width))
    fig.text(position='BC', offset='0/-2.4c', no_clip=True, font='9p,Helvetica,black',
             text='red: ridge   blue: trench (triangles on overriding plate)   black: other boundary   '
                  'arrows: plate velocity, 1 cm = 250 km/Myr')


def maps(paths, runs, labels, out, steps_wanted):
    """One figure per step, one map per run."""
    for step in steps_wanted:
        fig = pygmt.Figure()
        n = len(paths)
        with fig.subplot(nrows=n, ncols=1, figsize=('18c', '{}c'.format(10 * n)),
                         margins=['0.2c', '0.6c'], title='Step {}'.format(step)):
            for k, (path, (steps, plates), label) in enumerate(zip(paths, runs, labels)):
                with fig.set_panel(panel=k):
                    draw_map_panel(fig, path, steps, plates, step, label)
        _map_key(fig)
        fig.savefig(os.path.join(out, 'maps_{:04d}.png'.format(step)), dpi=200)


def sequence(paths, runs, labels, out, seq_steps, ncols=3):
    """The configuration through time as small multiples.  One run: a grid of time
    slices read left to right.  Several runs: one row per time slice, one column per run,
    so the runs can be compared at the same model time."""
    n_runs = len(paths)
    if n_runs == 1:
        nrows = int(np.ceil(len(seq_steps) / ncols))
        cells = [(k, seq_steps[k], 0) for k in range(len(seq_steps))]
        cols = ncols
    else:
        nrows, cols = len(seq_steps), n_runs
        cells = [(r * n_runs + c, st, c) for r, st in enumerate(seq_steps) for c in range(n_runs)]
    width = 8.5 * cols
    fig = pygmt.Figure()
    with pygmt.config(FONT_TITLE='10p', FONT_HEADING='14p', MAP_HEADING_OFFSET='1.2c'):
        with fig.subplot(nrows=nrows, ncols=cols, figsize=('{}c'.format(width), '{}c'.format(5.3 * nrows)),
                         margins=['0.15c', '0.35c'],
                         title='Plate configuration through time: ' + ' | '.join(labels)):
            for panel, step, run in cells:
                path, (steps, plates) = paths[run], runs[run]
                t = steps.loc[steps.step == step, 'time_myr']
                title = '{}t = {:.0f} Myr (step {})'.format(
                    labels[run] + ', ' if n_runs > 1 else '',
                    float(t.iloc[0]) if len(t) else float('nan'), step)
                with fig.set_panel(panel=panel):
                    draw_map_panel(fig, path, steps, plates, step, title, arrow_spacing=10.0,
                                   label_font='6p', annotate=True)
    _map_key(fig)
    fig.savefig(os.path.join(out, 'sequence.png'), dpi=200)


def saved_steps(dirs):
    """Steps saved in every one of the output folders."""
    return sorted(set.intersection(*[
        {int(f[9:13]) for f in os.listdir(d) if f.startswith('plate_id_') and f.endswith('.nc')}
        for d in dirs]))


def pick_sequence_steps(saved, n=12):
    """About n saved steps spread evenly over the run, always including the first and last."""
    if len(saved) <= n:
        return saved
    idx = np.unique(np.round(np.linspace(0, len(saved) - 1, n)).astype(int))
    return [saved[i] for i in idx]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('dirs', nargs='+')
    parser.add_argument('--labels', nargs='+')
    parser.add_argument('--out')
    parser.add_argument('--steps', nargs='+', type=int,
                        help='steps to map (default: first, middle and last saved step)')
    parser.add_argument('--sequence-steps', nargs='+', type=int,
                        help='steps for sequence.png (default: about --sequence-n saved steps spread over the run)')
    parser.add_argument('--sequence-n', type=int, default=12)
    parser.add_argument('--every', type=int, default=5, help='sample every Nth step in fu_plot')
    parser.add_argument('--min-area-frac', type=float, default=0.002,
                        help='leave plates smaller than this fraction of the sphere out of fu_plot')
    args = parser.parse_args()

    labels = args.labels or [os.path.basename(os.path.normpath(d)) for d in args.dirs]
    out = args.out or args.dirs[0]
    os.makedirs(out, exist_ok=True)
    runs = [load_run(d) for d in args.dirs]

    fu_plot(runs, labels, out, args.every, args.min_area_frac)
    timeseries(runs, labels, out)

    saved = saved_steps(args.dirs)
    seq_steps = args.sequence_steps or pick_sequence_steps(saved, args.sequence_n)
    if seq_steps:
        sequence(args.dirs, runs, labels, out, seq_steps)
    steps_wanted = args.steps
    if steps_wanted is None:
        steps_wanted = sorted({saved[0], saved[len(saved) // 2], saved[-1]}) if saved else []
    maps(args.dirs, runs, labels, out, steps_wanted)
    print('Figures written to', out)


if __name__ == '__main__':
    main()
