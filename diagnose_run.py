"""Diagnostic figures for one or more PlateDriver runs.

    python diagnose_run.py <output_dir> [<output_dir> ...] [--labels A B] [--out figs]
                           [--steps 0 50 100] [--every 5] [--min-area-frac 0.002]

Reads step_stats.csv and plate_stats.csv (written every step by run_SDP.py) and the
saved grids, and writes three PyGMT figures to --out (default: the first output_dir):

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
    lons = np.arange(-180 + spacing_deg / 2, 180, spacing_deg)
    lats = np.arange(-84 + spacing_deg / 2, 84, spacing_deg)
    lon2d, lat2d = np.meshgrid(lons, lats)
    lon, lat = lon2d.ravel(), lat2d.ravel()
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


def maps(paths, runs, labels, out, steps_wanted):
    for step in steps_wanted:
        fig = pygmt.Figure()
        n = len(paths)
        with fig.subplot(nrows=n, ncols=1, figsize=('18c', '{}c'.format(10 * n)),
                         margins=['0.2c', '0.6c'], title='Step {}'.format(step)):
            for k, (path, (steps, plates), label) in enumerate(zip(paths, runs, labels)):
                with fig.set_panel(panel=k):
                    fig.basemap(region='d', projection='N?', frame=['f', '+t' + label])
                    age_file = os.path.join(path, 'seafloor_age_{:04d}.nc'.format(step))
                    if os.path.exists(age_file):
                        with pygmt.config(COLOR_NAN='gray75'):   # continents
                            pygmt.makecpt(cmap='batlow', series=[0, 150], reverse=True)
                            fig.grdimage(age_file, cmap=True)
                    bfile = os.path.join(path, 'boundaries_{:04d}.geojson'.format(step))
                    if os.path.exists(bfile):
                        b = gpd.read_file(bfile)
                        for kind, pen in [('conservative', '0.8p,black'), ('ridge', '1.5p,red'),
                                          ('subduction', '1.5p,blue')]:
                            sel = b[b.boundary_type == kind]
                            if len(sel):
                                fig.plot(data=sel, pen=pen)
                    pid_file = os.path.join(path, 'plate_id_{:04d}.nc'.format(step))
                    rows = plates[plates.step == step]
                    if os.path.exists(pid_file) and len(rows):
                        plate_id = xr.open_dataarray(pid_file)
                        x, y, az, sp = _velocity_arrows(plate_id, rows)
                        fig.plot(x=x, y=y, direction=[az, sp / 250.0],
                                 style='V0.12c+e+a45+n', pen='0.6p,black', fill='black')
                        for lon, lat, text in _plate_labels(plate_id):
                            fig.text(x=lon, y=lat, text=text, font='9p,Helvetica-Bold,white',
                                     fill='black@40', clearance='0.05c')
                    srow = steps[steps.step == step]
                    if len(srow):
                        fig.text(position='BL', offset='0.1c/0.1c', font='8p,Helvetica,black',
                                 fill='white@20', justify='BL',
                                 text='RMS {:.0f} km/Myr, net rot {:.2f} deg/Myr; arrow 1 cm = 250 km/Myr'.format(
                                     float(srow.rms_speed_km_per_myr.iloc[0]),
                                     float(srow.net_rotation_deg_per_myr.iloc[0])))
        if n == 1:
            pygmt.makecpt(cmap='batlow', series=[0, 150], reverse=True)
            fig.colorbar(frame='xaf+lSeafloor age (Myr)', position='JBC+w8c/0.3c+o0/0.6c')
        fig.savefig(os.path.join(out, 'maps_{:04d}.png'.format(step)), dpi=200)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('dirs', nargs='+')
    parser.add_argument('--labels', nargs='+')
    parser.add_argument('--out')
    parser.add_argument('--steps', nargs='+', type=int,
                        help='steps to map (default: first, middle and last saved step)')
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

    steps_wanted = args.steps
    if steps_wanted is None:
        saved = sorted(set.intersection(*[
            {int(f[9:13]) for f in os.listdir(d) if f.startswith('plate_id_') and f.endswith('.nc')}
            for d in args.dirs]))
        steps_wanted = sorted({saved[0], saved[len(saved) // 2], saved[-1]}) if saved else []
    maps(args.dirs, runs, labels, out, steps_wanted)
    print('Figures written to', out)


if __name__ == '__main__':
    main()
