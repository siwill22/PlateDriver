import argparse
import os
import time

import pygplates
import numpy as np
import pandas as pd
import pygmt
import geopandas as gpd
import yaml

from gprm import ReconstructionModel
from gprm.utils.create_gpml import gpml2gdf

import plate_driver as pdriver
import self_driving_plates as sdp


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Run self-driving plates simulation.')
parser.add_argument('--config', default='config.yaml',
                    help='Path to YAML config file (default: config.yaml)')
args = parser.parse_args()

with open(args.config) as fh:
    cfg = yaml.safe_load(fh)

files       = cfg['input_files']
output_dir  = cfg.get('output_dir', './output')
run_cfg     = cfg['run']
forces      = cfg['forces']
params      = cfg['params']

reconstruction_time = run_cfg['reconstruction_time']
delta_time          = run_cfg['delta_time']
total_time          = run_cfg['total_time']

os.makedirs(output_dir, exist_ok=True)

def _out(filename):
    return os.path.join(output_dir, filename)


# ---------------------------------------------------------------------------
# Build reconstruction model
# ---------------------------------------------------------------------------

reconstruction_model = ReconstructionModel('')
reconstruction_model.add_rotation_model(files['rotation_model'])
for poly_file in files['dynamic_polygons']:
    reconstruction_model.add_dynamic_polygons(poly_file)
cratons = pygplates.FeatureCollection(files['cratons'])


# ---------------------------------------------------------------------------
# Initial timestep
# ---------------------------------------------------------------------------

snapshot = reconstruction_model.plate_snapshot(reconstruction_time)
result   = pdriver.get_torques_for_plate_snapshot(snapshot, reconstruction_time, delta_time, forces)

stage_pole_dict = {}
for res in result:
    stage_pole_dict[res[0]] = res[3]
    print(res[0], res[3])

cratons = sdp.nudge_polygons(cratons, stage_pole_dict)

rpolygons = sdp.nudge_topological_polygons(snapshot,
                                           forces=forces,
                                           reconstruction_time=reconstruction_time,
                                           delta_time=delta_time)

polarity_dict         = {}
cumulative_rotation_dict = {}
rotation_records      = []

(new_plate_polygons,
 stage_pole_dict,
 rasters,
 boundaries_gdf) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time,
                                         params=params, polarity_dict=polarity_dict)

tracers = sdp.generate_seafloor_age_tracers(
    new_plate_polygons, boundaries_gdf,
    spacing_deg=params['seafloor_spacing_deg'],
    spreading_rate_km_per_myr=params['spreading_rate_km_per_myr'],
    div_threshold=params['div_threshold_ridge'])

hotspot_locs = sdp.generate_hotspot_locations(n=params['n_hotspots'])
seamounts = gpd.GeoDataFrame(columns=['geometry', 'hotspot_id', 'age', 'active'],
                             geometry='geometry').set_crs(epsg=4326)


# ---------------------------------------------------------------------------
# Main time loop
# ---------------------------------------------------------------------------

for i in range(1, int(total_time / delta_time) + 1):
    t = i * delta_time

    print('t = {:.1f} Myr (step {})'.format(t, i))
    _step_times = {}
    _t0 = time.perf_counter()

    rpolygons = sdp.nudge_polygons(
        new_plate_polygons.explode(index_parts=False).reset_index(drop=True),
        stage_pole_dict)
    cratons = sdp.nudge_polygons(cratons, stage_pole_dict)
    _step_times['nudge_polygons'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    (new_plate_polygons,
     stage_pole_dict,
     rasters,
     boundaries_gdf) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time,
                                             cratons=cratons, params=params,
                                             polarity_dict=polarity_dict)
    _step_times['compute_timestep'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    for key, item in stage_pole_dict.items():
        print(key, item.get_lat_lon_euler_pole_and_angle_degrees())

    # Accumulate finite rotations: F(t) = S_t * F(t-dt)
    for pid, stage_rot in stage_pole_dict.items():
        prev = cumulative_rotation_dict.get(pid, pygplates.FiniteRotation())
        cumulative_rotation_dict[pid] = stage_rot * prev
        rotation_records.append((t, pid, cumulative_rotation_dict[pid]))
    _step_times['rotations'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    new_tracers = sdp.seed_ridge_tracers(boundaries_gdf,
                                         spacing_deg=params['tracer_spacing_deg'],
                                         threshold=params['div_threshold_ridge'],
                                         offset_deg=params['tracer_offset_deg'])
    tracers = pd.concat([tracers, new_tracers], ignore_index=True)
    tracers = sdp.advect_tracers(tracers, rasters['plate_id'], stage_pole_dict, delta_time,
                                  boundaries_gdf=boundaries_gdf,
                                  div_threshold=params['div_threshold_subduction'])
    tracers = sdp.deactivate_subducted_tracers(tracers, boundaries_gdf,
                                               threshold=params['div_threshold_subduction'])
    _step_times['tracers'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    seamounts = sdp.advect_tracers(seamounts, rasters['plate_id'], stage_pole_dict, delta_time,
                                   boundaries_gdf=boundaries_gdf,
                                   div_threshold=params['div_threshold_subduction'])
    seamounts = sdp.deactivate_subducted_tracers(seamounts, boundaries_gdf,
                                                  threshold=params['div_threshold_subduction'])
    new_seamounts = sdp.seed_hotspot_seamounts(hotspot_locs, prob=params['hotspot_prob'])
    if len(new_seamounts):
        seamounts = gpd.GeoDataFrame(
            pd.concat([seamounts, new_seamounts], ignore_index=True),
            geometry='geometry').set_crs(epsg=4326)
    _step_times['seamounts'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    active_tracers = tracers[tracers.active.astype(bool)]

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    _plot_interval = params.get('plot_interval')
    if _plot_interval is not None and i % _plot_interval == 0:
        _tp0 = time.perf_counter()
        fig = pygmt.Figure()
        fig.basemap(frame='f', projection='W10i', region='d')
        if len(active_tracers):
            vmax = 100
            pygmt.makecpt(cmap='cubhelix', series=[0, vmax], reverse=True)
            fig.plot(x=active_tracers.geometry.x.values, y=active_tracers.geometry.y.values,
                     style='c0.08c', fill=active_tracers.age.values.astype(float),
                     cmap=True, pen='0.1p,black')
        pygmt.makecpt(cmap='polar', series=[-0.001, 0.001], background='o')
        fig.plot(data=new_plate_polygons.iloc[:], pen='1p,black', transparency=50)
        if 'divergence' in boundaries_gdf.columns:
            ridge_bounds = boundaries_gdf[boundaries_gdf['divergence'] > params['div_threshold_ridge']]
            sub_bounds   = boundaries_gdf[boundaries_gdf['divergence'] < params['div_threshold_subduction']]
            cons_bounds  = boundaries_gdf[
                (boundaries_gdf['divergence'] >= params['div_threshold_subduction']) &
                (boundaries_gdf['divergence'] <= params['div_threshold_ridge'])
            ]
            if len(cons_bounds):
                fig.plot(data=cons_bounds, pen='1p,black')
            if len(ridge_bounds):
                fig.plot(data=ridge_bounds, pen='2p,red')
            if len(sub_bounds):
                fig.plot(data=sub_bounds, pen='2p,blue')
        glon, glat = np.meshgrid(rasters.lon, rasters.lat)
        stride = 10
        magnitude = np.sqrt(rasters['velocity_east'].data[::stride, ::stride].flatten()**2 +
                            rasters['velocity_north'].data[::stride, ::stride].flatten()**2)
        azimuth = np.degrees(np.arctan2(rasters['velocity_east'].data[::stride, ::stride].flatten(),
                                        rasters['velocity_north'].data[::stride, ::stride].flatten())) % 360
        fig.plot(data={'x': glon[::stride, ::stride].flatten(),
                       'y': glat[::stride, ::stride].flatten(),
                       'azimuth': azimuth,
                       'magnitude': magnitude / (2000. * delta_time)},
                 style='V0.15c+e+a45+ggray40+n', pen='0.1p,black', fill='black')
        fig.plot(data=gpml2gdf(cratons), pen='1p,black', fill='grey', transparency=50)
        label_pts = new_plate_polygons.representative_point()
        fig.text(x=label_pts.x, y=label_pts.y,
                 text=new_plate_polygons['PLATEID1'], font='10p,Helvetica-Bold,black', justify='CM')
        fig.savefig(_out('plot_{:04d}.png'.format(i)))
        _step_times['plot'] = time.perf_counter() - _tp0
    _t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # Periodic saves
    # ------------------------------------------------------------------
    if i % params['csv_save_interval'] == 0:
        active_tracers[['geometry', 'age']].assign(
            lon=active_tracers.geometry.x,
            lat=active_tracers.geometry.y
        ).drop(columns='geometry').to_csv(_out('tracers_{:04d}.csv'.format(i)), index=False)

        age_grid = None
        if len(active_tracers):
            age_grid = sdp.generate_seafloor_age_grid(
                _out('tracers_{:04d}.csv'.format(i)),
                craton_raster=rasters.get('craton'),
                spacing=params['raster_spacing'])
            age_grid.to_netcdf(_out('seafloor_age_{:04d}.nc'.format(i)))

        rasters['divergence'].to_netcdf(_out('divergence_{:04d}.nc'.format(i)))
        rasters['plate_id'].to_netcdf(_out('plate_id_{:04d}.nc'.format(i)))
        if 'craton' in rasters:
            rasters['craton'].to_netcdf(_out('craton_{:04d}.nc'.format(i)))

        craton_counts = sdp.craton_counts_from_raster(
            rasters['craton'] if 'craton' in rasters else None)

        if 'divergence' in boundaries_gdf.columns:
            bdf = boundaries_gdf[['geometry', 'plateid1_poly1', 'plateid1_poly2',
                                   'divergence']].copy()
            def _btype(d):
                if d > params['div_threshold_ridge']:      return 'ridge'
                if d < params['div_threshold_subduction']: return 'subduction'
                return 'conservative'
            bdf['boundary_type'] = bdf['divergence'].map(_btype)

            downgoing, overriding = [], []
            for _, br in bdf.iterrows():
                if br['boundary_type'] == 'subduction':
                    p1 = int(br['plateid1_poly1']) if pd.notna(br['plateid1_poly1']) else 0
                    p2 = int(br['plateid1_poly2']) if pd.notna(br['plateid1_poly2']) else 0
                    dg = sdp.resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts)
                    ov = p2 if dg == p1 else p1
                else:
                    dg, ov = 0, 0
                downgoing.append(dg)
                overriding.append(ov)
            bdf['downgoing_plate']  = downgoing
            bdf['overriding_plate'] = overriding
            bdf.to_file(_out('boundaries_{:04d}.geojson'.format(i)), driver='GeoJSON')

        earthquakes = sdp.generate_synthetic_earthquakes(
            boundaries_gdf, polarity_dict, new_plate_polygons,
            craton_counts=craton_counts,
            plate_id_raster=rasters['plate_id'],
            n_candidates=params['eq_n_candidates'],
            ridge_max_dist_deg=params['eq_ridge_max_dist_deg'],
            ridge_max_depth_km=params['eq_ridge_max_depth_km'],
            slab_dip_deg=params['eq_slab_dip_deg'],
            slab_max_depth_km=params['eq_slab_max_depth_km'],
            div_threshold_ridge=params['div_threshold_ridge'],
            div_threshold_subduction=params['div_threshold_subduction'],
            ridge_eq_prob=params['eq_ridge_prob'],
            subduction_eq_prob=params['eq_subduction_prob'])
        earthquakes.to_file(_out('earthquakes_{:04d}.geojson'.format(i)), driver='GeoJSON')

        topo = None
        if 'craton' in rasters:
            topo = sdp.compute_topography_raster(
                rasters['plate_id'], rasters['craton'], boundaries_gdf,
                polarity_dict=polarity_dict, craton_counts=craton_counts,
                continental_base_m=params['topo_continental_base_m'],
                noise_fraction=params['topo_noise_fraction'],
                spectral_beta=params['topo_spectral_beta'],
                andean_max_dist_deg=params['topo_andean_max_dist_deg'],
                andean_arc_dist_deg=params['topo_andean_arc_dist_deg'],
                andean_peak_m=params['topo_andean_peak_m'],
                himalayan_max_dist_deg=params['topo_himalayan_max_dist_deg'],
                himalayan_peak_m=params['topo_himalayan_peak_m'],
                div_threshold_subduction=params['div_threshold_subduction'])
            topo.to_netcdf(_out('topography_{:04d}.nc'.format(i)))

        if age_grid is not None:
            topo_bathy = sdp.generate_topo_bathy_grid(age_grid, topography_raster=topo)
            topo_bathy.to_netcdf(_out('topo_bathy_{:04d}.nc'.format(i)))

        sdp.write_rotation_file(rotation_records, _out('rotations_{:04d}.rot'.format(i)))

        _active_seamounts = seamounts[seamounts.active.astype(bool)]
        if len(_active_seamounts):
            _active_seamounts[['geometry', 'hotspot_id', 'age']].to_file(
                _out('seamounts_{:04d}.geojson'.format(i)), driver='GeoJSON')

        _step_times['save_outputs'] = time.perf_counter() - _t0

    total = sum(_step_times.values())
    print('  --- step timing (step {:d}, {:.1f}s total) ---'.format(i, total))
    for _name, _dt in sorted(_step_times.items(), key=lambda x: -x[1]):
        print('  {:20s}  {:5.1f}s  ({:4.1f}%)'.format(_name, _dt, 100 * _dt / total))


# ---------------------------------------------------------------------------
# Final saves (captures last step if not on a save interval)
# ---------------------------------------------------------------------------

active_seamounts = seamounts[seamounts.active.astype(bool)]
if len(active_seamounts):
    active_seamounts[['geometry', 'hotspot_id', 'age']].to_file(
        _out('seamounts_final.geojson'), driver='GeoJSON')
