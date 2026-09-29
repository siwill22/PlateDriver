import argparse
import os
import queue
import sys
import threading
import time
import traceback

import pygplates
import numpy as np
import pandas as pd
import pygmt
import geopandas as gpd
import yaml

from gprm import ReconstructionModel
from gprm.utils.create_gpml import gpml2gdf

import force_modes
import plate_driver as pdriver
import self_driving_plates as sdp
import step_stats


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


# Force model: a preset (simple | physical) plus per-component overrides
force_preset, modes = force_modes.resolve_force_modes(cfg)
print(force_modes.describe(force_preset, modes))
params['collision_mode'] = modes['collision']
params['drag_mode'] = modes['drag']
params['speed_normalisation'] = modes['speed_normalisation']
params['slab_pull_mode'] = modes['slab_pull']
params['slab_suction_mode'] = modes['slab_suction']
params['polarity_mode'] = modes['polarity']
params['ridge_push_mode'] = modes['ridge_push']
params['tracer_thinning_mode'] = modes['tracer_thinning']
params['plate_rebuild_mode'] = modes['plate_rebuild']
params['rifting_mode'] = modes['rifting']
params['collision_weld_mode'] = modes['collision_weld']
params['sutures_mode'] = modes['sutures']
params['rift_path_mode'] = modes['rift_path']
if modes['collision_weld'] == 'on' and modes['collision'] not in ('coupling', 'breakoff'):
    raise SystemExit('collision_weld: on needs collision: coupling or breakoff (they time the collisions)')
if modes['sutures'] == 'tracers' and modes['collision_weld'] != 'on':
    raise SystemExit('sutures: tracers are dropped at welds, so they need collision_weld: on')
_needs = sdp.force_fields_needed(params)     # per-step fields the force modes read


def _git_commit():
    try:
        import subprocess
        here = os.path.dirname(os.path.abspath(__file__))
        commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=here,
                                capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'],
                               cwd=here, capture_output=True, text=True).stdout.strip()
        return commit + ('-dirty' if dirty else '') if commit else 'unknown'
    except Exception:
        return 'unknown'


# Record exactly what this run used (rewritten after step 0 with the speed scale)
run_info = {'config_file': os.path.abspath(args.config),
            'git_commit': _git_commit(),
            'started': time.strftime('%Y-%m-%d %H:%M:%S'),
            'force_model': force_preset,
            'force_modes_resolved': modes,
            'config': cfg}


def _write_run_info():
    with open(_out('run_info.yaml'), 'w') as fh:
        yaml.safe_dump(run_info, fh, sort_keys=False)


_write_run_info()

# Per-step statistics are appended every step; start them afresh
_step_stats_csv  = _out('step_stats.csv')
_plate_stats_csv = _out('plate_stats.csv')
for _f in (_step_stats_csv, _plate_stats_csv):
    if os.path.exists(_f):
        os.remove(_f)


# All randomness (hotspots, seamounts, earthquakes, topography noise)
# comes from one generator, so a run is repeatable when random_seed is set.
rng = np.random.default_rng(params.get('random_seed'))
topo_noise_seed = int(rng.integers(2**31 - 1))   # same noise pattern for every saved frame

_missing = [f for f in [files['rotation_model'], files['cratons'], *files['dynamic_polygons']]
            if not os.path.exists(f)]
if _missing:
    raise SystemExit('Input file(s) listed in {} not found:\n  {}'.format(
        args.config, '\n  '.join(_missing)))


# ---------------------------------------------------------------------------
# Build reconstruction model
# ---------------------------------------------------------------------------

reconstruction_model = ReconstructionModel('')
reconstruction_model.add_rotation_model(files['rotation_model'])
for poly_file in files['dynamic_polygons']:
    reconstruction_model.add_dynamic_polygons(poly_file)
cratons = pygplates.FeatureCollection(files['cratons'])
if reconstruction_time != 0:
    # Cratons are stored at present-day positions; move them to the starting snapshot
    _rotation_model = pygplates.RotationModel(files['rotation_model'])
    cratons = sdp.nudge_polygons(cratons, {
        c.get_reconstruction_plate_id(): _rotation_model.get_rotation(
            reconstruction_time, c.get_reconstruction_plate_id())
        for c in cratons})


# ---------------------------------------------------------------------------
# Initial state (step 0, t = 0)
# ---------------------------------------------------------------------------
#
# Time convention: step i is the state at t = i * delta_time Myr after the starting
# snapshot.  Each step first moves everything (plates, cratons, tracers, seamounts)
# with the stage rotations computed at the previous step, then rebuilds the plate
# configuration and computes the stage rotations for the next move.

snapshot = reconstruction_model.plate_snapshot(reconstruction_time)
if not snapshot.resolved_topologies:
    raise SystemExit('No plate polygons resolve at reconstruction_time = {} Ma; check the '
                     'dynamic_polygons files and their valid time range.'.format(reconstruction_time))

# Torques from the GPlates boundary feature types (subduction zones, ridges) give a
# first set of stage rotations.  They are only used to classify the plate boundaries
# of the initial configuration; the plates are not moved with them.
result = pdriver.get_torques_for_plate_snapshot(snapshot, reconstruction_time, delta_time, forces)
seed_stage_poles = {res[0]: res[3] for res in result}

rpolygons = sdp.nudge_topological_polygons(snapshot)

polarity_dict            = {}
cumulative_rotation_dict = {}   # {pid: FiniteRotation from t=0 to the current step}
rotation_records         = []   # [(t, pid, FiniteRotation)] written to rotations_NNNN.rot
collision_state          = {}   # tracks collision age per plate pair; persists across timesteps
speed_state              = {}   # speed scale k (calibrate_once) and last step's scaling

# Trust the user's initial plate configuration — skip fragmentation check on step 0.
_params_step0 = {**params, 'fragmentation_enabled': False}

_slab_enabled = params.get('slab_tracers_enabled', False)   # mantle outputs (CMB heat, mantle_temp)
_slab_track = _slab_enabled or _needs['slab']                # slab tracers kept (outputs or forces)
_slab_columns = ['lon', 'lat', 'depth_km', 'age', 'active', 'sink_rate_km_per_myr', 'weight']
slab_tracers = pd.DataFrame(columns=_slab_columns) if _slab_track else None


def _new_seafloor_tracers(plate_polygons, bdry):
    return sdp.generate_seafloor_age_tracers(
        plate_polygons, bdry,
        spacing_deg=params['seafloor_spacing_deg'],
        spreading_rate_km_per_myr=params['spreading_rate_km_per_myr'],
        div_threshold=params['div_threshold_ridge'])


_step0_fields = None
if _needs['age'] or _needs['slab']:
    # The force modes read seafloor age and slab depth, which are seeded from the
    # step-0 boundaries.  A first, simple pass gives those boundaries (the plate
    # geometry does not depend on the force modes); the forces come from the second.
    _params_pass1 = {**_params_step0, 'drag_mode': 'none', 'speed_normalisation': 'none',
                     'speed_scale': None, 'slab_pull_mode': 'hard', 'slab_suction_mode': 'ratio',
                     'polarity_mode': 'id', 'ridge_push_mode': 'boundary'}
    _pass1 = sdp.compute_timestep(rpolygons, forces, seed_stage_poles, delta_time,
                                  cratons=cratons, params=_params_pass1, polarity_dict={},
                                  collision_state={}, speed_state={})
    tracers = _new_seafloor_tracers(_pass1[0], _pass1[3])
    _step0_fields = {'tracers': (tracers.geometry.x.values, tracers.geometry.y.values,
                                 tracers.age.values.astype(float), None)}
    if _needs['slab']:
        slab_tracers = sdp.seed_initial_slabs(
            _pass1[3], params.get('initial_slab_depth_km', 1320.0), delta_time,
            div_threshold_subduction=params['div_threshold_subduction'],
            slowing_factor=params.get('slab_slowing_factor', 4),
            max_sink_rate_km_per_myr=params.get('slab_max_sink_rate', 200.0),
            sink_rate_factor=params.get('slab_sink_rate_factor', 1.0))
        _step0_fields['slab_tracers'] = slab_tracers
        print('  [slab] {} initial slab tracers at the step-0 trenches, to {:.0f} km'.format(
            len(slab_tracers), params.get('initial_slab_depth_km', 1320.0)))
    del _pass1

force_diag = {}
(new_plate_polygons,
 stage_pole_dict,
 rasters,
 boundaries_gdf,
 _,
 collision_state,
 next_plate_id) = sdp.compute_timestep(rpolygons, forces, seed_stage_poles, delta_time,
                                       cratons=cratons, params=_params_step0,
                                       polarity_dict=polarity_dict,
                                       collision_state=collision_state,
                                       speed_state=speed_state,
                                       force_fields=_step0_fields,
                                       force_diag=force_diag)
if 'last' in speed_state:
    run_info['speed_scale_k'] = float(speed_state['last']['k_used'])
    run_info['rms_speed_unscaled_step0_km_per_myr'] = float(speed_state['last']['rms_unscaled'])
    _write_run_info()

if _step0_fields is None:
    tracers = _new_seafloor_tracers(new_plate_polygons, boundaries_gdf)

_hs_locs = params.get('hotspot_locations', None)
hotspot_locs = sdp.generate_hotspot_locations(n=params['n_hotspots'], locations=_hs_locs, seed=rng)
seamounts = gpd.GeoDataFrame(columns=['geometry', 'hotspot_id', 'age', 'active'],
                             geometry='geometry').set_crs(epsg=4326)

# Fragmentation / merging state (persists across timesteps)
_conservative_contact   = {}   # {(min_pid, max_pid): int consecutive near-zero steps}
_plate_cooldown         = {}   # {pid: float remaining_myr before eligible again}
_merge_pairs            = []   # [(pid_large, pid_small)] passed to next compute_timestep
_fragmentation_cooldown_pids = set()
_rift_state    = {}   # {pid: {'count', 'normal'}} plates whose best cut has passed the trial split
_rift_planned  = []   # [(pid, plane normal)] splits to make next step
_rift_cooldown = {}   # {pid: remaining_myr} plates just split, not tested again yet
_rift_events   = []
_welds         = []   # collision_weld: on, one row per weld
_sutures       = pd.DataFrame(columns=['lon', 'lat', 'suture_id', 'time_myr'])   # sutures: tracers
_next_suture_id = 1
_craton_plate_remap     = {}   # {craton_plateid: host_plateid} for absorbed cratons
if _slab_enabled:
    _mantle_spacing = params.get('mantle_grid_spacing_deg', 5.0)
    _cmb_nlat = len(np.arange(-90  + _mantle_spacing / 2, 90,  _mantle_spacing))
    _cmb_nlon = len(np.arange(-180 + _mantle_spacing / 2, 180, _mantle_spacing))
    cmb_heat = np.zeros((_cmb_nlat, _cmb_nlon))


# ---------------------------------------------------------------------------
# Async output functions
# ---------------------------------------------------------------------------

def _do_plot(i, active_tracers, new_plate_polygons, boundaries_gdf,
             rasters, seamounts, cratons_gdf, params, output_dir):
    fig = pygmt.Figure()
    fig.basemap(frame='f', projection='W10i', region='d')
    if len(active_tracers):
        vmax = 100
        pygmt.makecpt(cmap='cubhelix', series=[0, vmax], reverse=True)
        fig.plot(x=active_tracers.geometry.x.values, y=active_tracers.geometry.y.values,
                 style='c0.08c', fill=active_tracers.age.values.astype(float),
                 cmap=True, pen='0.1p,black')
    if 'divergence' in boundaries_gdf.columns:
        ridge_bounds = boundaries_gdf[boundaries_gdf['divergence'] > params['div_threshold_ridge']]
        sub_bounds   = boundaries_gdf[boundaries_gdf['divergence'] < params['div_threshold_subduction']]
        cons_bounds  = boundaries_gdf[
            (boundaries_gdf['divergence'] >= params['div_threshold_subduction']) &
            (boundaries_gdf['divergence'] <= params['div_threshold_ridge'])
        ]
        if len(cons_bounds):
            fig.plot(data=cons_bounds.copy(), pen='1p,black')
        if len(ridge_bounds):
            fig.plot(data=ridge_bounds.copy(), pen='2p,red')
        if len(sub_bounds):
            fig.plot(data=sub_bounds.copy(), pen='2p,blue')
    glon, glat = np.meshgrid(rasters.lon, rasters.lat)
    stride = 10
    magnitude = np.sqrt(rasters['velocity_east'].data[::stride, ::stride].flatten()**2 +
                        rasters['velocity_north'].data[::stride, ::stride].flatten()**2)
    azimuth = np.degrees(np.arctan2(rasters['velocity_east'].data[::stride, ::stride].flatten(),
                                    rasters['velocity_north'].data[::stride, ::stride].flatten())) % 360
    fig.plot(data={'x': glon[::stride, ::stride].flatten(),
                   'y': glat[::stride, ::stride].flatten(),
                   'azimuth': azimuth,
                   'magnitude': magnitude / 2000.},     # km/Myr → arrow length (cm)
             style='V0.15c+e+a45+ggray40+n', pen='0.1p,black', fill='black')
    _active_sm = seamounts[seamounts.active.astype(bool)]
    if len(_active_sm):
        fig.plot(x=_active_sm.geometry.x.values, y=_active_sm.geometry.y.values,
                 style='t0.2c', fill='red', pen='0.5p,black')
    fig.plot(data=cratons_gdf, pen='1p,black', fill='grey', transparency=50)
    label_pts = new_plate_polygons.representative_point()
    fig.text(x=label_pts.x, y=label_pts.y,
             text=new_plate_polygons['PLATEID1'], font='10p,Helvetica-Bold,black', justify='CM')
    fig.savefig(os.path.join(output_dir, 'plot_{:04d}.png'.format(i)))


def _do_save(i, active_tracers, rasters, boundaries_gdf, polarity_dict,
             new_plate_polygons, seamounts, rotation_records, params, output_dir,
             slab_tracers=None, cmb_heat=None, noise_seed=None):
    """Write all outputs for step i.  Each output is written independently, so one
    failure does not prevent the others; failures are reported with a traceback."""
    def _out(filename):
        return os.path.join(output_dir, filename.format(i))

    def _try(name, fn):
        try:
            return fn()
        except Exception:
            _report_error('save step {} ({})'.format(i, name))
            return None

    def _tracers():
        active_tracers[['geometry', 'age']].assign(
            lon=active_tracers.geometry.x,
            lat=active_tracers.geometry.y
        ).drop(columns='geometry').to_csv(_out('tracers_{:04d}.csv'), index=False)
        return True

    def _age_grid():
        grid = sdp.generate_seafloor_age_grid(
            _out('tracers_{:04d}.csv'),
            craton_raster=rasters.get('craton'),
            spacing=params['raster_spacing'])
        grid.to_netcdf(_out('seafloor_age_{:04d}.nc'))
        return grid

    wrote_tracers = _try('tracers', _tracers)
    age_grid = None
    if wrote_tracers and len(active_tracers):
        age_grid = _try('seafloor age grid', _age_grid)

    _try('divergence', lambda: rasters['divergence'].to_netcdf(_out('divergence_{:04d}.nc')))
    _try('plate_id', lambda: rasters['plate_id'].to_netcdf(_out('plate_id_{:04d}.nc')))
    if 'craton' in rasters:
        _try('craton', lambda: rasters['craton'].to_netcdf(_out('craton_{:04d}.nc')))

    craton_counts = sdp.craton_counts_from_raster(
        rasters['craton'] if 'craton' in rasters else None)
    _pol_counts = sdp.polarity_counts(craton_counts, params)   # counts that decide polarity

    def _boundaries():
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
                dg = sdp.resolve_downgoing_plate(p1, p2, polarity_dict, _pol_counts)
                ov = p2 if dg == p1 else p1
            else:
                dg, ov = 0, 0
            downgoing.append(dg)
            overriding.append(ov)
        bdf['downgoing_plate']  = downgoing
        bdf['overriding_plate'] = overriding
        bdf.to_file(_out('boundaries_{:04d}.geojson'), driver='GeoJSON')

    if 'divergence' in boundaries_gdf.columns and len(boundaries_gdf):
        _try('boundaries', _boundaries)

    def _earthquakes():
        earthquakes = sdp.generate_synthetic_earthquakes(
            boundaries_gdf, polarity_dict, new_plate_polygons,
            craton_counts=_pol_counts,
            plate_id_raster=rasters['plate_id'],
            n_candidates=params['eq_n_candidates'],
            ridge_max_dist_deg=params['eq_ridge_max_dist_deg'],
            ridge_max_depth_km=params['eq_ridge_max_depth_km'],
            slab_dip_deg=params['eq_slab_dip_deg'],
            slab_max_depth_km=params['eq_slab_max_depth_km'],
            div_threshold_ridge=params['div_threshold_ridge'],
            div_threshold_subduction=params['div_threshold_subduction'],
            ridge_eq_prob=params['eq_ridge_prob'],
            subduction_eq_prob=params['eq_subduction_prob'],
            rng=np.random.default_rng(None if noise_seed is None else noise_seed + i))
        earthquakes.to_file(_out('earthquakes_{:04d}.geojson'), driver='GeoJSON')

    _try('earthquakes', _earthquakes)

    def _topography():
        grid = sdp.compute_topography_raster(
            rasters['plate_id'], rasters['craton'], boundaries_gdf,
            polarity_dict=polarity_dict, craton_counts=craton_counts,
            polarity_craton_counts=_pol_counts,
            continental_base_m=params['topo_continental_base_m'],
            noise_fraction=params['topo_noise_fraction'],
            spectral_beta=params['topo_spectral_beta'],
            andean_max_dist_deg=params['topo_andean_max_dist_deg'],
            andean_arc_dist_deg=params['topo_andean_arc_dist_deg'],
            andean_peak_m=params['topo_andean_peak_m'],
            himalayan_max_dist_deg=params['topo_himalayan_max_dist_deg'],
            himalayan_peak_m=params['topo_himalayan_peak_m'],
            div_threshold_subduction=params['div_threshold_subduction'],
            noise_seed=noise_seed)
        grid.to_netcdf(_out('topography_{:04d}.nc'))
        return grid

    topo = _try('topography', _topography) if 'craton' in rasters else None

    if age_grid is not None:
        _try('topo_bathy', lambda: sdp.generate_topo_bathy_grid(
            age_grid, topography_raster=topo).to_netcdf(_out('topo_bathy_{:04d}.nc')))

    _try('rotations', lambda: sdp.write_rotation_file(rotation_records, _out('rotations_{:04d}.rot')))

    _active_seamounts = seamounts[seamounts.active.astype(bool)]
    if len(_active_seamounts):
        _try('seamounts', lambda: _active_seamounts[['geometry', 'hotspot_id', 'age']].to_file(
            _out('seamounts_{:04d}.geojson'), driver='GeoJSON'))

    if slab_tracers is not None and cmb_heat is not None:
        _try('mantle_temp', lambda: sdp.generate_mantle_temp_field(
            slab_tracers,
            spacing_deg=params.get('mantle_grid_spacing_deg', 5.0),
            depth_interval_km=params.get('mantle_depth_interval_km', 100.0),
            thermal_tau_myr=params.get('slab_thermal_tau_myr', 200.0),
            cmb_heat=cmb_heat,
            cmb_dVs_scale=params.get('cmb_dVs_scale', -0.5),
            cmb_decay_km=params.get('cmb_decay_km', 400.0),
            tomo_smooth_horiz_deg=params.get('tomo_smooth_horiz_deg', 0.0),
            tomo_smooth_depth_km=params.get('tomo_smooth_depth_km', 0.0)
        ).to_netcdf(_out('mantle_temp_{:04d}.nc')))


# ---------------------------------------------------------------------------
# Output worker thread
# ---------------------------------------------------------------------------

_output_queue  = queue.Queue(maxsize=3)
_output_errors = []


def _report_error(what):
    _output_errors.append(what)
    print('[output worker] ERROR in {}:\n{}'.format(what, traceback.format_exc()), flush=True)


def _output_worker():
    while True:
        item = _output_queue.get()
        if item is None:
            _output_queue.task_done()
            break
        task, payload = item
        try:
            if task == 'plot':
                _do_plot(**payload)
            elif task == 'save':
                _do_save(**payload)
        except Exception:
            _report_error('{} step {}'.format(task, payload.get('i')))
        finally:
            _output_queue.task_done()


_output_thread = threading.Thread(target=_output_worker, daemon=True)
_output_thread.start()


def _enqueue_outputs(i, force_save=False):
    """Queue the plot and/or save for step i (worker renders in the background)."""
    active_tracers = tracers[tracers.active.astype(bool)]

    _plot_interval = params.get('plot_interval')
    if _plot_interval is not None and i % _plot_interval == 0:
        _output_queue.put(('plot', dict(
            i=i,
            active_tracers=active_tracers,
            new_plate_polygons=new_plate_polygons,
            boundaries_gdf=boundaries_gdf,
            rasters=rasters,
            seamounts=seamounts,
            cratons_gdf=gpml2gdf(cratons),   # convert in main thread
            params=params,
            output_dir=output_dir,
        )))

    if force_save or i % params['csv_save_interval'] == 0:
        _save_payload = dict(
            i=i,
            active_tracers=active_tracers,
            rasters=rasters,
            boundaries_gdf=boundaries_gdf,
            polarity_dict=dict(polarity_dict),       # shallow copy — values not mutated
            new_plate_polygons=new_plate_polygons,
            seamounts=seamounts,
            rotation_records=rotation_records[:],    # list copy
            params=params,
            output_dir=output_dir,
            noise_seed=topo_noise_seed,
        )
        if _slab_enabled:
            _save_payload['slab_tracers'] = slab_tracers.copy()
            _save_payload['cmb_heat']     = cmb_heat.copy()
        _output_queue.put(('save', _save_payload))
        if force_diag.get('segments') is not None:
            force_diag['segments'].to_csv(_out('trench_forces_{:04d}.csv'.format(i)), index=False)
        if force_diag.get('ridge_push') is not None:
            force_diag['ridge_push'].to_csv(_out('ridge_push_torques_{:04d}.csv'.format(i)), index=False)
        if params['sutures_mode'] == 'tracers':
            _sutures.to_csv(_out('sutures_{:04d}.csv'.format(i)), index=False)
        return True
    return False


_step_rows = []
_plate_id_history = []   # plate-ID arrays of the last two steps, for boundary migration / jitter


def _record_step_stats(i, step_times=None, previous_poles=None):
    """Append this step's rows to step_stats.csv and plate_stats.csv."""
    step_row, plate_rows = step_stats.compute(
        i, i * delta_time, rasters, stage_pole_dict, delta_time, boundaries_gdf,
        polarity_dict, params, previous_stage_pole_dict=previous_poles)
    step_row['n_tracers'] = int(tracers.active.astype(bool).sum())
    step_row['n_slab_tracers'] = (int(slab_tracers['active'].astype(bool).sum())
                                  if _slab_track else 0)
    seg = force_diag.get('segments')
    if seg is not None and len(seg):
        # slab forces along the subduction segments (each boundary counted from both plates)
        down = seg[seg.downgoing & ~seg.collision]
        step_row['slab_pull_total'] = float((seg.pull_weight * seg.length_km).sum())
        step_row['slab_suction_total'] = float((seg.suction_weight * seg.length_km).sum())
        step_row['trench_slab_depth_mean_km'] = (float(np.average(down.slab_depth_km, weights=down.length_km))
                                                    if len(down) else np.nan)
        step_row['n_polarity_pairs'] = len(polarity_dict)
    # Boundary migration and jitter on the fixed grid: cells whose plate changed since the
    # last step, and cells that went A -> B -> A over the last three steps
    pid_now = rasters['plate_id'].values
    w = np.cos(np.radians(rasters['plate_id'].lat.values))[:, None] * np.ones_like(pid_now)
    if _plate_id_history:
        step_row['owner_change_area_frac'] = float(w[pid_now != _plate_id_history[-1]].sum() / w.sum())
    if len(_plate_id_history) == 2:
        flick = (pid_now == _plate_id_history[0]) & (pid_now != _plate_id_history[1])
        step_row['owner_flicker_area_frac'] = float(w[flick].sum() / w.sum())
    _plate_id_history.append(pid_now.copy())
    del _plate_id_history[:-2]
    last = speed_state.get('last', {})
    step_row['speed_scale_k'] = last.get('k_used', np.nan)
    step_row['rms_speed_unscaled_km_per_myr'] = last.get('rms_unscaled', np.nan)
    step_row['speed_capped'] = bool(last.get('capped', False))
    for name, dt in (step_times or {}).items():
        step_row['seconds_' + name] = dt
    step_row['seconds_step'] = sum((step_times or {}).values()) if step_times else np.nan
    # step_stats.csv is rewritten each step (columns can appear after step 0);
    # plate_stats.csv has fixed columns and is appended
    _step_rows.append(step_row)
    pd.DataFrame(_step_rows).to_csv(_step_stats_csv, index=False)
    step_stats.append_csv(_plate_stats_csv, plate_rows)
    return step_row


# ---------------------------------------------------------------------------
# Main time loop
# ---------------------------------------------------------------------------

n_steps = int(round(total_time / delta_time))
identity = pygplates.FiniteRotation()

try:
    _record_step_stats(0)
    _enqueue_outputs(0)
    for i in range(1, n_steps + 1):
        t = i * delta_time

        print('t = {:.1f} Myr (step {})'.format(t, i))
        _step_times = {}
        _t0 = time.perf_counter()

        # ------------------------------------------------------------------
        # Move plates and cratons with the stage rotations from the last step
        # ------------------------------------------------------------------
        moving_poles = stage_pole_dict
        previous_plate_id = rasters['plate_id']

        # Fields the force modes read, at this step's positions: the seafloor tracers
        # moved one step ahead (they are updated only after the forces), and the slab
        # tracers as they are (fixed in the mantle; depths lag one step)
        _fields = None
        if _needs['age'] or _needs['slab']:
            _fields = {'slab_tracers': slab_tracers}
            if _needs['age']:
                _fields['tracers'] = sdp.predict_tracers(
                    tracers.geometry.x.values, tracers.geometry.y.values,
                    tracers.age.values.astype(float), previous_plate_id, moving_poles, delta_time)
        force_diag = {}

        rpolygons = sdp.nudge_polygons(
            new_plate_polygons.explode(index_parts=False).reset_index(drop=True),
            moving_poles)
        # Absorbed cratons no longer have their own stage pole — supply the host's pole
        # so the craton feature moves with the plate it was absorbed into.
        _craton_stage_poles = dict(moving_poles)
        for _old_pid, _new_pid in _craton_plate_remap.items():
            if _new_pid in moving_poles:
                _craton_stage_poles[_old_pid] = moving_poles[_new_pid]
        cratons = sdp.nudge_polygons(cratons, _craton_stage_poles)
        _sutures = sdp.advect_sutures(_sutures, previous_plate_id, moving_poles)

        # Accumulate finite rotations: F(t) = S * F(t - dt), for the move just made
        for pid, stage_rot in moving_poles.items():
            prev = cumulative_rotation_dict.get(pid, identity)
            cumulative_rotation_dict[pid] = stage_rot * prev
            rotation_records.append((t, pid, cumulative_rotation_dict[pid]))
        _step_times['nudge_polygons'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

        # ------------------------------------------------------------------
        # Rebuild the plate configuration and compute the next stage rotations
        # ------------------------------------------------------------------
        # Rift splits decided last step, their cuts moved with the plate
        _rift_now = [(pid, sdp.rotate_vector(moving_poles.get(pid), n),
                      None if c is None else dict(xyz=sdp.rotate_xyz(moving_poles.get(pid), c['xyz']), side1=c['side1']))
                     for pid, n, c in _rift_planned]
        _rift_planned = []

        (new_plate_polygons,
         stage_pole_dict,
         rasters,
         boundaries_gdf,
         _new_plates,
         collision_state,
         next_plate_id) = sdp.compute_timestep(rpolygons, forces, moving_poles, delta_time,
                                               cratons=cratons, params=params,
                                               polarity_dict=polarity_dict,
                                               merge_pairs=_merge_pairs,
                                               fragmentation_cooldown_pids=_fragmentation_cooldown_pids,
                                               collision_state=collision_state,
                                               craton_id_remap=_craton_plate_remap,
                                               next_plate_id=next_plate_id,
                                               speed_state=speed_state,
                                               previous_plate_id=previous_plate_id,
                                               force_fields=_fields,
                                               force_diag=force_diag,
                                               rift_splits=_rift_now,
                                               rift_skip=set(_rift_cooldown),
                                               sutures=_sutures if len(_sutures) else None)
        _step_times['compute_timestep'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

        # Rifting (rifting: trial): a plate splits once its best cut has passed the trial
        # split on a similar great circle for rift_persist_myr in a row
        if params['rifting_mode'] == 'trial':
            for _child, _parent in force_diag.get('rift_applied', []):
                cratons = sdp.rebuild_split_cratons(cratons, rasters['craton'], _parent, _child,
                                                    _craton_plate_remap, params['raster_spacing'],
                                                    params.get('topo_simplify_tolerance', 0.4))
                _rift_cooldown[_child] = _rift_cooldown[_parent] = params.get('rift_cooldown_myr', 10.0)
                _rift_events.append(dict(time_myr=t, parent=_parent, child=_child))
                pd.DataFrame(_rift_events).to_csv(os.path.join(output_dir, 'rift_events.csv'), index=False)
            _persist = max(1, int(round(params.get('rift_persist_myr', 3.0) / delta_time)))
            _state = {}
            for _pid, _c in force_diag.get('rift', {}).items():
                _prev = _rift_state.get(_pid)
                _same = _prev is not None and abs(
                    sdp.rotate_vector(moving_poles.get(_pid), _prev['normal']) @ _c['normal']) > np.cos(np.radians(15.0))
                _count = _prev['count'] + 1 if _same else 1
                if _count >= _persist:
                    _rift_planned.append((_pid, _c['normal'], _c.get('cells')))
                    print('  [rift] plate {} passes the trial split for {} steps (tension / strength {:.2f}, '
                          'cut {:.0f} deg): splits next step'.format(_pid, _count, _c['ratio'], _c['length_deg']))
                else:
                    _state[_pid] = dict(count=_count, normal=_c['normal'])
            _rift_state = _state
            _rift_cooldown = {p: r - delta_time for p, r in _rift_cooldown.items() if r > delta_time}

        # New plates (fragments) inherit their parent's rotation history, so the
        # rotation file reconstructs them from t = 0 like any other plate.
        if _new_plates:
            print('  [fragment] {} new plate(s): {}'.format(
                len(_new_plates),
                '  '.join('plate {} → child {}'.format(parent, new)
                          for new, parent in _new_plates)))
            for _new_pid, _parent_pid in _new_plates:
                cumulative_rotation_dict[_new_pid] = cumulative_rotation_dict.get(_parent_pid, identity)
                rotation_records.extend([(tt, _new_pid, rot) for tt, p, rot in rotation_records
                                         if p == _parent_pid])

        # -----------------------------------------------------------------------
        # Fragmentation / merging bookkeeping
        # -----------------------------------------------------------------------
        _merge_pairs = []
        if params.get('merging_enabled', False) and 'divergence' in boundaries_gdf.columns:
            _merge_thresh   = params.get('merge_conservative_threshold', 0.05)  # km/Myr
            _merge_min_steps = int(params.get('merge_min_steps', 5))
            _merge_cooldown = params.get('merge_cooldown_myr', 10.0)
            _no_driver_req  = params.get('merge_require_no_driver', True)
            _div_sub        = params['div_threshold_subduction']
            _div_ridge      = params['div_threshold_ridge']
            _areas          = sdp.plate_areas(rasters['plate_id'])

            # Accumulate per-pair mean |divergence| and per-plate driver flags
            _pair_divs, _plate_driver = {}, {}
            for _, _br in boundaries_gdf.iterrows():
                _p1 = _br.get('plateid1_poly1'); _p2 = _br.get('plateid1_poly2')
                if not _p1 or not _p2 or _p1 == 0 or _p2 == 0:
                    continue
                _pair = (min(int(_p1), int(_p2)), max(int(_p1), int(_p2)))
                _div  = float(_br.get('divergence', 0))
                _pair_divs.setdefault(_pair, []).append(abs(_div))
                for _pid in (_pair[0], _pair[1]):
                    if _div > _div_ridge or _div < _div_sub:
                        _plate_driver[_pid] = True
                    else:
                        _plate_driver.setdefault(_pid, False)

            _new_contact = {}
            for _pair, _divs in _pair_divs.items():
                _mean_div = float(np.mean(_divs))
                _prev     = _conservative_contact.get(_pair, 0)
                _new_contact[_pair] = _prev + 1 if _mean_div < _merge_thresh else 0

                if _new_contact[_pair] >= _merge_min_steps:
                    _p1, _p2 = _pair
                    if _p1 in _plate_cooldown or _p2 in _plate_cooldown:
                        continue
                    _big, _small = (_p1, _p2) if _areas.get(_p1, 0) >= _areas.get(_p2, 0) else (_p2, _p1)
                    # only absorb a plate that is not itself driven by a ridge or trench
                    if _no_driver_req and _plate_driver.get(_small, False):
                        continue
                    _merge_pairs.append((_big, _small))
                    _plate_cooldown[_p1] = _merge_cooldown
                    _plate_cooldown[_p2] = _merge_cooldown
                    _new_contact[_pair]  = 0

            _conservative_contact = _new_contact

        # Welding (collision_weld: on): after weld_after_myr of continuous continental
        # collision the smaller plate is absorbed next step; sutures: tracers marks the join
        if params['collision_weld_mode'] == 'on':
            _weld_after = max(1, int(round(params.get('weld_after_myr', 30.0) / delta_time)))
            _areas_w = sdp.plate_areas(rasters['plate_id'])
            _used = {_p for _pr in _merge_pairs for _p in _pr}
            for _pair, _steps in sorted(collision_state.items()):
                _p1, _p2 = int(_pair[0]), int(_pair[1])
                if _steps < _weld_after or _p1 in _used or _p2 in _used or _p1 not in _areas_w or _p2 not in _areas_w:
                    continue
                if len(_areas_w) - len(_merge_pairs) <= 2:
                    continue            # never weld the world down to one plate
                _big, _small = (_p1, _p2) if _areas_w[_p1] >= _areas_w[_p2] else (_p2, _p1)
                _merge_pairs.append((_big, _small))
                _used |= {_p1, _p2}
                _rift_cooldown[_big] = params.get('rift_cooldown_myr', 10.0)
                _row = dict(time_myr=t, host=_big, absorbed=_small, collided_myr=_steps * delta_time, suture_id=0)
                _smp = force_diag.get('collision_samples', {}).get((_p1, _p2))
                if params['sutures_mode'] == 'tracers' and _smp:
                    _sutures = pd.concat([_sutures, sdp.suture_points(_smp, _next_suture_id, t)], ignore_index=True)
                    _row['suture_id'] = _next_suture_id
                    _next_suture_id += 1
                _welds.append(_row)
                pd.DataFrame(_welds).to_csv(os.path.join(output_dir, 'welds.csv'), index=False)
                print('  [weld] plates {} and {} collided for {:.0f} Myr: {} absorbed into {} next step{}'.format(
                    _p1, _p2, _steps * delta_time, _small, _big,
                    ', suture {} ({} points)'.format(_row['suture_id'], len(_smp)) if _row['suture_id'] else ''))

        if params.get('fragmentation_enabled', False) and _new_plates:
            _frag_cooldown = params.get('fragmentation_cooldown_myr', 5.0)
            for _new_pid, _parent_pid in _new_plates:
                _plate_cooldown[_new_pid]    = _frag_cooldown
                _plate_cooldown[_parent_pid] = _frag_cooldown

        # Decrement all cooldowns
        _plate_cooldown = {_pid: _rem - delta_time
                           for _pid, _rem in _plate_cooldown.items() if _rem > delta_time}
        _fragmentation_cooldown_pids = set(_plate_cooldown.keys())
        _step_times['frag_merge_bookkeeping'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

        # ------------------------------------------------------------------
        # Seafloor tracers: move with the same rotations as the plates, remove
        # subducted ones, then seed new crust at the ridges
        # ------------------------------------------------------------------
        _craton_counts = sdp.polarity_counts(
            sdp.craton_counts_from_raster(rasters['craton'] if 'craton' in rasters else None), params)
        tracers = sdp.advect_tracers(tracers, previous_plate_id, moving_poles, delta_time,
                                     boundaries_gdf=boundaries_gdf,
                                     div_threshold=params['div_threshold_subduction'],
                                     new_plate_id_raster=rasters['plate_id'])
        tracers = sdp.deactivate_subducted_tracers(tracers, boundaries_gdf,
                                                   threshold=params['div_threshold_subduction'],
                                                   plate_id_raster=rasters['plate_id'],
                                                   polarity_dict=polarity_dict,
                                                   craton_counts=_craton_counts)
        tracers = tracers[tracers.active.astype(bool)]
        new_tracers = sdp.seed_ridge_tracers(boundaries_gdf,
                                             spacing_deg=params['tracer_spacing_deg'],
                                             threshold=params['div_threshold_ridge'],
                                             offset_deg=params['tracer_offset_deg'])
        tracers = gpd.GeoDataFrame(pd.concat([tracers, new_tracers], ignore_index=True),
                                   geometry='geometry').set_crs(epsg=4326, allow_override=True)
        if params['tracer_thinning_mode'] == 'cell' and i % int(params.get('tracer_thin_interval', 10)) == 0:
            tracers = sdp.thin_tracers(tracers, rasters['plate_id'],
                                       cell_deg=params.get('tracer_thin_cell_deg', 0.5))
        _step_times['tracers'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

        if _slab_track:
            slab_tracers = sdp.update_slab_tracers(
                slab_tracers, boundaries_gdf, delta_time,
                div_threshold_subduction=params['div_threshold_subduction'],
                slowing_factor=params.get('slab_slowing_factor', 4),
                max_sink_rate_km_per_myr=params.get('slab_max_sink_rate', 200.0),
                sink_rate_factor=params.get('slab_sink_rate_factor', 1.0))
            _step_times['slab_tracers'] = time.perf_counter() - _t0; _t0 = time.perf_counter()
        if _slab_enabled:
            cmb_heat = sdp.update_cmb_heat(
                cmb_heat, slab_tracers, delta_time,
                spacing_deg=_mantle_spacing,
                heating_rate=params.get('cmb_heating_rate', 0.02),
                slab_suppression=params.get('cmb_slab_suppression', 0.05),
                surface_coupling=params.get('cmb_surface_coupling', 0.1),
                sigma_deg=params.get('cmb_sigma_deg', 10.0),
                max_heat=params.get('cmb_max_heat', 3.0))
            _step_times['cmb_heat'] = time.perf_counter() - _t0; _t0 = time.perf_counter()
        if _slab_track and 'divergence' in boundaries_gdf.columns:
            _sub = boundaries_gdf[boundaries_gdf['divergence'] < params['div_threshold_subduction']]
            _active_slabs = slab_tracers[slab_tracers['active'].astype(bool)]
            print('  [slab] sub-bdry div (km/Myr): min={:.3f} max={:.3f}  '
                  'sink_rate (km/Myr): max={:.1f}  '
                  'slablets: n={:d} depth={:.0f}-{:.0f} km'.format(
                      float(_sub['divergence'].min()) if len(_sub) else 0.0,
                      float(_sub['divergence'].max()) if len(_sub) else 0.0,
                      float(_active_slabs['sink_rate_km_per_myr'].max()) if len(_active_slabs) else 0.0,
                      len(_active_slabs),
                      float(_active_slabs['depth_km'].min()) if len(_active_slabs) else 0.0,
                      float(_active_slabs['depth_km'].max()) if len(_active_slabs) else 0.0))

        seamounts = sdp.advect_tracers(seamounts, previous_plate_id, moving_poles, delta_time,
                                       boundaries_gdf=boundaries_gdf,
                                       div_threshold=params['div_threshold_subduction'],
                                       new_plate_id_raster=rasters['plate_id'])
        seamounts = sdp.deactivate_subducted_tracers(seamounts, boundaries_gdf,
                                                      threshold=params['div_threshold_subduction'],
                                                      plate_id_raster=rasters['plate_id'],
                                                      polarity_dict=polarity_dict,
                                                      craton_counts=_craton_counts)
        new_seamounts = sdp.seed_hotspot_seamounts(hotspot_locs, prob=params['hotspot_prob'], rng=rng)
        if len(new_seamounts):
            _parts = [df for df in (seamounts, new_seamounts) if len(df)]
            seamounts = gpd.GeoDataFrame(
                pd.concat(_parts, ignore_index=True),
                geometry='geometry').set_crs(epsg=4326, allow_override=True)
        _step_times['seamounts'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

        # Plot / save (non-blocking; the worker writes files in the background, but
        # this waits when the queue is full).  The last step is always saved.
        _enqueue_outputs(i, force_save=(i == n_steps))
        _step_times['outputs_enqueue'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

        _stats_t0 = time.perf_counter()
        _row = _record_step_stats(i, {**_step_times, 'step_stats': 0.0}, previous_poles=moving_poles)
        _step_times['step_stats'] = time.perf_counter() - _stats_t0
        print('  [kinematics] plates={n_plates}  rms speed={rms_speed_km_per_myr:.1f} km/Myr  '
              'max plate rms={max_plate_rms_km_per_myr:.1f}  net rotation={net_rotation_deg_per_myr:.3f} deg/Myr'.format(**_row))

        total = sum(_step_times.values())
        print('  --- step timing (step {:d}, {:.1f}s total) ---'.format(i, total))
        for _name, _dt in sorted(_step_times.items(), key=lambda x: -x[1]):
            print('  {:20s}  {:5.1f}s  ({:4.1f}%)'.format(_name, _dt, 100 * _dt / total))

finally:
    # Drain output queue (also after an error or Ctrl-C in the main loop)
    _output_queue.join()      # wait for all enqueued plot/save tasks to finish
    _output_queue.put(None)   # stop worker
    _output_thread.join()

active_seamounts = seamounts[seamounts.active.astype(bool)]
if len(active_seamounts):
    active_seamounts[['geometry', 'hotspot_id', 'age']].to_file(
        _out('seamounts_final.geojson'), driver='GeoJSON')

if _output_errors:
    print('\n{} output task(s) failed:\n  {}'.format(len(_output_errors), '\n  '.join(_output_errors)))
    sys.exit(1)
