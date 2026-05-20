import pygplates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import pygmt
import geopandas as gpd

from gprm import ReconstructionModel
from gprm.utils.create_gpml import gpml2gdf, gdf2gpml
from ptt import separate_ridge_transform_segments
from gprm.utils.proximity import rasterize_polygons

import plate_driver as pdriver
import self_driving_plates as sdp



reconstruction_model = ReconstructionModel('')

#reconstruction_model.add_rotation_model('/Users/simon/GIT/octopus/Atlantis3/atlantis3_rotations.rot')
#reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis4_topologies.gpml')
#reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis4_topological_polygons.gpml')
#cratons = pygplates.FeatureCollection('/Users/simon/GIT/octopus/Atlantis3/atlantis4_cratons.geojson')

reconstruction_model.add_rotation_model('/Users/simon/GIT/octopus/Atlantis3/atlantis3_rotations.rot')
reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis6_topologies.gpml')
reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis6_topological_polygons.gpml')
cratons = pygplates.FeatureCollection('/Users/simon/GIT/octopus/Atlantis3/atlantis6_cratons.geojson')


reconstruction_time = 0
delta_time = 1.
total_time = 20.
forces = {'slab_pull': 8, 'slab_suction': 8, 'ridge_push': 3, 'continental_collision': 8}

params = {
    'raster_spacing':              0.5,        # degrees; controls velocity grid resolution and polygon transform
    'tessellation_deg':            0.2,      # polygon edge densification before rasterising
    'topo_simplify_tolerance':     0.5,        # polygon smoothing tolerance after polygonise
    'boundary_search_radius_deg':  2,        # plate ID neighbour search radius at boundaries
    'div_threshold_subduction':   -0.0001,   # cm/yr relative normal velocity below this → subduction
    'div_threshold_ridge':         0.2,   # cm/yr relative normal velocity above this → ridge
    'tracer_spacing_deg':          1.0,      # seed point spacing along ridge boundaries
    'tracer_offset_deg':           0.6,      # perpendicular offset from ridge axis at seeding
    'seafloor_spacing_deg':        2.0,      # grid spacing for initial seafloor age tracers
    'spreading_rate_km_per_myr':   50.0,     # assumed spreading rate for age estimation
    'csv_save_interval':           5,       # save active-tracer CSV every N timesteps
    'eq_n_candidates':             100000,   # random sphere points sampled for earthquake generation
    'eq_ridge_max_dist_deg':       1.0,      # max distance from ridge axis for ridge earthquakes
    'eq_ridge_max_depth_km':       50.0,     # max depth for ridge earthquakes (km)
    'eq_slab_dip_deg':             45.0,     # slab dip angle (degrees)
    'eq_slab_max_depth_km':        650.0,    # maximum slab depth (km)
    'eq_ridge_prob':               0.15,      # fraction of ridge candidates accepted (controls density)
    'eq_subduction_prob':          1.0,      # fraction of subduction candidates accepted
    'topo_continental_base_m':     400,      # background continental elevation (m)
    'topo_noise_fraction':         0.3,      # noise amplitude as fraction of local elevation
    'topo_spectral_beta':          2.0,      # SH power-law exponent (2 = brown / realistic terrain)
    'topo_andean_max_dist_deg':    11.0,      # cutoff distance from trench for Andean mountains
    'topo_andean_arc_dist_deg':    4.5,      # distance of arc-peak from trench (degrees)
    'topo_andean_peak_m':         3000,      # peak Andean elevation (m)
    'topo_himalayan_max_dist_deg': 8.0,      # half-width of Himalayan collision zone (degrees)
    'topo_himalayan_peak_m':      5000,      # peak Himalayan elevation (m)
    'n_hotspots':                   20,      # number of randomly placed hot spot sources
    'hotspot_prob':                0.3,      # probability of eruption per hot spot per timestep
    # To use fixed locations instead: set hotspot_locations = [(lon, lat), ...]
}


# The first time step is based on topologies (And plate boundary types)
snapshot = reconstruction_model.plate_snapshot(reconstruction_time)

result = pdriver.get_torques_for_plate_snapshot(snapshot, reconstruction_time, delta_time, forces)


stage_pole_dict = {}

for res in result:
    stage_pole_dict[res[0]] = res[3]
    print(res[0], res[3])
#print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])

cratons = sdp.nudge_polygons(cratons, stage_pole_dict)

rpolygons = sdp.nudge_topological_polygons(snapshot, 
                                           forces=forces, 
                                           reconstruction_time=0, 
                                           delta_time=delta_time)

polarity_dict = {}
cumulative_rotation_dict = {}
rotation_records = []

(new_plate_polygons,
 stage_pole_dict,
 rasters,
 boundaries_gdf) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time,
                                         params=params, polarity_dict=polarity_dict)

#print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])

tracers = sdp.generate_seafloor_age_tracers(
    new_plate_polygons, boundaries_gdf,
    spacing_deg=params['seafloor_spacing_deg'],
    spreading_rate_km_per_myr=params['spreading_rate_km_per_myr'],
    div_threshold=params['div_threshold_ridge'])

hotspot_locs = sdp.generate_hotspot_locations(n=params['n_hotspots'])
seamounts = gpd.GeoDataFrame(columns=['geometry', 'hotspot_id', 'age', 'active'],
                             geometry='geometry').set_crs(epsg=4326)

for i in range(1, int(total_time / delta_time) + 1):
    t = i * delta_time

    print('t = {:.1f} Myr (step {})'.format(t, i))
    
    rpolygons = sdp.nudge_polygons(
        new_plate_polygons.explode(index_parts=False).reset_index(drop=True), 
        stage_pole_dict)
    
    cratons = sdp.nudge_polygons(cratons, stage_pole_dict)


    (new_plate_polygons,
     stage_pole_dict,
     rasters,
     boundaries_gdf) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time,
                                             cratons=cratons, params=params,
                                             polarity_dict=polarity_dict)
    
    for key,item in stage_pole_dict.items():
        print(key,item.get_lat_lon_euler_pole_and_angle_degrees())

    # Accumulate finite rotations: F(t) = S_t * F(t-dt)
    t_record = t
    for pid, stage_rot in stage_pole_dict.items():
        prev = cumulative_rotation_dict.get(pid, pygplates.FiniteRotation())
        cumulative_rotation_dict[pid] = stage_rot * prev
        rotation_records.append((t_record, pid, cumulative_rotation_dict[pid]))

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

    fig = pygmt.Figure()
    fig.basemap(frame='f', projection='W10i', region='d')
    #fig.plot(data=clipped_polygons, fill='blue', transparency=50)
    pygmt.makecpt(cmap='polar', series=[-0.001,0.001], background='o')
    #fig.grdimage(rasters['divergence'], cmap=True)
    fig.plot(data=new_plate_polygons.iloc[:], pen='1p,black', transparency=50)
    # Overlay per-segment boundary type: red=ridge, blue=subduction, black=conservative
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
    glon,glat = np.meshgrid(rasters.lon, rasters.lat)
    stride=5
    magnitude = np.sqrt(rasters['velocity_east'].data[::stride,::stride].flatten()**2 + rasters['velocity_north'].data[::stride,::stride].flatten()**2)
    azimuth = np.degrees(np.arctan2(rasters['velocity_east'].data[::stride,::stride].flatten(), rasters['velocity_north'].data[::stride,::stride].flatten())) % 360
    fig.plot(data={'x': glon[::stride,::stride].flatten(),
                'y': glat[::stride,::stride].flatten(),
                'azimuth': azimuth,
                'magnitude': magnitude/(2000.*delta_time)},
            style='V0.15c+e+a45+ggray40+n', pen="0.1p,black", fill="black")
    fig.plot(data=gpml2gdf(cratons), pen='1p,black', fill='grey', transparency=50)
    label_pts = new_plate_polygons.representative_point()
    fig.text(x=label_pts.x, y=label_pts.y,
             text=new_plate_polygons['PLATEID1'], font='10p,Helvetica-Bold,black', justify='CM')
    active_tracers = tracers[tracers.active.astype(bool)]
    if len(active_tracers):
        vmax = 100 #max(float(active_tracers.age.max()), 1.0)
        pygmt.makecpt(cmap='cubhelix', series=[0, vmax], reverse=True)
        fig.plot(x=active_tracers.geometry.x.values, y=active_tracers.geometry.y.values,
                 style='c0.08c', fill=active_tracers.age.values.astype(float),
                 cmap=True, pen='0.1p,black')
    fig.savefig('./output/test_plot_{:04d}.png'.format(i))

    if i % params['csv_save_interval'] == 0:
        active_tracers[['geometry', 'age']].assign(
            lon=active_tracers.geometry.x,
            lat=active_tracers.geometry.y
        ).drop(columns='geometry').to_csv('./output/tracers_{:04d}.csv'.format(i), index=False)

        age_grid = None
        if len(active_tracers):
            age_grid = sdp.generate_seafloor_age_grid(
                './output/tracers_{:04d}.csv'.format(i),
                craton_raster=rasters.get('craton'),
                spacing=params['raster_spacing'])
            age_grid.to_netcdf('./output/seafloor_age_{:04d}.nc'.format(i))

        rasters['divergence'].to_netcdf('./output/divergence_{:04d}.nc'.format(i))
        rasters['plate_id'].to_netcdf('./output/plate_id_{:04d}.nc'.format(i))
        if 'craton' in rasters:
            rasters['craton'].to_netcdf('./output/craton_{:04d}.nc'.format(i))

        # Craton pixel counts for authoritative polarity determination
        craton_counts = sdp.craton_counts_from_raster(
            rasters['craton'] if 'craton' in rasters else None)

        # Save boundary segments with type and subduction polarity
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

            bdf.to_file('./output/boundaries_{:04d}.geojson'.format(i), driver='GeoJSON')

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
        earthquakes.to_file('./output/earthquakes_{:04d}.geojson'.format(i), driver='GeoJSON')

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
            topo.to_netcdf('./output/topography_{:04d}.nc'.format(i))

        if age_grid is not None:
            topo_bathy = sdp.generate_topo_bathy_grid(age_grid, topography_raster=topo)
            topo_bathy.to_netcdf('./output/topo_bathy_{:04d}.nc'.format(i))

sdp.write_rotation_file(rotation_records, './output/synthetic_rotations.rot')

active_seamounts = seamounts[seamounts.active.astype(bool)]
if len(active_seamounts):
    active_seamounts[['hotspot_id', 'age']].assign(
        lon=active_seamounts.geometry.x,
        lat=active_seamounts.geometry.y,
    ).to_csv('./output/hotspot_seamounts_final.csv', index=False)
