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
total_time = 60.
forces = {'slab_pull': 8, 'slab_suction': 8, 'ridge_push': 3, 'continental_collision': 8}

params = {
    'raster_spacing':              1,        # degrees; controls velocity grid resolution and polygon transform
    'tessellation_deg':            0.2,      # polygon edge densification before rasterising
    'topo_simplify_tolerance':     1,        # polygon smoothing tolerance after polygonise
    'boundary_search_radius_deg':  2,        # plate ID neighbour search radius at boundaries
    'div_threshold_subduction':   -1.0,   # cm/yr relative normal velocity below this → subduction
    'div_threshold_ridge':         1.0,   # cm/yr relative normal velocity above this → ridge
    'tracer_spacing_deg':          2.0,      # seed point spacing along ridge boundaries
    'tracer_offset_deg':           0.6,      # perpendicular offset from ridge axis at seeding
    'seafloor_spacing_deg':        2.0,      # grid spacing for initial seafloor age tracers
    'spreading_rate_km_per_myr':   50.0,     # assumed spreading rate for age estimation
    'csv_save_interval':           10,       # save active-tracer CSV every N timesteps
}


# The first time step is based on topologies (And plate boundary types)
snapshot = reconstruction_model.plate_snapshot(reconstruction_time)

result = pdriver.get_torques_for_plate_snapshot(snapshot, reconstruction_time, delta_time, forces)


stage_pole_dict = {}

for res in result:
    stage_pole_dict[res[0]] = res[3]
#    print(res[0], res[3])
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
    tracers = sdp.advect_tracers(tracers, rasters['plate_id'], stage_pole_dict, delta_time)
    tracers = sdp.deactivate_subducted_tracers(tracers, boundaries_gdf,
                                               threshold=params['div_threshold_subduction'])

    fig = pygmt.Figure()
    fig.coast(land='grey', projection='W10i', region='d')
    #fig.plot(data=clipped_polygons, fill='blue', transparency=50)
    pygmt.makecpt(cmap='polar', series=[-0.001,0.001], background='o')
    fig.grdimage(rasters['divergence'], cmap=True)
    fig.plot(data=new_plate_polygons.iloc[:], pen='1p,black', transparency=50)
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
        rasters['divergence'].to_netcdf('./output/divergence_{:04d}.nc'.format(i))
        rasters['plate_id'].to_netcdf('./output/plate_id_{:04d}.nc'.format(i))
        if 'craton' in rasters:
            rasters['craton'].to_netcdf('./output/craton_{:04d}.nc'.format(i))

sdp.write_rotation_file(rotation_records, './output/synthetic_rotations.rot')
