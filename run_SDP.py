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
reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis5_topologies.gpml')
reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis5_topological_polygons.gpml')
cratons = pygplates.FeatureCollection('/Users/simon/GIT/octopus/Atlantis3/atlantis5_cratons.geojson')


reconstruction_time = 0
delta_time = 1.
forces = {'slab_pull': 10,'slab_suction': 10, 'ridge_push':1}


# The first time step is based on topologies (And plate boundary types)
snapshot = reconstruction_model.plate_snapshot(reconstruction_time)

result = pdriver.get_torques_for_plate_snapshot(snapshot, reconstruction_time, delta_time, forces)


stage_pole_dict = {}

for res in result:
    stage_pole_dict[res[0]] = res[3]
    print(res[0], res[3])

print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])

cratons = sdp.nudge_polygons(cratons, stage_pole_dict)

rpolygons = sdp.nudge_topological_polygons(snapshot, 
                                           forces=forces, 
                                           reconstruction_time=0, 
                                           delta_time=delta_time)

(new_plate_polygons,
 stage_pole_dict,
 rasters, _) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time)

#print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])

tracers = gpd.GeoDataFrame(columns=['geometry', 'active', 'age']).set_crs(epsg=4326)

for t in range(100):

    print('Timestep {}'.format(t))
    
    rpolygons = sdp.nudge_polygons(
        new_plate_polygons.explode(index_parts=False).reset_index(drop=True), 
        stage_pole_dict)
    
    cratons = sdp.nudge_polygons(cratons, stage_pole_dict)


    (new_plate_polygons,
     stage_pole_dict,
     rasters,
     boundaries_gdf) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time, cratons=cratons)

    new_tracers = sdp.seed_ridge_tracers(boundaries_gdf, rasters['divergence'])
    tracers = pd.concat([tracers, new_tracers], ignore_index=True)
    tracers = sdp.advect_tracers(tracers, rasters['plate_id'], stage_pole_dict)
    tracers = sdp.deactivate_subducted_tracers(tracers, rasters['divergence'])

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
    fig.text(x=new_plate_polygons.centroid.x, y=new_plate_polygons.centroid.y,
             text=new_plate_polygons['PLATEID1'], font='10p,Helvetica-Bold,black', justify='CM')
    active_tracers = tracers[tracers.active.astype(bool)]
    inactive_tracers = tracers[~tracers.active.astype(bool)]
    if len(active_tracers):
        fig.plot(data=active_tracers, style='c0.04c', fill='darkmagenta', pen='0.2p,darkmagenta')
    if len(inactive_tracers):
        fig.plot(data=inactive_tracers, style='c0.04c', fill='orange', pen='0.2p,orange')
    fig.savefig('./output/test_plot_{:d}.png'.format(t))
    
    