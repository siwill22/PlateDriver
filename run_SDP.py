import pygplates
import numpy as np
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import pygmt

from gprm import ReconstructionModel
from gprm.utils.create_gpml import gpml2gdf, gdf2gpml
from ptt import separate_ridge_transform_segments
from gprm.utils.proximity import rasterize_polygons

import plate_driver as pdriver
import self_driving_plates as sdp



reconstruction_model = ReconstructionModel('')
reconstruction_model.add_rotation_model('/Users/simon/GIT/octopus/Atlantis3/atlantis3_rotations.rot')
reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis4_topologies.gpml')
reconstruction_model.add_dynamic_polygons('/Users/simon/GIT/octopus/Atlantis3/atlantis4_topological_polygons.gpml')



reconstruction_time = 0
delta_time = 10.
forces = {'slab_pull': 100,'slab_suction': 50}
#forces = {'ridge_push':1}


# The first time step is based on topologies (And plate boundary types)
snapshot = reconstruction_model.plate_snapshot(reconstruction_time)

result = pdriver.get_torques_for_plate_snapshot(snapshot, reconstruction_time, delta_time, forces)

stage_pole_dict = {}

for res in result:
    stage_pole_dict[res[0]] = res[3]

print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])


rpolygons = sdp.nudge_topological_polygons(snapshot, 
                                           forces=forces, 
                                           reconstruction_time=0, 
                                           delta_time=delta_time)

(new_plate_polygons, 
 stage_pole_dict) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time)

#print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])

for t in range(50):

    print('Timestep {}'.format(t))
    
    rpolygons = sdp.nudge_polygons(
        new_plate_polygons.explode(index_parts=False).reset_index(drop=True), 
        stage_pole_dict)

    (new_plate_polygons, 
     stage_pole_dict) = sdp.compute_timestep(rpolygons, forces, stage_pole_dict, delta_time)

    #print(stage_pole_dict)
    
    fig = pygmt.Figure()
    fig.coast(land='grey', projection='W10i', region='d')
    #fig.plot(data=clipped_polygons, fill='blue', transparency=50)
    fig.plot(data=new_plate_polygons.iloc[:], pen='red')
    #fig.show(width=1000)
    fig.savefig('./output/test_plot_{:d}.png'.format(t))
    
    