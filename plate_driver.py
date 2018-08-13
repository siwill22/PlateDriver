import pygplates
import numpy as np
import sys
sys.path.append('/Users/Simon/GIT/GPlatesClassStruggle/')
from reconstruction_classes import *

from ptt.subduction_convergence import find_overriding_and_subducting_plates




def get_velocity_vector_at_centroid(equivalent_stage_rotation, centroid_point, delta_time):
    
    velocity_vectors = pygplates.calculate_velocities(
                                        [centroid_point],
                                        equivalent_stage_rotation,
                                        delta_time)

    # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
    velocity_ne = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
        [centroid_point],
        velocity_vectors)[0]
    
    return velocity_ne

        
def get_torque_for_segment(segment, subducting_normal_reversal):
    
    # If segment is zero length, assume that the torque is zero
    if segment.is_zero_length():
        segment_torque = pygplates.Vector3D(0, 0, 0)
    
    else:
        normal = -subducting_normal_reversal * segment.get_great_circle_normal()

        segment_position_vector = segment.get_arc_point(0.5).to_xyz()

        segment_weight = segment.get_arc_length()
        #print segment_weight

        stage_pole_position = pygplates.Vector3D.cross(normal,segment_position_vector)
        rotation_rate_vector = pygplates.Vector3D(stage_pole_position.to_xyz()) * np.radians(0.1)

        v = pygplates.Vector3D.cross(rotation_rate_vector, segment_position_vector)

        segment_torque = pygplates.Vector3D.cross(segment_position_vector, v) * segment_weight

    return segment_torque


def get_torques_for_plate_snapshot(snapshot, reconstruction_time, delta_time,
                                   forces = ['slab_pull','slap_suction'], 
                                   slab_suction_to_slab_pull_ratio = 0.25):

    computed_plate_motions = []
    
    # loop over each plate, calling function to 
    for target_plate_id in snapshot.plate_ids:
        
        torque_for_plate = get_torque_for_plate_polygon(snapshot, target_plate_id, reconstruction_time, delta_time,
                                                        forces, slab_suction_to_slab_pull_ratio)
        
        computed_plate_motions.append(torque_for_plate)
        
    return computed_plate_motions
        

        
def get_torque_for_plate_polygon(snapshot, target_plate_id, reconstruction_time, delta_time,
                                 forces = ['slab_pull','slab_suction'], 
                                 slab_suction_to_slab_pull_ratio = 0.25):

    net_torque = pygplates.Vector3D(0, 0, 0)

    for topology in snapshot.resolved_topological_sections:

        shared_sub_segments = topology.get_shared_sub_segments()
        for shared_sub_segment in shared_sub_segments:

            # if the segment is not topologically correct, skip it 
            if len(shared_sub_segment.get_sharing_resolved_topologies()) != 2:
                #print 'skipping segment without exactly 2 shared polygons'
                continue
            if shared_sub_segment.get_resolved_feature().get_geometry() is None:
                continue

            # check if the segment is part of the boundary of the target plate
            if target_plate_id in [resolved_topology.get_feature().get_reconstruction_plate_id() for resolved_topology in shared_sub_segment.get_sharing_resolved_topologies()]:            

                # check if it is a subduction zone
                if shared_sub_segment.get_resolved_feature().get_feature_type() == pygplates.FeatureType.gpml_subduction_zone:

                    # find which of the plate pair is the overriding plate, and identify polarity
                    try:
                        (overriding_plate, 
                         subducting_plate, 
                         subduction_polarity) = find_overriding_and_subducting_plates(shared_sub_segment, reconstruction_time)
                    except:
                        #print shared_sub_segment.get_geometry()
                        break

                    # not necessary if the same topological error is caught previously? 
                    if overriding_plate is None:
                        continue

                    if subducting_plate.get_resolved_feature().get_reconstruction_plate_id() == target_plate_id:

                        if 'slab_pull' in forces:

                            if subduction_polarity == 'Left':
                                subducting_normal_reversal = 1
                            else:
                                subducting_normal_reversal = -1

                            for segment in shared_sub_segment.get_resolved_feature().get_geometry().get_segments():
                                segment_torque = get_torque_for_segment(segment, subducting_normal_reversal)

                                net_torque = net_torque + segment_torque

                    if overriding_plate.get_resolved_feature().get_reconstruction_plate_id() == target_plate_id:

                        if 'slab_suction' in forces:

                            if subduction_polarity == 'Left':
                                subducting_normal_reversal = -1
                            else:
                                subducting_normal_reversal = 1

                            for segment in shared_sub_segment.get_resolved_feature().get_geometry().get_segments():
                                segment_torque = get_torque_for_segment(segment, subducting_normal_reversal)

                                net_torque = net_torque + segment_torque * slab_suction_to_slab_pull_ratio

                    # TODO add more forces here, depending on plate boundary type


    # Now get the centroid of the plate polygon, and evaluate the computed rotation at the centroid   
    # also, evaluate the absolute rotation from the reconstruction model for comparison
    for topology in snapshot.resolved_topologies:
        if topology.get_resolved_feature().get_reconstruction_plate_id() == target_plate_id:
            centroid_point = topology.get_resolved_geometry().get_interior_centroid()

    if not net_torque.is_zero_magnitude():
        computed_stage_rotation = pygplates.FiniteRotation(pygplates.PointOnSphere(net_torque.to_normalised().to_xyz()).to_lat_lon(),
                                                           net_torque.get_magnitude())

        computed_velocity_ne = get_velocity_vector_at_centroid(computed_stage_rotation, centroid_point, delta_time)
    else:
        computed_stage_rotation = None
        computed_velocity_ne = None

    # get the reconstruction model vector at the centroid for comparison
    reconstruction_stage_rotation = snapshot.rotation_model.get_rotation(reconstruction_time,
                                                                         target_plate_id,
                                                                         from_time=reconstruction_time+delta_time,
                                                                         anchor_plate_id=0)
    recontruction_velocity_ne = get_velocity_vector_at_centroid(reconstruction_stage_rotation, centroid_point, delta_time)
    
    
    return [target_plate_id, centroid_point, net_torque, 
            computed_stage_rotation, computed_velocity_ne, 
            reconstruction_stage_rotation, recontruction_velocity_ne]

