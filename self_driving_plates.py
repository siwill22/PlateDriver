import pygplates
from gprm.utils.create_gpml import gpml2gdf, gdf2gpml
from gprm.utils.proximity import rasterize_polygons
from gprm.utils.geometry import wrap_polygon_feature, wrap_polygon_features
import plate_driver as pdriver
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from shapely.ops import unary_union
from shapely.geometry import box, LineString
import pygmt
import metpy.calc as mpcalc

import topojson as tp
from xrspatial.polygonize import polygonize



def nudge_topological_polygons(snapshot, forces=None, reconstruction_time=0, delta_time=1):
    '''
    Compute where polygons would end up at some point in the future
    '''
    rpolygons = []
    for topology in snapshot.resolved_topologies:

        plate_id = topology.get_resolved_feature().get_reconstruction_plate_id()

        if not forces:
            stage_rotation = snapshot.rotation_model.get_rotation(reconstruction_time,
                                                                  plate_id,
                                                                  from_time=reconstruction_time+delta_time,
                                                                  anchor_plate_id=0)
            nudged_geometry = stage_rotation * topology.get_resolved_geometry()
            
        else:
            result = pdriver.get_torque_for_plate_polygon(snapshot, 
                                                          plate_id, 
                                                          reconstruction_time, 
                                                          delta_time,
                                                          forces)

            nudged_geometry = result[3] * topology.get_resolved_geometry()

        nfeature = topology.get_feature()
        nfeature.set_geometry(nudged_geometry)
        rpolygons.append(nfeature)
    
    return rpolygons


def nudge_polygons(polygon_features, stage_pole_dict):
    
    if isinstance(polygon_features, gpd.GeoDataFrame):
        polygon_features = gdf2gpml(polygon_features.set_crs(epsg=4326))
    
    rpolygons = []
    for polygon_feature in polygon_features:
        plate_id = polygon_feature.get_reconstruction_plate_id()
        nudged_geometry = stage_pole_dict.get(plate_id, pygplates.FiniteRotation()) * polygon_feature.get_geometry()
        nfeature = polygon_feature.clone()
        nfeature.set_geometry(nudged_geometry)
        rpolygons.append(nfeature)
        
    return rpolygons


def tessellate_polygons(polygons, max_dist_deg=0.2):

    new_features = []
    for feature in polygons:
        new_feature = feature.clone()
        new_feature.set_geometry(feature.get_geometry().to_tessellated(np.radians(max_dist_deg)))
        new_features.append(new_feature)

    return new_features


def plate_polygons_to_wrapped_gdf(polygons):
    rp_gdf = gpml2gdf(polygons)

    date_line_wrapper = pygplates.DateLineWrapper(0.0)
    rp_gdf = wrap_polygon_features(rp_gdf)#.apply(lambda x: wrap_polygon_feature(x, date_line_wrapper), axis=1)

    #gdf_reversed = rp_gdf.iloc[::-1].reset_index(drop=True)
    return rp_gdf


def resolve_overlaps(gdf, priority_col):
    gdf = gdf.sort_values(priority_col).reset_index(drop=True)
    gdf['geometry'] = gdf.geometry.make_valid()
    
    result_geoms = []
    accumulated = None
    
    for i in range(len(gdf)):
        geom = gdf.geometry.iloc[i]
        if accumulated is not None:
            geom = geom.difference(accumulated)
        result_geoms.append(geom)
        # Build cumulative union of original geometries
        accumulated = gdf.geometry.iloc[i] if accumulated is None else accumulated.union(gdf.geometry.iloc[i])
    
    result = gdf.copy()
    result.geometry = result_geoms
    
    return result


def update_plate_id_raster(clipped_polygons, raster_spacing=1):
    valid = clipped_polygons[clipped_polygons.geometry.notna() & ~clipped_polygons.geometry.is_empty]
    plate_id_raster = rasterize_polygons(valid, zval_field='PLATEID1', sampling=raster_spacing)
    plate_id_raster = plate_id_raster.where(plate_id_raster>0, np.nan)

    plate_id_raster = fill_nan_spherical(plate_id_raster)
    
    return plate_id_raster


def update_plate_polygons(plate_id_raster, raster_spacing=1, topo_simplify=1):
    # Origin offset by half a cell so polygon edges align with raster cell boundaries
    transform = (raster_spacing, 0, -(180 + raster_spacing / 2),
                 0, raster_spacing, -(90 + raster_spacing / 2))
    new_plate_polygons = polygonize(plate_id_raster,
                                    return_type='geopandas', transform=transform)

    topo = tp.Topology(new_plate_polygons, prequantize=False)
    new_plate_polygons_smooth = topo.toposimplify(topo_simplify).to_gdf()

    # Clip the GeoDataFrame
    new_plate_polygons_smooth = new_plate_polygons_smooth.clip(box(-180, -90, 180, 90))

    new_plate_polygons_smooth = new_plate_polygons_smooth.rename(columns={'DN':'PLATEID1'})
    
    return new_plate_polygons_smooth.set_crs(epsg=4326)


def velocity_grid(plate_id_raster, plate_parameters, delta_time=1, raster_spacing=1):
    '''
    Given a raster defining the plateids around the margins of a deformation zone, extract rigid plate 
    velocities at these points
    '''
    
    velocity_type = 'MagAzim'

    # We decide to isolate points from the grid where the 
    velpts = pygmt.grd2xyz(plate_id_raster, skiprows=True)
    velpts = velpts.rename(columns={velpts.columns[2]:'pid'})

    # All domain points and associated (magnitude, azimuth, inclination) velocities for the current time.
    all_domain_points = []
    all_velocities = []
    plate_ids = []

    for i,velocity_point in velpts.iterrows():
        velocity_domain_point = pygplates.PointOnSphere(velocity_point.y, velocity_point.x)

        # Get the stage rotation of partitioning plate from 'time + delta_time' to 'time'.
        pid_int = int(velocity_point.pid)
        equivalent_stage_rotation = plate_parameters.get(pid_int, pygplates.FiniteRotation())

        # Calculate velocity at the velocity domain point.
        # This is from 'time + delta_time' to 'time' on the partitioning plate.
        velocity_vectors = pygplates.calculate_velocities(
            [velocity_domain_point],
            equivalent_stage_rotation,
            delta_time)

        # Convert global 3D velocity vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
        velocities = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                                [velocity_domain_point],
                                velocity_vectors)
        all_velocities.append(velocities[0])

    vel_east=[]
    vel_north=[]
    for velocity_vector in all_velocities:
        vel_east.append(velocity_vector.get_y())
        vel_north.append(velocity_vector.get_x())
        
    velpts['vel_east'] = vel_east
    velpts['vel_north'] = vel_north
    
    velocity_grid_east = pygmt.xyz2grd(x=velpts.x, y=velpts.y, z=velpts.vel_east, region='d', spacing=raster_spacing)
    velocity_grid_north = pygmt.xyz2grd(x=velpts.x, y=velpts.y, z=velpts.vel_north, region='d', spacing=raster_spacing)
    
    return velocity_grid_east, velocity_grid_north


def compute_boundary_divergence(boundaries_gdf, stage_pole_dict, delta_time, plate_id_raster,
                                plate_polygons_gdf=None):
    """Compute per-arc-segment relative normal velocity.

    Explodes each boundary LineString into individual 2-vertex segments and computes
    divergence independently for each.  Returns a GeoDataFrame where each row is one
    arc segment with its own 'divergence' value (cm/yr).  Positive = ridge, negative
    = subduction.  Different segments of the same plate-pair boundary can have
    different signs, correctly capturing oblique or mixed boundaries.

    Sign: dot(centroid_pid1_xyz, arc_normal) > 0 → normal points toward pid1 → sign=+1.
    Falls back to raster probe when plate_polygons_gdf is not provided.
    """
    identity = pygplates.FiniteRotation()
    OFFSET_DEG = 0.5

    pid_raster = plate_id_raster.rename({k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                                         if k in plate_id_raster.dims})

    # Build centroid lookup: pid → geocentric (x,y,z) tuple — stable sign reference
    pid_to_centroid_xyz = {}
    if plate_polygons_gdf is not None:
        for _, prow in plate_polygons_gdf.iterrows():
            pid = prow.get('PLATEID1')
            if pid and pid not in pid_to_centroid_xyz:
                rep = prow.geometry.representative_point()
                pid_to_centroid_xyz[pid] = pygplates.PointOnSphere(rep.y, rep.x).to_xyz()

    new_rows = []
    for _, row in boundaries_gdf.iterrows():
        pid1 = row.get('plateid1_poly1')
        pid2 = row.get('plateid1_poly2')
        seg   = row.geometry
        coords = list(seg.coords)

        if not pid1 or not pid2 or pid1 == 0 or pid2 == 0 or len(coords) < 2:
            r = row.to_dict()
            r['divergence'] = 0.0
            new_rows.append(r)
            continue

        for i in range(len(coords) - 1):
            try:
                p1 = pygplates.PointOnSphere(coords[i][1],   coords[i][0])
                p2 = pygplates.PointOnSphere(coords[i+1][1], coords[i+1][0])
                arc = pygplates.GreatCircleArc(p1, p2)
                seg_normal = arc.get_great_circle_normal()
            except pygplates.IndeterminateGreatCircleArcNormalError:
                continue

            # Midpoint: average lat/lon (accurate enough for short arc segments)
            mid_lat = (coords[i][1] + coords[i+1][1]) / 2.0
            mid_lon = (coords[i][0] + coords[i+1][0]) / 2.0
            pt = pygplates.PointOnSphere(mid_lat, mid_lon)

            seg_nl = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                [pt], [seg_normal])[0]
            n_north, n_east = seg_nl.get_x(), seg_nl.get_y()

            # Sign: centroid hemisphere test (preferred) or raster probe (fallback)
            if pid1 in pid_to_centroid_xyz:
                c = pid_to_centroid_xyz[pid1]
                dot_val = seg_normal.get_x()*c[0] + seg_normal.get_y()*c[1] + seg_normal.get_z()*c[2]
                sign = 1 if dot_val > 0 else -1
            elif pid2 in pid_to_centroid_xyz:
                c = pid_to_centroid_xyz[pid2]
                dot_val = seg_normal.get_x()*c[0] + seg_normal.get_y()*c[1] + seg_normal.get_z()*c[2]
                sign = -1 if dot_val > 0 else 1
            else:
                probe_lon = float(np.clip(mid_lon + OFFSET_DEG * n_east, -180, 180))
                probe_lat = float(np.clip(mid_lat + OFFSET_DEG * n_north,  -90,  90))
                try:
                    probe_val = int(pid_raster.sel(lon=probe_lon, lat=probe_lat,
                                                   method='nearest').item())
                except Exception:
                    probe_val = 0
                if probe_val == pid1:
                    sign = 1
                elif probe_val == pid2:
                    sign = -1
                else:
                    continue

            v1 = pygplates.calculate_velocities([pt], stage_pole_dict.get(pid1, identity), delta_time)[0]
            v2 = pygplates.calculate_velocities([pt], stage_pole_dict.get(pid2, identity), delta_time)[0]
            v1l = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down([pt], [v1])[0]
            v2l = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down([pt], [v2])[0]

            # divergence = (v_pid1 - v_pid2) · n_toward_pid1
            raw = (v1l.get_y() - v2l.get_y()) * n_east + (v1l.get_x() - v2l.get_x()) * n_north

            r = row.to_dict()
            r['geometry'] = LineString([coords[i], coords[i+1]])
            r['divergence'] = sign * raw
            new_rows.append(r)

    if not new_rows:
        result = boundaries_gdf.copy()
        result['divergence'] = 0.0
        return result

    return gpd.GeoDataFrame(new_rows, geometry='geometry', crs=boundaries_gdf.crs)


def polygons_to_stage_poles(polygons_gdf, plate_id_raster, forces,
                            craton_raster=None,
                            stage_pole_dict=None,
                            delta_time=1,
                            divergence_threshold_for_subduction=-1.0,
                            divergence_threshold_for_ridge_push=1.0,
                            boundary_search_radius_deg=2,
                            polarity_dict=None):

    if polarity_dict is None:
        polarity_dict = {}

    # Craton pixel count per plate ID — used to determine overriding plate at new boundaries
    craton_pixel_counts = {}
    if craton_raster is not None:
        vals = craton_raster.values.flatten()
        valid = vals[vals > 0]
        if len(valid):
            unique, counts = np.unique(valid, return_counts=True)
            craton_pixel_counts = {int(k): int(v) for k, v in zip(unique, counts)}

    new_stage_pole_dict = {}
    torque_dict = {}
    
    gpml_plates = gdf2gpml(polygons_gdf.explode(index_parts=False).reset_index(drop=True))

    for i,plate in enumerate(gpml_plates):
        
        # Dead variables?
        px = []
        py = []
        pe = []
        pn = []

        #if i!=1:
        #    continue
        if plate.get_reconstruction_plate_id() not in torque_dict.keys():
            net_torque = pygplates.Vector3D(0, 0, 0)
        else:
            net_torque = torque_dict[plate.get_reconstruction_plate_id()]

        for geometry in plate.get_geometries():
            #print(geometry.get_orientation())
            if geometry.get_orientation() == 'clockwise':
                clockwise_flag = 1
            else:
                clockwise_flag = -1
            for segment in geometry.get_segments():
                arc_point = segment.get_arc_point(0.5).to_lat_lon()
                #my,mx = arc_point
                id1, id2 = get_two_nearest_ids(arc_point,
                                               plate_id_raster.rename({k: v for k, v in {'x': 'lon', 'y': 'lat'}.items() if k in plate_id_raster.dims}),
                                               radius_deg=boundary_search_radius_deg)
                #print(id1, id2)
                # insert logic to determine overriding plate
                downgoing_plate = id1
                overriding_plate = id2

                # Per-segment relative normal velocity: positive = diverging, negative = converging.
                # outward normal of the current plate = -clockwise_flag * gca_normal
                # (gca_normal = A×B = left-hand normal; inward for CW polygons, outward for CCW)
                # divergence = (v_neighbor - v_current) · outward_normal
                pt_obj = pygplates.PointOnSphere(arc_point[0], arc_point[1])
                pid_current = plate.get_reconstruction_plate_id()
                pid_neighbor = id2 if id1 == pid_current else id1
                if stage_pole_dict is not None and pid_neighbor and pid_neighbor > 0:
                    try:
                        _v_cur = pygplates.calculate_velocities([pt_obj], stage_pole_dict.get(pid_current, pygplates.FiniteRotation()), delta_time)[0]
                        _v_nbr = pygplates.calculate_velocities([pt_obj], stage_pole_dict.get(pid_neighbor, pygplates.FiniteRotation()), delta_time)[0]
                        _vcl = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down([pt_obj], [_v_cur])[0]
                        _vnl = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down([pt_obj], [_v_nbr])[0]
                        _seg_n = segment.get_great_circle_normal()
                        _seg_nl = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down([pt_obj], [_seg_n])[0]
                        # multiply by -clockwise_flag to convert gca_normal to outward normal
                        divergence_value = -clockwise_flag * (
                            (_vnl.get_y() - _vcl.get_y()) * _seg_nl.get_y()
                            + (_vnl.get_x() - _vcl.get_x()) * _seg_nl.get_x())
                    except pygplates.IndeterminateGreatCircleArcNormalError:
                        divergence_value = 0.0
                else:
                    divergence_value = 0.0

                continental_collision = False
                if craton_raster is not None:
                    craton_id1, craton_id2 = get_two_nearest_ids(arc_point,
                                                craton_raster.rename({k: v for k, v in {'x': 'lon', 'y': 'lat'}.items() if k in craton_raster.dims}),
                                                radius_deg=boundary_search_radius_deg)
                    if craton_id1 and craton_id2 and craton_id1 > 0 and craton_id2 > 0:
                        print(craton_id1, craton_id2)
                        continental_collision = True


                segment_torque = None
                if divergence_value<divergence_threshold_for_subduction:

                    # Record overriding/downgoing polarity for this plate pair (sticky — only set once).
                    # Skip when craton_raster is None (bootstrap step) so the craton rule is
                    # not contaminated by an arbitrary oceanic fallback before craton info is available.
                    if (id1 and id2 and id1 > 0 and id2 > 0
                            and craton_raster is not None):
                        pair = (min(id1, id2), max(id1, id2))
                        if pair not in polarity_dict:
                            area1 = craton_pixel_counts.get(id1, 0)
                            area2 = craton_pixel_counts.get(id2, 0)
                            if area1 > 0 and area2 == 0:
                                polarity_dict[pair] = id2   # id1 continental → overriding; id2 downgoing
                            elif area2 > 0 and area1 == 0:
                                polarity_dict[pair] = id1   # id2 continental → overriding; id1 downgoing
                            elif area1 > 0 and area2 > 0:
                                polarity_dict[pair] = id2 if area1 >= area2 else id1  # larger craton = overriding
                            else:
                                polarity_dict[pair] = id1   # both oceanic: arbitrary

                    if continental_collision:

                        normal_reversal = -1 * clockwise_flag
                        segment_torque = pdriver.get_torque_for_segment(segment, normal_reversal)

                        net_torque = net_torque - segment_torque * forces.get('continental_collision', 0)

                    elif downgoing_plate:

                        if 'slab_pull' in forces:

                            normal_reversal = -1 * clockwise_flag
                            segment_torque = pdriver.get_torque_for_segment(segment, normal_reversal)

                            net_torque = net_torque + segment_torque * forces['slab_pull']

                    elif overriding_plate:

                        if 'slab_suction' in forces:

                            normal_reversal = -1 * clockwise_flag
                            segment_torque = pdriver.get_torque_for_segment(segment, normal_reversal)

                            net_torque = net_torque + segment_torque * forces['slab_suction']


                elif divergence_value>divergence_threshold_for_ridge_push:

                    if 'ridge_push' in forces:

                        normal_reversal = 1 * clockwise_flag
                        segment_torque = pdriver.get_torque_for_segment(segment, normal_reversal)

                        net_torque = net_torque + segment_torque * forces['ridge_push']
                    
                # Deprecated section, not used
                if segment_torque is not None:
                    if segment.is_zero_length():
                        continue
                    normal = segment.get_great_circle_normal() * -normal_reversal     
                    segment_velocity = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                        [arc_point],
                        [normal])[0]
                    px.append(arc_point[1])
                    py.append(arc_point[0])
                    pe.append(segment_velocity.get_y())
                    pn.append(segment_velocity.get_x())


        if not net_torque.is_zero_magnitude():
            computed_stage_rotation = pygplates.FiniteRotation(pygplates.PointOnSphere(net_torque.to_normalised().to_xyz()).to_lat_lon(),
                                                               net_torque.get_magnitude())
        else:
            computed_stage_rotation = pygplates.FiniteRotation()
        
        torque_dict[plate.get_reconstruction_plate_id()] = net_torque
        new_stage_pole_dict[plate.get_reconstruction_plate_id()] = computed_stage_rotation

    return new_stage_pole_dict


def compute_timestep(rpolygons, forces, stage_pole_dict, delta_time=1, cratons=None, params=None, polarity_dict=None):
    if params is None:
        params = {}

    tessellation_deg          = params.get('tessellation_deg',           0.2)
    raster_spacing            = params.get('raster_spacing',             1)
    topo_simplify             = params.get('topo_simplify_tolerance',    1)
    div_thresh_sub            = params.get('div_threshold_subduction',  -1.0)
    div_thresh_ridge          = params.get('div_threshold_ridge',        1.0)
    boundary_search_radius    = params.get('boundary_search_radius_deg', 2)

    new_features = tessellate_polygons(rpolygons, max_dist_deg=tessellation_deg)

    rp_gdf = plate_polygons_to_wrapped_gdf(new_features)

    clipped_polygons = resolve_overlaps(rp_gdf, priority_col='PLATEID1')
    clipped_polygons = clipped_polygons.explode(index_parts=False).reset_index(drop=True)
    
    plate_id_raster = update_plate_id_raster(clipped_polygons, raster_spacing=raster_spacing)
    
    #import matplotlib.pyplot as plt
    #fig,ax = plt.subplots()
    #plate_id_raster.plot(ax=ax)
    #clipped_polygons.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.5)
    #fig.savefig('./output/plate_id_raster_1{:d}.png'.format(np.random.randint(0,100000)))

    # Force values within cratons to overrule any plate polygon changes
    # (from a geological perspective, cratons will overide subduction zones).
    if cratons is not None:
        wrapped_cratons = plate_polygons_to_wrapped_gdf(cratons)
        valid_cratons = wrapped_cratons[wrapped_cratons.geometry.notna() & ~wrapped_cratons.geometry.is_empty]
        craton_raster = rasterize_polygons(valid_cratons, zval_field='PLATEID1', sampling=raster_spacing)
        plate_id_raster = craton_raster.where(craton_raster>0, plate_id_raster)

        #fig,ax = plt.subplots()
        #craton_raster.plot(ax=ax)
        ##new_plate_polygons.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.5)
        #fig.savefig('./output/craton_raster_{:d}.png'.format(np.random.randint(0,100000)))

    else:
        craton_raster = None

    
    new_plate_polygons = update_plate_polygons(plate_id_raster,
                                               raster_spacing=raster_spacing,
                                               topo_simplify=topo_simplify)
    new_plate_polygons = new_plate_polygons.explode(index_parts=False).reset_index(drop=True)
    new_plate_polygons = new_plate_polygons[new_plate_polygons.geometry.geom_type.isin(['Polygon'])]
    #print(new_plate_polygons)

    #import matplotlib.pyplot as plt
    #fig,ax = plt.subplots()
    #plate_id_raster.plot(ax=ax)
    #new_plate_polygons.plot(ax=ax, facecolor='none', edgecolor='black', linewidth=0.5)
    #fig.savefig('./output/plate_id_raster_2{:d}.png'.format(np.random.randint(0,100000)))

    # Handle polygon splitting
    # TODO make sure that polygon splits retain consistency with craton polygons
    new_plate_polygons = resolve_duplicate_plateids(new_plate_polygons, next_plateid=None)

    # Handle reduction in forces for plates with continents?
    # BUT we should do that in 'polygons_to_stage_poles'?

    boundaries_gdf = topo_polygons_to_boundaries(new_plate_polygons)
    boundaries_gdf = compute_boundary_divergence(boundaries_gdf, stage_pole_dict, delta_time, plate_id_raster,
                                                 plate_polygons_gdf=new_plate_polygons)

    velocity_grid_east, velocity_grid_north = velocity_grid(plate_id_raster, stage_pole_dict, delta_time,
                                                             raster_spacing=raster_spacing)
    # Visualisation only — Cartesian approximation, not used for boundary classification
    div = mpcalc.divergence(velocity_grid_east, velocity_grid_north)

    stage_pole_dict = polygons_to_stage_poles(new_plate_polygons, plate_id_raster, forces,
                                              craton_raster,
                                              stage_pole_dict=stage_pole_dict,
                                              delta_time=delta_time,
                                              divergence_threshold_for_subduction=div_thresh_sub,
                                              divergence_threshold_for_ridge_push=div_thresh_ridge,
                                              boundary_search_radius_deg=boundary_search_radius,
                                              polarity_dict=polarity_dict)
    #
    #print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])
    #print([(key,item) for key,item in stage_pole_dict.items()])
    
    rasters_dict = {
        'plate_id':        plate_id_raster.rename({'x': 'lon', 'y': 'lat'}),
        'velocity_east':   velocity_grid_east,
        'velocity_north':  velocity_grid_north,
        'divergence':      div,
    }
    if craton_raster is not None:
        rasters_dict['craton'] = craton_raster.rename({'x': 'lon', 'y': 'lat'})

    return new_plate_polygons, \
           stage_pole_dict, \
           xr.Dataset(rasters_dict), \
           boundaries_gdf


def seed_ridge_tracers(boundaries_gdf, spacing_deg=2.0,
                       threshold=1.0, offset_deg=0.6):
    """Seed tracer point pairs offset perpendicular to divergent boundary segments.

    Two points are seeded per position (one on each side of the ridge axis) so
    that after advection they are carried apart on their respective plates.
    offset_deg must exceed half the plate_id raster cell size (~0.5 deg) so each
    point clearly falls on one plate in the nearest-neighbour lookup.
    boundaries_gdf must have a 'divergence' column (from compute_boundary_divergence).
    Uses pygplates for all geometric operations.
    """
    new_points = []
    for _, row in boundaries_gdf.iterrows():
        div_val = row.get('divergence', 0.0)
        if div_val <= threshold:
            continue

        coords = list(row.geometry.coords)
        if len(coords) < 2:
            continue

        # Process each consecutive vertex pair as its own great-circle arc so
        # that curved boundaries are followed correctly, not approximated as a
        # single chord between the first and last coordinates.
        for i in range(len(coords) - 1):
            try:
                p1 = pygplates.PointOnSphere(coords[i][1],   coords[i][0])
                p2 = pygplates.PointOnSphere(coords[i+1][1], coords[i+1][0])
                arc = pygplates.GreatCircleArc(p1, p2)
                arc_normal = arc.get_great_circle_normal()
            except pygplates.IndeterminateGreatCircleArcNormalError:
                continue

            # Tessellate this arc segment at spacing_deg for even seed density
            try:
                seg_pts = list(
                    pygplates.PolylineOnSphere([p1, p2])
                    .to_tessellated(np.radians(spacing_deg))
                    .get_points())
            except Exception:
                seg_pts = [p1, p2]

            for pt in seg_pts:
                lat, lon = pt.to_lat_lon()
                seg_nl = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(
                    [pt], [arc_normal])[0]
                n_north, n_east = seg_nl.get_x(), seg_nl.get_y()
                for sign in [+1, -1]:
                    new_lon = np.clip(lon + sign * offset_deg * n_east,  -180, 180)
                    new_lat = np.clip(lat + sign * offset_deg * n_north,  -90,  90)
                    new_points.append({'geometry': Point(new_lon, new_lat),
                                       'active': True, 'age': 0})
    if not new_points:
        return gpd.GeoDataFrame(columns=['geometry', 'active', 'age']).set_crs(epsg=4326)
    return gpd.GeoDataFrame(new_points).set_crs(epsg=4326)


def generate_seafloor_age_tracers(plate_polygons, boundaries_gdf,
                                   spacing_deg=5.0, spreading_rate_km_per_myr=50.0,
                                   div_threshold=0.00025):
    """Seed initial tracers inside each plate polygon with age estimated from ridge distance.

    For each polygon, finds its own divergent boundary segments, builds a BallTree from
    densified ridge points, then assigns age = haversine_dist_km / spreading_rate to each
    interior grid point.  Points with no divergent boundary on their plate are skipped.
    boundaries_gdf must have a 'divergence' column (from compute_boundary_divergence).
    """
    from sklearn.neighbors import BallTree
    from shapely.geometry import Point

    EARTH_RADIUS_KM = 6371.0
    all_points = []

    for _, poly_row in plate_polygons.iterrows():
        pid = poly_row['PLATEID1']
        geom = poly_row.geometry
        if geom is None or geom.is_empty:
            continue

        # Divergent segments that border this plate — densify via pygplates per arc
        div_ridge_pts = []
        for _, brow in boundaries_gdf.iterrows():
            if brow['plateid1_poly1'] != pid and brow['plateid1_poly2'] != pid:
                continue
            div_val = brow.get('divergence', 0.0)
            if div_val <= div_threshold:
                continue
            coords = list(brow.geometry.coords)
            for i in range(len(coords) - 1):
                try:
                    p1 = pygplates.PointOnSphere(coords[i][1],   coords[i][0])
                    p2 = pygplates.PointOnSphere(coords[i+1][1], coords[i+1][0])
                    for pt in (pygplates.PolylineOnSphere([p1, p2])
                               .to_tessellated(np.radians(spacing_deg / 2))
                               .get_points()):
                        lat, lon = pt.to_lat_lon()
                        div_ridge_pts.append([lat, lon])
                except Exception:
                    div_ridge_pts.append([coords[i][1], coords[i][0]])

        if not div_ridge_pts:
            continue

        ridge_arr = np.deg2rad(np.array(div_ridge_pts))
        tree = BallTree(ridge_arr, metric='haversine')

        # Equal-area grid: uniform in (lon, sin(lat)) so density is ~constant on the sphere
        minx, miny, maxx, maxy = geom.bounds
        lons = np.arange(minx + spacing_deg / 2, maxx, spacing_deg)
        sin_step = np.sin(np.deg2rad(spacing_deg))
        sin_min = np.sin(np.deg2rad(miny))
        sin_max = np.sin(np.deg2rad(maxy))
        sins = np.arange(sin_min + sin_step / 2, sin_max, sin_step)
        lats = np.degrees(np.arcsin(np.clip(sins, -1.0, 1.0)))

        for lat in lats:
            for lon in lons:
                pt = Point(lon, lat)
                if not geom.contains(pt):
                    continue
                query = np.deg2rad([[lat, lon]])
                dist_rad, _ = tree.query(query, k=1)
                dist_km = dist_rad[0, 0] * EARTH_RADIUS_KM
                age = dist_km / spreading_rate_km_per_myr
                all_points.append({'geometry': pt, 'active': True, 'age': age})

    if not all_points:
        return gpd.GeoDataFrame(columns=['geometry', 'active', 'age']).set_crs(epsg=4326)
    return gpd.GeoDataFrame(all_points).set_crs(epsg=4326)


def advect_tracers(tracers, plate_id_raster, stage_pole_dict, delta_time=1,
                   boundaries_gdf=None, div_threshold=0.0):
    """Move active tracers one step using their plate's stage rotation.

    If boundaries_gdf is provided, tracers that cross a convergent boundary during
    this step (raster plate_id changes from old to new for a convergent pair) are
    immediately deactivated.  This catches the case where a tracer slips into the
    overriding plate's raster territory before the BallTree proximity check fires.
    """
    if tracers is None or len(tracers) == 0:
        return tracers
    active_mask = tracers.active.astype(bool)
    active = tracers[active_mask].copy()
    inactive = tracers[~active_mask]

    # Pre-build frozenset of convergent plate-pair keys for O(1) lookup
    convergent_pairs = set()
    if boundaries_gdf is not None and 'divergence' in boundaries_gdf.columns:
        conv = boundaries_gdf[boundaries_gdf['divergence'] < div_threshold]
        for _, brow in conv.iterrows():
            p1 = brow.get('plateid1_poly1', 0)
            p2 = brow.get('plateid1_poly2', 0)
            if p1 and p2 and p1 != 0 and p2 != 0:
                convergent_pairs.add(frozenset([p1, p2]))

    rows = []
    for _, row in active.iterrows():
        lon, lat = row.geometry.x, row.geometry.y
        old_pid = int(plate_id_raster.sel(lon=lon, lat=lat, method='nearest').item())
        rotation = stage_pole_dict.get(old_pid, pygplates.FiniteRotation())
        pt = rotation * pygplates.PointOnSphere(lat, lon)
        new_lat, new_lon = pt.to_lat_lon()
        r = row.to_dict()
        r['geometry'] = Point(new_lon, new_lat)
        r['age'] = row.age + delta_time
        # Deactivate if the tracer crossed a convergent plate boundary this step
        if convergent_pairs and old_pid > 0:
            new_pid = int(plate_id_raster.sel(lon=new_lon, lat=new_lat, method='nearest').item())
            if new_pid != old_pid and new_pid > 0 and frozenset([old_pid, new_pid]) in convergent_pairs:
                r['active'] = False
        rows.append(r)
    if not rows:
        return tracers
    updated = gpd.GeoDataFrame(rows, geometry='geometry').set_crs(epsg=4326)
    return pd.concat([updated, inactive], ignore_index=True)


def deactivate_subducted_tracers(tracers, boundaries_gdf, threshold=-1.0, max_dist_deg=1.0):
    """Mark active tracers as inactive when within max_dist_deg of a convergent boundary.

    Convergent boundaries are those with divergence < threshold in boundaries_gdf.
    Distance is computed to the full boundary geometry (not just midpoints) via
    pygplates tessellation and a haversine BallTree, matching the approach used in
    generate_seafloor_age_tracers.
    """
    from sklearn.neighbors import BallTree

    if tracers is None or len(tracers) == 0:
        return tracers
    if boundaries_gdf is None or len(boundaries_gdf) == 0:
        return tracers
    if 'divergence' not in boundaries_gdf.columns:
        return tracers

    tracers = tracers.copy()
    active_mask = tracers.active.astype(bool)
    if not active_mask.any():
        return tracers

    conv_boundaries = boundaries_gdf[boundaries_gdf['divergence'] < threshold]
    if len(conv_boundaries) == 0:
        return tracers

    # Tessellate convergent boundary segments into dense point arrays (pygplates spherical)
    SAMPLE_DEG = 0.5
    conv_pts = []
    for _, brow in conv_boundaries.iterrows():
        coords = list(brow.geometry.coords)
        if len(coords) < 2:
            continue
        try:
            pts = [pygplates.PointOnSphere(c[1], c[0]) for c in coords]
            polyline = pygplates.PolylineOnSphere(pts)
            for pt in polyline.to_tessellated(np.radians(SAMPLE_DEG)).get_points():
                lat, lon = pt.to_lat_lon()
                conv_pts.append([lat, lon])
        except Exception:
            for c in coords:
                conv_pts.append([c[1], c[0]])

    if not conv_pts:
        return tracers

    tree = BallTree(np.deg2rad(np.array(conv_pts)), metric='haversine')

    active_tracers = tracers[active_mask]
    coords_rad = np.deg2rad(np.column_stack([
        active_tracers.geometry.y.values,
        active_tracers.geometry.x.values,
    ]))
    dists_rad, _ = tree.query(coords_rad, k=1)
    dists_deg = np.degrees(dists_rad.flatten())

    deactivate_mask = dists_deg < max_dist_deg
    tracers.loc[active_tracers.index[deactivate_mask], 'active'] = False
    return tracers


def craton_counts_from_raster(craton_raster):
    """Return {plate_id: pixel_count} from a craton xarray DataArray.  Only counts pixels > 0."""
    if craton_raster is None:
        return {}
    vals = craton_raster.values.flatten()
    counts = {}
    for v in vals[vals > 0]:
        pid = int(v)
        counts[pid] = counts.get(pid, 0) + 1
    return counts


def resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts=None):
    """Return the downgoing plate ID for a convergent boundary between p1 and p2.

    Priority:
      1. Craton rule: plate with cratons is overriding; the other is downgoing.
         If both have cratons, the one with more pixels is overriding.
      2. polarity_dict (sticky result from previous steps).
      3. Fallback: p1 downgoing (arbitrary, consistent).
    """
    if p1 == 0 or p2 == 0:
        return p1
    if craton_counts:
        a1 = craton_counts.get(p1, 0)
        a2 = craton_counts.get(p2, 0)
        if a1 > 0 and a2 == 0:
            return p2           # p1 continental → overriding; p2 downgoing
        if a2 > 0 and a1 == 0:
            return p1           # p2 continental → overriding; p1 downgoing
        if a1 > 0 and a2 > 0:
            return p2 if a1 >= a2 else p1   # larger craton is overriding
    pair = (min(p1, p2), max(p1, p2))
    return polarity_dict.get(pair, p1)


def generate_synthetic_earthquakes(boundaries_gdf, polarity_dict=None, plate_polygons_gdf=None,
                                    craton_counts=None, plate_id_raster=None,
                                    n_candidates=100000,
                                    ridge_max_dist_deg=2.0,
                                    ridge_max_depth_km=50.0,
                                    slab_dip_deg=45.0,
                                    slab_max_depth_km=600.0,
                                    div_threshold_ridge=0.1,
                                    div_threshold_subduction=-0.1,
                                    ridge_eq_prob=1.0,
                                    subduction_eq_prob=1.0):
    """Generate a synthetic earthquake catalog from plate boundary geometry.

    Ridge earthquakes: random depth in [0, ridge_max_depth_km], symmetric about the
    ridge axis within ridge_max_dist_deg.

    Subduction earthquakes: overriding-plate side only (hemisphere test using
    polarity_dict + arc normal).  Depth increases linearly with horizontal distance
    from the boundary: depth = dist_km * tan(slab_dip_deg), capped at
    slab_max_depth_km.  Maximum horizontal extent = slab_max_depth_km / tan(dip).

    Returns a GeoDataFrame with columns:
        geometry, depth_km, boundary_type, plateid1, plateid2,
        downgoing_plate, overriding_plate
    """
    from sklearn.neighbors import BallTree

    if boundaries_gdf is None or 'divergence' not in boundaries_gdf.columns:
        return gpd.GeoDataFrame(
            columns=['geometry', 'depth_km', 'boundary_type',
                     'plateid1', 'plateid2', 'downgoing_plate', 'overriding_plate']
        ).set_crs(epsg=4326)

    if polarity_dict is None:
        polarity_dict = {}

    slab_dip_rad = np.radians(slab_dip_deg)
    slab_max_dist_km  = slab_max_depth_km / np.tan(slab_dip_rad)
    slab_max_dist_deg = slab_max_dist_km / 111.195
    max_query_deg = max(ridge_max_dist_deg, slab_max_dist_deg)

    # Centroid lookup for hemisphere test (which side is overriding plate)
    pid_to_centroid_xyz = {}
    if plate_polygons_gdf is not None:
        for _, prow in plate_polygons_gdf.iterrows():
            pid = prow.get('PLATEID1')
            if pid and pid not in pid_to_centroid_xyz:
                rep = prow.geometry.representative_point()
                pid_to_centroid_xyz[pid] = pygplates.PointOnSphere(rep.y, rep.x).to_xyz()

    # Tessellate active boundary segments; track iloc index into boundaries_gdf
    SAMPLE_DEG = 0.3
    pts_latlon = []
    pts_seg_iloc = []

    for iloc_i in range(len(boundaries_gdf)):
        row = boundaries_gdf.iloc[iloc_i]
        div = row.get('divergence', 0.0)
        if not (div > div_threshold_ridge or div < div_threshold_subduction):
            continue
        coords = list(row.geometry.coords)
        if len(coords) < 2:
            continue
        try:
            p1 = pygplates.PointOnSphere(coords[0][1], coords[0][0])
            p2 = pygplates.PointOnSphere(coords[1][1], coords[1][0])
            for pt in pygplates.PolylineOnSphere([p1, p2]).to_tessellated(
                    np.radians(SAMPLE_DEG)).get_points():
                lat, lon = pt.to_lat_lon()
                pts_latlon.append([lat, lon])
                pts_seg_iloc.append(iloc_i)
        except Exception:
            pass

    if not pts_latlon:
        return gpd.GeoDataFrame(
            columns=['geometry', 'depth_km', 'boundary_type',
                     'plateid1', 'plateid2', 'downgoing_plate', 'overriding_plate']
        ).set_crs(epsg=4326)

    # Uniform random points on sphere — Marsaglia's method, vectorised
    xyz_batches = []
    total = 0
    while total < n_candidates:
        batch = max(n_candidates * 4, 10000)
        x1 = np.random.uniform(-1, 1, batch)
        x2 = np.random.uniform(-1, 1, batch)
        s = x1**2 + x2**2
        ok = s < 1
        x1, x2, s = x1[ok], x2[ok], s[ok]
        sq = np.sqrt(1 - s)
        xyz_batches.append(np.column_stack([2*x1*sq, 2*x2*sq, 1 - 2*s]))
        total += int(ok.sum())
    xyz = np.vstack(xyz_batches)[:n_candidates]
    lats = np.degrees(np.arcsin(np.clip(xyz[:, 2], -1, 1)))
    lons = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))

    # Find nearest active boundary point for every candidate
    tree = BallTree(np.deg2rad(np.array(pts_latlon)), metric='haversine')
    dists_rad, nn_idx = tree.query(np.deg2rad(np.column_stack([lats, lons])), k=1)
    dists_deg = np.degrees(dists_rad.flatten())
    nn_idx    = nn_idx.flatten()

    earthquakes = []

    for ci in np.where(dists_deg < max_query_deg)[0]:
        lon_c    = float(lons[ci])
        lat_c    = float(lats[ci])
        dist_deg = float(dists_deg[ci])
        dist_km  = dist_deg * 111.195

        seg_row = boundaries_gdf.iloc[pts_seg_iloc[nn_idx[ci]]]
        div     = seg_row.get('divergence', 0.0)

        if div > div_threshold_ridge:
            # Ridge: shallow, symmetric about axis
            if dist_deg > ridge_max_dist_deg:
                continue
            if ridge_eq_prob < 1.0 and np.random.random() > ridge_eq_prob:
                continue
            earthquakes.append({
                'geometry':       Point(lon_c, lat_c),
                'depth_km':       round(float(np.random.uniform(0.0, ridge_max_depth_km)), 1),
                'boundary_type':  'ridge',
                'plateid1':       int(seg_row.get('plateid1_poly1', 0) or 0),
                'plateid2':       int(seg_row.get('plateid1_poly2', 0) or 0),
                'downgoing_plate':  0,
                'overriding_plate': 0,
            })

        elif div < div_threshold_subduction:
            # Subduction: overriding-plate side only, depth from slab dip
            p1 = int(seg_row.get('plateid1_poly1', 0) or 0)
            p2 = int(seg_row.get('plateid1_poly2', 0) or 0)
            if p1 == 0 or p2 == 0:
                continue

            dg        = resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts)
            overriding = p2 if dg == p1 else p1

            # Arc normal for hemisphere test
            coords = list(seg_row.geometry.coords)
            if len(coords) < 2:
                continue
            try:
                arc = pygplates.GreatCircleArc(
                    pygplates.PointOnSphere(coords[0][1], coords[0][0]),
                    pygplates.PointOnSphere(coords[1][1], coords[1][0]))
                seg_normal = arc.get_great_circle_normal()
            except pygplates.IndeterminateGreatCircleArcNormalError:
                continue

            # Is this candidate on the overriding-plate side?
            pt_xyz = pygplates.PointOnSphere(lat_c, lon_c).to_xyz()
            dot_pt = (seg_normal.get_x()*pt_xyz[0]
                    + seg_normal.get_y()*pt_xyz[1]
                    + seg_normal.get_z()*pt_xyz[2])

            if overriding in pid_to_centroid_xyz:
                c      = pid_to_centroid_xyz[overriding]
                dot_ov = (seg_normal.get_x()*c[0]
                        + seg_normal.get_y()*c[1]
                        + seg_normal.get_z()*c[2])
                if (dot_pt > 0) != (dot_ov > 0):
                    continue        # wrong side → downgoing plate, skip

            depth_km = dist_km * np.tan(slab_dip_rad)
            if depth_km > slab_max_depth_km:
                continue

            # Sanity check: reject if the point sits inside the downgoing plate's
            # surface expression — the slab dips away from that plate, not under it.
            if plate_id_raster is not None:
                try:
                    raster_pid = int(plate_id_raster.sel(
                        lon=lon_c, lat=lat_c, method='nearest').item())
                    if raster_pid == dg:
                        continue
                except Exception:
                    pass

            if subduction_eq_prob < 1.0 and np.random.random() > subduction_eq_prob:
                continue

            earthquakes.append({
                'geometry':        Point(lon_c, lat_c),
                'depth_km':        round(float(depth_km), 1),
                'boundary_type':   'subduction',
                'plateid1':        p1,
                'plateid2':        p2,
                'downgoing_plate':  dg,
                'overriding_plate': overriding,
            })

    if not earthquakes:
        return gpd.GeoDataFrame(
            columns=['geometry', 'depth_km', 'boundary_type',
                     'plateid1', 'plateid2', 'downgoing_plate', 'overriding_plate']
        ).set_crs(epsg=4326)

    return gpd.GeoDataFrame(earthquakes, geometry='geometry').set_crs(epsg=4326)


def write_rotation_file(rotation_records, path, fixed_plate_id=0):
    """Write accumulated finite rotations to a GPlates .rot file.

    rotation_records is a list of (time_myr, plate_id, FiniteRotation) tuples
    produced by accumulating stage poles across timesteps.  A t=0 identity
    entry is written automatically for every plate that appears in the records.
    """
    from collections import defaultdict

    by_plate = defaultdict(list)
    for time_myr, pid, rot in rotation_records:
        by_plate[pid].append((time_myr, rot))

    lines = []
    for pid in sorted(by_plate.keys()):
        # t=0 identity line (required by GPlates)
        lines.append('{:d}   {:8.4f}   {:9.4f}   {:10.4f}   {:9.4f}   {:d}'.format(
            pid, 0.0, 0.0, 0.0, 0.0, fixed_plate_id))
        for time_myr, rot in sorted(by_plate[pid], key=lambda x: x[0]):
            if rot.represents_identity_rotation():
                lat, lon, angle = 0.0, 0.0, 0.0
            else:
                lat, lon, angle = rot.get_lat_lon_euler_pole_and_angle_degrees()
            lines.append('{:d}   {:8.4f}   {:9.4f}   {:10.4f}   {:9.4f}   {:d}'.format(
                pid, time_myr, lat, lon, angle, fixed_plate_id))

    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def remove_boundary_segments(arc, lon_threshold=179.9, lat_threshold=89.9):
    """Remove segments along dateline/poles but keep connecting segments"""
    arc = np.asarray(arc)
    
    # Check which points are on boundaries
    on_dateline = np.abs(np.abs(arc[:, 0]) - 180) < (180 - lon_threshold)
    on_pole = np.abs(arc[:, 1]) > lat_threshold
    on_boundary = on_dateline | on_pole
    
    # Keep segments where at least one endpoint is NOT on boundary
    keep_segments = []
    for i in range(len(arc) - 1):
        if not (on_boundary[i] and on_boundary[i+1]):
            keep_segments.append(i)
    
    if not keep_segments:
        return []
    
    # Group consecutive segments into sub-arcs
    sub_arcs = []
    current = [arc[keep_segments[0]]]
    
    for i in range(len(keep_segments)):
        current.append(arc[keep_segments[i] + 1])
        if i < len(keep_segments) - 1 and keep_segments[i+1] != keep_segments[i] + 1:
            sub_arcs.append(np.array(current))
            current = [arc[keep_segments[i+1]]]
    
    sub_arcs.append(np.array(current))
    return sub_arcs


def topo_polygons_to_boundaries(plate_polygons):
    
    topo = tp.Topology(plate_polygons, prequantize=False, toposimplify=False)

    topo_dict = topo.to_dict()
    arcs = topo_dict['arcs']

    # Build arc-to-polygon mapping with correct decoding
    arc_to_polys = {}
    feat_to_plateid = {}  # Store PLATEID1 for each feature

    for obj_data in topo_dict['objects'].values():
        for feat_idx, feat in enumerate(obj_data['geometries']):
            # Store the PLATEID1 for this feature
            feat_to_plateid[feat_idx] = feat['properties'].get('PLATEID1')

            for ring in feat['arcs']:
                for arc_idx in ring:
                    # Decode TopoJSON arc index
                    if arc_idx < 0:
                        abs_idx = -arc_idx - 1
                    else:
                        abs_idx = arc_idx
                    arc_to_polys.setdefault(abs_idx, []).append(feat_idx)

    # Apply to all arcs with polygon IDs and PLATEID1
    boundaries = []
    for arc_idx, arc in enumerate(arcs):
        sub_arcs = remove_boundary_segments(arc)

        poly_list = arc_to_polys.get(arc_idx, [])
        unique_polys = list(set(poly_list))

        for sub_arc in sub_arcs:
            boundaries.append({
                'geometry': LineString(sub_arc),
                'poly1': unique_polys[0] if len(unique_polys) > 0 else None,
                'poly2': unique_polys[1] if len(unique_polys) > 1 else None,
                'plateid1_poly1': feat_to_plateid.get(unique_polys[0]) if len(unique_polys) > 0 else None,
                'plateid1_poly2': feat_to_plateid.get(unique_polys[1]) if len(unique_polys) > 1 else None
            })

    boundaries_gdf = gpd.GeoDataFrame(boundaries, geometry='geometry').set_crs(epsg=4326)

    return boundaries_gdf





def touches_dateline(geom, tol=1e-6):
    """Check if polygon has vertices at lon ≈ ±180"""
    coords = np.array(geom.exterior.coords)
    return np.any(np.abs(np.abs(coords[:, 0]) - 180) < tol)

def resolve_duplicate_plateids(gdf, next_plateid=None):
    if next_plateid is None:
        next_plateid = gdf['PLATEID1'].max() + 1
    
    gdf = gdf.copy()
    
    for plateid, group in gdf.groupby('PLATEID1'):
        if len(group) <= 1:
            continue
        
        # Identify which can be reassigned (not touching dateline)
        for idx in group.index:
            geom = gdf.loc[idx, 'geometry']
            if geom.geom_type == 'MultiPolygon':
                on_dateline = any(touches_dateline(p) for p in geom.geoms)
            else:
                on_dateline = touches_dateline(geom)
            
            gdf.loc[idx, '_on_dateline'] = on_dateline
        
        # Reassign all but one (prefer keeping dateline polygons unchanged)
        candidates = group.index[~gdf.loc[group.index, '_on_dateline']]
        keep_one = True
        
        for idx in candidates:
            if keep_one:
                keep_one = False
                continue
            gdf.loc[idx, 'PLATEID1'] = next_plateid
            next_plateid += 1
    
    gdf = gdf.drop(columns=['_on_dateline'], errors='ignore')
    return gdf



from pyproj import Geod
from shapely.geometry import Polygon, MultiPolygon
from shapely import get_coordinates, set_coordinates
import numpy as np

geod = Geod(ellps='WGS84')

def densify_great_circle(coords, max_dist_deg):
    if len(coords) < 2:
        return coords
    
    max_dist_m = max_dist_deg * 111_320
    
    new_coords = []
    for i in range(len(coords) - 1):
        lon1, lat1 = coords[i]
        lon2, lat2 = coords[i + 1]
        
        # Skip degenerate segments
        if (lon1, lat1) == (lon2, lat2):
            new_coords.append((lon1, lat1))
            continue
        
        _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
        
        if dist <= max_dist_m:
            new_coords.append((lon1, lat1))
        else:
            n_segments = int(np.ceil(dist / max_dist_m))
            points = geod.npts(lon1, lat1, lon2, lat2, n_segments - 1)
            new_coords.append((lon1, lat1))
            new_coords.extend(points)
    
    new_coords.append(coords[-1])
    return new_coords

def densify_polygon(geom, max_dist_deg=1.0):
    if geom is None or geom.is_empty:
        return geom
    
    if geom.geom_type == 'Polygon':
        exterior = densify_great_circle(list(geom.exterior.coords), max_dist_deg)
        interiors = []
        for ring in geom.interiors:
            if len(ring.coords) >= 4:  # valid ring needs at least 4 coords
                interiors.append(densify_great_circle(list(ring.coords), max_dist_deg))
        return Polygon(exterior, interiors)
    
    elif geom.geom_type == 'MultiPolygon':
        return MultiPolygon([densify_polygon(p, max_dist_deg) for p in geom.geoms])
    
    return geom  # return unchanged for non-polygon types



#####################################################################
from sklearn.neighbors import BallTree
import numpy as np

def fill_nan_spherical(da):
    vals = da.values
    lon, lat = np.meshgrid(da.x, da.y)
    
    valid = ~np.isnan(vals)
    invalid = np.isnan(vals)
    
    # BallTree with haversine needs (lat, lon) in radians
    valid_coords = np.deg2rad(np.column_stack([lat[valid], lon[valid]]))
    invalid_coords = np.deg2rad(np.column_stack([lat[invalid], lon[invalid]]))
    
    tree = BallTree(valid_coords, metric='haversine')
    _, idx = tree.query(invalid_coords, k=1)
    
    filled = vals.copy()
    filled[invalid] = vals[valid][idx.flatten()]
    
    return da.copy(data=filled)





#####################################################################
import numpy as np
from scipy.spatial import SphericalVoronoi
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from pyproj import Geod
import geopandas as gpd
from collections import defaultdict
from shapely.geometry import Point, MultiPoint
from shapely.ops import voronoi_diagram
from shapely.strtree import STRtree

def lonlat_to_xyz(lon, lat):
    lon_r, lat_r = np.radians(lon), np.radians(lat)
    return np.column_stack([
        np.cos(lat_r) * np.cos(lon_r),
        np.cos(lat_r) * np.sin(lon_r),
        np.sin(lat_r),
    ])

def xyz_to_lonlat(xyz):
    lon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    lat = np.degrees(np.arcsin(np.clip(xyz[:, 2], -1, 1)))
    return lon, lat


def densify_geodesic(geom, spacing_m=50_000):
    """Densify polygon boundaries using geodesic interpolation."""
    geod = Geod(ellps="WGS84")

    def densify_coords(coords):
        result = []
        for i in range(len(coords) - 1):
            lon1, lat1 = coords[i][:2]
            lon2, lat2 = coords[i + 1][:2]
            _, _, dist = geod.inv(lon1, lat1, lon2, lat2)
            n = max(int(dist / spacing_m), 1)
            result.append((lon1, lat1))
            result.extend(geod.npts(lon1, lat1, lon2, lat2, n))
        result.append(coords[-1][:2])
        return result

    if geom.geom_type == "Polygon":
        ext = densify_coords(list(geom.exterior.coords))
        ints = [densify_coords(list(r.coords)) for r in geom.interiors]
        return Polygon(ext, ints)
    elif geom.geom_type == "MultiPolygon":
        return MultiPolygon([densify_geodesic(p, spacing_m) for p in geom.geoms])
    return geom


def extract_boundary_points(geom):
    coords = []
    polys = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    for poly in polys:
        coords.extend(list(poly.exterior.coords)[:-1])
        for ring in poly.interiors:
            coords.extend(list(ring.coords)[:-1])
    return [(c[0], c[1]) for c in coords]


def voronoi_region_to_polygon(vertices_xyz):
    """Convert SphericalVoronoi region vertices (Nx3) to a lon/lat Polygon."""
    lon, lat = xyz_to_lonlat(vertices_xyz)

    # Detect dateline crossing and shift to [0, 360] if needed
    dlon = np.abs(np.diff(np.append(lon, lon[0])))
    if np.any(dlon > 180):
        lon = lon % 360

    ring = list(zip(lon, lat))
    ring.append(ring[0])
    poly = Polygon(ring)
    return poly.buffer(0) if not poly.is_valid else poly


from shapely.geometry import box

VALID_BOUNDS = box(-180, -90, 180, 90)

def expand_to_fill_gaps_spherical(gdf, spacing_m=50_000):
    all_union = unary_union(gdf.geometry)
    envelope = VALID_BOUNDS  # strict geographic bounds
    gap = envelope.difference(all_union)

    if gap.is_empty:
        return gdf.copy()

    # 1. Densify boundaries, deduplicate shared-boundary points
    unique_points = []
    unique_labels = []
    seen = {}
    BOUNDARY_LABEL = -1  # sentinel for edge-of-world points

    for i, geom in enumerate(gdf.geometry):
        dense = densify_geodesic(geom, spacing_m)
        pts = extract_boundary_points(dense)
        for lon, lat in pts:
            key = (round(lon, 8), round(lat, 8))
            if key not in seen:
                seen[key] = len(unique_points)
                unique_points.append(Point(lon, lat))
                unique_labels.append({i})
            else:
                unique_labels[seen[key]].add(i)

    # 2. Add dense boundary generators along dateline and poles
    #    These act as hard walls — their Voronoi cells won't be assigned to any polygon
    spacing_deg = spacing_m / 111_000  # rough conversion
    for lat in np.arange(-90, 90, spacing_deg):
        for lon in [-180.0, 180.0]:
            key = (round(lon, 8), round(lat, 8))
            if key not in seen:
                seen[key] = len(unique_points)
                unique_points.append(Point(lon, lat))
                unique_labels.append({BOUNDARY_LABEL})
    for lon in np.arange(-180, 180, spacing_deg):
        for lat in [-90.0, 90.0]:
            key = (round(lon, 8), round(lat, 8))
            if key not in seen:
                seen[key] = len(unique_points)
                unique_points.append(Point(lon, lat))
                unique_labels.append({BOUNDARY_LABEL})

    # 3. Planar Voronoi
    mp = MultiPoint(unique_points)
    vor = voronoi_diagram(mp, envelope=envelope)

    # 4. Assign every cell to nearest generator
    point_tree = STRtree(unique_points)
    region_polys = defaultdict(list)

    for cell in vor.geoms:
        idx = point_tree.nearest(cell.representative_point())
        for lbl in unique_labels[idx]:
            if lbl != BOUNDARY_LABEL:  # skip edge-of-world cells
                region_polys[lbl].append(cell)

    # 5. Clip Voronoi territory to gap only, union with original, clip to valid bounds
    result = gdf.copy()
    for i in range(len(gdf)):
        original = gdf.geometry.iloc[i]
        if region_polys[i]:
            territory = unary_union(region_polys[i])
            gap_portion = territory.intersection(gap)
            expanded = unary_union([original, gap_portion]).intersection(VALID_BOUNDS)
            result.loc[result.index[i], "geometry"] = expanded

    # 6. Validate and repair geometries
    result["geometry"] = result.geometry.make_valid()
    # Remove any zero-width slivers from make_valid producing GeometryCollections
    from shapely.geometry import GeometryCollection, MultiPolygon
    def extract_polygons(geom):
        if geom.geom_type in ("Polygon", "MultiPolygon"):
            return geom
        elif geom.geom_type == "GeometryCollection":
            polys = [g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
            if polys:
                return unary_union(polys)
        return geom
    
    result["geometry"] = result.geometry.apply(extract_polygons)

    return result


import numpy as np
import geopandas as gpd
from shapely.geometry import LineString

EARTH_RADIUS_KM = 6371.0

def haversine(lon1, lat1, lon2, lat2):
    """Haversine distance in degrees (approximate, for comparison purposes)."""
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))  # returns great-circle distance in degrees


def is_dateline_segment(p1, p2, tol=1e-6):
    """Check if both endpoints lie on the dateline (lon ≈ ±180)."""
    return (abs(abs(p1[0]) - 180) < tol) and (abs(abs(p2[0]) - 180) < tol)


def get_two_nearest_ids(coordinate, da, radius_deg=2.0):
    """
    For a line segment, find the two nearest unique integers
    from an xarray DataArray within a radius of the segment midpoint.
    
    da: xarray.DataArray with dims (lat, lon) containing integers.
    radius_deg: search radius in degrees (great-circle).
    """
    #mid = segment.interpolate(0.5, normalized=True)
    #mx, my = mid.x, mid.y
    my, mx = coordinate

    # Generous bounding box (haversine will do the real filtering)
    bbox_buf = radius_deg * 1.5
    subset = da.sel(
        lon=slice(mx - bbox_buf, mx + bbox_buf),
        lat=slice(my - bbox_buf, my + bbox_buf)  # flip if lat descending
    )

    if subset.size == 0:
        return None, None

    vals = subset.values
    lons, lats = np.meshgrid(subset.lon.values, subset.lat.values)

    # Compute haversine distances from midpoint to all grid nodes
    dists_deg = haversine(mx, my, lons, lats)

    # Mask to within radius
    within = dists_deg <= radius_deg
    if not within.any():
        return None, None

    vals_within = vals[within]
    dists_within = dists_deg[within]

    # Remove NaNs
    valid = ~np.isnan(vals_within)
    vals_within = vals_within[valid].astype(int)
    dists_within = dists_within[valid]

    if len(vals_within) == 0:
        return None, None

    # For each unique integer, find minimum distance
    unique_vals = np.unique(vals_within)
    distances = {}
    for uv in unique_vals:
        distances[uv] = dists_within[vals_within == uv].min()

    sorted_ids = sorted(distances, key=distances.get)
    if len(sorted_ids) >= 2:
        return sorted_ids[0], sorted_ids[1]
    elif len(sorted_ids) == 1:
        return sorted_ids[0], None
    return None, None


def process_geodataframe(gdf, da, radius_deg=2.0, skip_dateline=True):
    results = []
    for idx, row in gdf.iterrows():
        coords = list(row.geometry.exterior.coords)
        for i in range(len(coords) - 1):
            p1, p2 = coords[i], coords[i + 1]

            if is_dateline_segment(p1, p2):
                if skip_dateline:
                    continue
                # If not skipping, flag it
                results.append({
                    'polygon_idx': idx,
                    'segment_idx': i,
                    'segment': LineString([p1, p2]),
                    'dateline': True,
                    'id_1': None,
                    'id_2': None
                })
                continue

            seg = LineString([p1, p2])
            mid = seg.interpolate(0.5, normalized=True)
            id1, id2 = get_two_nearest_ids((mid.y,mid.x), da, radius_deg)
            results.append({
                'polygon_idx': idx,
                'segment_idx': i,
                'segment': seg,
                'dateline': False,
                'id_1': id1,
                'id_2': id2
            })
            
    return results




