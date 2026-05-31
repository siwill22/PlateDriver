import time
import math
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
    valid = clipped_polygons.dropna(subset=['geometry'])
    valid = valid[~valid.geometry.is_empty]
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

    # Fix any invalid geometries produced by toposimplify before clipping
    new_plate_polygons_smooth.geometry = new_plate_polygons_smooth.geometry.make_valid()
    new_plate_polygons_smooth = new_plate_polygons_smooth.clip(box(-180, -90, 180, 90))

    new_plate_polygons_smooth = new_plate_polygons_smooth.rename(columns={'DN':'PLATEID1'})
    
    return new_plate_polygons_smooth.set_crs(epsg=4326)


def velocity_grid(plate_id_raster, plate_parameters, delta_time=1, raster_spacing=1):
    '''
    Given a raster defining the plateids around the margins of a deformation zone, extract rigid plate
    velocities at these points.

    Grouped per plate: one pygplates call per unique plate ID instead of one per pixel,
    reducing ~259,200 calls to ~15-20.
    '''
    lons = plate_id_raster.x.values   # (nlon,)
    lats = plate_id_raster.y.values   # (nlat,)
    pids = plate_id_raster.values     # (nlat, nlon)

    nlat, nlon = pids.shape
    vel_east  = np.zeros((nlat, nlon), dtype=float)
    vel_north = np.zeros((nlat, nlon), dtype=float)

    lon2d, lat2d = np.meshgrid(lons, lats)   # (nlat, nlon)

    for pid_val in np.unique(pids):
        if np.isnan(pid_val):
            continue
        pid_int  = int(pid_val)
        rotation = plate_parameters.get(pid_int, pygplates.FiniteRotation())

        rows, cols = np.where(pids == pid_val)
        lats_m = lat2d[rows, cols]
        lons_m = lon2d[rows, cols]

        pts = [pygplates.PointOnSphere(float(la), float(lo))
               for la, lo in zip(lats_m, lons_m)]

        vecs = pygplates.calculate_velocities(pts, rotation, delta_time)
        ned  = pygplates.LocalCartesian.convert_from_geocentric_to_north_east_down(pts, vecs)

        for k, nv in enumerate(ned):
            vel_east[rows[k], cols[k]]  = nv.get_y()
            vel_north[rows[k], cols[k]] = nv.get_x()

    coords = {'y': lats, 'x': lons}
    da_east  = xr.DataArray(vel_east,  dims=['y', 'x'], coords=coords)
    da_north = xr.DataArray(vel_north, dims=['y', 'x'], coords=coords)
    return da_east, da_north


def compute_boundary_divergence(boundaries_gdf, stage_pole_dict, delta_time, plate_id_raster,
                                plate_polygons_gdf=None):
    """Compute per-arc-segment relative normal velocity.

    Explodes each boundary LineString into individual 2-vertex segments and computes
    divergence independently for each.  Returns a GeoDataFrame where each row is one
    arc segment with its own 'divergence' value (cm/yr).  Positive = ridge, negative
    = subduction.  Different segments of the same plate-pair boundary can have
    different signs, correctly capturing oblique or mixed boundaries.

    Sign is determined by probing the plate-ID raster on both sides of each segment
    midpoint — more robust than a centroid hemisphere test for large or concave plates.
    """
    identity = pygplates.FiniteRotation()
    OFFSET_DEG = 0.5

    pid_raster = plate_id_raster.rename({k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                                         if k in plate_id_raster.dims})

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

            # Sign: probe raster on both sides of the segment midpoint.
            # Probing both directions handles triple-junction edges where one probe
            # may land on a third plate, and avoids the centroid hemisphere test which
            # fails for large or concave plate geometries.
            try:
                pos_lon = float(np.clip(mid_lon + OFFSET_DEG * n_east,  -180, 180))
                pos_lat = float(np.clip(mid_lat + OFFSET_DEG * n_north,  -90,  90))
                neg_lon = float(np.clip(mid_lon - OFFSET_DEG * n_east,  -180, 180))
                neg_lat = float(np.clip(mid_lat - OFFSET_DEG * n_north,  -90,  90))
                pos_val = int(pid_raster.sel(lon=pos_lon, lat=pos_lat, method='nearest').item())
                neg_val = int(pid_raster.sel(lon=neg_lon, lat=neg_lat, method='nearest').item())
            except Exception:
                continue

            if pos_val == pid1 or neg_val == pid2:
                sign = 1    # normal points toward pid1
            elif pos_val == pid2 or neg_val == pid1:
                sign = -1   # normal points toward pid2
            else:
                continue    # neither probe hit pid1 or pid2 (e.g. triple junction)

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
                            polarity_dict=None,
                            collision_pairs=None,
                            collision_mode='resistance'):

    if polarity_dict is None:
        polarity_dict = {}

    if collision_pairs is None:
        collision_pairs = set()

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
                        # print(craton_id1, craton_id2)
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
                        if collision_mode == 'coupling':
                            # New behaviour: contribute zero torque; record pair for post-step
                            # velocity blending.  The convergence driving force is gone because
                            # there is no longer a subducting slab to pull the plate forward.
                            if pid_current and pid_neighbor and pid_current > 0 and pid_neighbor > 0:
                                collision_pairs.add((min(pid_current, pid_neighbor),
                                                     max(pid_current, pid_neighbor)))
                        else:
                            # Original behaviour: subtract a resistance torque scaled by the
                            # continental_collision force parameter.
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

    return new_stage_pole_dict, collision_pairs


def apply_collision_coupling(stage_pole_dict, collision_pairs,
                              collision_state, plate_pixel_counts, params):
    """Blend stage poles of actively colliding plate pairs toward each other.

    Simulates the progressive suturing of continents (India–Asia style): the
    longer two cratons have been in contact, the more their velocities are
    coupled.  Coupling is area-weighted so the smaller plate adjusts more.

    collision_state is mutated in-place: ages for active pairs are incremented,
    pairs no longer colliding are pruned.
    """
    coupling_max  = params.get('collision_coupling_max', 0.0)
    timescale     = params.get('collision_coupling_timescale', 30)
    if coupling_max <= 0:
        return

    def rot_to_vec(rot):
        if rot.represents_identity_rotation():
            return np.zeros(3)
        lat, lon, ang = rot.get_lat_lon_euler_pole_and_angle_degrees()
        ang_r = math.radians(ang)
        lat_r = math.radians(lat)
        lon_r = math.radians(lon)
        return np.array([math.cos(lat_r) * math.cos(lon_r),
                         math.cos(lat_r) * math.sin(lon_r),
                         math.sin(lat_r)]) * ang_r

    def vec_to_rot(v):
        ang = float(np.linalg.norm(v))
        if ang < 1e-12:
            return pygplates.FiniteRotation()
        axis = v / ang
        lat = math.degrees(math.asin(float(np.clip(axis[2], -1.0, 1.0))))
        lon = math.degrees(math.atan2(float(axis[1]), float(axis[0])))
        return pygplates.FiniteRotation((lat, lon), math.degrees(ang))

    # Update collision ages: increment active pairs, prune inactive ones
    for pair in list(collision_state.keys()):
        if pair not in collision_pairs:
            del collision_state[pair]
    for pair in collision_pairs:
        collision_state[pair] = collision_state.get(pair, 0) + 1

    for (pid_lo, pid_hi), age in collision_state.items():
        if pid_lo not in stage_pole_dict or pid_hi not in stage_pole_dict:
            continue
        alpha    = min(coupling_max, age / timescale * coupling_max)
        area_lo  = plate_pixel_counts.get(pid_lo, 1)
        area_hi  = plate_pixel_counts.get(pid_hi, 1)
        total    = area_lo + area_hi
        vA       = rot_to_vec(stage_pole_dict[pid_lo])
        vB       = rot_to_vec(stage_pole_dict[pid_hi])
        v_target = (area_lo * vA + area_hi * vB) / total
        stage_pole_dict[pid_lo] = vec_to_rot(vA + alpha * (v_target - vA))
        stage_pole_dict[pid_hi] = vec_to_rot(vB + alpha * (v_target - vB))
        print(f'  [collision] pair ({pid_lo},{pid_hi}) age={age} coupling={alpha:.3f}')


def apply_plate_merges(plate_id_raster, merge_pairs):
    """Replace pixels of smaller plate IDs with larger plate IDs.
    merge_pairs: iterable of (pid_large, pid_small).
    """
    if not merge_pairs:
        return plate_id_raster
    vals = plate_id_raster.values.copy()
    for pid_large, pid_small in merge_pairs:
        mask = vals == pid_small
        if mask.any():
            vals[mask] = pid_large
            print(f'  [merge] Plate {pid_small} absorbed into {pid_large} ({int(mask.sum())} px)')
    return xr.DataArray(vals, dims=plate_id_raster.dims, coords=plate_id_raster.coords)


def detect_and_split_plates(plate_id_raster, new_plate_polygons, params,
                             craton_raster=None, cooldown_pids=None):
    """Detect geometric necks and split implausible plates in the raster.

    Screens plates by concavity ratio and aspect ratio, then uses
    skimage.morphology.skeletonize + distance transform to locate the
    narrowest neck and executes a raster cut perpendicular to the skeleton
    at that point.

    Returns (updated_raster, split_info) where split_info is a list of
    (new_pid, parent_pid) for stage-pole inheritance by the caller.
    """
    try:
        from skimage.morphology import skeletonize
        from skimage.measure import label as sk_label
        from scipy.ndimage import distance_transform_edt
    except ImportError:
        print('  [fragment] skimage/scipy not available — skipping')
        return plate_id_raster, []

    if cooldown_pids is None:
        cooldown_pids = set()

    min_concavity = params.get('min_concavity_ratio', 0.85)
    max_aspect    = params.get('max_aspect_ratio', 4.0)
    rsp           = params.get('raster_spacing', 1.0)
    min_child_px  = max(4, int(params.get('min_child_area_deg2', 25.0) / rsp ** 2))

    vals = plate_id_raster.values.copy()
    lats = plate_id_raster.y.values   # latitude of each raster row
    split_info = []
    next_pid = int(new_plate_polygons['PLATEID1'].max()) + 1

    # Detect dateline-spanning plates: any plate with raster cells in BOTH the
    # first (westernmost) and last (easternmost) columns wraps around ±180° and
    # must never be fragmented — toposimplify can move polygon vertices enough
    # that touches_dateline() misses them, but the raster check is always robust.
    _first_col = vals[:, 0]
    _last_col  = vals[:, -1]
    _dateline_pids = (
        set(int(v) for v in _first_col if not np.isnan(v) and v > 0) &
        set(int(v) for v in _last_col  if not np.isnan(v) and v > 0)
    )

    for _, row in new_plate_polygons.iterrows():
        pid = int(row['PLATEID1'])
        if pid in cooldown_pids:
            continue
        if pid in _dateline_pids:
            continue
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        hull_area = geom.convex_hull.area
        if hull_area < 1e-10:
            continue
        concavity = geom.area / hull_area
        b = geom.bounds
        dx, dy = b[2] - b[0], b[3] - b[1]
        if min(dx, dy) < 1e-10:
            continue

        # FIX 1: Physical aspect ratio — scale longitude extent by cos(lat).
        # Without this, circumpolar plates (large dx in degrees, small dy) always
        # trigger fragmentation even when they are physically compact.
        lat_c   = (b[1] + b[3]) / 2.0
        dx_phys = dx * np.cos(np.radians(lat_c))
        aspect  = max(dx_phys, dy) / max(min(dx_phys, dy), 1e-10)

        if concavity >= min_concavity and aspect <= max_aspect:
            continue

        mask = (vals == pid).astype(np.uint8)
        if int(mask.sum()) < max(min_child_px * 2, 20):
            continue

        # FIX 2: Equal-area skeleton and distance transform.
        # The raster is equirectangular — pixels at high latitudes are physically
        # much narrower in the E-W direction. Skeletonize and EDT in raw pixel
        # space give wrong medial-axis paths and wrong "radius" values near poles.
        # Solution: extract a tight sub-array and compress each row by cos(lat)
        # so all pixels represent equal physical area before skeletonizing.

        pr, pc = np.where(mask)
        r0, r1 = int(pr.min()), int(pr.max())
        c0, c1 = int(pc.min()), int(pc.max())
        sub      = mask[r0:r1+1, c0:c1+1]
        sub_lats = lats[r0:r1+1]
        cos_lats = np.abs(np.cos(np.radians(sub_lats)))
        max_cos  = float(cos_lats.max())
        nr, nc   = sub.shape

        # Guard: if a plate's bounding box spans > 75% of the raster longitude
        # it is either dateline-wrapping or genuinely circumpolar; skip it.
        # (Dateline-spanning plates are already caught by _dateline_pids above,
        # but this handles circumpolar plates whose bbox is wide in pixel space.)
        full_width_frac = (c1 - c0 + 1) / mask.shape[1]
        if full_width_frac > 0.75:
            continue

        # Build equal-area mask: for each row resample columns by cos(lat).
        # Each column in ea-space represents equal physical longitude width.
        ea = np.zeros((nr, nc), dtype=np.uint8)
        for ri in range(nr):
            scale = cos_lats[ri] / max_cos   # 1.0 at most-equatorial row
            if scale >= 0.99:
                ea[ri] = sub[ri]
            else:
                ea_w = max(1, int(round(nc * scale)))
                ea_c = np.arange(ea_w)
                # Inverse sample: which source column feeds each ea column?
                src_c = np.clip(np.round(ea_c / scale).astype(int), 0, nc - 1)
                ea[ri, :ea_w] = sub[ri, src_c]

        def _ea_col_to_orig(ri, c_ea):
            scale = cos_lats[ri] / max_cos
            c_sub = int(round(c_ea / max(scale, 1e-3)))
            return c0 + min(c_sub, nc - 1)

        try:
            skeleton = skeletonize(ea)
        except Exception:
            continue

        skel_pts = np.argwhere(skeleton)
        if len(skel_pts) < 3:
            continue

        # EDT on the equal-area mask: pixel distances now approximate physical
        # distances, so the minimum gives the physically narrowest neck.
        dist       = distance_transform_edt(ea)
        skel_radii = dist[skeleton > 0]
        neck_idx   = int(np.argmin(skel_radii))
        neck_row_ea, neck_col_ea = skel_pts[neck_idx]
        neck_radius = float(skel_radii[neck_idx])

        # Map neck back to full-raster coordinates
        neck_row  = r0 + int(neck_row_ea)
        neck_col  = _ea_col_to_orig(int(neck_row_ea), int(neck_col_ea))

        # Skeleton principal direction from nearby points in ea space
        win = max(3, int(neck_radius * 2))
        nearby = skel_pts[
            (np.abs(skel_pts[:, 0] - neck_row_ea) <= win) &
            (np.abs(skel_pts[:, 1] - neck_col_ea) <= win)
        ]
        if len(nearby) < 2:
            continue
        centered = nearby - nearby.mean(axis=0)
        _, eigvecs = np.linalg.eigh(centered.T @ centered)
        main_dir = eigvecs[:, -1]
        cut_dir  = np.array([-main_dir[1], main_dir[0]])   # perpendicular

        # Draw cutting line in original raster space, mapping each ea point back
        cut_mask = mask.copy()
        cut_len  = int(np.ceil(neck_radius * 2.5)) + 2
        for t in range(-cut_len, cut_len + 1):
            r_ea = int(neck_row_ea) + int(round(t * cut_dir[0]))
            c_ea = int(neck_col_ea) + int(round(t * cut_dir[1]))
            if not (0 <= r_ea < nr):
                continue
            r_orig = r0 + r_ea
            c_orig = _ea_col_to_orig(r_ea, c_ea)
            if 0 <= r_orig < mask.shape[0] and 0 <= c_orig < mask.shape[1]:
                cut_mask[r_orig, c_orig] = 0

        # FIX 3: 4-connectivity (face-adjacent only).
        # With 8-connectivity, two regions that meet only diagonally across a
        # 1-pixel cut still count as connected, so the split never fires.
        labeled = sk_label(cut_mask, connectivity=1)
        labels  = [l for l in np.unique(labeled) if l > 0]
        if len(labels) < 2:
            continue

        areas = {l: int((labeled == l).sum()) for l in labels}
        valid = [l for l, a in areas.items() if a >= min_child_px]
        if len(valid) < 2:
            continue
        valid.sort(key=lambda l: areas[l], reverse=True)

        # Larger child keeps parent PID; smaller child gets new PID
        child_label = valid[1]
        child_mask  = (labeled == child_label) & (mask == 1)
        new_pid = next_pid; next_pid += 1
        vals[child_mask] = new_pid
        split_info.append((new_pid, pid))
        print(f'  [fragment] Plate {pid} → child {new_pid} '
              f'(concavity={concavity:.2f} aspect={aspect:.1f} '
              f'neck_r={neck_radius:.1f}px child_px={areas[child_label]})')

    return xr.DataArray(vals, dims=plate_id_raster.dims, coords=plate_id_raster.coords), split_info


def compute_timestep(rpolygons, forces, stage_pole_dict, delta_time=1, cratons=None, params=None,
                     polarity_dict=None, merge_pairs=None, fragmentation_cooldown_pids=None,
                     collision_state=None):
    if params is None:
        params = {}

    tessellation_deg          = params.get('tessellation_deg',           0.2)
    raster_spacing            = params.get('raster_spacing',             1)
    topo_simplify             = params.get('topo_simplify_tolerance',    1)
    div_thresh_sub            = params.get('div_threshold_subduction',  -1.0)
    div_thresh_ridge          = params.get('div_threshold_ridge',        1.0)
    boundary_search_radius    = params.get('boundary_search_radius_deg', 2)
    collision_mode            = params.get('collision_mode',             'resistance')

    _ct = {}
    _t0 = time.perf_counter()

    new_features = tessellate_polygons(rpolygons, max_dist_deg=tessellation_deg)
    _ct['tessellate'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    rp_gdf = plate_polygons_to_wrapped_gdf(new_features)
    _ct['wrap_polygons'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    clipped_polygons = resolve_overlaps(rp_gdf, priority_col='PLATEID1')
    clipped_polygons = clipped_polygons.explode(index_parts=False).reset_index(drop=True)
    _ct['resolve_overlaps'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    plate_id_raster = update_plate_id_raster(clipped_polygons, raster_spacing=raster_spacing)
    _ct['plate_id_raster'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    # Force values within cratons to overrule any plate polygon changes
    if cratons is not None:
        wrapped_cratons = plate_polygons_to_wrapped_gdf(cratons)
        valid_cratons = wrapped_cratons.dropna(subset=['geometry'])
        valid_cratons = valid_cratons[~valid_cratons.geometry.is_empty]
        craton_raster = rasterize_polygons(valid_cratons, zval_field='PLATEID1', sampling=raster_spacing)
        plate_id_raster = craton_raster.where(craton_raster>0, plate_id_raster)
    else:
        craton_raster = None
    _ct['craton_raster'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    new_plate_polygons = update_plate_polygons(plate_id_raster,
                                               raster_spacing=raster_spacing,
                                               topo_simplify=topo_simplify)
    new_plate_polygons = new_plate_polygons.explode(index_parts=False).reset_index(drop=True)
    new_plate_polygons = new_plate_polygons[new_plate_polygons.geometry.geom_type.isin(['Polygon'])]
    _ct['update_polygons'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    new_plate_polygons = resolve_duplicate_plateids(new_plate_polygons, next_plateid=None,
                                                    cratons_gdf=valid_cratons if cratons is not None else None,
                                                    plate_id_raster=plate_id_raster)
    _ct['resolve_duplicates'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    def _rebuild_polygons(raster):
        polys = update_plate_polygons(raster, raster_spacing=raster_spacing,
                                      topo_simplify=topo_simplify)
        polys = polys.explode(index_parts=False).reset_index(drop=True)
        polys = polys[polys.geometry.geom_type.isin(['Polygon'])]
        return resolve_duplicate_plateids(polys, next_plateid=None,
                                          cratons_gdf=valid_cratons if cratons is not None else None,
                                          plate_id_raster=raster)

    # Apply merges requested by caller (plates whose relative motion dropped to zero)
    if merge_pairs:
        plate_id_raster    = apply_plate_merges(plate_id_raster, merge_pairs)
        new_plate_polygons = _rebuild_polygons(plate_id_raster)
        _ct['merge'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    # Plate fragmentation: split geometrically implausible plates at their neck
    _split_info = []
    if params.get('fragmentation_enabled', False):
        plate_id_raster, _split_info = detect_and_split_plates(
            plate_id_raster, new_plate_polygons, params,
            craton_raster=craton_raster,
            cooldown_pids=fragmentation_cooldown_pids or set())
        if _split_info:
            for new_pid, parent_pid in _split_info:
                if parent_pid in stage_pole_dict:
                    stage_pole_dict[new_pid] = stage_pole_dict[parent_pid]
            new_plate_polygons = _rebuild_polygons(plate_id_raster)
        _ct['fragment'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    boundaries_gdf = topo_polygons_to_boundaries(new_plate_polygons,
                                                  raster_spacing=raster_spacing)
    _ct['topo_boundaries'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    boundaries_gdf = compute_boundary_divergence(boundaries_gdf, stage_pole_dict, delta_time, plate_id_raster,
                                                 plate_polygons_gdf=new_plate_polygons)
    _ct['boundary_divergence'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    velocity_grid_east, velocity_grid_north = velocity_grid(plate_id_raster, stage_pole_dict, delta_time,
                                                             raster_spacing=raster_spacing)
    du_dx = np.gradient(velocity_grid_east.values, velocity_grid_east.x.values, axis=1)
    dv_dy = np.gradient(velocity_grid_north.values, velocity_grid_north.y.values, axis=0)
    div = xr.DataArray(du_dx + dv_dy,
                       dims=velocity_grid_east.dims,
                       coords=velocity_grid_east.coords)
    _ct['velocity_grid'] = time.perf_counter() - _t0; _t0 = time.perf_counter()

    stage_pole_dict, _collision_pairs = polygons_to_stage_poles(
                                              new_plate_polygons, plate_id_raster, forces,
                                              craton_raster,
                                              stage_pole_dict=stage_pole_dict,
                                              delta_time=delta_time,
                                              divergence_threshold_for_subduction=div_thresh_sub,
                                              divergence_threshold_for_ridge_push=div_thresh_ridge,
                                              boundary_search_radius_deg=boundary_search_radius,
                                              polarity_dict=polarity_dict,
                                              collision_mode=collision_mode)

    if collision_mode == 'coupling' and collision_state is not None and _collision_pairs:
        _pvals = plate_id_raster.values.flatten()
        _pvalid = _pvals[~np.isnan(_pvals) & (_pvals > 0)]
        _pu, _pc = np.unique(_pvalid, return_counts=True)
        _plate_pixel_counts = {int(k): int(v) for k, v in zip(_pu, _pc)}
        apply_collision_coupling(stage_pole_dict, _collision_pairs,
                                 collision_state, _plate_pixel_counts, params)

    _ct['stage_poles'] = time.perf_counter() - _t0

    _ct_total = sum(_ct.values())
    print('  [compute_timestep {:.1f}s]'.format(_ct_total),
          '  '.join('{}: {:.1f}s ({:.0f}%)'.format(k, v, 100*v/_ct_total)
                    for k, v in sorted(_ct.items(), key=lambda x: -x[1])))
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
           boundaries_gdf, \
           _split_info, \
           collision_state if collision_state is not None else {}


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

    Vectorised: raster lookups are done as single numpy index arrays; rotations are
    grouped per plate so pygplates is called once per plate, not once per tracer.
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
            if pd.notna(p1) and pd.notna(p2) and p1 != 0 and p2 != 0:
                convergent_pairs.add(frozenset([int(p1), int(p2)]))

    lons = active.geometry.x.values.copy()
    lats = active.geometry.y.values.copy()
    ages = active.age.values.copy().astype(float) + delta_time
    active_flags = active.active.values.copy().astype(bool)

    # Vectorised nearest-neighbour raster lookup using searchsorted
    # Coordinate names may be 'x'/'y' or 'lon'/'lat' depending on caller
    rast_lons = (plate_id_raster.lon if 'lon' in plate_id_raster.coords
                 else plate_id_raster.x).values
    rast_lats = (plate_id_raster.lat if 'lat' in plate_id_raster.coords
                 else plate_id_raster.y).values
    rast_vals = plate_id_raster.values          # (nlat, nlon)
    nlon, nlat = len(rast_lons), len(rast_lats)

    def _raster_pids(lo_arr, la_arr):
        ci = np.clip(np.searchsorted(rast_lons, lo_arr), 0, nlon - 1)
        ri = np.clip(np.searchsorted(rast_lats, la_arr), 0, nlat - 1)
        return rast_vals[ri, ci]

    old_pids = _raster_pids(lons, lats)

    new_lons = lons.copy()
    new_lats = lats.copy()

    # Group by plate — one pygplates call per unique plate ID
    for pid_val in np.unique(old_pids):
        if np.isnan(pid_val):
            continue
        pid_int = int(pid_val)
        rotation = stage_pole_dict.get(pid_int, pygplates.FiniteRotation())
        if rotation.represents_identity_rotation():
            continue
        mask = old_pids == pid_val
        idxs = np.where(mask)[0]
        for k in idxs:
            pt = rotation * pygplates.PointOnSphere(float(lats[k]), float(lons[k]))
            new_lats[k], new_lons[k] = pt.to_lat_lon()

    # Vectorised convergent-crossing check
    if convergent_pairs:
        new_pids = _raster_pids(new_lons, new_lats)
        moved = (new_pids != old_pids) & (old_pids > 0) & (new_pids > 0)
        for k in np.where(moved)[0]:
            if frozenset([int(old_pids[k]), int(new_pids[k])]) in convergent_pairs:
                active_flags[k] = False

    updated = active.copy()
    updated['geometry'] = gpd.points_from_xy(new_lons, new_lats)
    updated['age'] = ages
    updated['active'] = active_flags
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


def update_slab_tracers(slab_tracers, boundaries_gdf, delta_time,
                         div_threshold_subduction=-0.0001,
                         slowing_factor=4,
                         sample_deg=1.0,
                         upper_mantle_km=670.0,
                         cmb_depth_km=2890.0,
                         max_sink_rate_km_per_myr=200.0):
    """Seed new slablets at subduction zones and sink existing ones downward.

    slab_tracers is a plain DataFrame with columns:
        lon, lat, depth_km, age, active, sink_rate_km_per_myr
    All sinking is vertical — lon/lat stay fixed after seeding.
    Sinking rate: upper mantle at full rate, lower mantle at rate/slowing_factor.
    1 cm/yr convergence = 10 km/Myr sinking.
    max_sink_rate_km_per_myr caps the rate; synthetic plates can have unrealistically
    large convergence velocities, so this prevents slablets from racing to the CMB
    in just a few steps (real maximum is ~100-200 km/Myr).
    """
    slab_tracers = slab_tracers.copy()

    if len(slab_tracers) > 0:
        active = slab_tracers['active'].astype(bool)
        upper = active & (slab_tracers['depth_km'] < upper_mantle_km)
        lower = active & (slab_tracers['depth_km'] >= upper_mantle_km)
        slab_tracers.loc[upper, 'depth_km'] += (
            slab_tracers.loc[upper, 'sink_rate_km_per_myr'] * delta_time)
        slab_tracers.loc[lower, 'depth_km'] += (
            slab_tracers.loc[lower, 'sink_rate_km_per_myr'] / slowing_factor * delta_time)
        slab_tracers.loc[active, 'age'] += delta_time
        slab_tracers.loc[slab_tracers['depth_km'] >= cmb_depth_km, 'active'] = False

    if boundaries_gdf is None or 'divergence' not in boundaries_gdf.columns:
        return slab_tracers

    sub = boundaries_gdf[boundaries_gdf['divergence'] < div_threshold_subduction]
    new_rows = []
    for _, brow in sub.iterrows():
        sink_rate = min(abs(float(brow['divergence'])) * 10.0,  # cm/yr → km/Myr
                        max_sink_rate_km_per_myr)
        coords = list(brow.geometry.coords)
        for i in range(len(coords) - 1):
            lon1, lat1 = coords[i]
            lon2, lat2 = coords[i + 1]
            dist = np.hypot(lon2 - lon1, lat2 - lat1)
            n = max(int(dist / sample_deg), 1)
            for lo, la in zip(np.linspace(lon1, lon2, n, endpoint=False),
                              np.linspace(lat1, lat2, n, endpoint=False)):
                new_rows.append({'lon': lo, 'lat': la,
                                 'depth_km': 0.0, 'age': 0.0,
                                 'active': True,
                                 'sink_rate_km_per_myr': sink_rate})
    if new_rows:
        slab_tracers = pd.concat([slab_tracers, pd.DataFrame(new_rows)],
                                  ignore_index=True)
    return slab_tracers


def update_cmb_heat(cmb_heat, slab_tracers, delta_time,
                    spacing_deg=5.0,
                    cmb_depth_km=2890.0,
                    cmb_layer_km=400.0,
                    heating_rate=0.02,
                    slab_suppression=0.05,
                    surface_coupling=0.1,
                    sigma_deg=10.0,
                    max_heat=3.0):
    """Accumulate CMB heat each timestep, building LLSVP-like warm regions.

    Core heating increments the 2D field (nlat × nlon) everywhere each step;
    slablets suppress it with depth-weighted coupling.  A wide horizontal
    Gaussian keeps the pattern large-scale.  Positive values = warm regions
    (slow dVs); passed to generate_mantle_temp_field as the upwelling anomaly.

    Two physical mechanisms are combined:
      - Thermal (material): slablets near the CMB suppress with full weight (1.0).
      - Dynamic pressure (incompressibility): surface slablets suppress with
        surface_coupling weight, representing instantaneous Stokes-flow coupling.
    Coupling varies linearly with depth from surface_coupling at 0 km to 1.0 at
    cmb_depth_km, so the LLSVP anti-correlation with subduction appears from the
    first timestep rather than waiting ~50 Myr for material to arrive physically.

    Parameters
    ----------
    heating_rate : float
        Heat units added per Myr everywhere (core heating).
    slab_suppression : float
        Heat units removed per slablet per Myr (at full coupling = CMB depth).
    surface_coupling : float
        Fractional coupling at depth=0 (0 = no instantaneous signal, 1 = full).
        Default 0.1 represents ~10% of the CMB suppression felt immediately from
        surface subduction via dynamic pressure.
    sigma_deg : float
        Horizontal Gaussian smoothing (degrees).  LLSVPs are large-scale features.
    max_heat : float
        Safety cap on accumulated heat units.
    """
    from scipy.ndimage import gaussian_filter

    nlat, nlon = cmb_heat.shape

    # Core heating everywhere
    cmb_heat = cmb_heat + heating_rate * delta_time

    # Suppress from ALL active slablets with depth-weighted coupling.
    # Coupling linearly interpolates from surface_coupling (dynamic pressure, depth=0)
    # to 1.0 (physical material contact, depth=cmb_depth_km).
    if len(slab_tracers) > 0:
        all_active = slab_tracers[slab_tracers['active'].astype(bool)]
        if len(all_active) > 0:
            lons_f   = all_active['lon'].values.astype(float)
            lats_f   = all_active['lat'].values.astype(float)
            depths_f = all_active['depth_km'].values.astype(float)

            coupling = surface_coupling + (1.0 - surface_coupling) * (depths_f / cmb_depth_km)
            weights  = slab_suppression * delta_time * coupling

            lon_idx = np.clip(np.floor((lons_f + 180) / spacing_deg).astype(int), 0, nlon - 1)
            lat_idx = np.clip(np.floor((lats_f + 90)  / spacing_deg).astype(int), 0, nlat - 1)
            suppress = np.zeros_like(cmb_heat)
            np.add.at(suppress, (lat_idx, lon_idx), weights)
            cmb_heat = np.maximum(cmb_heat - suppress, 0.0)

    # Wide horizontal Gaussian — LLSVPs are large-scale features
    sig = sigma_deg / spacing_deg
    cmb_heat = gaussian_filter(cmb_heat, sigma=[sig, sig], mode=['reflect', 'wrap'])

    return np.minimum(cmb_heat, max_heat)


def generate_mantle_temp_field(slab_tracers,
                                spacing_deg=5.0,
                                depth_interval_km=100.0,
                                max_depth_km=2890.0,
                                temp_anomaly_K=-500.0,
                                thermal_tau_myr=200.0,
                                dVs_per_K=-0.001,
                                sigma_deg=3.0,
                                sigma_depth_km=150.0,
                                cmb_heat=None,
                                cmb_dVs_scale=-0.5,
                                cmb_decay_km=400.0,
                                tomo_smooth_horiz_deg=0.0,
                                tomo_smooth_depth_km=0.0):
    """Bin slab tracers onto a 3D grid and return a dVs/Vs anomaly Dataset.

    Each slablet contributes a cold temperature anomaly (temp_anomaly_K, default -500 K)
    that decays exponentially with age (thermal_tau_myr).  The summed ΔT per grid cell
    is Gaussian-smoothed and converted to dVs/Vs % (dVs_per_K, default -0.001 %/K,
    so a fresh slab gives -0.5% per slablet per cell).

    If cmb_heat (a 2D nlat×nlon array from update_cmb_heat) is provided, a warm
    negative-dVs anomaly is added in the deepest layers, decaying upward with a
    Gaussian envelope of half-width cmb_decay_km.

    Returns xr.Dataset with variables 'dVs' and optionally 'cmb_heat'
    (dims: depth, lat, lon).
    """
    from scipy.ndimage import gaussian_filter

    lons_g  = np.arange(-180 + spacing_deg / 2, 180, spacing_deg)
    lats_g  = np.arange( -90 + spacing_deg / 2,  90, spacing_deg)
    depths_g = np.arange(depth_interval_km / 2, max_depth_km, depth_interval_km)
    nlon, nlat, ndepth = len(lons_g), len(lats_g), len(depths_g)
    field = np.zeros((ndepth, nlat, nlon), dtype=float)

    if len(slab_tracers) > 0:
        active = slab_tracers[slab_tracers['active'].astype(bool)]
        if len(active) > 0:
            ages_f   = active['age'].values.astype(float)
            lons_f   = active['lon'].values.astype(float)
            lats_f   = active['lat'].values.astype(float)
            depths_f = active['depth_km'].values.astype(float)

            dT  = temp_anomaly_K * np.exp(-ages_f / thermal_tau_myr)
            dVs = dT * dVs_per_K

            lon_idx = np.clip(
                np.floor((lons_f + 180) / spacing_deg).astype(int), 0, nlon - 1)
            lat_idx = np.clip(
                np.floor((lats_f + 90)  / spacing_deg).astype(int), 0, nlat - 1)
            dep_idx = np.clip(
                np.floor(depths_f / depth_interval_km).astype(int), 0, ndepth - 1)
            np.add.at(field, (dep_idx, lat_idx, lon_idx), dVs)

            sig_sp  = sigma_deg       / spacing_deg
            sig_dep = sigma_depth_km  / depth_interval_km
            field = gaussian_filter(field, sigma=[sig_dep, sig_sp, sig_sp],
                                    mode=['reflect', 'reflect', 'wrap'])

    # Inject LLSVP-like warm anomaly from CMB heat accumulation.
    # Iterate deepest→shallowest; break once decay becomes negligible.
    if cmb_heat is not None and cmb_heat.shape == (nlat, nlon):
        _n_injected = 0
        for dep_i in range(ndepth - 1, -1, -1):
            dist_from_cmb = max_depth_km - depths_g[dep_i]
            decay = np.exp(-(dist_from_cmb / cmb_decay_km) ** 2)
            if decay < 1e-4:
                break
            field[dep_i] += cmb_heat * cmb_dVs_scale * decay
            _n_injected += 1
        print(f'  [mantle] LLSVP injected into {_n_injected} depth layers; '
              f'cmb_heat range {cmb_heat.min():.3f}–{cmb_heat.max():.3f}; '
              f'dVs at CMB: {field[-1].min():.4f}–{field[-1].max():.4f} %')
    elif cmb_heat is not None:
        print(f'  [mantle] WARNING: cmb_heat shape {cmb_heat.shape} != expected ({nlat},{nlon}) — LLSVP skipped')

    for dep_i in range(ndepth):
        field[dep_i] -= field[dep_i].mean()

    if tomo_smooth_horiz_deg > 0 or tomo_smooth_depth_km > 0:
        sig_h = tomo_smooth_horiz_deg / spacing_deg
        sig_d = tomo_smooth_depth_km  / depth_interval_km
        field = gaussian_filter(field,
                                sigma=[max(sig_d, 1e-9), max(sig_h, 1e-9), max(sig_h, 1e-9)],
                                mode=['reflect', 'reflect', 'wrap'])
        for dep_i in range(ndepth):
            field[dep_i] -= field[dep_i].mean()

    data_vars = {
        'dVs': xr.DataArray(
            field, dims=['depth', 'lat', 'lon'],
            coords={'depth': depths_g, 'lat': lats_g, 'lon': lons_g},
            attrs={'units': '%',
                   'long_name': 'dVs/Vs anomaly (positive = cold fast slab; negative = warm slow region)',
                   'temp_anomaly_K': temp_anomaly_K,
                   'thermal_tau_myr': thermal_tau_myr,
                   'dVs_per_K': dVs_per_K}),
    }
    if cmb_heat is not None and cmb_heat.shape == (nlat, nlon):
        data_vars['cmb_heat'] = xr.DataArray(
            cmb_heat, dims=['lat', 'lon'],
            coords={'lat': lats_g, 'lon': lons_g},
            attrs={'units': 'heat units',
                   'long_name': 'CMB heat accumulation (proxy for LLSVP intensity)'})
    return xr.Dataset(data_vars)


def generate_hotspot_locations(n=20, locations=None, seed=None):
    """Return a GeoDataFrame of fixed hot spot source points.

    If locations is a list of (lon, lat) tuples those are used directly.
    Otherwise n points are drawn uniformly on the sphere via the normal-vector
    method (generate (x,y,z) ~ N(0,1), normalise, convert to lon/lat).
    """
    if locations is not None:
        rows = [{'hotspot_id': i, 'geometry': Point(lon, lat)}
                for i, (lon, lat) in enumerate(locations)]
    else:
        rng = np.random.default_rng(seed)
        xyz = rng.standard_normal((n, 3))
        xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)
        lons = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
        lats = np.degrees(np.arcsin(xyz[:, 2]))
        rows = [{'hotspot_id': i, 'geometry': Point(float(lons[i]), float(lats[i]))}
                for i in range(n)]
    return gpd.GeoDataFrame(rows, geometry='geometry').set_crs(epsg=4326)


def seed_hotspot_seamounts(hotspot_locs, prob=0.3, rng=None):
    """Probabilistically create new seamounts at fixed hot spot locations.

    For each hot spot, a seamount is created with probability prob.
    Returns a GeoDataFrame (possibly empty) with columns matching the tracers
    GeoDataFrame so that advect_tracers and deactivate_subducted_tracers work
    without modification.
    """
    if rng is None:
        rng = np.random.default_rng()
    rows = []
    for _, hs in hotspot_locs.iterrows():
        if rng.uniform() < prob:
            rows.append({
                'geometry':    hs.geometry,
                'hotspot_id':  int(hs['hotspot_id']),
                'age':         0.0,
                'active':      True,
            })
    if not rows:
        return gpd.GeoDataFrame(columns=['geometry', 'hotspot_id', 'age', 'active'],
                                geometry='geometry').set_crs(epsg=4326)
    return gpd.GeoDataFrame(rows, geometry='geometry').set_crs(epsg=4326)


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
                'plateid1':       int(v) if pd.notna(v := seg_row.get('plateid1_poly1', 0)) else 0,
                'plateid2':       int(v) if pd.notna(v := seg_row.get('plateid1_poly2', 0)) else 0,
                'downgoing_plate':  0,
                'overriding_plate': 0,
            })

        elif div < div_threshold_subduction:
            # Subduction: overriding-plate side only, depth from slab dip
            _r1, _r2 = seg_row.get('plateid1_poly1', 0), seg_row.get('plateid1_poly2', 0)
            p1 = int(_r1) if pd.notna(_r1) else 0
            p2 = int(_r2) if pd.notna(_r2) else 0
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


def compute_topography_raster(plate_id_raster, craton_raster, boundaries_gdf,
                               polarity_dict=None, craton_counts=None,
                               continental_base_m=500,
                               noise_fraction=0.15,
                               spectral_beta=2.0,
                               andean_max_dist_deg=5.0,
                               andean_arc_dist_deg=1.5,
                               andean_peak_m=5000,
                               himalayan_max_dist_deg=4.0,
                               himalayan_peak_m=8000,
                               div_threshold_subduction=-0.0001,
                               coastal_margin_deg=2.0):
    """Generate a synthetic elevation raster from tectonic setting.

    Continental pixels (craton_raster > 0) receive:
      - Background: continental_base_m + bandlimited Gaussian noise
      - Andean zone: Gaussian elevation peak at andean_arc_dist_deg from the trench,
        only on the overriding (continental) plate side of ocean-continent subduction
      - Himalayan zone: linear ramp peaked at the collision boundary for
        continent-continent convergent boundaries (both plates have cratons)
    Ocean pixels remain at 0.  Returns an xr.DataArray with lon/lat dims.
    """
    from sklearn.neighbors import BallTree
    import pyshtools as pysh
    from scipy.interpolate import RegularGridInterpolator

    if polarity_dict is None:
        polarity_dict = {}
    if craton_counts is None:
        craton_counts = {}

    def _ensure_lonlat(da):
        return da.rename({k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                          if k in da.dims})

    pid_r  = _ensure_lonlat(plate_id_raster)
    crat_r = _ensure_lonlat(craton_raster)

    lons = pid_r.lon.values
    lats = pid_r.lat.values
    raster_spacing = float(abs(lons[1] - lons[0]))

    lon2d, lat2d = np.meshgrid(lons, lats)          # (nlat, nlon)
    nlat, nlon   = lon2d.shape

    craton_vals = np.nan_to_num(crat_r.values, nan=0.0).astype(float)
    is_craton   = craton_vals > 0
    elev        = np.zeros((nlat, nlon), dtype=float)

    # --- Background continental elevation (flat; noise added after all zones) ---
    elev[is_craton] = float(continental_base_m)

    # --- Classify convergent segments: Andean vs Himalayan ---
    andean_pts   = []   # (lat, lon)
    andean_ovpid = []   # overriding plate ID for each Andean BallTree point
    himalayan_pts = []

    TESS_DEG = 0.5
    if boundaries_gdf is not None and 'divergence' in boundaries_gdf.columns:
        for iloc_i in range(len(boundaries_gdf)):
            row = boundaries_gdf.iloc[iloc_i]
            if row.get('divergence', 0.0) >= div_threshold_subduction:
                continue
            _r1, _r2 = row.get('plateid1_poly1', 0), row.get('plateid1_poly2', 0)
            p1 = int(_r1) if pd.notna(_r1) else 0
            p2 = int(_r2) if pd.notna(_r2) else 0
            if p1 == 0 or p2 == 0:
                continue

            dg = resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts)
            ov = p2 if dg == p1 else p1
            a_ov = craton_counts.get(ov, 0)
            a_dg = craton_counts.get(dg, 0)
            is_andean    = a_ov > 0
            is_himalayan = a_ov > 0 and a_dg > 0

            if not (is_andean or is_himalayan):
                continue

            coords = list(row.geometry.coords)
            if len(coords) < 2:
                continue
            try:
                pp1 = pygplates.PointOnSphere(coords[0][1], coords[0][0])
                pp2 = pygplates.PointOnSphere(coords[1][1], coords[1][0])
                seg_pts = (pygplates.PolylineOnSphere([pp1, pp2])
                           .to_tessellated(np.radians(TESS_DEG)).get_points())
                for pt in seg_pts:
                    lat, lon = pt.to_lat_lon()
                    if is_andean:
                        andean_pts.append([lat, lon])
                        andean_ovpid.append(ov)
                    if is_himalayan:
                        himalayan_pts.append([lat, lon])
            except Exception:
                pass

    # --- Andean elevation ---
    if andean_pts:
        tree = BallTree(np.deg2rad(andean_pts), metric='haversine')
        ri, ci = np.where(is_craton)
        pts_rad = np.deg2rad(np.column_stack([lat2d[ri, ci], lon2d[ri, ci]]))
        dists_deg, nn_idx = tree.query(pts_rad, k=1)
        dists_deg = np.degrees(dists_deg.flatten())
        nn_idx    = nn_idx.flatten()

        for k in np.where(dists_deg < andean_max_dist_deg)[0]:
            if int(craton_vals[ri[k], ci[k]]) != andean_ovpid[nn_idx[k]]:
                continue    # pixel is in the downgoing plate — skip
            d = dists_deg[k]
            arc_elev = andean_peak_m * np.exp(
                -((d - andean_arc_dist_deg) / andean_arc_dist_deg) ** 2)
            elev[ri[k], ci[k]] = max(elev[ri[k], ci[k]], arc_elev)

    # --- Himalayan elevation ---
    if himalayan_pts:
        tree = BallTree(np.deg2rad(himalayan_pts), metric='haversine')
        pts_rad = np.deg2rad(np.column_stack([lat2d.flatten(), lon2d.flatten()]))
        dists_deg, _ = tree.query(pts_rad, k=1)
        dists_deg = np.degrees(dists_deg.flatten()).reshape(nlat, nlon)

        him_zone = (dists_deg < himalayan_max_dist_deg) & is_craton
        him_elev = himalayan_peak_m * np.maximum(
            0.0, 1.0 - dists_deg / himalayan_max_dist_deg)
        elev[him_zone] = np.maximum(elev[him_zone], him_elev[him_zone])

    # --- Elevation-scaled spherical noise via spherical harmonic synthesis ---
    # lmax chosen so DH grid spacing ≈ raster_spacing: spacing = 90/(lmax+1)
    lmax = max(int(90.0 / raster_spacing) - 1, 10)
    l = np.arange(lmax + 1, dtype=float)
    power = np.where(l > 0, (l + 1.0) ** (-spectral_beta), 0.0)

    clm = pysh.SHCoeffs.from_random(power)
    sh_grid = clm.expand(grid='DH')
    noise_arr = sh_grid.to_array()
    noise_arr /= noise_arr.std()

    dh_lats = sh_grid.lats()   # descending degrees
    dh_lons = sh_grid.lons()   # 0 → ~360

    interp = RegularGridInterpolator(
        (dh_lats[::-1], dh_lons), noise_arr[::-1, :],
        method='linear', bounds_error=False, fill_value=0.0)

    lon2d_360 = lon2d % 360.0
    pts = np.column_stack([lat2d.flatten(), lon2d_360.flatten()])
    noise = interp(pts).reshape(nlat, nlon)
    std = noise.std()
    if std > 0:
        noise /= std

    elev[is_craton] += noise[is_craton] * elev[is_craton] * noise_fraction
    elev[is_craton]  = np.maximum(elev[is_craton], 100.0)

    # --- Coastal taper: smooth elevation to zero within coastal_margin_deg of the
    #     craton edge, eliminating the cliff where continent meets ocean.
    #     distance_transform_edt gives the Euclidean distance (in pixels) to the
    #     nearest non-craton pixel; multiply by raster_spacing to get degrees.
    from scipy.ndimage import distance_transform_edt
    dist_deg = distance_transform_edt(is_craton) * raster_spacing
    t = np.clip(dist_deg / coastal_margin_deg, 0.0, 1.0)
    taper = t * t * (3.0 - 2.0 * t)   # smoothstep: 0 at coast, 1 inland
    elev *= taper                       # ocean pixels are already 0; no-op there

    return xr.DataArray(elev, dims=['lat', 'lon'],
                        coords={'lat': lats, 'lon': lons},
                        attrs={'units': 'meters', 'long_name': 'synthetic elevation'})


def generate_seafloor_age_grid(tracer_csv, craton_raster=None, region='d', spacing='1d'):
    """Interpolate seafloor age tracers to a global grid and optionally mask continents.

    The tracer CSV is expected to have columns [age, lon, lat] (the format
    written by run_SDP.py).  blockmedian reduces duplicate points per cell
    before sphinterpolate produces the final grid.

    craton_raster : xr.DataArray with plate-ID values > 0 on continental pixels.
        Where it is positive the returned age grid is set to NaN.
    """
    # pygmt spacing must be a string like '0.5d'; convert float/int if needed
    if not isinstance(spacing, str):
        spacing = f'{spacing}d'

    res = pygmt.blockmedian(tracer_csv, region=region, incols='1,2,0',
                            spacing=spacing, output_type='numpy')
    age_grid = pygmt.sphinterpolate(res, region=region, spacing=spacing)

    if craton_raster is not None:
        # Normalise dim names so interp works regardless of x/y vs lon/lat
        rename = {k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                  if k in craton_raster.dims}
        cr = craton_raster.rename(rename) if rename else craton_raster
        cr_on_grid = cr.interp(lat=age_grid.lat, lon=age_grid.lon, method='nearest')
        age_grid = age_grid.where(cr_on_grid <= 0)

    age_grid.attrs.update({'units': 'Myr', 'long_name': 'seafloor age'})
    return age_grid


def generate_topo_bathy_grid(seafloor_age_grid, topography_raster=None,
                              default_ocean_depth_m=-3500):
    """Combine GDH1 bathymetry with synthetic topography into a seamless grid.

    GDH1 model (Stein & Stein 1992):
      depth = -(2600 + 365 * sqrt(age))          for age < 20 Myr
      depth = -(5651 - 2473 * exp(-0.0278*age))  for age >= 20 Myr

    topography_raster : xr.DataArray, optional.  Continental elevations (> 0 m).
        Where present these overwrite the bathymetry.  If None, only bathymetry
        is returned.
    default_ocean_depth_m : fill value for ocean pixels with no age data.
        Set to None to leave those pixels as NaN.
    """
    age = seafloor_age_grid.values.copy().astype(float)

    # GDH1 — work in-place on a pre-NaN'd array
    bathy = np.full_like(age, np.nan)
    young = np.isfinite(age) & (age < 20.)
    old   = np.isfinite(age) & (age >= 20.)
    bathy[young] = -(2600. + 365. * np.sqrt(age[young]))
    bathy[old]   = -(5651. - 2473. * np.exp(-0.0278 * age[old]))

    combined = bathy.copy()

    # Fill uncovered ocean pixels with a default depth
    if default_ocean_depth_m is not None:
        combined[~np.isfinite(combined)] = float(default_ocean_depth_m)

    # Overwrite continental pixels with topography
    if topography_raster is not None:
        rename = {k: v for k, v in {'x': 'lon', 'y': 'lat'}.items()
                  if k in topography_raster.dims}
        topo = topography_raster.rename(rename) if rename else topography_raster
        topo_on_grid = topo.interp(lat=seafloor_age_grid.lat,
                                   lon=seafloor_age_grid.lon, method='nearest')
        topo_vals = topo_on_grid.values
        continental = topo_vals > 0
        combined[continental] = topo_vals[continental]

    return xr.DataArray(combined,
                        dims=seafloor_age_grid.dims,
                        coords=seafloor_age_grid.coords,
                        attrs={'units': 'meters', 'long_name': 'topography-bathymetry'})


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


def remove_boundary_segments(arc, lon_threshold=None, lat_threshold=None,
                             raster_spacing=0.25):
    """Remove arc segments that run along the dateline or poles.

    Only called for arcs shared between two different plates (outer-boundary
    arcs are excluded upstream).  Uses a narrow threshold (raster_spacing/4)
    so only vertices at exactly ±180°/±90° (from clip) are filtered, while
    real plate-boundary vertices at ±179.875° (0.25° grid cell edges) are kept.
    """
    if lon_threshold is None:
        lon_threshold = 180.0 - raster_spacing / 4   # e.g. 179.9375 for 0.25° grid — catches
    if lat_threshold is None:                        # exactly ±180° but not 179.875° cell edges
        lat_threshold = 90.0 - raster_spacing / 4
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


def topo_polygons_to_boundaries(plate_polygons, raster_spacing=0.25):

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
        poly_list = arc_to_polys.get(arc_idx, [])
        unique_polys = list(set(poly_list))

        # Only process arcs shared between exactly two polygons with different PLATEIDs.
        # Outer-boundary arcs (world edge at ±180°/±90°) have len==1 and are never real
        # plate boundaries; same-plate arcs (duplicate ID fragments) are internal artifacts.
        if len(unique_polys) < 2:
            continue
        pid1 = feat_to_plateid.get(unique_polys[0])
        pid2 = feat_to_plateid.get(unique_polys[1])
        if pid1 == pid2:
            continue

        sub_arcs = remove_boundary_segments(arc, raster_spacing=raster_spacing)

        for sub_arc in sub_arcs:
            boundaries.append({
                'geometry': LineString(sub_arc),
                'poly1': unique_polys[0],
                'poly2': unique_polys[1],
                'plateid1_poly1': pid1,
                'plateid1_poly2': pid2
            })

    boundaries_gdf = gpd.GeoDataFrame(boundaries, geometry='geometry').set_crs(epsg=4326)

    return boundaries_gdf





def touches_dateline(geom, tol=1e-6):
    """Check if polygon has vertices at lon ≈ ±180"""
    coords = np.array(geom.exterior.coords)
    return np.any(np.abs(np.abs(coords[:, 0]) - 180) < tol)

def resolve_duplicate_plateids(gdf, next_plateid=None, cratons_gdf=None,
                                plate_id_raster=None):
    if next_plateid is None:
        next_plateid = gdf['PLATEID1'].max() + 1

    # Pre-compute which plate IDs span the dateline using the raster.
    # A plate that has cells in both the first and last raster column wraps around
    # ±180° and must NOT be split into fragments with different IDs — toposimplify
    # can move polygon vertices far enough from ±180° that touches_dateline() misses them.
    dateline_spanning_ids = set()
    if plate_id_raster is not None:
        first_col = plate_id_raster.values[:, 0]
        last_col  = plate_id_raster.values[:, -1]
        ids_west = set(int(v) for v in first_col if not np.isnan(v) and v > 0)
        ids_east = set(int(v) for v in last_col  if not np.isnan(v) and v > 0)
        dateline_spanning_ids = ids_west & ids_east

    gdf = gdf.copy()

    # Pre-compute a union geometry per craton plate ID so we can check intersection cheaply
    craton_union_by_id = {}
    if cratons_gdf is not None and len(cratons_gdf):
        for cid, cgroup in cratons_gdf.groupby('PLATEID1'):
            craton_union_by_id[int(cid)] = cgroup.union_all()

    for plateid, group in gdf.groupby('PLATEID1'):
        if len(group) <= 1:
            continue

        # All fragments of a dateline-spanning plate are the same physical plate —
        # keep them all with the original ID regardless of vertex positions.
        if int(plateid) in dateline_spanning_ids:
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
        candidates = list(group.index[~gdf.loc[group.index, '_on_dateline']])

        # If this plate ID belongs to a craton, keep the fragment that contains the
        # craton geometry — prevents the craton-override raster from creating an orphan
        # island the following timestep when the larger plate retains the original ID
        if int(plateid) in craton_union_by_id:
            craton_geom = craton_union_by_id[int(plateid)]
            preferred = next(
                (idx for idx in candidates
                 if gdf.loc[idx, 'geometry'].intersects(craton_geom)),
                None)
            if preferred is not None:
                candidates = [preferred] + [c for c in candidates if c != preferred]

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




