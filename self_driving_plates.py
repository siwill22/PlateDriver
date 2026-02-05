import pygplates
from gprm.utils.create_gpml import gpml2gdf, gdf2gpml
from gprm.utils.proximity import rasterize_polygons
from gprm.utils.geometry import wrap_polygon_feature, wrap_polygon_features
import plate_driver as pdriver
import numpy as np
from shapely.ops import unary_union
from shapely.geometry import box
import pygmt
import metpy.calc as mpcalc

import topojson as tp
from xrspatial.experimental.polygonize import polygonize



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


def nudge_polygons(polygon_gdf, stage_pole_dict):
    
    polygon_features = gdf2gpml(polygon_gdf.set_crs(epsg=4326))
    
    rpolygons = []
    for polygon_feature in polygon_features:
        plate_id = polygon_feature.get_reconstruction_plate_id()
        nudged_geometry = stage_pole_dict[plate_id] * polygon_feature.get_geometry()
        nfeature = polygon_feature.clone()
        nfeature.set_geometry(nudged_geometry)
        rpolygons.append(nfeature)
        
    return rpolygons


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


def update_plate_id_raster(clipped_polygons):
    plate_id_raster = rasterize_polygons(clipped_polygons, zval_field='PLATEID1')
    plate_id_raster = plate_id_raster.where(plate_id_raster>0, np.nan)

    plate_id_raster = fill_nan_spherical(plate_id_raster)
    
    return plate_id_raster


def update_plate_polygons(plate_id_raster):
    # Build affine from your raster metadata
    transform = (1, 0, -180.5,
                 0, 1, -90.5)
    new_plate_polygons = polygonize(plate_id_raster,#[1:,1:], 
                                    return_type='geopandas', transform=transform)

    topo = tp.Topology(new_plate_polygons, prequantize=False)
    new_plate_polygons_smooth = topo.toposimplify(1).to_gdf()  # adjust tolerance

    # Clip the GeoDataFrame
    new_plate_polygons_smooth = new_plate_polygons_smooth.clip(box(-180, -90, 180, 90))

    new_plate_polygons_smooth = new_plate_polygons_smooth.rename(columns={'DN':'PLATEID1'})
    
    return new_plate_polygons_smooth.set_crs(epsg=4326)


def velocity_grid(plate_id_raster, plate_parameters, delta_time=1):
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
        if int(velocity_point.pid)==0:
            equivalent_stage_rotation = pygplates.FiniteRotation()
        else:
            equivalent_stage_rotation = plate_parameters[int(velocity_point.pid)] #rotation_model.get_rotation(reconstruction_time, int(velocity_point.pid), reconstruction_time + delta_time)

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
    
    velocity_grid_east = pygmt.xyz2grd(x=velpts.x, y=velpts.y, z=velpts.vel_east, region='d', spacing=1)
    velocity_grid_north = pygmt.xyz2grd(x=velpts.x, y=velpts.y, z=velpts.vel_north, region='d', spacing=1) 
    
    return velocity_grid_east, velocity_grid_north


def polygons_to_stage_poles(polygons_gdf, divergence_raster, plate_id_raster, forces,
                            divergence_threshold_for_subduction = -0.00025,
                            divergence_threshold_for_ridge_push = 0.00025):
    
    stage_pole_dict = {}
    torque_dict = {}
    
    gpml_plates = gdf2gpml(polygons_gdf.explode(index_parts=False).reset_index(drop=True))

    for i,plate in enumerate(gpml_plates):
        
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
                                               radius_deg=2)
                #print(id1, id2)
                # insert logic to determine overriding plate
                downgoing_plate = id1
                overriding_plate = id2

                # This tells us if boundary is divergent or 
                divergence_value = divergence_raster.interp(lon=arc_point[1], lat=arc_point[0]).item()

                segment_torque = None
                if divergence_value<divergence_threshold_for_subduction:

                    if downgoing_plate:

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
            computed_stage_rotation = None
        
        torque_dict[plate.get_reconstruction_plate_id()] = net_torque
        stage_pole_dict[plate.get_reconstruction_plate_id()] = computed_stage_rotation
        
    return stage_pole_dict


def compute_timestep(rpolygons, forces, stage_pole_dict, delta_time=1):
    
    rp_gdf = plate_polygons_to_wrapped_gdf(rpolygons)

    clipped_polygons = resolve_overlaps(rp_gdf, priority_col='PLATEID1')

    plate_id_raster = update_plate_id_raster(clipped_polygons)

    new_plate_polygons = update_plate_polygons(plate_id_raster)
    new_plate_polygons = new_plate_polygons.explode(index_parts=False).reset_index(drop=True)
    new_plate_polygons = new_plate_polygons[new_plate_polygons.geometry.geom_type.isin(['Polygon'])]
    #print(new_plate_polygons)

    velocity_grid_east, velocity_grid_north = velocity_grid(plate_id_raster, stage_pole_dict, delta_time)
    #
    div = mpcalc.divergence(velocity_grid_east, velocity_grid_north)
    #
    stage_pole_dict = polygons_to_stage_poles(new_plate_polygons, div, plate_id_raster, forces)
    #
    #print([(key,item.get_lat_lon_euler_pole_and_angle_degrees()) for key,item in stage_pole_dict.items()])

    
    return new_plate_polygons, stage_pole_dict








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
            mid = segment.interpolate(0.5, normalized=True)
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




