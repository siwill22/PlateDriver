# PlateDriver — Self-Driving Plate Tectonics Simulation

## Inspiration

The core approach in PlateDriver — computing forces on plate boundaries and integrating them forward to derive stage rotations — is heavily inspired by the [Tectonic Explorer](https://github.com/concord-consortium/tectonic-explorer) project by the Concord Consortium. That project demonstrates that a small set of physically motivated boundary forces (slab pull, ridge push, collision resistance) is sufficient to produce broadly realistic plate motion on an interactive timescale. PlateDriver extends the same idea into a full simulation pipeline built on [GPlates](https://www.gplates.org) / pygplates, with raster-based plate tracking, tracer advection, and output compatible with GPlates workflows.

## What this is

PlateDriver is a forward model of plate tectonic evolution. Starting from an initial set of tectonic plates (defined in GPlates format), the model computes the forces acting on each plate at every timestep, derives stage rotations, and evolves the plate geometry forward in time. The plates "drive themselves" — no prescribed rotation model is needed after the initial snapshot.

The model tracks:
- **Plate polygons** — boundaries evolve through each timestep
- **Seafloor tracers** — massless particles seeded at spreading ridges and advected with the plates, deactivated when subducted
- **Seafloor age** — estimated from tracer distance to nearest ridge, incremented each step
- **Finite rotations** — accumulated from stage poles, written to a GPlates-compatible `.rot` file
- **Synthetic earthquakes** — generated near active boundaries with physically motivated depth geometry
- **Synthetic topography** — continental elevations with Andean/Himalayan mountain belts and elevation-scaled terrain noise
- **Seafloor bathymetry** — derived from seafloor age via the GDH1 age-depth relation

## Inputs

All input files are GPlates-format and currently point to files in `/Users/simon/GIT/octopus/Atlantis3/`:

| File | Purpose |
|---|---|
| `atlantis3_rotations.rot` | Rotation model for the initial plate snapshot |
| `atlantis6_topologies.gpml` | Topological plate polygons (dynamic) |
| `atlantis6_topological_polygons.gpml` | Additional topological features |
| `atlantis6_cratons.geojson` | Craton (continental block) polygons, used for collision detection and upper-plate determination |

These are a synthetic ("Atlantis") plate model, not Earth. Paths are hardcoded in `run_SDP.py` and should be updated for any new model.

## Outputs

All outputs go to `./output/` (must exist before running):

| File pattern | Interval | Contents |
|---|---|---|
| `test_plot_NNNN.png` | Every step | Map: plate polygons, boundary types (red=ridge, blue=subduction, black=conservative), velocity vectors, cratons, tracers coloured by age |
| `tracers_NNNN.csv` | Every `csv_save_interval` steps | Active tracer lon/lat/age |
| `divergence_NNNN.nc` | Every `csv_save_interval` steps | Divergence field (xarray, lon/lat) |
| `plate_id_NNNN.nc` | Every `csv_save_interval` steps | Plate ID raster (xarray, lon/lat) |
| `craton_NNNN.nc` | Every `csv_save_interval` steps | Craton ID raster (xarray, lon/lat) |
| `boundaries_NNNN.geojson` | Every `csv_save_interval` steps | Boundary segments with plate IDs, divergence, boundary type, and downgoing/overriding plate for subduction zones |
| `earthquakes_NNNN.geojson` | Every `csv_save_interval` steps | Synthetic earthquake locations with depth, boundary type, and plate IDs |
| `topography_NNNN.nc` | Every `csv_save_interval` steps | Synthetic continental elevation (m); ocean pixels = 0 |
| `seafloor_age_NNNN.nc` | Every `csv_save_interval` steps | Spherically interpolated seafloor age grid (Myr); continental pixels = NaN |
| `topo_bathy_NNNN.nc` | Every `csv_save_interval` steps | Combined topography-bathymetry (m); GDH1 ocean depths merged with continental topography |
| `synthetic_rotations.rot` | End of run | GPlates rotation file with finite rotations for all plates at all timesteps |

## How to run

```bash
conda activate pygmt17
python run_SDP.py
```

All configuration is in the top section of `run_SDP.py`. No command-line arguments.

## Key parameters

### Timing
```python
reconstruction_time = 0     # starting snapshot time (Ma)
delta_time          = 1.    # timestep size (Myr)
total_time          = 100.  # run duration (Myr); loop runs steps 1..int(total_time/delta_time)
```

### Forces (`forces` dict)
Magnitudes are dimensionless scale factors applied to the torque from each boundary segment type.
```python
forces = {
    'slab_pull':            8,   # attractive force on downgoing plate at convergent boundaries
    'slab_suction':         8,   # attractive force on overriding plate at convergent boundaries
    'ridge_push':           3,   # repulsive force at divergent boundaries
    'continental_collision': 8,  # resistance (subtracted) at continent-continent boundaries
}
```

### Numerical/algorithmic (`params` dict)
```python
params = {
    # Core rasterisation
    'raster_spacing':              0.5,     # degrees; controls velocity grid and rasterisation resolution
    'tessellation_deg':            0.2,     # polygon edge densification before rasterising
    'topo_simplify_tolerance':     0.5,     # TopoJSON simplification tolerance after polygonise
    'boundary_search_radius_deg':  2,       # plate ID neighbour search radius at segment midpoints

    # Boundary type thresholds
    'div_threshold_subduction':   -0.0001,  # relative normal velocity (cm/yr) below this → convergent
    'div_threshold_ridge':         0.2,     # relative normal velocity above this → divergent

    # Seafloor tracers
    'tracer_spacing_deg':          1.0,     # along-ridge spacing for newly seeded tracer pairs
    'tracer_offset_deg':           0.6,     # perpendicular offset from ridge axis when seeding
    'seafloor_spacing_deg':        2.0,     # grid spacing for initial seafloor age tracers
    'spreading_rate_km_per_myr':   50.0,    # assumed spreading rate for initial age estimation

    # Output interval
    'csv_save_interval':           5,       # save NetCDF/CSV/GeoJSON outputs every N timesteps

    # Synthetic earthquakes
    'eq_n_candidates':             100000,  # random sphere points sampled per step
    'eq_ridge_max_dist_deg':       1.0,     # max distance from ridge axis for ridge earthquakes
    'eq_ridge_max_depth_km':       50.0,    # max depth for ridge earthquakes (km)
    'eq_slab_dip_deg':             45.0,    # assumed slab dip angle (degrees)
    'eq_slab_max_depth_km':        650.0,   # maximum slab depth (km)
    'eq_ridge_prob':               0.15,    # fraction of ridge candidates accepted (density control)
    'eq_subduction_prob':          1.0,     # fraction of subduction candidates accepted

    # Synthetic topography
    'topo_continental_base_m':     400,     # background continental elevation (m)
    'topo_noise_fraction':         0.15,    # noise amplitude as fraction of local elevation
    'topo_spectral_beta':          2.0,     # SH power-law exponent (2 = brown / realistic terrain)
    'topo_andean_max_dist_deg':    11.0,    # cutoff distance from trench for Andean mountains
    'topo_andean_arc_dist_deg':    4.5,     # distance of arc-peak from trench (degrees)
    'topo_andean_peak_m':          3000,    # peak Andean elevation (m)
    'topo_himalayan_max_dist_deg': 8.0,     # half-width of Himalayan collision zone (degrees)
    'topo_himalayan_peak_m':       5000,    # peak Himalayan elevation (m)
}
```

---

## Codebase reference (for further development)

### File structure

| File | Role |
|---|---|
| `self_driving_plates.py` | Library — all simulation logic as importable functions |
| `run_SDP.py` | Driver script — configuration, bootstrapping, main loop, plotting, output |
| `plate_driver.py` | External dependency — torque computation from boundary segments |

### Per-timestep data flow

```
nudge_polygons(new_plate_polygons, stage_pole_dict)
    → rpolygons (list of pygplates features, geometries advected by previous stage poles)

compute_timestep(rpolygons, forces, stage_pole_dict, delta_time, cratons, params, polarity_dict)
    → new_plate_polygons (GeoDataFrame, PLATEID1 col)
    → stage_pole_dict    (Dict[int, FiniteRotation])
    → rasters            (xr.Dataset: plate_id, velocity_east, velocity_north, divergence, craton)
    → boundaries_gdf     (GeoDataFrame: geometry(LineString), plateid1_poly1, plateid1_poly2,
                          divergence — one row per arc segment)

[accumulate finite rotations]
    cumulative_rotation_dict[pid] = stage_rot * cumulative_rotation_dict.get(pid, identity)
    rotation_records.append((t, pid, cumulative_rotation_dict[pid]))

seed_ridge_tracers(boundaries_gdf, ...)
    → new_tracers (GeoDataFrame: geometry, active=True, age=0)

pd.concat([tracers, new_tracers])

advect_tracers(tracers, rasters['plate_id'], stage_pole_dict, delta_time, boundaries_gdf)
    → tracers (active points moved by plate stage rotation; age += delta_time;
               deactivated if plate-crossing detected at a convergent boundary)

deactivate_subducted_tracers(tracers, boundaries_gdf, threshold)
    → tracers (active=False for remaining tracers near convergent segments)

[at csv_save_interval steps]
    generate_seafloor_age_grid(tracer_csv, craton_raster, spacing)   → seafloor_age_NNNN.nc
    craton_counts_from_raster(rasters['craton'])                      → craton_counts dict
    [boundary GeoJSON with resolve_downgoing_plate]                   → boundaries_NNNN.geojson
    generate_synthetic_earthquakes(boundaries_gdf, ...)               → earthquakes_NNNN.geojson
    compute_topography_raster(plate_id_raster, craton_raster, ...)    → topography_NNNN.nc
    generate_topo_bathy_grid(age_grid, topo)                          → topo_bathy_NNNN.nc

write_rotation_file(rotation_records, path)   # called once after loop
```

### `compute_timestep` internal pipeline

```
tessellate_polygons           densify polygon edges to max_dist_deg
plate_polygons_to_wrapped_gdf wrap at dateline
resolve_overlaps              clip overlapping polygons by PLATEID1 priority
update_plate_id_raster        rasterize to plate_id grid; fill NaN via BallTree haversine
[craton raster]               rasterize cratons, overwrite plate_id where craton > 0
update_plate_polygons         polygonize raster → toposimplify → clip to globe
resolve_duplicate_plateids    rename duplicate IDs from polygon splits
topo_polygons_to_boundaries   extract shared boundary segments via TopoJSON
compute_boundary_divergence   per-arc-segment relative normal velocity (not raster divergence);
                              sign determined by hemisphere test against plate centroids
velocity_grid                 compute east/north velocity at each raster cell via pygplates
mpcalc.divergence             divergence of velocity field (raster only; for visualisation)
polygons_to_stage_poles       integrate boundary torques → stage rotation per plate
```

### Key data structures

**`stage_pole_dict`**: `Dict[int, pygplates.FiniteRotation]`  
Maps plate ID → stage rotation for one `delta_time` step. Used for advection and velocity computation.

**`polarity_dict`**: `Dict[Tuple[int,int], int]`  
Maps `(min_pid, max_pid)` → downgoing plate ID at each convergent boundary pair. Populated when craton data is available; uses craton-area heuristic (plate with larger craton pixel count = overriding). Not yet used to modify force directions in `polygons_to_stage_poles` — see Known Issues.

**`craton_counts`**: `Dict[int, int]`  
Maps plate ID → number of craton raster pixels on that plate. Recomputed at each save interval from `craton_counts_from_raster`. Used by `resolve_downgoing_plate` to authoritatively assign subduction polarity for output files.

**`cumulative_rotation_dict`**: `Dict[int, pygplates.FiniteRotation]`  
Running total finite rotation per plate. Updated each step as `F_new = S * F_old` (left-multiply preserves temporal order).

**`rotation_records`**: `List[Tuple[float, int, pygplates.FiniteRotation]]`  
Flat list of `(time_myr, plate_id, finite_rotation)` for the full run. Written to `.rot` at end.

**`tracers`**: `GeoDataFrame` with columns `geometry` (Point), `active` (bool), `age` (float, Myr).  
Always use `.active.astype(bool)` when masking — the column may have `dtype=object` after `pd.concat`.

**`rasters`**: `xr.Dataset` with variables `plate_id`, `velocity_east`, `velocity_north`, `divergence`, and optionally `craton`. All on a regular lon/lat grid at `raster_spacing` degree resolution.

**`boundaries_gdf`**: `GeoDataFrame` produced by `topo_polygons_to_boundaries` then exploded by `compute_boundary_divergence`. Each row is a single 2-vertex arc segment with columns `plateid1_poly1`, `plateid1_poly2`, and `divergence` (relative normal velocity, cm/yr; positive = diverging). Boundary type is derived on-the-fly by comparing `divergence` against `div_threshold_ridge` and `div_threshold_subduction`.

### Boundary divergence computation

Each boundary segment is classified independently. `compute_boundary_divergence` iterates over every consecutive vertex pair in every boundary line and computes the relative normal velocity for that arc:

1. For the arc from vertex A to vertex B, the great-circle normal is computed via `GreatCircleArc.get_great_circle_normal()`
2. The sign is determined by a hemisphere test: `dot(centroid_of_plate1, arc_normal) > 0` → normal points toward plate 1 → positive = plate 1 diverging away
3. Velocities at the arc midpoint are computed for both plates using `pygplates.calculate_velocities`
4. The relative velocity is projected onto the arc normal to give the divergence value

This is spherically correct and has no polar singularity. The raster-based `mpcalc.divergence` is retained only for visualisation.

### Synthetic earthquake generation

`generate_synthetic_earthquakes` uses Marsaglia's method to generate uniformly distributed candidate points on the sphere, then:

- **Ridge earthquakes**: accept candidates within `eq_ridge_max_dist_deg` of a divergent segment; assign random shallow depth up to `eq_ridge_max_depth_km`; thinned by `eq_ridge_prob`
- **Subduction earthquakes**: accept candidates on the overriding-plate side of convergent segments (hemisphere test); depth = `distance_from_trench_km × tan(slab_dip)`; capped at `eq_slab_max_depth_km`; sanity-checked against plate ID raster to reject points that fall within the downgoing plate's surface expression

### Synthetic topography

`compute_topography_raster` builds elevation in three layers (all ocean pixels = 0):

1. **Background**: all craton pixels set to `topo_continental_base_m`
2. **Andean belt**: craton pixels on the overriding plate within `topo_andean_max_dist_deg` of an ocean-continent convergent boundary receive a Gaussian elevation profile peaked at `topo_andean_arc_dist_deg` from the trench
3. **Himalayan belt**: pixels near continent-continent convergent boundaries (both plates have cratons) receive a linear ramp from `topo_himalayan_peak_m` at the suture to 0 at `topo_himalayan_max_dist_deg`
4. **Terrain noise**: after all smooth zones are assembled, spherically-correct bandlimited noise is added via spherical harmonic synthesis (`pyshtools`). Power spectrum ∝ `(l+1)^{-β}`. Noise amplitude scales as `elevation × noise_fraction`, so mountains are rough and lowlands are gentle.

### Seafloor bathymetry

`generate_seafloor_age_grid` uses `pygmt.blockmedian` + `pygmt.sphinterpolate` to interpolate tracer ages onto a regular grid (masked where craton > 0).

`generate_topo_bathy_grid` applies the GDH1 age-depth relation:
- Age < 20 Myr: `depth = -(2600 + 365 × √age)`
- Age ≥ 20 Myr: `depth = -(5651 − 2473 × exp(−0.0278 × age))`

Ocean pixels with no tracer coverage receive a default depth (`default_ocean_depth_m`, default −3500 m). Continental pixels (topography > 0) overwrite the bathymetry.

### Bootstrap sequence (before the main loop)

1. `reconstruction_model.plate_snapshot(0)` → initial GPlates topology snapshot
2. `pdriver.get_torques_for_plate_snapshot(...)` → initial stage poles from topology (used only once)
3. `nudge_polygons(cratons, stage_pole_dict)` → move cratons to position after first step
4. `nudge_topological_polygons(snapshot, forces, ...)` → first polygon geometries from topology
5. `compute_timestep(...)` → establishes `new_plate_polygons`, `rasters`, `boundaries_gdf`
6. `generate_seafloor_age_tracers(...)` → seed initial tracer field with estimated ages

After bootstrap, the topological snapshot is no longer used; the simulation is self-contained.

### Subduction polarity determination

For each convergent boundary pair, the downgoing plate is resolved by `resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts)`:

1. If one plate has craton pixels and the other does not → continental plate is overriding
2. If both have craton pixels → larger craton pixel count is overriding
3. If neither has craton pixels → fall back to `polarity_dict` sticky assignment (arbitrary on first encounter; consistent thereafter)

`polarity_dict` is only populated when craton data is available (guarded by `craton_raster is not None`) to prevent contamination of the dict during the bootstrap step.

### Rotation file format

GPlates `.rot` format, one line per plate per timestep:
```
plate_id   time_Ma   euler_lat   euler_lon   angle_deg   fixed_plate_id
```
Fixed plate is `000` (absolute frame). Each plate gets a `t=0` identity line automatically. Identity rotations (zero angle, undefined pole) are written as `0.0 0.0 0.0`. Lines are sorted by plate ID then time.

---

## Known issues and TODO

### MEDIUM PRIORITY — Subduction polarity not yet used in force calculations

**Location**: `polygons_to_stage_poles` in `self_driving_plates.py`.

**Problem**: `polarity_dict` is populated correctly but the slab pull and slab suction forces are applied to whichever plate happens to be `id1` or `id2` from the raster lookup, not the physically correct downgoing/overriding plate.

**Fix**: Replace the hardcoded assignment with:
```python
pair = (min(id1, id2), max(id1, id2))
downgoing_plate  = polarity_dict.get(pair, id1)
overriding_plate = id2 if downgoing_plate == id1 else id1
```

---

### LOW PRIORITY — Torque accumulation between plates is order-dependent

**Location**: `polygons_to_stage_poles`, the `torque_dict` internal to that function.

**Problem**: `torque_dict` is initialised empty at the start of each call and the check `if plate_id not in torque_dict` will always be True for the first encounter of each plate, so net_torque always starts from zero. This is correct behaviour — but the variable serves no purpose and should be removed to avoid confusion.

---

### LOW PRIORITY — `plate_driver.py` dependency is opaque

The `pdriver` module (`plate_driver.py`) is used for torque computation but is not documented here. Key functions used:
- `pdriver.get_torques_for_plate_snapshot(snapshot, time, delta_time, forces)` — bootstrap only
- `pdriver.get_torque_for_plate_polygon(snapshot, plate_id, time, delta_time, forces)` — in `nudge_topological_polygons`
- `pdriver.get_torque_for_segment(segment, normal_reversal)` — in `polygons_to_stage_poles`

---

### LOW PRIORITY — Unused code in `self_driving_plates.py`

Several functions defined in the lower half of the file appear to be unused by the current simulation pipeline (e.g. `densify_geodesic`, `expand_to_fill_gaps_spherical`, `voronoi_region_to_polygon`). These appear to be experimental utilities from earlier development. They should either be documented as utilities or moved to a separate file.
