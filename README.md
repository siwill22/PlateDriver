# PlateDriver — Self-Driving Plate Tectonics Simulation

## What this is

PlateDriver is a forward model of plate tectonic evolution. Starting from an initial set of tectonic plates (defined in GPlates format), the model computes the forces acting on each plate at every timestep, derives stage rotations, and evolves the plate geometry forward in time. The plates "drive themselves" — no prescribed rotation model is needed after the initial snapshot.

The model tracks:
- **Plate polygons** — boundaries evolve through each timestep
- **Seafloor tracers** — massless particles seeded at spreading ridges and advected with the plates, deactivated when subducted
- **Seafloor age** — estimated from tracer distance to nearest ridge, incremented each step
- **Finite rotations** — accumulated from stage poles, written to a GPlates-compatible `.rot` file

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
| `test_plot_NNNN.png` | Every step | Map: divergence field, plate polygons, velocity vectors, cratons, tracers coloured by age |
| `tracers_NNNN.csv` | Every `csv_save_interval` steps | Active tracer lon/lat/age |
| `divergence_NNNN.nc` | Every `csv_save_interval` steps | Divergence field (xarray, lon/lat) |
| `plate_id_NNNN.nc` | Every `csv_save_interval` steps | Plate ID raster (xarray, lon/lat) |
| `craton_NNNN.nc` | Every `csv_save_interval` steps | Craton ID raster (xarray, lon/lat) |
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
total_time          = 60.   # run duration (Myr); loop runs steps 1..int(total_time/delta_time)
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
    'raster_spacing':           1,        # degrees; controls velocity grid and rasterisation resolution
    'tessellation_deg':         0.2,      # polygon edge densification before rasterising
    'topo_simplify_tolerance':  1,        # TopoJSON simplification tolerance after polygonise
    'boundary_search_radius_deg': 2,      # plate ID neighbour search radius at segment midpoints
    'div_threshold_subduction': -0.00025, # divergence below this → convergent boundary
    'div_threshold_ridge':       0.0001,  # divergence above this → divergent boundary
    'tracer_spacing_deg':        2.0,     # along-ridge spacing for newly seeded tracer pairs
    'tracer_offset_deg':         0.6,     # perpendicular offset from ridge axis when seeding
    'seafloor_spacing_deg':      2.0,     # grid spacing for initial seafloor age tracers
    'spreading_rate_km_per_myr': 50.0,    # assumed spreading rate for initial age estimation
    'csv_save_interval':         10,      # save NetCDF/CSV outputs every N timesteps
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
    → boundaries_gdf     (GeoDataFrame: geometry(LineString), poly1, poly2,
                          plateid1_poly1, plateid1_poly2)

[accumulate finite rotations]
    cumulative_rotation_dict[pid] = stage_rot * cumulative_rotation_dict.get(pid, identity)
    rotation_records.append((t, pid, cumulative_rotation_dict[pid]))

seed_ridge_tracers(boundaries_gdf, rasters['divergence'], ...)
    → new_tracers (GeoDataFrame: geometry, active=True, age=0)

pd.concat([tracers, new_tracers])

advect_tracers(tracers, rasters['plate_id'], stage_pole_dict, delta_time)
    → tracers (active points moved by their plate's stage rotation; age += delta_time)

deactivate_subducted_tracers(tracers, rasters['divergence'], threshold)
    → tracers (active=False where divergence < threshold)

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
velocity_grid                 compute east/north velocity at each raster cell via pygplates
mpcalc.divergence             divergence of velocity field (see Known Issues)
polygons_to_stage_poles       integrate boundary torques → stage rotation per plate
```

### Key data structures

**`stage_pole_dict`**: `Dict[int, pygplates.FiniteRotation]`  
Maps plate ID → stage rotation for one `delta_time` step. Used for advection and velocity computation.

**`polarity_dict`**: `Dict[Tuple[int,int], int]`  
Maps `(min_pid, max_pid)` → downgoing plate ID at each convergent boundary pair. Populated on first encounter using craton-area heuristic (larger craton = overriding); sticky thereafter. Currently built but not yet used to modify force calculations.

**`cumulative_rotation_dict`**: `Dict[int, pygplates.FiniteRotation]`  
Running total finite rotation per plate. Updated each step as `F_new = S * F_old` (left-multiply preserves temporal order).

**`rotation_records`**: `List[Tuple[float, int, pygplates.FiniteRotation]]`  
Flat list of `(time_myr, plate_id, finite_rotation)` for the full run. Written to `.rot` at end.

**`tracers`**: `GeoDataFrame` with columns `geometry` (Point), `active` (bool), `age` (float, Myr).  
Always use `.active.astype(bool)` when masking — the column may have `dtype=object` after `pd.concat`.

**`rasters`**: `xr.Dataset` with variables `plate_id`, `velocity_east`, `velocity_north`, `divergence`, and optionally `craton`. All on a regular lon/lat grid at `raster_spacing` degree resolution. `plate_id` and `craton` have `lon`/`lat` dims (renamed from `x`/`y` after `rasterize_polygons`).

**`boundaries_gdf`**: `GeoDataFrame` produced by `topo_polygons_to_boundaries`. Each row is a shared boundary arc (LineString) with columns `poly1`, `poly2` (polygon indices) and `plateid1_poly1`, `plateid1_poly2` (plate IDs on each side). Dateline/pole-crossing arcs are split and cleaned.

### Bootstrap sequence (before the main loop)

1. `reconstruction_model.plate_snapshot(0)` → initial GPlates topology snapshot
2. `pdriver.get_torques_for_plate_snapshot(...)` → initial stage poles from topology (used only once)
3. `nudge_polygons(cratons, stage_pole_dict)` → move cratons to position after first step
4. `nudge_topological_polygons(snapshot, forces, ...)` → first polygon geometries from topology
5. `compute_timestep(...)` → establishes `new_plate_polygons`, `rasters`, `boundaries_gdf`
6. `generate_seafloor_age_tracers(...)` → seed initial tracer field with estimated ages

After bootstrap, the topological snapshot is no longer used; the simulation is self-contained.

### Subduction polarity determination (`polarity_dict`)

For each convergent boundary pair encountered for the first time, the downgoing plate is assigned by:
1. If one plate has cratons and the other does not → continental plate is overriding
2. If both have cratons → larger craton pixel count is overriding
3. If neither has cratons → `id1` is assigned as downgoing (arbitrary; consistent thereafter)

Craton pixel counts are precomputed once per call to `polygons_to_stage_poles` by counting pixels in `craton_raster` per plate ID.

**Current status**: `polarity_dict` is populated correctly each timestep but the stored polarity is not yet used to modify force directions in `polygons_to_stage_poles`. The `downgoing_plate = id1` / `overriding_plate = id2` assignment remains the same arbitrary ordering as before. This is the next intended use of `polarity_dict`.

### Rotation file format

GPlates `.rot` format, one line per plate per timestep:
```
plate_id   time_Ma   euler_lat   euler_lon   angle_deg   fixed_plate_id
```
Fixed plate is `000` (absolute frame). Each plate gets a `t=0` identity line automatically. Identity rotations (zero angle, undefined pole) are written as `0.0 0.0 0.0`. Lines are sorted by plate ID then time.

---

## Known issues and TODO

### HIGH PRIORITY — Divergence calculation is incorrect near the poles

**Location**: `self_driving_plates.py`, line `div = mpcalc.divergence(velocity_grid_east, velocity_grid_north)` inside `compute_timestep`.

**Problem**: `mpcalc.divergence` falls back to a Cartesian finite-difference calculation (treating degrees as a distance unit) because the input DataArrays have no MetPy CRS metadata. The correct spherical formula is:

```
div = 1/(R·cos φ) · [∂u/∂λ + ∂(v·cos φ)/∂φ]
```

Even implementing this formula directly would not fully solve the problem, because the east/north velocity decomposition itself is ill-conditioned near the poles: the east unit vector rotates continuously with longitude at any fixed latitude, so a uniformly moving plate produces rapidly varying `u` components near the pole from purely geometric causes (not physics). The `1/cos φ` denominator then amplifies this noise without bound.

**Correct solution**: Replace the raster-based divergence entirely with a **per-segment relative normal velocity** approach:
- For each boundary segment, both plate IDs are already known from `get_two_nearest_ids`
- Compute the velocity of each plate at the segment midpoint using `pygplates.calculate_velocities` and the respective stage rotations
- The dot product of relative velocity `(v₁ − v₂)` with the outward boundary normal gives convergence (negative) or divergence (positive)
- This is correct everywhere on the sphere, has no coordinate singularity, and is more physically natural

**Implications of this refactor**: The divergence raster is currently consumed by:
1. `polygons_to_stage_poles` — boundary segment classification (ridge vs subduction)
2. `seed_ridge_tracers` — locate ridge positions for tracer seeding
3. `deactivate_subducted_tracers` — detect subduction zones via raster lookup

After the refactor, item 1 would be computed inline in `polygons_to_stage_poles`. Items 2 and 3 would query `boundaries_gdf` (which carries per-segment divergence flags) rather than sampling the raster.

The divergence raster can still be produced for visualisation (it is plotted and saved as NetCDF), but should not be used as the primary signal for boundary classification.

---

### MEDIUM PRIORITY — Subduction polarity not yet used in force calculations

**Location**: `polygons_to_stage_poles` in `self_driving_plates.py`, lines ~231–233.

**Problem**: `polarity_dict` is populated correctly but `downgoing_plate = id1` / `overriding_plate = id2` remains an arbitrary assignment. The slab pull and slab suction forces are applied to whichever plate happens to be `id1` or `id2` from the raster lookup, not the physically correct plate.

**Fix**: Replace the hardcoded assignment with:
```python
pair = (min(id1, id2), max(id1, id2))
downgoing_plate  = polarity_dict.get(pair, id1)
overriding_plate = id2 if downgoing_plate == id1 else id1
```

---

### LOW PRIORITY — Torque accumulation between plates is order-dependent

**Location**: `polygons_to_stage_poles`, the `torque_dict` internal to that function.

**Problem**: `torque_dict` is initialised empty at the start of each call and populated as plates are processed. The check `if plate_id not in torque_dict.keys()` will always be True for the first encounter of each plate, so net_torque always starts from zero. This is correct behaviour — but the variable `torque_dict` serves no purpose and should be removed to avoid confusion.

---

### LOW PRIORITY — `plate_driver.py` dependency is opaque

The `pdriver` module (`plate_driver.py`) is used for torque computation but is not documented here. Key functions used:
- `pdriver.get_torques_for_plate_snapshot(snapshot, time, delta_time, forces)` — bootstrap only
- `pdriver.get_torque_for_plate_polygon(snapshot, plate_id, time, delta_time, forces)` — in `nudge_topological_polygons`
- `pdriver.get_torque_for_segment(segment, normal_reversal)` — in `polygons_to_stage_poles`

---

### LOW PRIORITY — `gpd` import unused after tracer initialisation change

`run_SDP.py` imports `geopandas as gpd` but no longer directly references it after the initial empty `GeoDataFrame` was replaced by `generate_seafloor_age_tracers`. The import is harmless but can be removed if desired.

---

### LOW PRIORITY — Unused code in `self_driving_plates.py`

Several functions defined below line 694 appear to be unused by the current simulation pipeline (e.g. `densify_geodesic`, `expand_to_fill_gaps_spherical`, `voronoi_region_to_polygon`). These appear to be experimental utilities from earlier development. They should either be documented as utilities or moved to a separate file.
