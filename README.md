# PlateDriver — Self-Driving Plate Tectonics Simulation

## Inspiration

The core approach — computing forces on plate boundaries and integrating them forward to derive stage rotations — is inspired by the [Tectonic Explorer](https://github.com/concord-consortium/tectonic-explorer) project by the Concord Consortium. PlateDriver extends that idea into a full simulation pipeline built on [GPlates](https://www.gplates.org) / pygplates, with raster-based plate tracking, tracer advection, mantle structure synthesis, and output compatible with GPlates workflows.

## What this is

PlateDriver is a forward model of plate tectonic evolution. Starting from an initial set of tectonic plates (defined in GPlates format), the model computes the forces acting on each plate at every timestep, derives stage rotations, and evolves the plate geometry forward in time. The plates "drive themselves" — no prescribed rotation model is needed after the initial snapshot.

The model tracks:
- **Plate polygons** — boundaries evolve through each timestep driven by slab pull, ridge push, and collision forces
- **Seafloor tracers** — massless particles seeded at spreading ridges and advected with the plates, deactivated when subducted
- **Seafloor age** — estimated from tracer positions, interpolated to a regular grid each save interval
- **Hotspots and seamounts** — fixed mantle hotspot sources that erupt seamounts advected by plate motion and deactivated at subduction zones
- **Synthetic earthquakes** — generated near active boundaries with physically motivated depth geometry
- **Synthetic topography** — continental elevations with Andean/Himalayan mountain belts and spectrally-coloured terrain noise
- **Seafloor bathymetry** — derived from seafloor age via the GDH1 age-depth relation
- **Slab tracers and mantle structure** — subducted material tracked vertically through the mantle; combined with a CMB heat accumulation proxy to produce a 3D dVs/Vs anomaly field analogous to seismic tomography
- **Finite rotations** — accumulated from stage poles, written to a GPlates-compatible `.rot` file
- **Plate fragmentation** — geometrically implausible plates (U-shapes, thin necks) can be split at their narrowest point; disabled by default
- **Plate merging** — plates with persistently near-zero relative motion can be absorbed into their neighbour; disabled by default

## Inputs

Input files are in GPlates format, specified in `config.yaml`:

| Key | Purpose |
|---|---|
| `rotation_model` | Rotation model for the initial plate snapshot (`.rot`) |
| `dynamic_polygons` | Topological plate polygon files (`.gpml`), one or more |
| `cratons` | Craton (continental block) polygons (`.geojson`), used for collision detection and subduction polarity |

## How to run

```bash
conda activate pygmt17
python run_SDP.py --config config.yaml
```

All configuration is in `config.yaml`. The `--config` argument is required.

## Outputs

All outputs go to the directory specified by `output_dir` in `config.yaml`.

| File pattern | Interval | Contents |
|---|---|---|
| `plot_NNNN.png` | Every `plot_interval` steps | Map: plate polygons, boundary types, velocity vectors, cratons, seafloor tracers coloured by age, seamounts (red triangles) |
| `tracers_NNNN.csv` | Every `csv_save_interval` steps | Active seafloor tracer lon/lat/age |
| `boundaries_NNNN.geojson` | Every `csv_save_interval` steps | Boundary arc segments with divergence (cm/yr), boundary type, plate IDs, and downgoing/overriding plates for subduction zones |
| `earthquakes_NNNN.geojson` | Every `csv_save_interval` steps | Synthetic earthquake locations with depth, boundary type, plate IDs |
| `seamounts_NNNN.geojson` | Every `csv_save_interval` steps | Active seamount locations with hotspot ID and age |
| `seafloor_age_NNNN.nc` | Every `csv_save_interval` steps | Spherically interpolated seafloor age grid (Myr); continental pixels = NaN |
| `topography_NNNN.nc` | Every `csv_save_interval` steps | Synthetic continental elevation (m); ocean pixels = 0 |
| `topo_bathy_NNNN.nc` | Every `csv_save_interval` steps | Combined topography-bathymetry (m); GDH1 ocean depths merged with continental topography |
| `mantle_temp_NNNN.nc` | Every `csv_save_interval` steps | 3D dVs/Vs anomaly field (depth × lat × lon); positive = cold fast slab, negative = warm slow region. Also contains `cmb_heat` 2D layer. Only written when `slab_tracers_enabled: true` |
| `plate_id_NNNN.nc` | Every `csv_save_interval` steps | Plate ID raster (lon/lat) |
| `craton_NNNN.nc` | Every `csv_save_interval` steps | Craton ID raster (lon/lat) |
| `seamounts_final.geojson` | End of run | All active seamounts at final timestep |
| `rotations_NNNN.rot` | Every `csv_save_interval` steps | GPlates rotation file with finite rotations for all plates up to step NNNN |

## Key parameters (`config.yaml`)

### Run timing

```yaml
run:
  reconstruction_time: 0      # starting snapshot time (Ma)
  delta_time: 1.0             # timestep size (Myr)
  total_time: 200.0           # run duration (Myr)
```

### Forces

Dimensionless scale factors applied to the torque from each boundary type.

```yaml
forces:
  slab_pull:            6     # attractive force on downgoing plate at convergent boundaries
  slab_suction:         6     # attractive force on overriding plate at convergent boundaries
  ridge_push:           3     # repulsive force at divergent boundaries
  continental_collision: 18   # resistance at continent-continent boundaries
```

### Core numerical parameters

```yaml
params:
  raster_spacing:              0.25    # degrees; velocity grid and rasterisation resolution
  tessellation_deg:            0.2     # polygon edge densification before rasterising
  topo_simplify_tolerance:     0.4     # TopoJSON simplification tolerance after polygonise
  boundary_search_radius_deg:  2       # plate ID search radius at boundary segment midpoints

  div_threshold_subduction:   -0.0001  # relative normal velocity (cm/yr) below this → convergent
  div_threshold_ridge:         0.2     # relative normal velocity above this → divergent

  tracer_spacing_deg:          1.0     # along-ridge spacing for newly seeded tracer pairs
  tracer_offset_deg:           0.6     # perpendicular offset from ridge axis when seeding
  seafloor_spacing_deg:        1.0     # grid spacing for initial seafloor age tracers
  spreading_rate_km_per_myr:   50.0    # assumed spreading rate for initial age estimation

  csv_save_interval:           5       # save outputs every N timesteps
  plot_interval:               1       # make plots every N timesteps; set to null to disable
```

### Synthetic earthquakes

```yaml
params:
  eq_n_candidates:             100000  # random sphere points sampled per step
  eq_ridge_max_dist_deg:       1.0     # max distance from ridge for ridge earthquakes
  eq_ridge_max_depth_km:       50.0    # max depth for ridge earthquakes (km)
  eq_slab_dip_deg:             45.0    # assumed slab dip angle (degrees)
  eq_slab_max_depth_km:        650.0   # maximum slab depth (km)
  eq_ridge_prob:               0.15    # fraction of ridge candidates accepted (density control)
  eq_subduction_prob:          1.0     # fraction of subduction candidates accepted
```

### Synthetic topography

```yaml
params:
  topo_continental_base_m:     400     # background continental elevation (m)
  topo_noise_fraction:         0.3     # noise amplitude as fraction of local elevation
  topo_spectral_beta:          2.0     # SH power-law exponent (2 = brown / realistic terrain)
  topo_andean_max_dist_deg:    11.0    # cutoff distance from trench for Andean mountains
  topo_andean_arc_dist_deg:    4.5     # distance of arc-peak from trench (degrees)
  topo_andean_peak_m:          3000    # peak Andean elevation (m)
  topo_himalayan_max_dist_deg: 8.0     # half-width of Himalayan collision zone (degrees)
  topo_himalayan_peak_m:       5000    # peak Himalayan elevation (m)
```

### Hotspots and seamounts

```yaml
params:
  n_hotspots:    30     # number of randomly placed hotspot sources (ignored if hotspot_locations set)
  hotspot_prob:  1.0    # probability of eruption per hotspot per timestep

  # Optional: fix hotspot positions as [[lon, lat], ...]; overrides n_hotspots
  # hotspot_locations:
  #   - [25.0, -30.0]
  #   - [-150.0, 20.0]
```

Hotspot positions are fixed in the mantle frame (not reconstructed). Seamounts are advected with their host plate and deactivated when subducted.

### Slab tracer system and mantle tomography

Tracks subducted material vertically through the mantle and synthesises a 3D seismic velocity anomaly field.

```yaml
params:
  slab_tracers_enabled:    true   # enable slab tracking; produces mantle_temp_NNNN.nc
  slab_slowing_factor:     4      # lower-mantle sinking slowdown (Lithgow-Bertelloni & Richards 1998)
  slab_max_sink_rate:      400    # km/Myr; cap on sinking rate (real Earth max ~100–200 km/Myr)
  slab_thermal_tau_myr:    4      # Myr; thermal decay timescale (lower = faster fade)
  mantle_grid_spacing_deg: 5.0    # degrees; horizontal resolution of 3D mantle field
  mantle_depth_interval_km: 100   # km; depth layer thickness
```

Sinking rate: upper mantle (0–670 km) at `|convergence_vel| × 10 km/Myr`; lower mantle at `rate / slab_slowing_factor`; capped at `slab_max_sink_rate`. Slablets warm up exponentially with e-folding time `slab_thermal_tau_myr`. The dVs field is zero-mean at each depth level (matching real tomography convention) and optionally smoothed to simulate finite-resolution imaging.

### CMB heat accumulation (LLSVP proxy)

A 2D CMB heat field accumulates where subduction is absent, producing warm low-velocity anomalies in the deepest mantle layers.

```yaml
params:
  cmb_heating_rate:     0.02   # heat units added per Myr everywhere (core heating)
  cmb_slab_suppression: 0.05   # heat units removed per slablet per Myr (at full coupling)
  cmb_surface_coupling: 0.1    # fractional coupling at surface (0 = thermal only, 1 = instantaneous)
  cmb_sigma_deg:        10.0   # horizontal smoothing of CMB heat field (degrees)
  cmb_max_heat:         3.0    # cap on heat accumulation
  cmb_dVs_scale:        -0.5   # dVs % per heat unit (negative = warm/slow)
  cmb_decay_km:         400.0  # Gaussian decay length upward from CMB (km)
```

The `cmb_surface_coupling` parameter implements an instantaneous dynamic pressure coupling: surface slablets suppress CMB heat with weight `surface_coupling`, deep slablets with weight up to 1.0 (linearly interpolated with depth). This allows the CMB field to respond to surface subduction from the first timestep, without waiting for material to physically reach the base of the mantle.

### Tomographic smoothing

Blurs the final dVs field to simulate the resolution limits of real seismic tomography.

```yaml
params:
  tomo_smooth_horiz_deg:  5.0    # horizontal Gaussian sigma (degrees; 0 = off)
  tomo_smooth_depth_km:   150.0  # vertical Gaussian sigma (km; 0 = off)
```

Typical global tomography models have horizontal resolution ~500 km (~5°) and vertical resolution ~150 km.

### Plate fragmentation

Geometrically implausible plates are detected by shape screening and optionally split at their narrowest point. Both mechanisms are **disabled by default**.

```yaml
params:
  fragmentation_enabled:       false  # split implausible plates at their neck
  min_concavity_ratio:         0.85   # concavity (area/convex_hull_area) below this → check neck
  max_aspect_ratio:            4.0    # bounding-box aspect ratio above this → check neck
  min_child_area_deg2:         25.0   # minimum child plate area; prevents slivers
  fragmentation_cooldown_myr:  5.0    # Myr before a plate can fragment again
```

Each candidate plate is rasterized and skeletonized (`skimage.morphology.skeletonize`). The skeleton point with the smallest inscribed-circle radius is the neck. A great-circle cut is drawn through the neck, perpendicular to the local skeleton. The raster is flood-filled to assign the two children: the larger child keeps the parent PLATEID; the smaller child receives a new PLATEID. Both children inherit the parent stage pole and diverge naturally in subsequent timesteps as the new boundary is classified by `compute_boundary_divergence`.

### Plate merging

Plates with persistently near-zero relative motion are absorbed into their largest neighbour.

```yaml
params:
  merging_enabled:              false  # absorb plates with near-zero relative motion
  merge_conservative_threshold: 0.05   # cm/yr; |divergence| below this = near-zero
  merge_min_steps:              5      # consecutive near-zero steps before merge fires
  merge_require_no_driver:      true   # only merge if smaller plate has no ridge/subduction
  merge_cooldown_myr:           10.0   # Myr before merged plate is eligible again
```

A `_conservative_contact` counter is maintained across timesteps for each plate pair. When the mean absolute divergence stays below `merge_conservative_threshold` for `merge_min_steps` consecutive steps, the smaller plate's raster pixels are relabelled with the larger plate's ID. The merged polygon is re-derived from the updated raster; no special stage-pole initialisation is needed. A cooldown suppresses re-merging or re-fragmentation shortly after an event.

---

## Analysis scripts

### `analyse_boundaries.py`

Reads boundary GeoJSON time series and produces a CSV of boundary lengths and crust flux:

```bash
python analyse_boundaries.py --output_dir ./output
```

Output: `boundary_stats.csv` with columns `timestep`, `ridge_length_km`, `subduction_length_km`, `conservative_length_km`, `ridge_flux_km2_per_myr`, `subduction_flux_km2_per_myr`. Arc lengths computed using pygplates spherical geometry.

---

## Codebase reference

### File structure

| File | Role |
|---|---|
| `self_driving_plates.py` | Library — all simulation logic as importable functions |
| `run_SDP.py` | Driver script — config loading, bootstrapping, main loop, plotting, output |
| `plate_driver.py` | External dependency — torque computation from boundary segments |
| `analyse_boundaries.py` | Standalone post-processing — boundary length and flux time series |
| `generate_topography.py` | Standalone script — topography generation utilities |
| `generate_trench_raster.py` | Standalone script — trench raster generation |
| `blend_topobathy.py` | Standalone script — topography/bathymetry blending utilities |

### Async output worker

Plotting and file I/O are handled by a background daemon thread so the main simulation loop is never blocked waiting for a PNG to render or a NetCDF to be written. The worker (`_output_worker` in `run_SDP.py`) consumes from `_output_queue` (a `queue.Queue(maxsize=3)`). Each iteration the main loop enqueues a `('plot', payload)` task and, at save intervals, a `('save', payload)` task, then immediately proceeds to the next timestep.

The `maxsize=3` limit applies backpressure: if the worker falls more than 3 tasks behind (e.g. when saving very large mantle temperature grids) the main loop pauses briefly rather than accumulating unbounded state. At the end of the run, `_output_queue.join()` ensures all queued tasks are fully written before the process exits.

Mutable state is snapshotted before enqueueing: `polarity_dict` and `rotation_records` are shallow-copied; `cmb_heat` and `slab_tracers` are copied with `.copy()`; cratons are converted from pygplates FeatureCollection to GeoDataFrame in the main thread before being passed to the worker.

### Per-timestep data flow

```
nudge_polygons(new_plate_polygons, stage_pole_dict)
    → rpolygons (list of pygplates features, geometries advected by previous stage poles)

compute_timestep(rpolygons, forces, stage_pole_dict, delta_time, cratons, params, polarity_dict,
                 merge_pairs, fragmentation_cooldown_pids)
    → new_plate_polygons (GeoDataFrame, PLATEID1 col)
    → stage_pole_dict    (Dict[int, FiniteRotation])
    → rasters            (xr.Dataset: plate_id, velocity_east, velocity_north, divergence, craton)
    → boundaries_gdf     (GeoDataFrame: geometry(LineString), plateid1_poly1, plateid1_poly2,
                          divergence — one row per 2-vertex arc segment)
    → split_info         (list of (new_pid, parent_pid) from fragmentation, or [])

[accumulate finite rotations]
    cumulative_rotation_dict[pid] = stage_rot * cumulative_rotation_dict.get(pid, identity)

[fragmentation / merging bookkeeping]
    _conservative_contact updated per boundary pair
    merge candidates fired if counter ≥ merge_min_steps → _merge_pairs for next compute_timestep
    _plate_cooldown decremented; _fragmentation_cooldown_pids updated

seed_ridge_tracers(boundaries_gdf, ...)   → new seafloor tracers
advect_tracers(tracers, ...)              → tracers moved by stage rotation; age += delta_time
deactivate_subducted_tracers(tracers, boundaries_gdf)  → tracers near convergent segments deactivated

advect_tracers(seamounts, ...)            → seamounts moved by stage rotation; age += delta_time
deactivate_subducted_tracers(seamounts, boundaries_gdf)
seed_hotspot_seamounts(hotspot_locs, prob)  → new seamounts at hotspot locations

[if slab_tracers_enabled]
    update_slab_tracers(slab_tracers, boundaries_gdf, delta_time, ...)
        → sink existing slablets; seed new ones at subduction zones
    update_cmb_heat(cmb_heat, slab_tracers, delta_time, ...)
        → increment CMB heat everywhere; suppress where slablets present

[at csv_save_interval steps]
    generate_seafloor_age_grid(tracers, craton_raster)  → seafloor_age_NNNN.nc
    [boundary GeoJSON with resolve_downgoing_plate]      → boundaries_NNNN.geojson
    generate_synthetic_earthquakes(boundaries_gdf, ...)  → earthquakes_NNNN.geojson
    compute_topography_raster(...)                       → topography_NNNN.nc
    generate_topo_bathy_grid(age_grid, topo)             → topo_bathy_NNNN.nc
    generate_mantle_temp_field(slab_tracers, cmb_heat)  → mantle_temp_NNNN.nc

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
resolve_duplicate_plateids    rename duplicate PLATEID1s from polygon splits;
                              prefers the fragment containing craton geometry
[apply_plate_merges]          if merge_pairs non-empty: replace smaller-plate pixels in raster,
                              re-polygonize (merging_enabled only)
[detect_and_split_plates]     shape-screen plates → skeletonize → neck cut → flood-fill;
                              assign new PLATEID to smaller child (fragmentation_enabled only)
topo_polygons_to_boundaries   extract shared boundary segments via TopoJSON
compute_boundary_divergence   per-arc-segment relative normal velocity; sign determined
                              by probing the plate-ID raster on both sides of each segment
velocity_grid                 compute east/north velocity at each raster cell via pygplates
polygons_to_stage_poles       integrate boundary torques → stage rotation per plate
```

### Key data structures

**`stage_pole_dict`**: `Dict[int, pygplates.FiniteRotation]`  
Maps plate ID → stage rotation for one `delta_time` step.

**`polarity_dict`**: `Dict[Tuple[int,int], int]`  
Maps `(min_pid, max_pid)` → downgoing plate ID at each convergent boundary pair. Uses craton-area heuristic: plate with larger craton pixel count is overriding.

**`craton_counts`**: `Dict[int, int]`  
Maps plate ID → craton raster pixel count. Used by `resolve_downgoing_plate`.

**`tracers`** / **`seamounts`**: `GeoDataFrame` with columns `geometry` (Point), `active` (bool), `age` (float, Myr). Always use `.active.astype(bool)` when masking.

**`slab_tracers`**: `DataFrame` with columns `lon`, `lat`, `depth_km`, `age` (Myr since subduction), `active` (bool), `sink_rate_km_per_myr`.

**`cmb_heat`**: `ndarray` (nlat × nlon) of heat accumulation units. Built up over time by `update_cmb_heat`; injected into the deep mantle by `generate_mantle_temp_field`.

**`rasters`**: `xr.Dataset` with variables `plate_id`, `velocity_east`, `velocity_north`, `divergence`, and optionally `craton`. All on a regular lon/lat grid at `raster_spacing` degree resolution.

**`boundaries_gdf`**: `GeoDataFrame` where each row is a single 2-vertex arc segment with columns `plateid1_poly1`, `plateid1_poly2`, and `divergence` (cm/yr; positive = diverging).

### Boundary divergence computation

Each boundary segment is classified independently. `compute_boundary_divergence` iterates over every consecutive vertex pair in every boundary line:

1. The great-circle arc normal is computed via `GreatCircleArc.get_great_circle_normal()`
2. The sign is determined by probing the plate-ID raster on **both sides** of the segment midpoint at `OFFSET_DEG = 0.5°`. If the positive-direction probe returns pid1 (or the negative-direction probe returns pid2), sign = +1. This is robust for large, concave, or dateline-spanning plates where a centroid hemisphere test can give wrong signs
3. Velocities at the midpoint are computed for both plates using `pygplates.calculate_velocities`
4. The relative velocity is projected onto the arc normal to give the divergence value

This is spherically correct and has no polar singularity.

### Slab tracer physics

- **Seeding**: new slablets created along subduction segments at `sample_deg` spacing, placed at depth = 0 km
- **Sinking**: upper mantle (0–670 km): rate = `min(|divergence| × 10, max_sink_rate)` km/Myr; lower mantle: rate / `slab_slowing_factor`; deactivated on reaching CMB (~2890 km)
- **Thermal fading**: `dT = T₀ × exp(-age / slab_thermal_tau_myr)` where `T₀ = -500 K`
- **dVs conversion**: `dVs = dT × dVs_per_K` (-0.001 %/K) → cold slab = positive dVs (fast waves)
- **Final field**: Gaussian-smoothed in 3D, zero-meaned at each depth level, optionally blurred by `tomo_smooth_horiz_deg` / `tomo_smooth_depth_km`

### Synthetic earthquake generation

`generate_synthetic_earthquakes` uses Marsaglia's method to generate uniformly distributed candidate points on the sphere:

- **Ridge earthquakes**: within `eq_ridge_max_dist_deg` of a divergent segment; random shallow depth ≤ `eq_ridge_max_depth_km`; thinned by `eq_ridge_prob`
- **Subduction earthquakes**: on the overriding-plate side of convergent segments; depth = `distance_from_trench × tan(slab_dip_deg)`; capped at `eq_slab_max_depth_km`

### Synthetic topography

`compute_topography_raster` builds elevation in layers (ocean pixels = 0):

1. **Background**: craton pixels → `topo_continental_base_m`
2. **Andean belt**: overriding-plate cratons near ocean-continent convergence → Gaussian profile peaked at `topo_andean_arc_dist_deg` from trench
3. **Himalayan belt**: pixels near continent-continent boundaries → linear ramp from `topo_himalayan_peak_m` at suture to 0 at `topo_himalayan_max_dist_deg`
4. **Terrain noise**: bandlimited spherical harmonic noise (`pyshtools`), power ∝ `(l+1)^{-β}`, amplitude scales with local elevation

### Seafloor bathymetry

`generate_topo_bathy_grid` applies the GDH1 age-depth relation:
- Age < 20 Myr: `depth = -(2600 + 365 × √age)`
- Age ≥ 20 Myr: `depth = -(5651 − 2473 × exp(−0.0278 × age))`

Ocean pixels with no tracer coverage receive a default depth of −3500 m. Continental pixels overwrite the bathymetry.

### Subduction polarity determination

`resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts)`:
1. One plate has craton pixels, other does not → continental plate is overriding
2. Both have craton pixels → larger craton pixel count is overriding
3. Neither has craton pixels → fall back to `polarity_dict` sticky assignment

### Rotation file format

GPlates `.rot` format:
```
plate_id   time_Ma   euler_lat   euler_lon   angle_deg   000
```
Fixed plate is `000` (absolute frame). Each plate gets a `t=0` identity line. Sorted by plate ID then time.

---

## Known issues and TODO

### MEDIUM PRIORITY — Subduction polarity not used in force calculations

`polarity_dict` is populated correctly but slab pull and slab suction forces in `polygons_to_stage_poles` are applied based on arbitrary `id1`/`id2` ordering from the raster lookup, not the physically resolved downgoing/overriding plate. Fix: look up `polarity_dict` when assigning forces at each convergent boundary pair.

### LOW PRIORITY — `plate_driver.py` dependency is opaque

The `pdriver` module is used for torque computation but not documented here. Key functions:
- `pdriver.get_torques_for_plate_snapshot(snapshot, time, delta_time, forces)` — bootstrap only
- `pdriver.get_torque_for_segment(segment, normal_reversal)` — inside `polygons_to_stage_poles`
