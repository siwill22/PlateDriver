# OldDriver — Self-Driving Plate Tectonics Simulation

## Inspiration

The core approach — computing forces on plate boundaries and integrating them forward to derive stage rotations — is inspired by the [Tectonic Explorer](https://github.com/concord-consortium/tectonic-explorer) project by the Concord Consortium. OldDriver extends that idea into a full simulation pipeline built on [GPlates](https://www.gplates.org) / pygplates, with raster-based plate tracking, tracer advection, mantle structure synthesis, and output compatible with GPlates workflows.

The code is only intended for teaching, specifically to generate learning materials (synthetic datasets for Earth-like planets) used during the 'Southern Ocean Marine Geosciences' course at OUC Qingdao.

## What this is

OldDriver is a forward model of plate tectonic evolution. Starting from an initial set of tectonic plates (defined in GPlates format), the model computes the forces acting on each plate at every timestep, derives stage rotations, and evolves the plate geometry forward in time. The plates "drive themselves" — no prescribed rotation model is needed after the initial snapshot.

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
- **Plate fragmentation** — plates with a thin neck can be split at their narrowest point; off unless `fragmentation_enabled: true` (the shipped `config.yaml` turns it on)
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

All configuration is in the YAML file given by `--config` (default `config.yaml`). Missing input files are reported before the run starts.

### Requirements

Python ≥ 3.10 with pygplates (≥ 0.36; developed with 1.0), [gprm](https://github.com/siwill22/GPlatesReconstructionModel), PyGMT, geopandas, shapely ≥ 2, xarray, numpy, scipy, scikit-learn, scikit-image, pyproj, topojson, xarray-spatial, pyshtools and PyYAML. The `pygmt17` conda environment has all of these.

### Units and time

- **Velocities and boundary divergence** are in **km/Myr**, which is numerically the same as **mm/yr** (1 cm/yr = 10 km/Myr).
- **Time**: step `i` is the state `t = i × delta_time` Myr after the starting snapshot. Step 0 is the starting configuration itself and is saved and plotted like any other step. The final step is always saved.
- **Each step** moves the plates, cratons, seafloor tracers and seamounts together, using the stage rotations computed at the previous step. It then rebuilds the plate configuration and computes the stage rotations for the next move.

## Outputs

All outputs go to the directory specified by `output_dir` in `config.yaml`.

| File pattern | Interval | Contents |
|---|---|---|
| `plot_NNNN.png` | Every `plot_interval` steps | Map: plate polygons, boundary types, velocity vectors, cratons, seafloor tracers coloured by age, seamounts (red triangles) |
| `tracers_NNNN.csv` | Every `csv_save_interval` steps | Active seafloor tracer lon/lat/age |
| `boundaries_NNNN.geojson` | Every `csv_save_interval` steps | Boundary arc segments with divergence (km/Myr = mm/yr), boundary type, plate IDs, and downgoing/overriding plates for subduction zones |
| `earthquakes_NNNN.geojson` | Every `csv_save_interval` steps | Synthetic earthquake locations with depth, boundary type, plate IDs |
| `seamounts_NNNN.geojson` | Every `csv_save_interval` steps | Active seamount locations with hotspot ID and age |
| `seafloor_age_NNNN.nc` | Every `csv_save_interval` steps | Spherically interpolated seafloor age grid (Myr); continental pixels = NaN |
| `topography_NNNN.nc` | Every `csv_save_interval` steps | Synthetic continental elevation (m); ocean pixels = 0 |
| `topo_bathy_NNNN.nc` | Every `csv_save_interval` steps | Combined topography-bathymetry (m); GDH1 ocean depths merged with continental topography |
| `mantle_temp_NNNN.nc` | Every `csv_save_interval` steps | 3D dVs/Vs anomaly field (depth × lat × lon); positive = cold fast slab, negative = warm slow region. Also contains `cmb_heat` 2D layer. Only written when `slab_tracers_enabled: true` |
| `plate_id_NNNN.nc` | Every `csv_save_interval` steps | Plate ID raster (lon/lat) |
| `craton_NNNN.nc` | Every `csv_save_interval` steps | Craton ID raster (lon/lat) |
| `divergence_NNNN.nc` | Every `csv_save_interval` steps | Horizontal divergence of the plate velocity field (1/Myr); ≈0 inside rigid plates, positive at spreading and negative at converging boundaries |
| `seamounts_final.geojson` | End of run | All active seamounts at final timestep |
| `step_stats.csv` | Every step | One row per step: global RMS plate speed, fastest plate, net rotation of the lithosphere (rate and pole), trench and ridge length, tracer counts, slab-force totals and mean trench slab depth (slab modes), seconds per step by phase |
| `trench_forces_NNNN.csv` | Every `csv_save_interval` steps | One row per trench segment of each plate: position, plates, downgoing or not, convergence, slab depth, seafloor age, pull and suction weights. Only with `slab_pull: ramp/slab_length`, `slab_suction: depth_split` or `polarity: age` |
| `plate_stats.csv` | Every step | One row per plate per step: angular velocity (pole, rate), RMS speed, area, continental fraction, perimeter and trench (downgoing/overriding), ridge and conservative lengths |
| `run_info.yaml` | Start of run | The config used, the resolved force modes and the git commit |
| `rotations_NNNN.rot` | Every `csv_save_interval` steps | GPlates rotation file with finite rotations for all plates up to step NNNN. The time column is elapsed model time (see [Rotation file format](#rotation-file-format)) |

## Key parameters (`config.yaml`)

### Run timing

```yaml
run:
  reconstruction_time: 0      # starting snapshot of the input rotation model (Ma)
  delta_time: 1.0             # timestep size (Myr)
  total_time: 200.0           # run duration (Myr)
```

### Force model

```yaml
force_model: simple      # simple | physical
force_modes: {}          # per-component overrides, e.g. {collision: coupling}
```

`force_model` picks a preset for every component of the force model; an entry under
`force_modes` overrides one component. `simple` is the original boundary-force model
(and the fastest). `physical` is the v1 force model (basal drag, age-dependent slab
pull and ridge push, and so on), which is being built in stages; see
`force_model_v1_plan.md`. Until it is complete, choosing a mode that does not exist
yet stops the run with a message. The resolved modes are printed at the start of the
run and written to `run_info.yaml`. The collision mode (`resistance` or `coupling`)
is `force_modes.collision`; the older `params.collision_mode` is still read.

Implemented so far:

- **`drag`:**
  - `none`: torque = angular velocity, as in the original model.
  - `area`: basal drag over the plate area, ω = D⁻¹τ, where D is the area tensor
    Σ w(𝟙 − r rᵀ)dA; w = `continent_drag_factor` under continents and 1 under ocean.
  - `area_edge` (`physical`): the same, plus boundary resistance equal to basal drag over
    a strip `edge_drag_width_km` wide along every boundary.
- **`speed_normalisation`:**
  - `none`: velocities are multiplied by `speed_scale` (default 1).
  - `calibrate_once` (`physical`): the scale factor k is set at step 0 so that the
    area-weighted RMS plate speed equals `target_rms_speed`, then kept for the whole run,
    so RMS speed still varies with the boundary configuration.
  - `per_step`: k is recomputed every step, which fixes the RMS speed.

  In every mode, a step whose RMS speed would exceed `max_rms_speed_factor` ×
  `target_rms_speed` is scaled down to that limit and logged. The k used is recorded in
  `run_info.yaml` and in each row of `step_stats.csv`. Setting `speed_scale` fixes k,
  which is useful for comparing two configs at the same scale.
- **`slab_pull`:**
  - `hard`: `forces.slab_pull` on every convergent segment of the downgoing plate, as
    in the original model.
  - `ramp`: the same, × min(1, convergence / `pull_ramp_km_per_myr`).
  - `slab_length` (`physical`): × √(min(age, A)/A) × min(1, z / `slab_ref_depth_km`).
    Here age is the seafloor age 1° inside the downgoing plate, A =
    `slab_age_saturation_myr`, and z is the depth of the deepest slab tracer within
    `slab_search_radius_deg` on the overriding side. The pull lasts as long as the slab
    is there, whatever the convergence rate, and stops when the segment turns divergent.
- **`slab_suction`:**
  - `ratio`: `forces.slab_suction` on the overriding plate, as in the original model.
  - `depth_split` (`physical`): the part of the slab below `slab_ref_depth_km` pulls
    *both* plates toward the trench, with weight `forces.slab_suction_deep` ×
    min(1, (z − z_ref) / z_ref).
- **`polarity`** (which plate subducts when a pair first converges; the decision is
  kept):
  - `id`: the lower plate ID.
  - `age` (`physical`): the plate with the older seafloor along the contact.
  - The continent rule comes first. With `id` it looks at whole plates: a plate with
    continent anywhere overrides one without, on every call. With `age` it looks at the
    contact: a plate overrides if its own continent lies within
    `boundary_search_radius_deg` along at least half the contact. The decision is then
    kept, except that it flips if continent reaches the trench on the downgoing side
    against an oceanic overriding side.

- **`ridge_push`:**
  - `boundary`: `forces.ridge_push` on every divergent segment, as in the original model.
  - `body` (`physical`): a body force over young seafloor, (weight / cutoff) ×
    ∇min(age, `ridge_push_age_cutoff_myr`), pointing toward older seafloor. The gradient
    is taken within each plate. One full flank from the ridge out to the cutoff gives the
    same force per unit ridge length as `boundary`. A plate whose seafloor never reaches the
    cutoff gets less push. Each saved step writes `ridge_push_torques_NNNN.csv` with the
    body and boundary-mode torque of every plate.

The modes that read seafloor age or slab depth build them only when selected. The
age grid is on the `physics_grid_deg` grid: each cell takes the nearest seafloor
tracer on its own plate, moved one step ahead. With `slab_length` or `depth_split`,
slab tracers are kept even if `slab_tracers_enabled` is false (that flag then only
controls the mantle outputs). The step-0 trenches start with slab columns down to
`initial_slab_depth_km`, as the model would have built them. Each saved step writes
`trench_forces_NNNN.csv` with the slab depth, seafloor age and weights of every
trench segment.

### Forces

Dimensionless scale factors applied to the torque from each boundary type.

```yaml
forces:
  slab_pull:            6     # attractive force on downgoing plate at convergent boundaries
  slab_suction:         6     # attractive force on overriding plate at convergent boundaries (slab_suction: ratio)
  slab_suction_deep:    3     # suction from the lower-mantle slab, on both plates (slab_suction: depth_split)
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

  div_threshold_subduction:   -0.0001  # relative normal velocity (km/Myr = mm/yr) below this → convergent
  div_threshold_ridge:         0.2     # relative normal velocity (km/Myr) above this → divergent

  tracer_spacing_deg:          1.0     # along-ridge spacing for newly seeded tracer pairs
  tracer_offset_deg:           0.6     # perpendicular offset from ridge axis when seeding
  seafloor_spacing_deg:        1.0     # grid spacing for initial seafloor age tracers
  spreading_rate_km_per_myr:   50.0    # assumed spreading rate for initial age estimation

  csv_save_interval:           5       # save outputs every N timesteps
  plot_interval:               1       # make plots every N timesteps; set to null to disable
  random_seed:                 null    # integer → repeatable run (hotspots, seamounts, earthquakes, terrain noise)
```

Simplification works in lon/lat, but pygplates reads every edge as a great-circle arc. So polygon edges are densified to ≤ 1° wherever polygons are handed to pygplates, and the extracted boundary lines get the same treatment. Without this, long edges would bow toward the pole. Boundaries are extracted from the undensified polygons, so shared edges are always recognised.

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
  slab_sink_rate_factor:   1.0    # upper-mantle sinking rate = factor × convergence rate
  slab_max_sink_rate:      400    # km/Myr; cap on sinking rate (real Earth max ~100–200 km/Myr)
  slab_thermal_tau_myr:    4      # Myr; thermal decay timescale (lower = faster fade)
  mantle_grid_spacing_deg: 5.0    # degrees; horizontal resolution of 3D mantle field
  mantle_depth_interval_km: 100   # km; depth layer thickness
```

Sinking rate: upper mantle (0–670 km) at `slab_sink_rate_factor × |convergence|` km/Myr; lower mantle at `rate / slab_slowing_factor`; capped at `slab_max_sink_rate`.

> Runs made before the unit fix sank slabs 10× faster: convergence in km/Myr was treated as cm/yr and multiplied by 10. Set `slab_sink_rate_factor: 10` to reproduce them. With the physical factor of 1, you may want a longer `slab_thermal_tau_myr` so slabs stay visible in the lower mantle.

Each slablet has a weight of `delta_time × min(1, |convergence| / 10 km/Myr)`. This keeps the mantle signal independent of the time step and stops near-stalled boundaries from building up slablets. Slablets warm up exponentially with e-folding time `slab_thermal_tau_myr` and are removed when they reach the CMB. Horizontal smoothing uses a spherical Gaussian (great-circle distance, correct across ±180° and the poles). The dVs field has zero area-weighted mean at each depth level, matching the real tomography convention, and can optionally be smoothed to simulate finite-resolution imaging.

### CMB heat accumulation (LLSVP proxy)

A 2D CMB heat field accumulates where subduction is absent, producing warm low-velocity anomalies in the deepest mantle layers.

```yaml
params:
  cmb_heating_rate:     0.02   # heat units added per Myr everywhere (core heating)
  cmb_slab_suppression: 0.05   # heat units removed per full-weight slablet per Myr (at full coupling)
  cmb_surface_coupling: 0.1    # fractional coupling at surface (0 = thermal only, 1 = instantaneous)
  cmb_sigma_deg:        10.0   # horizontal smoothing of CMB heat field (degrees per √Myr; spherical Gaussian)
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

Plates with a thin neck can be split at their narrowest point. Both fragmentation and merging are off unless enabled; the shipped `config.yaml` turns fragmentation on.

```yaml
params:
  fragmentation_enabled:       true   # split plates at a thin neck
  min_neck_ratio:              0.35   # cut if neck half-width / inscribed radius of the smaller part < this
  min_concavity_ratio:         0.85   # screening: concavity (area / convex-hull area) below this → search for a neck
  max_aspect_ratio:            2.0    # screening: aspect ratio above this → search for a neck
  min_child_area_deg2:         25.0   # minimum child plate area; prevents slivers
  min_peninsula_area_deg2:     5.0    # smaller minimum for very thin necks (neck_ratio < 0.2)
  fragmentation_cooldown_myr:  5.0    # Myr before a plate can fragment again
```

Areas are spherical, given in equatorial-equivalent deg².

1. **Screening.** Concavity and aspect ratio are measured in an equal-area projection centred on the plate. Only plates that fail either test are examined further.
2. **Neck search.** The plate is eroded by increasing great-circle distances from its edge. The first distance at which the eroded core breaks into two substantial pieces gives the neck half-width. Each plate pixel is then assigned to the nearer piece, which places the cut through the neck.
3. **Split rule.** A U- or C-shaped plate of uniform width, or any convex plate, has no neck and is never cut.
4. **Cratons.** Cratons stay with the parent plate ID. A cut that would separate cratons of the same plate is skipped.
5. **New plates.** The child plate gets a new, never-reused PLATEID. It inherits the parent's stage pole and full rotation history, then diverges naturally as its new boundary is classified.

Disconnected pieces of one plate ID (for example, after overlaps are resolved) are handled the same way at every step:
- Pieces smaller than `min_peninsula_area_deg2` are absorbed by their neighbour.
- Larger pieces become child plates.

### Plate merging

Plates with persistently near-zero relative motion are absorbed into their largest neighbour.

```yaml
params:
  merging_enabled:              false  # absorb plates with near-zero relative motion
  merge_conservative_threshold: 0.05   # km/Myr (= mm/yr); |divergence| below this = near-zero
  merge_min_steps:              5      # consecutive near-zero steps before merge fires
  merge_require_no_driver:      true   # only merge if smaller plate has no ridge/subduction
  merge_cooldown_myr:           10.0   # Myr before merged plate is eligible again
```

A `_conservative_contact` counter is maintained across timesteps for each plate pair. When the mean absolute divergence stays below `merge_conservative_threshold` for `merge_min_steps` consecutive steps, the smaller plate's raster pixels are relabelled with the larger plate's ID. Plate size is spherical area. With `merge_require_no_driver`, a merge only happens if the smaller plate has no ridge or subduction boundary. Cratons of an absorbed plate move with the host from then on. The merged polygon is re-derived from the updated raster; no special stage-pole initialisation is needed. A cooldown suppresses re-merging or re-fragmentation shortly after an event.

---

## Analysis scripts

### `analyse_boundaries.py`

Reads boundary GeoJSON time series and produces a CSV of boundary lengths and crust flux:

```bash
python analyse_boundaries.py --output_dir ./output [--delta_time 1.0]
```

Output: `boundary_stats.csv` with columns `timestep` (plus `time_myr` when `--delta_time` is given), `ridge_length_km`, `subduction_length_km`, `conservative_length_km`, `ridge_flux_km2_per_myr`, `subduction_flux_km2_per_myr`. Arc lengths computed using pygplates spherical geometry; flux = length × |divergence| (km × km/Myr).

### Ties on the grid (`tie_break`)

On a regular grid, a cell is often exactly as near to two cells of different plates, for
example when a gap between plate polygons is filled. `tie_break: hash` settles
each tie with a coin flip seeded by `tie_break_seed` and the cell's position. It is unbiased
and the same every step, and it does not depend on how the nearest cells are found.
`tie_break: previous_owner` (the default) gives the cell to the plate that held it in the
previous step if that plate is one of the tied ones, and otherwise uses the coin flip. This roughly halves A → B → A jitter of boundary cells
without slowing boundary migration (see `step_stats.csv`: `owner_change_area_frac`,
`owner_flicker_area_frac`). The model is sensitive to such details: two runs that differ
only in how ties are broken diverge visibly after a few tens of Myr, while their statistics
stay alike.

### `diagnose_run.py`

Diagnostic figures (PyGMT) for one or more runs, overlaid for comparison:

```bash
python diagnose_run.py ./output_a ./output_b --labels simple physical [--out figs] [--steps 0 50 100]
```

- `fu_plot.png`: plate RMS speed against trench fraction and against continental fraction (Forsyth & Uyeda 1975 found speed rises with the first and falls with the second).
- `timeseries.png`: RMS speed and fastest plate, net rotation, number of plates, seconds per step.
- `maps_NNNN.png`: seafloor age, boundaries by type and plate velocity arrows.

### `regression_check.py`

Compares two output folders file by file (grids, CSVs, GeoJSONs, rotation files) and exits with status 1 if anything differs:

```bash
python regression_check.py ./output_reference ./output_new [--tol 0]
```

Runs with the same `random_seed` are repeatable, so `force_model: simple` is checked against a reference run this way.

---

## Codebase reference

### File structure

| File | Role |
|---|---|
| `self_driving_plates.py` | Library — all simulation logic as importable functions |
| `run_SDP.py` | Driver script — config loading, bootstrapping, main loop, plotting, output |
| `plate_driver.py` | Torque helpers: segment torque, torque → stage rotation, and the bootstrap torques from the GPlates boundary feature types |
| `force_modes.py` | Force-model presets (`simple`, `physical`), mode choices and their validation |
| `step_stats.py` | Per-step kinematic statistics (`step_stats.csv`, `plate_stats.csv`) |
| `diagnose_run.py` | Diagnostic figures for one or more runs |
| `regression_check.py` | File-by-file comparison of two output folders |
| `analyse_boundaries.py` | Standalone post-processing — boundary length and flux time series |
| `generate_topography.py` | Standalone script: blends synthetic topography with a reference bathymetry and adds trenches (from `boundaries_NNNN.geojson`, ocean pixels only) and optional seamounts. Accepts −180..180 or 0..360 grids. Holds the shared helpers used by the two scripts below |
| `generate_trench_raster.py` | Standalone script: trench-depth anomaly raster from `boundaries_NNNN.geojson` (writes `trench_NNNN.nc` next to the input by default) |
| `blend_topobathy.py` | Standalone script: the same blend as `generate_topography.py`, but takes a pre-computed trench NetCDF. `--topo`, `--bathy`, `--trench` and `--output` are required |
| `generate_ocean_currents.py` | Standalone script: Stommel (1948) wind-driven barotropic currents from a `topo_bathy_NNNN.nc` (or Stommel–Arons abyssal circulation with `--deep`). `--topo` is required. Idealised winds are easterly trades from the equator to 30° and westerlies peaking at 60° (`--phi-max`). Default drag `--R-drag 1e-6` s⁻¹. Currents are depth-averaged, and ψ = 0 on all land, so there is no net flow between landmasses (e.g. no circumpolar current) |

### Async output worker

Plotting and file I/O are handled by a background daemon thread so the main simulation loop is never blocked waiting for a PNG to render or a NetCDF to be written. The worker (`_output_worker` in `run_SDP.py`) consumes from `_output_queue` (a `queue.Queue(maxsize=3)`). At plot intervals the main loop enqueues a `('plot', payload)` task and at save intervals a `('save', payload)` task, then immediately proceeds to the next timestep.

Each output file is written independently: if one fails (for example the seafloor age grid), the error is printed with a full traceback and the other outputs for that step are still written. The run ends with a list of failed outputs and a non-zero exit code. The queue is drained even if the main loop stops with an error or Ctrl-C.

The `maxsize=3` limit applies backpressure: if the worker falls more than 3 tasks behind (e.g. when saving very large mantle temperature grids) the main loop pauses briefly rather than accumulating unbounded state. At the end of the run, `_output_queue.join()` ensures all queued tasks are fully written before the process exits.

Mutable state is snapshotted before enqueueing: `polarity_dict` and `rotation_records` are shallow-copied; `cmb_heat` and `slab_tracers` are copied with `.copy()`; cratons are converted from pygplates FeatureCollection to GeoDataFrame in the main thread before being passed to the worker.

### Per-timestep data flow

```
[step 0]
    pdriver.get_torques_for_plate_snapshot(...)   → seed stage poles from GPlates boundary types
                                                    (only used to classify the initial boundaries)
    compute_timestep(initial polygons, seed poles, ...) → configuration at t = 0, stage poles S0
    generate_seafloor_age_tracers, generate_hotspot_locations; outputs for step 0

[step i ≥ 1, t = i × delta_time]
    moving_poles = stage poles from step i-1
    nudge_polygons(plates, moving_poles); nudge_polygons(cratons, moving_poles + craton remap)
    cumulative_rotation_dict[pid] = moving_pole * cumulative_rotation_dict[pid]   → rotation_records(t)

    compute_timestep(rpolygons, forces, moving_poles, delta_time, cratons, params, polarity_dict,
                     merge_pairs, fragmentation_cooldown_pids, collision_state, craton_id_remap,
                     next_plate_id)
        → new_plate_polygons (GeoDataFrame, PLATEID1 col)
        → stage_pole_dict    (Dict[int, FiniteRotation] — for the next move)
        → rasters            (xr.Dataset: plate_id, velocity_east, velocity_north, divergence, craton)
        → boundaries_gdf     (one row per 2-vertex arc segment, divergence in km/Myr)
        → new_plates         ([(new_pid, parent_pid)] — children copy the parent's rotation history)
        → collision_state, next_plate_id

    [fragmentation / merging bookkeeping]  merge candidates → _merge_pairs for the next step; cooldowns

    advect_tracers(tracers, previous plate_id raster, moving_poles)   → same move as the plates
    deactivate_subducted_tracers(tracers, boundaries_gdf, ...)        → only on the downgoing plate
    seed_ridge_tracers(boundaries_gdf, ...)                           → new crust, age 0
    [same for seamounts; then seed_hotspot_seamounts]
    [if slab_tracers_enabled] update_slab_tracers, update_cmb_heat

    [at plot_interval / csv_save_interval steps, and always at the last step]
        plot_NNNN.png, tracers, seafloor age, boundaries, earthquakes, topography,
        topo_bathy, mantle_temp, plate_id, craton, divergence, rotations_NNNN.rot
```

### `compute_timestep` internal pipeline

```
tessellate_polygons           densify polygon edges to tessellation_deg (great circles)
plate_polygons_to_wrapped_gdf wrap at dateline, keeping interior rings (holes)
resolve_overlaps              clip overlapping polygons by PLATEID1 priority
update_plate_id_raster        rasterize to plate_id grid; fill NaN via BallTree haversine
[apply_plate_merges]          merge_pairs from the caller (merging_enabled only)
[craton raster]               rasterize cratons (with remap of absorbed plates), overwrite plate_id
split_disconnected_plates     disconnected pieces of one ID → new plates (or absorbed if tiny)
update_plate_polygons         polygonize raster → toposimplify → clip → densify edges to ≤ 1°
[detect_and_split_plates]     shape screening → spherical erosion neck search → cut
                              (fragmentation_enabled only)
detect_enclosed_plates        absorb plates enclosed by one larger plate
topo_polygons_to_boundaries   extract shared boundary segments via TopoJSON
compute_boundary_divergence   per-arc-segment relative normal velocity; sign from probing the
                              plate-ID raster on both sides of each segment
velocity_grid                 east/north velocity at each raster cell (ω × r)
spherical_divergence          divergence raster (1/Myr)
polygons_to_stage_poles       integrate boundary torques → stage rotation per plate
[apply_collision_coupling]    force_modes.collision: coupling only
```

### Key data structures

**`stage_pole_dict`**: `Dict[int, pygplates.FiniteRotation]`  
Maps plate ID → stage rotation for one `delta_time` step. The rotation angle is the net torque magnitude × `delta_time`, so the torque acts as an angular rate and plate speeds do not depend on the step size.

**`polarity_dict`**: `Dict[Tuple[int,int], int]`  
Maps `(min_pid, max_pid)` → downgoing plate ID at each convergent boundary pair. Uses craton-area heuristic: plate with larger craton pixel count is overriding.

**`craton_counts`**: `Dict[int, int]`  
Maps plate ID → craton raster pixel count. Used by `resolve_downgoing_plate`.

**`tracers`** / **`seamounts`**: `GeoDataFrame` with columns `geometry` (Point), `active` (bool), `age` (float, Myr). Always use `.active.astype(bool)` when masking.

**`slab_tracers`**: `DataFrame` with columns `lon`, `lat`, `depth_km`, `age` (Myr since subduction), `active` (bool), `sink_rate_km_per_myr`, `weight`.

**`cmb_heat`**: `ndarray` (nlat × nlon) of heat accumulation units. Built up over time by `update_cmb_heat`; injected into the deep mantle by `generate_mantle_temp_field`.

**`rasters`**: `xr.Dataset` with variables `plate_id`, `velocity_east`, `velocity_north`, `divergence`, and optionally `craton`. All on a regular lon/lat grid at `raster_spacing` degree resolution.

**`boundaries_gdf`**: `GeoDataFrame` where each row is a single 2-vertex arc segment with columns `plateid1_poly1`, `plateid1_poly2`, and `divergence` (km/Myr = mm/yr; positive = diverging).

### Boundary divergence computation

Each boundary segment is classified independently. `compute_boundary_divergence` iterates over every consecutive vertex pair in every boundary line:

1. The great-circle arc normal is computed via `GreatCircleArc.get_great_circle_normal()`
2. The sign is determined by probing the plate-ID raster on **both sides** of the segment midpoint, 0.5° along the great circle through the arc normal (lookups wrap across ±180°). Probes that hit pid1 on the positive side or pid2 on the negative side vote for sign = +1; the opposite hits vote for −1. This is robust for large, concave, or dateline-spanning plates, where a centroid hemisphere test can give the wrong sign. Segments whose side cannot be determined (e.g. at triple junctions) are kept with divergence 0.
3. Velocities at the midpoint are computed for both plates as ω × r, which is numerically identical to `pygplates.calculate_velocities`
4. The relative velocity is projected onto the arc normal to give the divergence value

In `polygons_to_stage_poles` the same relative velocity classifies each polygon segment:
- **Convergent:** slab pull is applied if the plate is the downgoing plate and slab suction if it is the overriding plate (polarity from `resolve_downgoing_plate`). Continent–continent segments get collision resistance instead.
- **Divergent:** ridge push.

Outward normals are derived per ring from its orientation, so clockwise and anticlockwise polygons and interior rings all get the correct sign.

This is spherically correct and has no polar singularity.

### Slab tracer physics

- **Seeding**: new slablets are placed every `sample_deg` (1°, great-circle) along subduction segments at depth 0 km. Each has weight `delta_time × min(1, |divergence| / 10 km/Myr)`
- **Sinking**: upper mantle (0–670 km): rate = `min(slab_sink_rate_factor × |divergence|, max_sink_rate)` km/Myr; lower mantle: rate / `slab_slowing_factor`. Slablets are removed on reaching the CMB (~2890 km)
- **Thermal fading**: `dT = T₀ × weight × exp(-age / slab_thermal_tau_myr)` where `T₀ = -500 K`
- **dVs conversion**: `dVs = dT × dVs_per_K` (-0.001 %/K) → cold slab = positive dVs (fast waves)
- **Final field**: spherical Gaussian horizontally, Gaussian in depth. Each depth level has zero area-weighted mean. It can optionally be blurred by `tomo_smooth_horiz_deg` / `tomo_smooth_depth_km`

### Synthetic earthquake generation

`generate_synthetic_earthquakes` draws uniformly distributed candidate points on the sphere and assigns each to the nearest active boundary (great-circle distance):

- **Ridge earthquakes**: within `eq_ridge_max_dist_deg` of a divergent segment; random shallow depth ≤ `eq_ridge_max_depth_km`; thinned by `eq_ridge_prob`
- **Subduction earthquakes**: only where the plate-ID raster shows the overriding plate; depth = `distance_from_trench × tan(slab_dip_deg)`; capped at `eq_slab_max_depth_km`

### Synthetic topography

`compute_topography_raster` builds elevation in layers (ocean pixels = 0):

1. **Background**: craton pixels → `topo_continental_base_m`
2. **Andean belt**: overriding-plate cratons near ocean-continent convergence → Gaussian profile peaked at `topo_andean_arc_dist_deg` from trench
3. **Himalayan belt**: convergent segments with continental crust on *both* sides nearby (probed up to `topo_himalayan_max_dist_deg` from that segment, not anywhere on the plate) → linear ramp from `topo_himalayan_peak_m` at the suture to 0 at `topo_himalayan_max_dist_deg`
4. **Terrain noise**: bandlimited spherical harmonic noise (`pyshtools`), power ∝ `(l+1)^{-β}`, amplitude scales with local elevation. The noise pattern is the same in every saved frame of a run (seeded from `random_seed`)
5. **Coastal taper**: elevations fall smoothly to 0 within 2° (great-circle distance) of the coast

### Seafloor bathymetry

`generate_topo_bathy_grid` applies the GDH1 age-depth relation:
- Age < 20 Myr: `depth = -(2600 + 365 × √age)`
- Age ≥ 20 Myr: `depth = -(5651 − 2473 × exp(−0.0278 × age))`

Ocean pixels with no tracer coverage receive a default depth of −3500 m. Continental pixels overwrite the bathymetry.

### Subduction polarity determination

`resolve_downgoing_plate(p1, p2, polarity_dict, craton_counts)`:
1. One plate has craton pixels, other does not → continental plate is overriding
2. Both have craton pixels → larger craton pixel count is overriding
3. Neither has craton pixels → fall back to `polarity_dict` sticky assignment, set
   when the pair first converges (lower plate ID).

With `polarity: age`, steps 1–2 are replaced by the contact rule above and the
`polarity_dict` entry is final: callers pass `sdp.polarity_counts(...)`, which is empty
in that mode.

### Rotation file format

GPlates `.rot` format:
```
plate_id   time   euler_lat   euler_lon   angle_deg   000
```
Fixed plate is `000` (absolute frame). Each plate gets a `t=0` identity line. Sorted by plate ID then time.

The `time` column is **elapsed model time** in Myr since the starting snapshot, not "Ma before present". The finite rotation at time `t` moves the step-0 plate geometry to its position at step `t / delta_time`. To use the file in GPlates, treat the starting configuration as "present day" (0 Ma) and the simulated future as increasing "age". Plates created by fragmentation carry their parent's rotations up to the split.

---

## Known issues and TODO

- **Overlap priority**: where moved plates overlap, the lower PLATEID wins (`resolve_overlaps`). A physically based rule, for example letting the overriding plate win at trenches, is not implemented.
- **Boundaries on the antimeridian**: a boundary that runs exactly along ±180° is not extracted, because the polygons are split at the dateline and the two sides are not recognised as a shared edge.
- **Terrain noise is fixed to the grid**: continents move over a stationary noise pattern rather than carrying their texture with them.
