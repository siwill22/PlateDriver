# Force model v1 — plan

Draft prepared with Claude Code. Decisions come from a design interview (Q1–Q18); the
question number is given in brackets.

## Goal

Plausible synthetic worlds first; speed matters more than strict realism (Q1). Forces
remain dimensionless, tunable weights, but take the physically correct *form*.
Every change is a mode flag. A `simple` preset reproduces today's model exactly: it
is the fastest option and useful for teaching.

## Decisions

| # | Topic | `physical` (default) | Also available | `simple` = today |
|---|---|---|---|---|
| Q2 | Plate mobility | Drag tensor: ω = k · I⁻¹ τ | — | ω = τ (no drag) |
| Q3 | Continental keels | w = `continent_drag_factor` (3) on craton cells, 1 on ocean | optional `craton_drag_factor` | — |
| Q4 | Seafloor age | Existing age tracers, queried every step | — | age used for output only |
| Q5 | Slab pull magnitude | √min(age, 80)/√80 × min(1, L_slab/660 km) | age only | uniform |
| Q6 | Ridge push | Body force over ocean < 80 Myr | `boundary` | `boundary` |
| Q7 | Slab suction | Deep slab (> 660 km) → suction on *both* plates | `ratio` | overriding plate only, same weight as pull |
| Q8 | Pull on/off | Pull ∝ slab present (L_slab); convergence only seeds slab | `ramp`, `hard` | `hard` |
| Q9 | Polarity | Older ocean subducts; sticky per pair; continent rule first | `id` | lower plate ID |
| Q10 | Net rotation | Left free; logged per step | — | not logged |
| Q11 | Speed scale | `calibrate_once` at step 0 to `target_rms_speed`, then frozen; cap at 3× | `per_step`, `none`; `speed_scale` override | `none` |
| Q12 | Collisions | `breakoff`: slab detaches after `breakoff_delay_myr`, then coupling | `viscous` (replaces `resistance`), `coupling` | `resistance` |
| — | Transform resistance | Not in v1 (Forsyth & Uyeda: small) | — | — |
| Q13 | Physics grid | Separate `physics_grid_deg: 1.0`, built lazily | — | not built |

Rule (Q13): a field (physics grid, age raster, slab-length query) is built only if a
selected mode reads it. The `simple` preset must cost exactly what a step costs today.

## Config schema (Q14)

```yaml
force_model: physical            # preset: simple | physical
force_modes:                     # optional; any key here overrides the preset
  drag: tensor                   # none | tensor
  ridge_push: body               # boundary | body
  slab_pull: slab_length         # hard | ramp | slab_length
  slab_suction: depth_split      # ratio | depth_split
  polarity: age                  # id | age
  collision: breakoff            # resistance | viscous | coupling | breakoff
  speed_normalisation: calibrate_once   # none | calibrate_once | per_step
params:
  continent_drag_factor: 3.0
  craton_drag_factor: null
  target_rms_speed: 50           # km/Myr (= mm/yr)
  max_rms_speed_factor: 3.0      # cap = factor × target; each capped step is logged
  speed_scale: null              # a number fixes k and skips calibration (compare runs)
  slab_ref_depth_km: 660
  slab_age_saturation_myr: 80
  ridge_push_age_cutoff_myr: 80
  pull_ramp_km_per_myr: 10       # ramp mode only
  breakoff_delay_myr: 10
  physics_grid_deg: 1.0
```

When the run starts, the log prints the resolved modes and the k it used.

## Formulation notes

**Spherical geometry.** Everything here is on a sphere, and nothing uses a flat
approximation:

- **Tensors:** for each plate p, the drag tensor is I_p = Σ w (𝟙 − r rᵀ) dA over its
  cells; J_p is the same sum with w = 1. r is the unit position vector and dA uses
  cos(lat) weights on the physics grid.
- **Drag:** a rigid plate rotating at ω against a viscous mantle feels torque −c I ω,
  so the balance is ω = k I⁻¹ τ.
- **RMS speed:** found exactly from the tensors, with no velocity grid:
  RMS² = R² Σ_p ω_pᵀ J_p ω_p / A_total.
- **Net rotation:** ω_nr = (Σ_p J_p)⁻¹ Σ_p J_p ω_p. Log it to `net_rotation.csv`, with
  one row per step giving the pole, the rate and the RMS speed.
- **Frame:** drag acts relative to the mantle, so the frame is fixed to the mantle. That
  is the same frame as the hotspots, and net rotation is a real output (Q10).

**Ridge push as a body force (Q6).** The force per unit area is f = −C ∇age for ocean
cells younger than the age cutoff, and 0 elsewhere. The plate's torque is Σ r × f dA.

- Plan to add a gradient only after the age raster has been smoothed spherically; the
  nearest-tracer raster is blocky.
- The gradient is wrap-aware, with ∂/∂x = (1/(R cos φ)) ∂/∂λ.

*Correction to what I said in Q6:* the magnitude is ∝ |∇age|, not constant per unit
area. The force F(t) grows linearly with age, so the force per unit area is dF/dx. Its
integral from the ridge out to 80 Myr is therefore C·80 per unit length of ridge,
whatever the spreading rate. That matches Parsons & Richter (1980). A constant force
per area would give fast-spreading plates too much push.

**Age raster (Q4, Q13).** Assign each physics-grid cell the age of the nearest active
seafloor tracer, using one BallTree per step. Craton cells, and cells with no tracer
within 2 cells, get NaN.

**Slab length (Q5, Q7, Q8).**

- For each trench segment, take the deepest slab tracer within
  `boundary_search_radius_deg` of the segment, on its overriding side. Call its depth z.
  One BallTree on the slab tracers serves all segments.
- Pull = w_sp · f_age · min(1, z/660), where f_age is the age term from Q5.
- Suction = w_ss · min(1, max(0, z − 660)/660). It acts on both plates, toward the
  trench.
- The final definition will be fixed in S3, with a figure.

**Trench age (Q5).** This is the age of the nearest seafloor tracer on the downgoing
side, probed at an offset along the segment normal.

**Polarity (Q9).** Keep the existing continent rule in `resolve_downgoing_plate` as the
first test. For ocean–ocean pairs, at first convergence the plate whose ocean is older
along the contact goes down. The result is stored in `polarity_dict`, and its format
stays the same.

**Collisions (Q12).**

- `viscous`: push-back ∝ w_col · max(0, −div) / v_ref, which can stop plates but never
  reverse them.
- `breakoff`: once `collision_state` age exceeds `breakoff_delay_myr`, deactivate the
  slab tracers attached to the pair so the pull goes to zero, then apply coupling.

## Speed work (Q15)

1. Store tracers as numpy arrays instead of GeoDataFrames, in `advect_tracers`,
   `deactivate_subducted_tracers` and `seed_ridge_tracers`. They currently take about
   0.3–0.4 s per step, and the target is under 0.05 s (an estimate). The CSV output
   stays the same.
2. Share lookups within each step: one tree on the seafloor tracers and one on the slab
   tracers, reused by the age raster, trench ages and slab length.
3. After S5, spend a fixed amount of time profiling `resolve_overlaps`, `wrap_polygons`,
   `plate_id_raster` and `fragment`, at 0.1–0.3 s each.
4. No GPU work in v1: the work is many small geometry operations, and PyTorch in pygmt17
   is broken.

## Validation (Q16)

**Test case:** atlantis66 (the current `config.yaml`), 100 steps, dt = 1, fixed
`random_seed`.

`diagnose_run.py <output_dir>` makes PyGMT figures:

1. **Forsyth & Uyeda plot:** plate speed against trench fraction (the trend should be
   positive) and against continental fraction (the trend should be negative).
2. **Time series:** RMS speed, net rotation, maximum plate speed and step time.
3. **Map snapshots:** velocity arrows over the age grid, with trenches coloured by slab
   length.
4. **Regression test:** `force_model: simple` must reproduce the `review-fixes` outputs
   exactly (rotations, plate IDs, tracers) over 10 steps.

Each gate shows the new mode side by side with `simple`.

## Stages (Q17)

Each stage ends at a review gate, and nothing moves to the next stage without approval.

| Stage | Contents | Gate |
|---|---|---|
| S0 | `diagnose_run.py`; preset and mode plumbing; mode log; save regression baseline; clean up warnings (`SettingWithCopyWarning`, the concat `FutureWarning` at run_SDP.py:603) | Figures for the current model; regression passes |
| S1 | Drag tensor with continent factor; speed calibration, cap and `speed_scale`; `net_rotation.csv` | F&U plot and time series compared with S0 |
| S2 | Numpy tracers; shared lookups; lazy physics grid and age raster | Step time; age grid unchanged; regression passes |
| S3 | Slab-length pull, depth-split suction, age polarity | F&U plot; slab-length maps |
| S4 | Ridge-push body force | Maps and time series |
| S5 | `breakoff` and `viscous` collisions | Collision sequence maps |
| S6 | Profiling pass | Step-time table |
| S7 | Update README and config documentation, with the sources below | Read-through |

## Branches (Q18)

- `review-fixes` holds the earlier review fixes and speed-ups, committed as 31832cb.
- `force-model-v1` branches from it, and this work goes there.
- Nothing is merged to `master` without approval.

## Known input issue

`config.yaml` pairs `atlantis3_rotations.rot`, which defines plates 101–901, with the
atlantis66 topologies and cratons, which use plates 111–777. This does no harm at
`reconstruction_time: 0`, because every rotation is then the identity. At any start
time above 0, those plates would silently get no rotation.

## Sources

- **Textbooks:**
  - Turcotte & Schubert, *Geodynamics* (3rd ed., 2014)
  - Kearey, Klepeis & Vine, *Global Tectonics* (3rd ed., 2009)
- **Drag and keels:** Solomon & Sleep (1974); Forsyth & Uyeda (1975); Chapple & Tullis
  (1977); Zhong (2001).
- **Slab pull and suction:** Carlson, Hilde & Uyeda (1983); Schellart (2004); Conrad &
  Lithgow-Bertelloni (2002, 2004).
- **Ridge push:** Lister (1975); Parsons & Richter (1980).
- **Coupled solutions:** Hager & O'Connell (1981); Ricard & Vigny (1989).
- **Net rotation:** Argus & Gordon (1991); Gripp & Gordon (2002).
