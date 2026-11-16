# Public release plan: triage of the framing problems

Draft prepared with Claude Code, 2026-11-16. Status: **in progress** (below). Revised the same day: slablet can be
made public too (Simon). That replaces decision 8 of the cross-repo register (slablet private, and optional in
PlateDriver).

## Progress

**A1, done (2026-11-16):**
- **Code:** `overlap: plate_id | polarity` (`f114119`, `0bfb67a`). The first version used one drawing order and
  failed on real cycles of polarity (777 beneath 666 beneath 555 beneath 777). The version kept settles each overlap
  by its own pair.
- **Gate:** `show-me/2026-11-16-overlap_polarity/`. Today's code with `plate_id` reproduces the 5 Oct runs bit for
  bit.
- **Largest plate:** with `polarity` it peaks at 32–34 % against 35–37 %, and breaks up sooner.
- **Old ocean:** not reduced. Ocean older than 200 Myr is 9–21 % at 500 Myr, against 9–30 %. So the guess below,
  that the overlap rule also causes the old ocean, is not supported, and v3's case stands.
- **Defaults:** `physical` and `simple` still use `plate_id` until Simon decides (decision 1).

**A2, done:** the `physical` preset now rebuilds plates in tiles, so the dateline seam is no longer a boundary
(`c15aeb6`).

**A3, done:**
- **Gate:** `show-me/2026-11-16-no_neck_split/`.
- **How often:** neck splitting fires rarely, 4–9 times per 500 Myr.
- **Without it:** the largest plate stays larger for longer in two of three seeds, peaking at 36–38 % against
  32–34 %. Plate counts and old ocean are alike.
- **Decision:** the rift test does not replace neck splitting, so it stays on, as a shape rule.

**Decision 1, tried and assessed (Simon: "let's try your recommendation"):**
- **`physical`:** `overlap: polarity` (`f150b35`), as in A1.
- **`simple`:** tried with `polarity` (`show-me/2026-11-16-simple_overlap/`) and returned to `plate_id`.
  - `simple` has no drag and no speed calibration (RMS 140 mm/yr at the start). In 500-Myr runs it ended as one plate
    in 1 of 3 seeds with `plate_id` and 2 of 3 with `polarity`.
  - With `polarity` its speeds ran away, up to 490 mm/yr.
  - A one-plate world hung the run; runs now end at that step with a message.
  - The README says `simple` is for runs of up to about 200 Myr.

**C, done (C1–C4):**
- rotation files count time back from the saved step (`42a14ed`);
- the README states the purpose, which plate-map changes are force-driven and which are rules, and how Earth values
  are used (`ad800a0`);
- the names users see say "continents", not "cratons", and old names are still read (`407a908`);
- the v1 plan's claim about flat geometry is corrected (`c31d640`);
- four unused modes are removed (`1c63261`), and old configs that name the remaining behaviour still run;
- the name is PlateDriver (`103d7cd`). The GitHub repository is renamed at the release, after asking Simon.

**B2, done:** the Atlantis starting map is in `examples/atlantis/`, which Simon agreed to publish (`48e8da7`). Input
paths are relative to the PlateDriver folder.

**B3, done:** smoke tests and a test of the drag tensor on a spherical cap (`2b412ba`, `0bfb67a`).

**S, done except S5 (the tag):**
- **Licence and tests:** slablet has an MIT licence (`16dcc83`), and its contract test passes on a fresh clone
  (`2ce084f`).
- **S1:** the S40RTS stage-0 work is committed at Simon's request (slablet `3d98c18`, `5479451`; PlateDriver
  `207d67c`). slablet's `main` is fast-forwarded to its working branch locally; nothing is pushed.
- **S3:** config paths use `$SLABLET_DATA` and `$GPRM_CACHE` (`76dabf5`). The configs resolve to the same files as
  before, and 151 tests pass.

**S4, done:** nothing is held back (decision 7), and the status page is brought up to date (slablet `3049e15`).

**Tests, night of 16-17 November** (Simon: "extensive and careful tests of the existing functionality"):
- **Fresh install.** `environment.yml` builds an environment from conda-forge and GitHub (pandas 3.1, numpy 2.5,
  OpenBLAS 0.3.34). It found two faults no test had: breakoff wrote into a read-only pandas array (`0a3a062`), and the
  drag solve's SVD failed on one OpenBLAS (`a796f16`, now falls back to QR). slablet's own tests pass there after a
  test-helper fix (slablet `b11e8f4`). The viewer installs, builds and passes its three browser checks from a clean
  clone, on runs exported from the clone.
- **Mode matrix** (`show-me/2026-11-16-mode_matrix/`): every allowed mode value runs in both environments, and the two
  environments agree to 1e-12. Fixed on the way: output folders named with `=` or `+` lost their age grids
  (`96abdcb`, `93873fa`); `rifting: trial` silently did nothing in `simple` and now stops with a message (`e2b272f`);
  a null `random_seed` is now recorded so the run can be repeated (`a6ca02c`).
- **Seed ensemble** (`show-me/2026-11-16-seeds_current/`): eight more seeds of the current model to 500 Myr; none
  failed, every invariant holds, and the new rotation files reconstruct each run.
- **Reproduction:** today's code reproduces the 16 November polarity run bit for bit over 100 Myr, from its old config
  (`cratons` and the `craton_*` names). slablet's mantle field on Muller 2022 is bit-identical to a run from `fc4d706`.
- **1000 Myr** (`show-me/2026-11-17-long_current/`): the current model on four seeds to its full length. All finish,
  and the first 500 Myr of seeds 0-2 are bit-identical to the earlier runs. But in three of four seeds the
  continents shatter from about 400 Myr (35-51 plates over 750-1000 Myr, 13-23 of them slivers): merging holds
  the small-plate rifting for 500 Myr, not for 1000. This belongs with v3's questions about plate size.
- **Numerical sensitivity** (`show-me/2026-11-17-numerical_sensitivity/`): the force grid is converged at 1 deg. The time
  step and the plate grid are not: a finer step or plate grid leaves more starting ocean unsubducted, so part of the
  defaults' consumption of old ocean is numerical. Two gate-tool faults with time steps other than 1 Myr are fixed
  (`741e72b`), and the least-cost rift cut now scales with the force grid (`a815951`).

**Still open:**
- **`show-me/` references.** The README cites six `show-me/` folders as evidence, and `show-me/` is not in the
  repository. Before the release, either publish the folders it cites (for example under `docs/results/`) or drop
  the references.
- **B1:** PlateDriver requires a tagged slablet, once slablet is public.
- **Release:** tag slablet; merge `mantle-extract` into PlateDriver's `master`; rename the GitHub repository. Each
  outward step waits for Simon.

## Why

PlateDriver is public on GitHub as OldDriver, but its default branch is still at force model v1, and the work since
(57 commits on `mantle-extract`) is to be released properly soon. A review on 16 November found nine places where the
code's description of itself disagrees with what it does. Here we triage them against the aims Simon set for the
release:

- PlateDriver stays a teaching tool. Research use is welcome if the code stays simple, cheap to run and readable.
- How to make a synthetic Earth driven by plate forces has high priority.
- Synthetic and real plate histories should work side by side.
- slablet can be made public as well.

Each problem was checked against the code and, where possible, measured on the saved current-model runs
(`runs/current_model/tb0`–`tb2`, 500 Myr each).

## What the checks found

### At many trenches the wrong plate is consumed

Where moved plates overlap, the plate with the lower ID keeps the area: the plates are drawn into the plate-ID raster
highest ID first (`compute_timestep`, `self_driving_plates.py:2638`). Continents are drawn over the plates afterwards,
so this decides the outcome wherever the overlap is ocean. Which plate subducts is decided elsewhere, by the force
model's polarity rule (older ocean, continents override). The two agree only when the downgoing plate happens to have
the higher ID.

We measured the share of subduction flux (convergence × trench length, trenches with no continent within 1°) where
the downgoing plate keeps the overlap, using the `trench_forces_NNNN.csv` files. In the three current-model runs:

- **Overall:** it is 32–41 %.
- **Through time:** it rises from 11 % in the first 100 Myr to 62 % by 500 Myr. The old plates carry the old ocean
  that subducts, and they also have the low IDs.
- **In `simple`:** for pairs the continent rule does not settle, the lower ID is also the plate that subducts. So the
  rule is wrong at every such trench.

When the downgoing plate keeps the overlap, the trench travels with that plate. Its seafloor never reaches the trench,
and the overriding plate's edge is consumed instead.

Plate 111, the lowest ID, shows the consequence.

- **Size:** in all three seeds it is the largest plate at 500 Myr, at 22–34 % of the sphere, grown from 6 %.
- **Trenches:** from 200 Myr it has 5,000–25,000 km of trench where it subducts, and it keeps the overlap along all
  of it.
- **Net effect:** it gains area at its ridges and loses none at its trenches.

This is a plausible cause of the two symptoms that motivated force model v3: a giant ocean plate (37–49 % at
1000 Myr) and much ocean older than 200 Myr. Those may come from this rule rather than from the lack of subduction
initiation. The scripts for these numbers are not in the repository yet; stage A puts the measure into the gate
table. (The A1 gate supports this for the giant plate but not for the old ocean; see Progress.)

### Other measurements

| Item | Finding |
|---|---|
| Ridge push, `body` against `boundary` torque (1,887 plate-steps) | Median angle 14°. 19 % of the torque (weighted) is more than 30° apart. Median magnitude ratio 0.61. The two modes differ by design where flanks are incomplete, so this bounds the lon/lat-gradient error rather than measuring it. |
| Step time, current model | Median 3.7 s, about 1 h for 1000 Myr. Of this, 3.1 s is in `compute_timestep` (rebuild, forces, rift test), 0.17 s in seafloor tracers and 0.04 s in slab tracers. The mantle is cheap; the rebuild and the rift test are the cost. |
| What the default forces need from slablet | `MantleParams`, `initial_slab_columns`, `seed_slabs`, `sink_slabs`, `slab_depth_near`, all in slablet's core. |
| PlateDriver's input paths | `config.yaml` and `configs/current_model.yaml` point to `/Users/simon/GIT/octopus/Atlantis3` (616 kB). octopus is public and MIT-licensed, but `Atlantis3/` is untracked there, so the starting map is not public yet. Four archived notebooks have absolute paths too. |
| PlateDriver's tests | Three files: the slablet adapter, the sphere helpers and `plate_step`. None runs the model. |
| slablet's state | Private, with no licence. `main` is 28 commits behind `asthenosphere-layer`. 33 files hold another session's uncommitted work (stage 0 of the S40RTS plan). 40 tracked files name local paths, mostly tomography models under `~/Data/SeismicTomography`. The scoring tools need these models, and not all of them can be redistributed. Its docs discuss IPV and defamation, which stay private. |

## Triage

P0 must be done before the release; P1 comes soon after; "later" is research.

| # | Problem | What it does | Priority | Stage |
|---|---|---|---|---|
| 5 | Overlaps go to the lower plate ID, not the overriding plate | The wrong plate is consumed in 32–41 % of ocean subduction (62 % by 500 Myr); a likely cause of the giant plate and old ocean | P0, correctness | A |
| 7 | The `physical` preset counts the dateline seam as boundary in edge drag | The default force balance is wrong by grid geometry (speed scale about 5 %) | P0 | A |
| 2 | The default `physical` preset needs slablet, which is private | A public clone cannot run its own default or the teaching model | P0, blocker | S, B |
| 13 | slablet is not ready to be public (licence, paths, uncommitted work, docs) | It cannot be published as it is | P0, blocker | S |
| 10 | Hard-coded input paths in PlateDriver | A public clone cannot find the starting map | P0, blocker | B |
| 12 | No test runs the model | Nothing catches a broken release | P0 | B |
| 9 | Rotation files put elapsed time in GPlates' age column | Run histories play backwards in GPlates, and against the age grids | P0 | C |
| 8 | "Craton" means all continent, strong cores, or deep keels | Misleading names. Renaming after the release breaks users | P0 | C |
| 1 | The README says teaching only; the repo is a research hub | Confuses readers; no effect on results | P0, docs | C |
| 3 | "Earth is not a target", yet the clock is set by Earth's plate speed | Framing; the calibration is hidden | P0 docs, P1 physics | C, E |
| 4 | Forces set speeds, but rules set the plate layout | Framing; shape-based splitting duplicates the stress-based rift test | P0 docs, P1 test | C, A |
| 11 | A 4,200-line library and a 1,000-line script with module-level state | Hard to read; the synthetic and real uses cannot share the force code | P1 | D |
| 9b | The contract's time sign is opposite to `plate_step`'s | A trap for any shared physics that turns two times into a duration | P1 | D |
| 6 | Slabs are vertical for the mantle but dip 45° for earthquakes; continents are polygons but rebuilt from the raster | Earthquakes and synthetic tomography disagree; hidden continent is lost at rifts | P1 | F |
| 7b | Ridge push gradient taken along lon/lat | Ridge push is biased towards N–S and E–W | P1 | F |

## Organising synthetic and real side by side

The two uses differ in one thing: where plate motions come from. In a synthetic world the force model computes them
and the plates are moved. In a real reconstruction the rotation file gives them. Everything else is the same: the
plate state, the boundaries, seafloor age, the datasets and the mantle.

With slablet public, the plate state already exists: slablet's plate-history contract (`slablet/plate_history.py`,
`PlateStep`). PlateDriver emits it (`plate_step`), and slablet's real driver (`RealPlateHistory`) emits it for any
GPlates reconstruction. We therefore propose that PlateDriver's force model and its datasets read a `PlateStep`,
whichever driver made it:

- On a synthetic step, the predicted motions move the plates.
- On a real step, they are compared with the reconstruction's own motions. This tests the forces and can calibrate
  their weights against Earth's plate motions.

The second use also resolves problem 3. Earth's plate motions may calibrate the physics, while Earth's present-day
outcomes (ocean ages, plate sizes) remain a reference.

The split between the two repositories:

- **slablet (public):**
  - the contract;
  - the real-reconstruction driver;
  - the mantle model, dynamic topography and flow;
  - the asthenosphere and volcanism models;
  - the observation and scoring tools.
- **PlateDriver (public):**
  - the synthetic driver: one forward step (move, rebuild, resolve overlaps), and the rules that change the plate
    layout;
  - the force model, which reads a `PlateStep`;
  - the teaching datasets: earthquakes, topography, bathymetry, seamounts and rotation files;
  - two entry points, `make_synthetic_earth.py` (with `run_SDP.py` kept as an alias) and
    `forces_on_reconstruction.py`.

PlateDriver depends on a tagged slablet release. slablet never imports PlateDriver, as its import test already
enforces.

**What the contract lacks for the forces.** It needs three extensions to the register's decision 4, which go into the
register when they are built:

- plate IDs on trench samples (downgoing and overriding) and on ridge samples (either side);
- all boundary segments, transforms included, with their plate pair, length and normal, because edge drag and
  collision contacts need them;
- one stated time direction for both drivers (problem 9b).

**Like for like.** The force model gets one input route, `force_inputs(step)`. It turns any `PlateStep` into what the
model uses: the raster, the boundary segments, polarity, ages and continents. Before the real driver uses this route,
a test runs it on saved synthetic steps and shows it gives the same torques as the current internal route, to a
stated tolerance. The internal route can then be retired in a gated stage.

This replaces the earlier proposal, in which PlateDriver owned its own state and reconstruction reader and slablet
stayed private behind an adapter. Several things go:

- the second reconstruction reader;
- the sphere helpers copied into PlateDriver (register decision 7);
- the optional-slablet machinery (`require_slablet`).

The cost is that PlateDriver must be shielded from slablet's research churn by pinning a release.

**Spherical geometry:**

- **Plate IDs:** real plate IDs come from the contract's exact point queries on a cell-centred grid (`on_grid`), not
  from planar rasterising.
- **Comparing motions:** predicted and actual motions are compared as 3-D velocity vectors, weighted by cell area
  (cos lat).
- **Frame:** the force balance assumes a mantle at rest, so a reconstruction must be in a mantle frame. Spin-axis
  frames are refused, as slablet already does.

## Stages

Every stage starts plain and ends with a gate.

- **Physics stages:** three tie seeds × 500 Myr of the current model (about 35 min with three workers), led by the
  map sequence, with the gate table and figures in `show-me/`.
- **Refactoring stages:** checked with `regression_check.py --tol 0`.

### A. Overlaps by polarity (problems 5, 7, 4)

1. **Draw the overriding plate last.** Plates go into the raster downgoing first and overriding last, using each
   pair's polarity from the previous step: an ordering of the "subducts beneath" pairs, with ties and cycles broken by
   ID. This becomes a new mode, `overlap: plate_id | polarity`.
   - **Gate:** the current model with `overlap: polarity`.
   - **Measures, added to the gate table:**
     - the share of subduction flux where the downgoing plate keeps the overlap (0 at ocean trenches is the check);
     - the largest plate;
     - the ocean-age distribution;
     - the plate count.
2. **Tiles on by default.** Set `rebuild_tiles: on` in the `physical` preset. The current model has it already.
3. **Without shape splitting.** Run the current model without neck fragmentation, which the stress-based rift test
   may make unnecessary.

Then re-read v3's motivation against the A1 gate. If the giant plate and old ocean are gone or much reduced, v3's
subduction initiation is no longer the priority; otherwise v3 resumes at S1.

### S. slablet goes public (problems 2, 13), in parallel with A

1. **Committed and merged.** The S40RTS stage-0 work (33 files) is committed or set aside by the session that owns
   it, and `asthenosphere-layer` is merged into `main`.
2. **Licence.** MIT, as for PlateDriver and octopus.
3. **Data paths.**
   - Local data sit under `$SLABLET_DATA` (which `cooling.py` already reads), with a documented layout, and configs
     and analysis scripts read it instead of `/Users/simon/Data`.
   - The model and the real driver run with no local data beyond what gprm downloads.
   - Each scoring tool says which data it needs and where to get it.
4. **Docs.** All of `docs/plans/` and `docs/status.md` go public (decision 7). The AAD paper draft
   stays out; it is local to PlateDriver's `show-me/`.
5. **Clean clone.** Tests pass from a clean clone, and a release is tagged (for example v0.1).

### B. A public clone runs (problems 2, 10, 12)

1. **slablet required.** PlateDriver requires slablet at the tagged release, installed with pip from GitHub. Out go
   `require_slablet`, the optional imports, and the copied sphere helpers with their test. The gate tools' check
   against an uncommitted slablet stays, for development.
2. **Bundled starting map.** `examples/atlantis/` gets the Atlantis starting map, if Simon agrees to publish it,
   and the configs use relative paths. The archived notebooks get relative paths, or the folder is dropped.
3. **Smoke tests.** Each runs 20 steps of `simple`, and 20 of `physical`, from the bundled map. One more checks the
   torque on a spherical-cap plate against its analytic value.

### C. Conventions and framing before the release (problems 1, 3, 4, 8, 9)

1. **Rotation time counted back.** Rotation files count time back from the saved step, which becomes "present day"
   (0 Ma). GPlates then plays a run forward, and the files agree with that step's age grid. `rotation_time: elapsed`
   keeps the old form.
2. **Names:**
   - The input `cratons` becomes `continents`, and the old key is accepted with a warning.
   - The craton raster becomes the continent raster, with the viewer updated.
   - "Craton cores" become "cores". These could later serve as the synthetic cratons that the cross-repo register
     defers.
3. **README:**
   - **What this is:** teaching first, with research use through the same code and slablet.
   - **Forces and rules:** which motions are force-driven and which layout changes are rules, with each rule's
     trigger.
   - **Earth values, restated:** Earth's present-day outcomes are a reference. Earth's plate motions may calibrate the
     physics. The speed scale is a calibration to Earth's mean plate speed, and the README says so.
   - **Corrections:** known issues updated, and the v1 plan's claim that nothing uses a flat approximation corrected.
4. **Pruning.** Modes that are in no preset, and are not kept for a stated reason, are removed (decision 5).

**Release:**
- tag slablet;
- merge `mantle-extract` into PlateDriver's `master`, pinned to that tag;
- tag PlateDriver and update its GitHub description.

### D. Code layout around the contract (problems 11, 9b)

1. **Moves only.** Modules are moved one at a time, each bit-identical for the current model and for `simple`, in
   this order: forces, seafloor, datasets, events, the synthetic step.
2. **The driver loop.** It becomes a function with explicit state instead of a module-level script. `run_SDP.py
   --config` keeps working.
3. **Contract extensions,** built in slablet with tests:
   - plate IDs on trench and ridge samples;
   - the boundary segments;
   - one time direction. `plate_step` labels its step time = −(t + Δt), so both drivers run towards smaller times.

   PlateDriver's `plate_step` fills all of these.
4. **One input route.** `force_inputs(step)`, with the test that it matches the internal route on saved synthetic
   steps.

### E. Forces on a real reconstruction (problem 3; the priority use)

1. **Reader.** `RealPlateHistory` supplies the extended contract for Müller et al. (2022). It needs its seafloor
   tracker cache, which takes about 6 min to build once. Slab depth is first assumed full at every trench, as at a
   synthetic run's step 0.
2. **Evaluate, without fitting.** Compute each force's torque and the predicted plate motions with the current
   weights. Compare them with the reconstruction's motions at 0 Ma, then every 10 Myr to 100 Ma. Lead with maps of
   predicted and actual velocity; per-plate angle and speed ratio follow.
3. **Fit the weights.**
   - **Method:** predicted rotation is linear in the force weights for a fixed drag shape, so this is a least-squares
     fit. It is done per time and over 0–100 Ma, with a grid search over the continent drag factor and the edge
     width.
   - **Report:** which weights Earth's plate motions constrain.
   - **Use:** a preset option takes the fitted weights. `target_rms_speed` stays as the simple calibration.
4. **Datasets for Earth.** The same datasets (earthquakes, topography, ages) are made from real steps, to set beside a
   synthetic planet's in teaching.

### F. Consistent slabs and forces (problems 6, 7b)

1. **Dipping slabs.** Upper-mantle slabs dip towards the overriding plate. slablet already has dipping slabs, and
   PlateDriver supplies the down-dip azimuth from polarity. Earthquakes come from the slab itself, so seismicity and
   synthetic tomography agree.
2. **Ridge push on the sphere.** It comes from a gradient taken on the sphere, per plate and smoothed. The direction
   error is measured before and after.
3. **Continents cut at rifts.** Continents are cut as polygons at rifts, so no hidden continent is lost.

### Later (research, each with its own plan)

- Synthetic forces read their own `PlateStep`, and the internal route is retired (gated, not bit-identical).
- v3's subduction initiation, if stage A leaves it needed.
- Forward runs from a real snapshot (for example 100 Ma), compared with the reconstruction over the first 20–50 Myr,
  before seeds diverge.
- A generator of starting maps (random continents and plates), so a new planet needs no hand-drawn map.
- Per-segment polarity.

## Decisions for Simon (recommendation first)

1. **Overlap rule in `simple`:**
   - **Recommended:** switch `simple` to `overlap: polarity` too, and renew its regression reference once, keeping
     `plate_id` for old runs.
   - **Alternative:** keep `simple` exactly as the original model.
2. **slablet in PlateDriver:**
   - **Recommended:** a required dependency, pinned to a tagged release.
   - **Alternative:** optional, as now, with `simple` runnable without it.
3. **Where the plate state lives:**
   - **Recommended:** slablet's contract, extended, read by PlateDriver's forces and datasets from either driver.
   - **Alternative:** PlateDriver's own state and reader, the earlier proposal.
4. **Order:**
   - **Recommended:** A and S in parallel, then B and C; release slablet, then PlateDriver; then D, E and F. The
     release then fixes wrong results and blockers with small changes, and the large refactor runs on a tested base.
   - **Alternative:** refactor before the release.
5. **Modes to remove.** Candidates:
   - `plate_rebuild: rotated`, kept only for comparison;
   - `collision_weld_scope: plate`, which collapses the world;
   - `collision_tear`, which never passes;
   - `slab_pull: ramp`.
6. **Name:** PlateDriver (Simon, 16 Nov). The GitHub repository OldDriver is renamed at the release.
7. **slablet's docs in public:** everything is published as it is, with the full git history (Simon, 16 Nov). The
   recommendation had been to keep the cross-repo register and the status page private until reviewed.

## Budgets

- **Step time:** `simple` and the current model get no slower than now (median 3.7 s per step for the current model).
  A new feature states its cost.
- **Rules:** every rule that changes the plate layout is listed in the README with its trigger.
- **Modes:** new physics is a mode, off until its gate passes. This is unchanged.
- **slablet changes:** they reach PlateDriver only through a new tagged release.
