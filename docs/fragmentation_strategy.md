# Stress-Based Plate Fragmentation Strategy

> **Superseded** by the rifting trial split of `force_model_v2_plan.md` (2026-09). Kept for the
> record of the original neck-fragmentation idea, which still runs as `fragmentation_enabled`.

## Problem

Plates in the simulation can evolve into implausible shapes (U-shapes, thin isthmuses) because
the model has no intraplate rupture logic. In reality, such stress concentrations would cause the
plate to split along a predictable line. This document describes a strategy to detect and execute
physically motivated plate fragmentation.

---

## Why Plates Get U-Shaped

The current pipeline:
1. Rotates plates rigidly each step
2. Resolves overlaps via priority (lower PLATEID wins)
3. Re-polygonizes the raster

Because there is no rupture criterion, a plate that loses material asymmetrically (e.g. one arm
subducted faster than the other) can develop non-convex, geometrically implausible shapes.
The narrowest cross-section of such a plate concentrates all transmitted force but is never checked.

---

## Recommended Strategy: Medial-Axis Stress Proxy

### Core idea

At the "neck" of a plate, the full net boundary force is transmitted through a narrow cross-section.
Stress ∝ F_perp / width. When this ratio exceeds a threshold, split the plate there.

### Step 1 — Shape screening (cheap filter)

After re-polygonization each timestep, compute per plate:

- **Concavity ratio**: `area / convex_hull_area` — a fully convex plate = 1.0; U-shapes ≪ 1.0
- **Aspect ratio**: `bounding_box_long / bounding_box_short`

Only run the full stress calculation for plates below a concavity threshold (e.g. < 0.85) or
above an aspect ratio threshold (e.g. > 4). This keeps the common case cheap.

### Step 2 — Medial axis and neck detection

For candidate plates, rasterize the plate polygon at the existing grid spacing and apply
`skimage.morphology.skeletonize` to get the medial axis. At each skeleton point, the
**inscribed circle radius** = distance to nearest boundary pixel. The **neck** = skeleton
point with minimum inscribed radius.

The raster skeleton can be used to identify the neck region, then snapped to a great-circle
arc perpendicular to the local skeleton direction for vector accuracy on the sphere.

### Step 3 — Force component at the neck

The net torque per plate is already computed in `get_torque_for_segment()` in `plate_driver.py`.
Convert the net torque vector to a net surface force vector at the plate centroid. Project that
force onto the direction **perpendicular to the skeleton** at the neck point. This gives `F_perp`.

Stress proxy: `σ_proxy = |F_perp| / neck_width`

### Step 4 — Fragmentation criterion and split line

If `σ_proxy > σ_threshold` (a tunable config parameter, e.g. `fragmentation_threshold: 5e-3`):

1. Draw a **great-circle arc** through the neck point, perpendicular to the local skeleton
   direction (i.e. perpendicular to the long axis of the plate at that location).
2. Use Shapely `split()` (or raster flood-fill from the neck row/column) to divide the plate
   polygon into two children.
3. Assign the larger child the parent PLATEID; give the smaller child a new PLATEID (extending
   the `resolve_duplicate_plateids()` logic in `self_driving_plates.py`).
4. Initialize the child stage poles by inheriting the parent's stage pole — they will diverge
   over subsequent steps as their boundary forces differ.

### Step 5 — New boundary seeding

The fracture line becomes a **ridge** (divergent boundary) if the two fragments are being pulled
apart, or a **transform** if motion is parallel. The existing `compute_boundary_divergence()`
function handles this classification automatically once the two plates have separate stage poles.
No special-casing is needed.

---

## Integration Points in Existing Code

| What | Where |
|---|---|
| Per-step polygon list | `self_driving_plates.py` — after `resolve_duplicate_plateids()` |
| Net torque per plate | `plate_driver.py` → `get_torque_for_segment()` accumulation loop |
| Config parameters | `config.yaml` — add `fragmentation_threshold`, `min_concavity_ratio` |
| New PLATEID allocation | `resolve_duplicate_plateids()` in `self_driving_plates.py` |
| Raster skeleton | New helper using `skimage.morphology.skeletonize` on the plate ID raster |

The natural insertion point is **after `resolve_duplicate_plateids()` and before computing
stage poles for the next step** — so the new child plates get fresh torque calculations
immediately.

---

## Geological Plausibility Considerations

- **Rift zones follow pre-existing weaknesses** — the medial axis of a U-shape naturally
  points to the geometrically weakest location, which is a reasonable first approximation.
- **Continental vs oceanic lithosphere**: oceanic lithosphere is weaker in tension; the
  `cratons` raster already flags continental crust. A lower `σ_threshold` for oceanic plates
  would be geologically reasonable.
- **Minimum plate size**: prevent trivially small slivers from forming by enforcing a minimum
  child area (e.g. > 5° × 5° equivalent).
- **Rate limiting**: fragmentation should not happen every step — add a per-plate cooldown
  (e.g. no re-fragmentation within 5 Myr of last split) to avoid cascades.

---

## Verification

1. Run a simulation that historically produces a U-shaped plate; confirm fragmentation fires
   at the neck and the two children diverge to form a new ridge.
2. Check that stage poles for child plates are physically reasonable (close to parent's at t=0,
   then diverging naturally).
3. Verify the new boundary is classified correctly (ridge vs transform) by `compute_boundary_divergence()`.
4. Inspect output `boundaries_NNNN.geojson` to confirm the fracture line appears as a new
   plate boundary segment.
5. Confirm no duplicate PLATEIDs appear in the raster after the split step.
