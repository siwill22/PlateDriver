import type { Camera, PerspectiveCamera } from 'three';

import { VelocityField } from '../../vendor/petrify/js/index.js';
import { ThreeProjector } from '../core/boundaries';
import { FlatProjector } from '../core/flatProjector';
import type { ProjectionMode } from '../core/projection';

/** One plate at one saved step, from plates/NNNN.json (plate_stats.csv). */
export interface PlateInfo {
  id: number;
  pole_lat: number | null;
  pole_lon: number | null;
  rate_deg_per_myr: number | null;
  rms_speed_km_per_myr: number | null;
  area_frac: number | null;
  continental_frac: number | null;
  label_lon: number | null;
  label_lat: number | null;
}

export interface PlatesFrame {
  step: number;
  time: number;
  plates: PlateInfo[];
}

const DEG = Math.PI / 180;

type Projected = [number, number, number] | null;

/**
 * A flat projector for VelocityField that drops arrows crossing the map's seam.
 *
 * VelocityField projects each arrow's base and then its tip, and draws a straight
 * screen line between them, with no seam handling of its own. On a flat map an arrow
 * whose tip is across the antimeridian would be drawn right across the map. This wrapper
 * relies on that base-then-tip call order: after a base projects, the next call is its
 * tip, and the tip is refused (the arrow skipped) when the two straddle the seam.
 */
class SeamSafeArrowProjector {
  private base: number[] | null = null;
  constructor(private readonly inner: FlatProjector) {}
  update(w: number, h: number): void { this.inner.update(w, h); this.base = null; }
  project(v: ArrayLike<number>): Projected {
    if (this.base === null) {
      const p = this.inner.project(v);
      if (p) this.base = [v[0], v[1], v[2]];
      return p;
    }
    const base = this.base;
    this.base = null;
    return this.inner.seamSplit(base, v) ? null : this.inner.project(v);
  }
}
export const ARROW_DEG_PER_KM_MYR = 0.12;

/**
 * Plate velocity arrows and plate-ID labels, drawn on a 2D canvas over the globe.
 *
 * Arrows are petrify's VelocityField: omega x r at fixed points, computed in the export
 * from each plate's pole and rate, with each tip found along the great circle the plate
 * moves on (see velocities.js). Labels sit at the plate's area-weighted centroid, as in
 * diagnose_run.py's maps. Both are projected with the same projector as the boundaries,
 * so they hide at the horizon the same way.
 */
export class PlateOverlay {
  readonly canvas: HTMLCanvasElement;
  private readonly ctx: CanvasRenderingContext2D;
  private projector!: ThreeProjector | FlatProjector;
  private arrowProjector!: ThreeProjector | SeamSafeArrowProjector;
  private velocity: VelocityField | null = null;
  private plates: PlatesFrame | null = null;
  showArrows = true;
  showLabels = true;

  constructor(camera: Camera, mode: ProjectionMode = 'globe') {
    this.setCamera(camera, mode);
    this.canvas = document.createElement('canvas');
    Object.assign(this.canvas.style, { position: 'fixed', left: '0', top: '0', pointerEvents: 'none' });
    document.body.appendChild(this.canvas);
    this.ctx = this.canvas.getContext('2d')!;
    this.resize();
  }

  /** A new camera is built for every projection, so the projectors are rebuilt too. */
  setCamera(camera: Camera, mode: ProjectionMode): void {
    if (mode === 'globe') {
      this.projector = new ThreeProjector(camera as PerspectiveCamera);
      this.arrowProjector = this.projector;
    } else {
      const flat = new FlatProjector(camera);
      flat.setFlatMode(mode);
      this.projector = flat;
      this.arrowProjector = new SeamSafeArrowProjector(flat);
    }
  }

  resize(): void {
    const dpr = Math.min(devicePixelRatio, 2);
    Object.assign(this.canvas.style, { width: `${innerWidth}px`, height: `${innerHeight}px` });
    this.canvas.width = Math.round(innerWidth * dpr);
    this.canvas.height = Math.round(innerHeight * dpr);
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  async loadVelocities(url: string): Promise<void> {
    // Black arrows and a 0.5 km/Myr floor, as diagnose_run.py draws them. Length: 0.12 deg
    // of arc per km/Myr, i.e. the distance the plate moves in about 13 Myr, so a
    // 30 km/Myr plate gets a 3.6 deg arrow at the 10 deg arrow spacing.
    this.velocity = await VelocityField.load(url, {
      colour: 'rgba(0,0,0,0.85)', minSpeed: 0.5, scale: ARROW_DEG_PER_KM_MYR,
    });
  }

  setFrame(time: number, plates: PlatesFrame | null): void {
    this.velocity?.setTime(time);
    this.plates = plates;
  }

  draw(): void {
    const ctx = this.ctx;
    ctx.save();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.restore();
    this.projector.update(innerWidth, innerHeight);
    if (this.showArrows && this.velocity) {
      this.arrowProjector.update(innerWidth, innerHeight);
      this.velocity.draw(ctx, this.arrowProjector);
    }
    if (this.showLabels && this.plates) this.drawLabels();
  }

  private drawLabels(): void {
    const ctx = this.ctx;
    ctx.save();
    ctx.font = 'bold 11px ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    const v = [0, 0, 0];
    for (const p of this.plates!.plates) {
      if (p.label_lon === null || p.label_lat === null) continue;
      const lon = p.label_lon * DEG, lat = p.label_lat * DEG;
      v[0] = Math.cos(lat) * Math.cos(lon); v[1] = Math.cos(lat) * Math.sin(lon); v[2] = Math.sin(lat);
      const s = this.projector.project(v);
      if (!s) continue;
      const text = String(p.id);
      const w = ctx.measureText(text).width + 6;
      ctx.fillStyle = 'rgba(0,0,0,0.6)';
      ctx.fillRect(s[0] - w / 2, s[1] - 8, w, 16);
      ctx.fillStyle = '#ffffff';
      ctx.fillText(text, s[0], s[1] + 0.5);
    }
    ctx.restore();
  }
}
