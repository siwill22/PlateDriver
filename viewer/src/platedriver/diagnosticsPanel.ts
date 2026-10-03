/**
 * Whole-run diagnostics above the time bar (export_run.py's diagnostics.json): one strip
 * per series, every model step, with the rift and weld times as ticks along the top.
 *
 * The x axis is the time slider's own: each strip spans exactly the slider's track, so a
 * point sits straight above the thumb position for that time. A vertical line marks the
 * current time; hovering reads off every series at the pointer, and a click moves the
 * slider there.
 */

export interface DiagLine { name: string; values: (number | null)[] }
export interface DiagSeries { id: string; label: string; units: string; lines: DiagLine[] }
export interface Diagnostics {
  time: number[];
  series: DiagSeries[];
  events: Record<string, number[]>;
}

const UI_FONT = 'ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
const LINE_COLORS = ['#6fe3ff', '#ffb35c', '#b6f27a', '#ff7aa8'];
const EVENT_COLORS: Record<string, string> = { rift: '#ff5a5a', weld: '#6fe3ff' };
const STRIP_H = 40;
const EVENT_H = 14;
const THUMB_HALF = 7.5;   // the slider thumb is 15 px wide: its centre runs 7.5 px in from each end

export class DiagnosticsPanel {
  readonly el: HTMLDivElement;
  private readonly canvas: HTMLCanvasElement;
  private readonly readout: HTMLDivElement;
  private data: Diagnostics | null = null;
  private shown = new Set<string>();
  private current = 0;
  private hoverT: number | null = null;
  private enabled = true;

  /** track(): the slider track's left and right edges in viewport px. */
  constructor(private readonly track: () => { left: number; right: number },
              private readonly onPick: (t: number) => void) {
    this.el = document.createElement('div');
    Object.assign(this.el.style, {
      position: 'fixed', left: '16px', right: '16px', bottom: '62px', zIndex: '10',
      padding: '6px 0 4px', borderRadius: '6px', display: 'none',
      background: 'rgba(12,16,22,0.82)', border: '1px solid rgba(255,255,255,0.10)',
      color: 'rgba(235,240,250,0.92)', font: `11px/1.3 ${UI_FONT}`,
    });
    this.canvas = document.createElement('canvas');
    Object.assign(this.canvas.style, { display: 'block', width: '100%', cursor: 'crosshair' });
    this.readout = document.createElement('div');
    Object.assign(this.readout.style, {
      position: 'fixed', zIndex: '11', pointerEvents: 'none', display: 'none', whiteSpace: 'pre',
      padding: '4px 8px', borderRadius: '4px', background: 'rgba(12,16,22,0.92)',
      border: '1px solid rgba(255,255,255,0.15)', color: 'rgba(235,240,250,0.95)',
      font: `11px/1.4 ${UI_FONT}`, fontVariantNumeric: 'tabular-nums',
    });
    this.el.appendChild(this.canvas);
    document.body.append(this.el, this.readout);

    this.canvas.addEventListener('mousemove', (e) => { this.hoverT = this.timeAt(e.clientX); this.draw(e); });
    this.canvas.addEventListener('mouseleave', () => {
      this.hoverT = null; this.readout.style.display = 'none'; this.draw();
    });
    this.canvas.addEventListener('click', (e) => {
      const t = this.timeAt(e.clientX);
      if (t !== null) this.onPick(t);
    });
    addEventListener('resize', () => this.draw());
  }

  setData(d: Diagnostics | null): void {
    this.data = d;
    this.updateVisibility();
    this.draw();
  }

  /** Which series ids are drawn, in the data's order. */
  setShown(ids: Iterable<string>): void {
    this.shown = new Set(ids);
    this.updateVisibility();
    this.draw();
  }

  set visible(v: boolean) { this.enabled = v; this.updateVisibility(); this.draw(); }

  setTime(t: number): void { this.current = t; this.draw(); }

  /** Height in px while shown (0 when hidden), so other bottom-anchored boxes can clear it. */
  get height(): number { return this.el.style.display === 'none' ? 0 : this.el.offsetHeight; }

  private series(): DiagSeries[] {
    return this.data ? this.data.series.filter((s) => this.shown.has(s.id)) : [];
  }

  private updateVisibility(): void {
    this.el.style.display = this.enabled && this.data && this.series().length ? 'block' : 'none';
  }

  private range(): [number, number] {
    const t = this.data!.time;
    return [t[0], t[t.length - 1]];
  }

  /** x in canvas px of time t, on the slider track's scale. */
  private xOf(t: number): number {
    const box = this.canvas.getBoundingClientRect();
    const { left, right } = this.track();
    const [t0, t1] = this.range();
    const a = left + THUMB_HALF - box.left;
    const b = right - THUMB_HALF - box.left;
    return a + (t1 > t0 ? (t - t0) / (t1 - t0) : 0) * (b - a);
  }

  private timeAt(clientX: number): number | null {
    if (!this.data) return null;
    const box = this.canvas.getBoundingClientRect();
    const { left, right } = this.track();
    const [t0, t1] = this.range();
    const f = (clientX - (left + THUMB_HALF)) / (right - left - 2 * THUMB_HALF);
    if (f < -0.02 || f > 1.02 || box.width <= 0) return null;
    return Math.round(t0 + Math.min(1, Math.max(0, f)) * (t1 - t0));
  }

  private draw(e?: MouseEvent): void {
    if (this.el.style.display === 'none' || !this.data) return;
    const rows = this.series();
    const d = this.data;
    const cssW = this.canvas.clientWidth;
    const cssH = EVENT_H + rows.length * STRIP_H;
    const dpr = Math.min(devicePixelRatio, 2);
    this.canvas.style.height = `${cssH}px`;
    this.canvas.width = Math.round(cssW * dpr);
    this.canvas.height = Math.round(cssH * dpr);
    const ctx = this.canvas.getContext('2d')!;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);
    ctx.font = `11px ${UI_FONT}`;
    const [t0, t1] = this.range();
    const xa = this.xOf(t0);
    const xb = this.xOf(t1);

    // events: rift and weld ticks
    ctx.textBaseline = 'middle';
    ctx.fillStyle = 'rgba(235,240,250,0.62)';
    ctx.fillText('rifts / welds', 14, EVENT_H / 2);
    for (const [kind, times] of Object.entries(d.events)) {
      ctx.strokeStyle = EVENT_COLORS[kind] ?? '#ffffff';
      ctx.lineWidth = 1;
      const y0 = kind === 'rift' ? 2 : EVENT_H / 2;
      ctx.beginPath();
      for (const t of times) {
        const x = this.xOf(t);
        ctx.moveTo(x + 0.5, y0);
        ctx.lineTo(x + 0.5, y0 + EVENT_H / 2 - 2);
      }
      ctx.stroke();
    }

    const readout: string[] = [];
    const hi = this.hoverT === null ? -1 : nearest(d.time, this.hoverT);
    if (hi >= 0) readout.push(`t = ${d.time[hi]} Myr`);

    rows.forEach((s, r) => {
      const top = EVENT_H + r * STRIP_H;
      const plotTop = top + 4;
      const plotH = STRIP_H - 10;
      let lo = Infinity;
      let up = -Infinity;
      for (const l of s.lines) for (const v of l.values) if (v !== null) { lo = Math.min(lo, v); up = Math.max(up, v); }
      if (!Number.isFinite(lo)) { lo = 0; up = 1; }
      if (lo >= 0 && lo < 0.25 * up) lo = 0;     // start at zero unless the series sits far from it
      if (up - lo < 1e-9) up = lo + 1;
      const yOf = (v: number) => plotTop + plotH - ((v - lo) / (up - lo)) * plotH;

      ctx.strokeStyle = 'rgba(255,255,255,0.08)';
      ctx.beginPath();
      ctx.moveTo(xa, top + 0.5);
      ctx.lineTo(xb, top + 0.5);
      ctx.stroke();

      // labels in the left gutter, the range at the right
      ctx.fillStyle = 'rgba(235,240,250,0.92)';
      ctx.textAlign = 'left';
      ctx.fillText(s.label, 14, plotTop + 6);
      ctx.fillStyle = 'rgba(235,240,250,0.55)';
      ctx.fillText(s.units, 14, plotTop + 18);
      if (s.lines.length > 1) {
        // the line names, side by side, under the units
        let x = 14;
        s.lines.forEach((l, k) => {
          ctx.fillStyle = LINE_COLORS[k % LINE_COLORS.length];
          ctx.fillText(l.name, x, plotTop + 30);
          x += ctx.measureText(l.name).width + 8;
        });
      }
      ctx.fillStyle = 'rgba(235,240,250,0.55)';
      ctx.textAlign = 'left';
      ctx.fillText(fmt(up), xb + 6, plotTop + 4);
      ctx.fillText(fmt(lo), xb + 6, plotTop + plotH - 2);

      s.lines.forEach((l, k) => {
        ctx.strokeStyle = LINE_COLORS[k % LINE_COLORS.length];
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        let pen = false;
        l.values.forEach((v, i) => {
          if (v === null) { pen = false; return; }
          const x = this.xOf(d.time[i]);
          const y = yOf(v);
          if (pen) ctx.lineTo(x, y); else ctx.moveTo(x, y);
          pen = true;
        });
        ctx.stroke();
        if (hi >= 0) {
          const v = l.values[hi];
          readout.push(`${s.label}${s.lines.length > 1 ? ` (${l.name})` : ''}: ${v === null ? '–' : fmt(v)} ${s.units}`);
        }
      });
    });

    // current time and hover lines
    const mark = (t: number, color: string) => {
      const x = Math.round(this.xOf(t)) + 0.5;
      ctx.strokeStyle = color;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x, 0);
      ctx.lineTo(x, cssH);
      ctx.stroke();
    };
    mark(this.current, 'rgba(111,227,255,0.9)');
    if (hi >= 0) mark(d.time[hi], 'rgba(255,255,255,0.45)');

    if (e && hi >= 0) {
      this.readout.textContent = readout.join('\n');
      this.readout.style.display = 'block';
      const w = this.readout.offsetWidth;
      const h = this.readout.offsetHeight;
      this.readout.style.left = `${Math.min(innerWidth - w - 8, e.clientX + 12)}px`;
      this.readout.style.top = `${Math.max(8, e.clientY - h - 12)}px`;
    }
  }
}

function nearest(times: number[], t: number): number {
  let best = 0;
  for (let i = 1; i < times.length; i++) if (Math.abs(times[i] - t) < Math.abs(times[best] - t)) best = i;
  return best;
}

function fmt(v: number): string {
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 100) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  if (a >= 1) return v.toFixed(2);
  return v.toFixed(3);
}
