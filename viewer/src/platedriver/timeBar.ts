/**
 * The model-time slider: one bar along the bottom, 0 -> end, left to right.
 *
 * Labelled in forward model time t (Myr), never Ma -- a PlateDriver run is not
 * real Earth history, and an age readout would invite students to read it as
 * one (viewer_plan.md, "Time").
 *
 * The slider only offers the saved steps: `setSteps` takes the list of times
 * a run saved, and the thumb moves by index through them, so every value it
 * can report is one that has data.
 */

const UI_FONT = 'ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", sans-serif';

function ensureStyles(): void {
  if (document.getElementById('pd-timebar-style')) return;
  const style = document.createElement('style');
  style.id = 'pd-timebar-style';
  // Slider thumb/track pseudo-elements can't be reached from inline styles.
  style.textContent = `
    .pd-scrub { -webkit-appearance: none; appearance: none; flex: 1; height: 16px;
      background: transparent; cursor: pointer; }
    .pd-scrub::-webkit-slider-runnable-track { height: 3px; border-radius: 2px;
      background: rgba(255,255,255,0.2); }
    .pd-scrub::-moz-range-track { height: 3px; border-radius: 2px;
      background: rgba(255,255,255,0.2); }
    .pd-scrub::-webkit-slider-thumb { -webkit-appearance: none; appearance: none;
      width: 15px; height: 15px; margin-top: -6px; border: none; border-radius: 50%;
      background: #6fe3ff; box-shadow: 0 0 0 3px rgba(111,227,255,0.22); }
    .pd-scrub::-moz-range-thumb { width: 15px; height: 15px; border: none;
      border-radius: 50%; background: #6fe3ff; box-shadow: 0 0 0 3px rgba(111,227,255,0.22); }
  `;
  document.head.appendChild(style);
}

export class TimeBar {
  private readonly slider: HTMLInputElement;
  private readonly value: HTMLSpanElement;
  private readonly endLabel: HTMLSpanElement;
  private steps: number[] = [0];

  constructor(private readonly onTime: (t: number) => void) {
    ensureStyles();
    const bar = document.createElement('div');
    Object.assign(bar.style, {
      position: 'fixed', left: '16px', right: '16px', bottom: '16px', zIndex: '10',
      display: 'flex', alignItems: 'center', gap: '12px',
      padding: '8px 14px', borderRadius: '6px',
      background: 'rgba(12,16,22,0.82)', border: '1px solid rgba(255,255,255,0.10)',
      color: 'rgba(235,240,250,0.92)', font: `12px/1.4 ${UI_FONT}`,
    });

    const name = document.createElement('span');
    name.textContent = 'model time t (Myr)';
    name.style.color = 'rgba(235,240,250,0.62)';

    this.value = document.createElement('span');
    Object.assign(this.value.style, { minWidth: '4.5em', fontWeight: '600', fontVariantNumeric: 'tabular-nums' });

    this.endLabel = document.createElement('span');
    this.endLabel.style.color = 'rgba(235,240,250,0.62)';

    this.slider = document.createElement('input');
    this.slider.type = 'range';
    this.slider.className = 'pd-scrub';
    this.slider.min = '0';
    this.slider.step = '1';
    this.slider.addEventListener('input', () => { this.onTime(this.time); this.refresh(); });

    bar.append(name, this.value, this.slider, this.endLabel);
    document.body.appendChild(bar);
  }

  /** The saved step times, in Myr, ascending. Keeps the nearest saved time to
   *  the current one, so switching run holds the time step. */
  setSteps(times: number[]): void {
    const keep = this.time;
    this.steps = times.length ? [...times] : [0];
    this.slider.max = String(this.steps.length - 1);
    this.endLabel.textContent = `${this.steps[this.steps.length - 1]}`;
    this.slider.value = String(nearestIndex(this.steps, keep));
    this.refresh();
  }

  /** Move to the saved step nearest `t` (does not fire onTime). */
  setTime(t: number): void {
    this.slider.value = String(nearestIndex(this.steps, t));
    this.refresh();
  }

  /** The slider track's left and right edges, in viewport px. */
  trackRect(): { left: number; right: number } {
    const r = this.slider.getBoundingClientRect();
    return { left: r.left, right: r.right };
  }

  get time(): number {
    return this.steps[Number(this.slider.value)] ?? 0;
  }

  private refresh(): void {
    this.value.textContent = `t = ${this.time}`;
  }
}

/** Index of the saved step nearest `t`. */
export function nearestIndex(steps: readonly number[], t: number): number {
  let best = 0;
  for (let i = 1; i < steps.length; i++) {
    if (Math.abs(steps[i] - t) < Math.abs(steps[best] - t)) best = i;
  }
  return best;
}
