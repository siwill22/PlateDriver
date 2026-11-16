import GUI from 'lil-gui';
import { Color, Scene, WebGLRenderer, type Camera, type PerspectiveCamera } from 'three';
import type { OrbitControls } from 'three/addons/controls/OrbitControls.js';

import '../core/lilgui-theme.css';
import { BoundaryOverlay } from '../core/boundaries';
import { lonLatToVec3 } from '../core/constants';
import { Graticule } from '../core/graticule';
import { OceanSurface } from '../core/oceanSurface';
import {
  createProjectionCamera, createProjectionControls, PROJECTION_LABEL, PROJECTION_ORDER,
  updateProjectionCameraAspect, type ProjectionMode,
} from '../core/projection';
import { applyChromeLightness, resolveTheme } from '../core/theme';
import { ContinentLayer } from './continentLayer';
import { DiagnosticsPanel, type Diagnostics } from './diagnosticsPanel';
import { PlateOverlay, type PlatesFrame } from './plateOverlay';
import { ColorBar, ScalarLayer } from './scalarLayer';
import {
  DATA, fetchJson, fetchRun, fetchRuns, stepTag, type RunEntry, type RunInfo,
} from './runData';
import { nearestIndex, TimeBar } from './timeBar';

// A light, neutral globe, so boundaries can keep diagnose_run.py's map colours:
// ridge red, trench blue, other black.
const theme = resolveTheme('newsprint');
applyChromeLightness(theme.lightness);

/** Full entries, not partial ones: BoundaryOverlay.load's override replaces each type's
 *  style wholesale. PlateDriver's 'conservative' is exported as 'other'. */
const BOUNDARY_STYLE = {
  subduction: { stroke: '#0000ff', width: 1.9, label: 'Trench' },
  ridge: { stroke: '#ff0000', width: 1.6, label: 'Ridge' },
  transform: { stroke: '#000000', width: 1.0, label: 'Transform' },
  other: { stroke: '#000000', width: 1.0, label: 'Other boundary' },
};

const params = new URLSearchParams(location.search);

/** What the panel controls. `projection` starts from ?proj= (globe by default). */
const view = {
  projection: (PROJECTION_ORDER.find((m) => m === params.get('proj')) ?? 'globe') as ProjectionMode,
  continents: true,
  belts: true,
  strength: params.get('layer') === 'strength',
  strengthOpacity: 0.9,
  boundaries: true,
  arrows: true,
  plateIds: true,
  diagnostics: true,
};

/** Diagnostics series drawn above the time bar, by id (export_run.py's diagnostics()). */
const diagShown: Record<string, boolean> = {
  plates: true, rms: true, boundaries: true, largest: true, netrot: false, hidden: false,
};
const DIAG_LABEL: Record<string, string> = {
  plates: 'Plates', rms: 'RMS speed', boundaries: 'Boundary length', largest: 'Largest plate',
  netrot: 'Net rotation', hidden: 'Continent hidden',
};

const scene = new Scene();
const renderer = new WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
renderer.setClearColor(new Color(theme.page));
document.body.appendChild(renderer.domElement);

// A projection switch builds a new camera and new controls (core/projection.ts):
// the globe and the flat maps need different camera types.
let mode: ProjectionMode = view.projection;
let camera: Camera = createProjectionCamera(mode, innerWidth / innerHeight);
let controls: OrbitControls = createProjectionControls(mode, camera, renderer.domElement);

// ?view=<lon>,<lat> points the globe camera at that spot (used by the screenshot checks).
const viewParam = params.get('view');
if (viewParam && mode === 'globe') {
  const [lon, lat] = viewParam.split(',').map(Number);
  const [x, y, z] = lonLatToVec3(lon, lat, camera.position.length());
  camera.position.set(x, y, z);
  camera.lookAt(0, 0, 0);
}

const ocean = new OceanSurface(mode);
ocean.applyTheme(theme);
scene.add(ocean.mesh);

const continents = new ContinentLayer(mode);
scene.add(continents.mesh);

const strength = new ScalarLayer(mode);
scene.add(strength.mesh);
const colorBar = new ColorBar();

const graticule = new Graticule(30, 30);
graticule.setProjection(mode);
graticule.applyTheme(theme);
scene.add(graticule.lines);

// Overlay canvases, in paint order: boundaries, then arrows and labels on top.
const boundaries = new BoundaryOverlay(camera as PerspectiveCamera);
boundaries.setCamera(camera, mode);
const plateOverlay = new PlateOverlay(camera, mode);

function applyProjection(next: ProjectionMode): void {
  if (next === mode) return;
  mode = next;
  controls.dispose();
  camera = createProjectionCamera(mode, innerWidth / innerHeight);
  controls = createProjectionControls(mode, camera, renderer.domElement);
  ocean.setProjection(mode);
  continents.setProjection(mode);
  strength.setProjection(mode);
  graticule.setProjection(mode);
  boundaries.setCamera(camera, mode);
  plateOverlay.setCamera(camera, mode);
}

const status = document.createElement('div');
Object.assign(status.style, {
  position: 'fixed', left: '16px', top: '16px', zIndex: '10',
  padding: '8px 12px', borderRadius: '6px',
  background: 'rgba(12,16,22,0.82)', border: '1px solid rgba(255,255,255,0.10)',
  color: 'rgba(235,240,250,0.92)', whiteSpace: 'pre',
  font: '12px/1.5 ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif',
});
document.body.appendChild(status);

// --- panel ------------------------------------------------------------------
// lil-gui writes the new value into `view` before onChange fires, so each handler
// acts on the value it is given and never compares it against `view`.

const gui = new GUI({ title: 'PlateDriver' });
gui.domElement.style.zIndex = '10';
gui.add(view, 'projection', Object.fromEntries(PROJECTION_ORDER.map((m) => [PROJECTION_LABEL[m], m])))
  .name('Projection').onChange((m: ProjectionMode) => { applyProjection(m); });
const layers = gui.addFolder('Layers');
layers.add(view, 'continents').name('Continents').onChange((v: boolean) => { continents.visible = v; });
const beltsCtl = layers.add(view, 'belts').name('Belts (orange)')
  .onChange((v: boolean) => { continents.belts = v; });
const strengthCtl = layers.add(view, 'strength').name('Rift strength').onChange((v: boolean) => {
  strength.visible = v;
  colorBar.visible = v && strengthGrid !== null;
  void applyTime(timeBar.time);
});
const strengthOpacityCtl = layers.add(view, 'strengthOpacity', 0, 1, 0.05).name('Strength opacity')
  .onChange((v: number) => { strength.opacity = v; });
layers.add(view, 'boundaries').name('Plate boundaries').onChange((v: boolean) => { boundaries.visible = v; });
layers.add(view, 'arrows').name('Velocity arrows').onChange((v: boolean) => { plateOverlay.showArrows = v; });
layers.add(view, 'plateIds').name('Plate IDs').onChange((v: boolean) => { plateOverlay.showLabels = v; });

const timeBar = new TimeBar((t) => { void applyTime(t); });
const diag = new DiagnosticsPanel(() => timeBar.trackRect(), (t) => {
  timeBar.setTime(t);
  void applyTime(timeBar.time);
});
const applyDiagShown = () => diag.setShown(Object.keys(diagShown).filter((k) => diagShown[k]));
applyDiagShown();
const diagFolder = gui.addFolder('Diagnostics');
diagFolder.add(view, 'diagnostics').name('Show plot').onChange((v: boolean) => { diag.visible = v; });
for (const id of Object.keys(diagShown)) {
  diagFolder.add(diagShown, id).name(DIAG_LABEL[id]).onChange(() => applyDiagShown());
}

let runEntry: RunEntry | null = null;
let run: RunInfo | null = null;
let continentGrid: { nlon: number; nlat: number } | null = null;
let strengthGrid: { nlon: number; nlat: number; path: string } | null = null;
const strengthCache = new Map<number, Promise<Uint8Array | null>>();
const platesCache = new Map<number, Promise<PlatesFrame | null>>();
const continentCache = new Map<number, Promise<Uint8Array | null>>();
/** Bumped on every applyTime, so a slow fetch can't overwrite a later step. */
let request = 0;

function cached<T>(cache: Map<number, Promise<T | null>>, step: number, load: () => Promise<T>): Promise<T | null> {
  let p = cache.get(step);
  if (!p) {
    p = load().catch(() => null);
    cache.set(step, p);
  }
  return p;
}

const platesAt = (step: number) => cached(platesCache, step,
  () => fetchJson<PlatesFrame>(`${DATA}/${runEntry!.path}/plates/${stepTag(step)}.json`));

const continentsAt = (step: number) => cached(continentCache, step, async () => {
  const r = await fetch(`${DATA}/${runEntry!.path}/grids/continents/${stepTag(step)}.bin`);
  if (!r.ok) throw new Error(String(r.status));
  return new Uint8Array(await r.arrayBuffer());
});

const strengthAt = (step: number) => cached(strengthCache, step, async () => {
  const r = await fetch(`${DATA}/${runEntry!.path}/${strengthGrid!.path}/${stepTag(step)}.bin`);
  if (!r.ok) throw new Error(String(r.status));
  return new Uint8Array(await r.arrayBuffer());
});

async function applyTime(t: number): Promise<void> {
  if (!run) { status.textContent = `PlateDriver viewer\nno run loaded   t = ${t} Myr`; return; }
  const info = run.steps[nearestIndex(run.steps.map((s) => s.time), t)];
  const mine = ++request;
  void boundaries.setAge(info.time);
  // The strength grid is only fetched while its layer is on.
  const [plates, grid, sgrid] = await Promise.all([
    platesAt(info.step), continentGrid ? continentsAt(info.step) : Promise.resolve(null),
    strengthGrid && view.strength ? strengthAt(info.step) : Promise.resolve(null),
  ]);
  if (mine !== request) return;
  diag.setTime(info.time);
  plateOverlay.setFrame(info.time, plates);
  if (continentGrid) continents.setGrid(grid, continentGrid.nlon, continentGrid.nlat);
  if (strengthGrid) strength.setGrid(sgrid, strengthGrid.nlon, strengthGrid.nlat);
  status.textContent = `${run.label}\n`
    + `${run.force_model}   t = ${info.time} Myr   (step ${info.step})   ${info.n_plates} plates`;
}

async function loadRun(entry: RunEntry): Promise<void> {
  runEntry = entry;
  run = await fetchRun(entry);
  platesCache.clear();
  continentCache.clear();
  strengthCache.clear();
  const sl = run.layers?.find((l) => l.id === 'strength' && l.kind === 'scalar');
  strengthGrid = sl && run.grid ? { nlon: run.grid.nlon, nlat: run.grid.nlat, path: sl.path } : null;
  if (sl) {
    strength.setStyle({ scale: sl.scale ?? 0.01, max: sl.max ?? 1, nodata: sl.nodata ?? 255 });
    colorBar.set(`${sl.label ?? 'Rift strength'} (${sl.units ?? ''})`, sl.max ?? 1);
  } else {
    strength.setGrid(null, 1, 1);
  }
  strengthCtl.enable(Boolean(sl));
  strengthOpacityCtl.enable(Boolean(sl));
  strength.visible = view.strength;
  colorBar.visible = view.strength && Boolean(sl);
  const layer = run.layers?.find((l) => l.id === 'continents');
  continentGrid = layer && run.grid ? { nlon: run.grid.nlon, nlat: run.grid.nlat } : null;
  // Belts only exist with continent_strength: cores.
  beltsCtl.enable(Boolean(layer?.has_cores));
  const [, , diagData] = await Promise.all([
    boundaries.load(`${DATA}/${entry.path}/boundaries/manifest.json`, BOUNDARY_STYLE),
    plateOverlay.loadVelocities(`${DATA}/${entry.path}/velocities.json`),
    fetchJson<Diagnostics>(`${DATA}/${entry.path}/diagnostics.json`).catch(() => null),
  ]);
  diag.setData(diagData);
  timeBar.setSteps(run.steps.map((s) => s.time));
  await applyTime(timeBar.time);
}

async function boot(): Promise<void> {
  timeBar.setSteps([0]);
  let runs: RunEntry[] = [];
  try { runs = await fetchRuns(); } catch { /* no runs exported yet */ }
  // ?run=<name> picks a run; otherwise the most recently exported one.
  const want = params.get('run');
  const entry = runs.find((r) => r.name === want) ?? runs[runs.length - 1];
  if (entry) {
    const pick = { run: entry.name };
    gui.add(pick, 'run', Object.fromEntries(runs.map((r) => [r.label, r.name]))).name('Run')
      .onChange((name: string) => {
        const next = runs.find((r) => r.name === name);
        if (next) void loadRun(next);
      });
    await loadRun(entry);
  }
  else await applyTime(0);
  // ?t=<Myr> opens at a model time (used by the screenshot checks).
  const t = params.get('t');
  if (t !== null && run) { timeBar.setTime(Number(t)); await applyTime(timeBar.time); }
}

addEventListener('resize', () => {
  updateProjectionCameraAspect(camera, innerWidth / innerHeight);
  renderer.setSize(innerWidth, innerHeight);
  boundaries.setRect({ x: 0, y: 0, width: innerWidth, height: innerHeight });
  plateOverlay.resize();
});

renderer.setAnimationLoop(() => {
  controls.update();
  renderer.render(scene, camera);
  boundaries.draw();
  plateOverlay.draw();
});

void boot();
