/** The exported run files written by prep/export_run.py. */

export const DATA = `${import.meta.env.BASE_URL}data`;

export interface RunEntry {
  name: string;
  label: string;
  path: string;
}

export interface StepInfo {
  step: number;
  time: number;
  n_plates: number;
}

export interface LayerInfo {
  id: string;
  kind: string;
  path: string;
  has_cores?: boolean;
  /** 'scalar' layers: value = byte * scale, `nodata` not drawn, colour map from 0 to max. */
  scale?: number;
  nodata?: number;
  max?: number;
  label?: string;
  units?: string;
}

export interface RunInfo {
  name: string;
  label: string;
  force_model: string;
  force_modes: Record<string, string>;
  steps: StepInfo[];
  /** Grid layout of every exported grid (Archive layout: nodes, poles included). */
  grid?: { deg: number; nlon: number; nlat: number };
  layers?: LayerInfo[];
}

async function fetchJson<T>(url: string): Promise<T> {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url}: ${r.status}`);
  return r.json() as Promise<T>;
}

export function fetchRuns(): Promise<RunEntry[]> {
  return fetchJson<RunEntry[]>(`${DATA}/runs.json`);
}

export function fetchRun(entry: RunEntry): Promise<RunInfo> {
  return fetchJson<RunInfo>(`${DATA}/${entry.path}/run.json`);
}

export const stepTag = (step: number): string => String(step).padStart(4, '0');

export { fetchJson };
