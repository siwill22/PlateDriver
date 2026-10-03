import {
  ClampToEdgeWrapping, DataTexture, LinearFilter, Mesh, NearestFilter, RedFormat, RepeatWrapping,
  RGBAFormat, ShaderMaterial, UnsignedByteType, type BufferGeometry,
} from 'three';

import { LIGHT_DIR, R_SURFACE } from '../core/constants';
import { GEOGRAPHIC_GLSL } from '../core/glsl/geographic';
import {
  createSurfaceGeometry, isFlat, PROJECTION_UNIFORM, type ProjectionMode,
} from '../core/projection';

/**
 * A scalar grid draped on the surface through a colour map (export_run.py's 'scalar'
 * layers): one byte per node, value = byte * scale, `nodata` not drawn.
 *
 * Same grid layout and nearest-node lookup as ContinentLayer (Geode Archive layout), so
 * each node's value is shown as it is, never blended with its neighbour. The colour map
 * runs from 0 to `max`.
 */

const VERT = /* glsl */ `
varying vec3 vWorldPos;
void main() {
  vWorldPos = (modelMatrix * vec4(position, 1.0)).xyz;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

const FRAG = /* glsl */ `
precision highp float;
${GEOGRAPHIC_GLSL}
uniform sampler2D uValues;
uniform sampler2D uColormap;
uniform vec2 uGridSize;          // nlon, nlat
uniform float uScale;            // value per byte
uniform float uMax;              // value at the top of the colour map
uniform float uNoData;           // byte not drawn
uniform float uOpacity;
uniform float uProjectionMode;   // 0 globe, 1 Plate Carree, 2 Robinson
uniform vec3 uLightDir;
uniform float uShadeStrength;
varying vec3 vWorldPos;

void main() {
  vec2 ll;
  if (uProjectionMode < 0.5) ll = worldToGeographic(vWorldPos);
  else if (uProjectionMode < 1.5) ll = worldToGeographicFlat(vWorldPos);
  else {
    ll = worldToGeographicRobinson(vWorldPos);
    if (robinsonOffMap(ll)) discard;
  }
  float i = (ll.x + PI) / (2.0 * PI) * uGridSize.x;
  float j = (ll.y + 0.5 * PI) / PI * (uGridSize.y - 1.0);
  vec2 uv = vec2((floor(i + 0.5) + 0.5) / uGridSize.x, (floor(j + 0.5) + 0.5) / uGridSize.y);
  float b = floor(texture2D(uValues, uv).r * 255.0 + 0.5);
  if (abs(b - uNoData) < 0.5) discard;
  float x = clamp(b * uScale / uMax, 0.0, 1.0);
  vec3 col = texture2D(uColormap, vec2(x, 0.5)).rgb;
  vec3 n = normalize(vWorldPos);
  float ndl = dot(n, normalize(uLightDir)) * 0.5 + 0.5;
  float shade = mix(1.0, 0.5 + 0.5 * ndl * ndl, uShadeStrength);
  gl_FragColor = vec4(col * shade, uOpacity);
}
`;

/** viridis at eleven stops (matplotlib), interpolated to 256 entries. */
const VIRIDIS = [
  '#440154', '#482475', '#414487', '#355f8d', '#2a788e', '#21918c',
  '#22a884', '#44bf70', '#7ad151', '#bddf26', '#fde725',
];

function hexRgb(hex: string): [number, number, number] {
  const n = parseInt(hex.slice(1), 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

/** 256 RGBA entries of the colour map. */
export function colormapBytes(stops: string[] = VIRIDIS): Uint8Array {
  const rgb = stops.map(hexRgb);
  const out = new Uint8Array(256 * 4);
  for (let k = 0; k < 256; k++) {
    const x = (k / 255) * (rgb.length - 1);
    const a = Math.min(Math.floor(x), rgb.length - 2);
    const f = x - a;
    for (let c = 0; c < 3; c++) out[4 * k + c] = Math.round(rgb[a][c] * (1 - f) + rgb[a + 1][c] * f);
    out[4 * k + 3] = 255;
  }
  return out;
}

export interface ScalarStyle {
  scale: number;
  max: number;
  nodata: number;
}

export class ScalarLayer {
  readonly mesh: Mesh;
  private readonly material: ShaderMaterial;
  private texture: DataTexture | null = null;
  private mode: ProjectionMode;
  private hasData = false;
  private enabled = false;

  constructor(mode: ProjectionMode = 'globe') {
    this.mode = mode;
    const cmap = new DataTexture(colormapBytes(), 256, 1, RGBAFormat, UnsignedByteType);
    cmap.magFilter = LinearFilter;
    cmap.minFilter = LinearFilter;
    cmap.wrapS = ClampToEdgeWrapping;
    cmap.needsUpdate = true;
    this.material = new ShaderMaterial({
      vertexShader: VERT,
      fragmentShader: FRAG,
      transparent: true,
      uniforms: {
        uValues: { value: null },
        uColormap: { value: cmap },
        uGridSize: { value: [720, 361] },
        uScale: { value: 0.01 },
        uMax: { value: 1 },
        uNoData: { value: 255 },
        uOpacity: { value: 0.9 },
        uProjectionMode: { value: PROJECTION_UNIFORM[mode] },
        uLightDir: { value: LIGHT_DIR.clone() },
        uShadeStrength: { value: isFlat(mode) ? 0 : 1 },
      },
    });
    this.mesh = new Mesh(this.buildGeometry(mode), this.material);
    // Above the continents (2), below the graticule and the overlay canvases.
    this.mesh.renderOrder = 3;
    this.mesh.frustumCulled = false;
    this.mesh.visible = false;
  }

  private buildGeometry(mode: ProjectionMode): BufferGeometry {
    return createSurfaceGeometry(mode, R_SURFACE * 1.0009);
  }

  setStyle(style: ScalarStyle): void {
    const u = this.material.uniforms;
    u.uScale.value = style.scale;
    u.uMax.value = style.max > 0 ? style.max : 1;
    u.uNoData.value = style.nodata;
  }

  set opacity(v: number) { this.material.uniforms.uOpacity.value = v; }
  set visible(v: boolean) { this.enabled = v; this.mesh.visible = v && this.hasData; }

  /** One step's bytes, nlon x nlat, longitude fastest. */
  setGrid(bytes: Uint8Array | null, nlon: number, nlat: number): void {
    this.hasData = bytes !== null;
    this.mesh.visible = this.enabled && this.hasData;
    if (!bytes) return;
    if (!this.texture || this.texture.image.width !== nlon || this.texture.image.height !== nlat) {
      this.texture?.dispose();
      this.texture = new DataTexture(bytes, nlon, nlat, RedFormat, UnsignedByteType);
      this.texture.magFilter = NearestFilter;
      this.texture.minFilter = NearestFilter;
      this.texture.wrapS = RepeatWrapping;
      this.texture.wrapT = ClampToEdgeWrapping;
      this.texture.generateMipmaps = false;
      this.material.uniforms.uValues.value = this.texture;
      this.material.uniforms.uGridSize.value = [nlon, nlat];
    } else {
      (this.texture.image.data as Uint8Array).set(bytes);
    }
    this.texture.needsUpdate = true;
  }

  setProjection(mode: ProjectionMode): void {
    if (mode === this.mode) return;
    this.mode = mode;
    this.mesh.geometry.dispose();
    this.mesh.geometry = this.buildGeometry(mode);
    this.material.uniforms.uProjectionMode.value = PROJECTION_UNIFORM[mode];
    this.material.uniforms.uShadeStrength.value = isFlat(mode) ? 0 : 1;
  }
}

/** A horizontal colour bar with a title and end labels, for a ScalarLayer. */
export class ColorBar {
  readonly el: HTMLDivElement;
  private readonly title: HTMLDivElement;
  private readonly lo: HTMLSpanElement;
  private readonly hi: HTMLSpanElement;

  constructor() {
    this.el = document.createElement('div');
    Object.assign(this.el.style, {
      position: 'fixed', left: '16px', top: '74px', zIndex: '10', width: '240px',
      padding: '8px 12px', borderRadius: '6px', display: 'none',
      background: 'rgba(12,16,22,0.82)', border: '1px solid rgba(255,255,255,0.10)',
      color: 'rgba(235,240,250,0.92)',
      font: '12px/1.4 ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif',
    });
    this.title = document.createElement('div');
    const canvas = document.createElement('canvas');
    canvas.width = 256;
    canvas.height = 1;
    Object.assign(canvas.style, { width: '100%', height: '12px', display: 'block', margin: '4px 0 2px' });
    const ctx = canvas.getContext('2d')!;
    const img = ctx.createImageData(256, 1);
    img.data.set(colormapBytes());
    ctx.putImageData(img, 0, 0);
    const ends = document.createElement('div');
    Object.assign(ends.style, { display: 'flex', justifyContent: 'space-between' });
    this.lo = document.createElement('span');
    this.hi = document.createElement('span');
    ends.append(this.lo, this.hi);
    this.el.append(this.title, canvas, ends);
    document.body.appendChild(this.el);
  }

  set(title: string, max: number): void {
    this.title.textContent = title;
    this.lo.textContent = '0';
    this.hi.textContent = max.toFixed(2);
  }

  set visible(v: boolean) { this.el.style.display = v ? 'block' : 'none'; }
}
