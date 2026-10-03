import {
  DataTexture, Mesh, NearestFilter, RedFormat, RepeatWrapping, ClampToEdgeWrapping,
  ShaderMaterial, UnsignedByteType, Vector4, type BufferGeometry,
} from 'three';

import { LIGHT_DIR, R_SURFACE } from '../core/constants';
import { GEOGRAPHIC_GLSL } from '../core/glsl/geographic';
import {
  createSurfaceGeometry, isFlat, PROJECTION_UNIFORM, type ProjectionMode,
} from '../core/projection';

/**
 * The continents, as a class grid draped on the surface: 0 ocean (not drawn),
 * 1 continent, 2 craton core, 3 belt (export_run.py's continent_classes).
 *
 * The grid is in Geode's Archive layout -- nodes at lon -180 + i * 360 / nlon and
 * lat -90 + j * 180 / (nlat - 1) -- and is looked up nearest-node, so a class is never
 * blended with its neighbour. lon/lat come from the fragment's own position through
 * GEOGRAPHIC_GLSL's inverse for the active projection, so the same layer draws on the
 * globe and on both flat maps.
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
uniform sampler2D uClasses;
uniform vec2 uGridSize;          // nlon, nlat
uniform vec4 uColors[4];         // rgba per class; alpha 0 = not drawn
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
  // Node i sits at u = i / nlon of the way round; +0.5 texel lands on its texel centre.
  // Latitude nodes include both poles, so nlat - 1 intervals.
  float i = (ll.x + PI) / (2.0 * PI) * uGridSize.x;
  float j = (ll.y + 0.5 * PI) / PI * (uGridSize.y - 1.0);
  vec2 uv = vec2((floor(i + 0.5) + 0.5) / uGridSize.x, (floor(j + 0.5) + 0.5) / uGridSize.y);
  int c = int(texture2D(uClasses, uv).r * 255.0 + 0.5);
  vec4 col = c == 1 ? uColors[1] : c == 2 ? uColors[2] : c == 3 ? uColors[3] : uColors[0];
  if (col.a < 0.01) discard;
  vec3 n = normalize(vWorldPos);
  float ndl = dot(n, normalize(uLightDir)) * 0.5 + 0.5;
  float shade = mix(1.0, 0.5 + 0.5 * ndl * ndl, uShadeStrength);
  gl_FragColor = vec4(col.rgb * shade, col.a);
}
`;

/** Class colours: continent grey as diagnose_run's maps (gray75), belts orange as the
 *  R5 figures. Cores are the continent colour; the belts between them are what shows. */
export const CONTINENT_COLORS = {
  continent: '#bfbfbf',
  core: '#bfbfbf',
  belt: '#ff9a1f',
};

function rgba(hex: string, a = 1): Vector4 {
  const n = parseInt(hex.slice(1), 16);
  return new Vector4(((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255, a);
}

export class ContinentLayer {
  readonly mesh: Mesh;
  private readonly material: ShaderMaterial;
  private texture: DataTexture | null = null;
  private mode: ProjectionMode;
  private showBelts = true;
  private hasData = false;
  private enabled = true;

  constructor(mode: ProjectionMode = 'globe') {
    this.mode = mode;
    this.material = new ShaderMaterial({
      vertexShader: VERT,
      fragmentShader: FRAG,
      transparent: true,
      uniforms: {
        uClasses: { value: null },
        uGridSize: { value: [720, 361] },
        uColors: { value: [rgba('#000000', 0), rgba('#000000', 0), rgba('#000000', 0), rgba('#000000', 0)] },
        uProjectionMode: { value: PROJECTION_UNIFORM[mode] },
        uLightDir: { value: LIGHT_DIR.clone() },
        uShadeStrength: { value: isFlat(mode) ? 0 : 1 },
      },
    });
    this.applyColors();
    this.mesh = new Mesh(this.buildGeometry(mode), this.material);
    // Above the ocean (1), below the graticule and the overlay canvases.
    this.mesh.renderOrder = 2;
    this.mesh.frustumCulled = false;
    this.mesh.visible = false;
  }

  private buildGeometry(mode: ProjectionMode): BufferGeometry {
    return createSurfaceGeometry(mode, R_SURFACE * 1.0006);
  }

  private applyColors(): void {
    const c = this.material.uniforms.uColors.value as Vector4[];
    c[1] = rgba(CONTINENT_COLORS.continent);
    c[2] = rgba(CONTINENT_COLORS.core);
    c[3] = rgba(this.showBelts ? CONTINENT_COLORS.belt : CONTINENT_COLORS.core);
  }

  set belts(v: boolean) { this.showBelts = v; this.applyColors(); }
  set visible(v: boolean) { this.enabled = v; this.mesh.visible = v && this.hasData; }

  /** One step's class bytes, nlon x nlat, longitude fastest. */
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
      this.material.uniforms.uClasses.value = this.texture;
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
