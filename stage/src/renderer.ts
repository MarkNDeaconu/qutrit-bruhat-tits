/**
 * renderer.ts — pixi.js v8 Poincaré-disk scene.
 *
 * The ball to depth 8 (13,121 vertices) is generated procedurally from
 * addresses (tree.ts) and laid out by the contract's pure layout function
 * (hyperbolic.ts). The camera is a Möbius transform: panning to a vertex
 * animates the parameter c along the hyperbolic segment to its disk point.
 * Implements the StagePort consumed by walk.ts.
 */

import {
  Application,
  Container,
  Graphics,
  Particle,
  ParticleContainer,
  Sprite,
  Texture,
} from 'pixi.js';
import {
  type Complex,
  type Mobius,
  addressToDisk,
  mobiusApply,
  poincareSegmentPoints,
} from './hyperbolic.js';
import { generateBall, parentOf } from './tree.js';
import { MAX_DEPTH, type StagePort } from './walk.js';

// ---------------------------------------------------------------------------

export const BALL_DEPTH = 8;
const ARC_DEPTH = 4; // hyperbolic arcs up to this child depth, lines beyond
const ARC_POINTS = 12;
const ANIM_EDGE_DETAIL = 6; // edge LOD while the camera animates
const HIT_RADIUS_PX = 12;

const COLOR_PURE = 0xd55e00; // vermilion (Okabe–Ito)
const COLOR_ALT = 0x0072b2; // azure
const COLOR_GOLD = 0xe5b80b; // origin halo

interface Theme {
  bg: number;
  edge: number;
  edgeAlpha: number;
  boundary: number;
  boundaryAlpha: number;
  trail: number;
  pulse: number;
  flash: number;
  endpoint: number; // green marker for the current walk endpoint
  cutoff: number; // red marker where a path leaves the renderable range
  additive: boolean;
}

const DARK: Theme = {
  bg: 0x0b0e1a,
  edge: 0xffffff,
  edgeAlpha: 0.18,
  boundary: 0xffffff,
  boundaryAlpha: 0.22,
  trail: 0xffd166,
  pulse: 0xfff3c4,
  flash: 0xffffff,
  endpoint: 0x2ee36e, // vivid green
  cutoff: 0xff3b3b, // red
  additive: true,
};

const LIGHT: Theme = {
  bg: 0xf5f3ee,
  edge: 0x223047,
  edgeAlpha: 0.28,
  boundary: 0x223047,
  boundaryAlpha: 0.35,
  trail: 0xc9920a,
  pulse: 0xa86f00,
  flash: 0x7a5800,
  endpoint: 0x0a8f44, // darker green for light bg
  cutoff: 0xc81e1e, // darker red for light bg
  additive: false,
};

interface VRec {
  addr: string;
  depth: number;
  z0: Complex; // layout disk coordinate (pure function of the address)
  zx: number; // camera-transformed disk coordinate
  zy: number;
  particle: Particle;
}

interface Tween {
  start: number;
  dur: number;
  update: (t: number) => void;
  resolve: () => void;
  done: boolean;
}

export const easeInOut = (t: number): number =>
  t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;

function vertexRadius(depth: number): number {
  return Math.max(7 * Math.pow(0.82, depth), 1.2);
}

// ---------------------------------------------------------------------------

export class Renderer implements StagePort {
  app!: Application;

  private world!: Container;
  private boundaryG!: Graphics;
  private edgesG!: Graphics;
  private trailG!: Graphics;
  private liveTrailG!: Graphics; // partial segment revealed under a moving pulse
  private verticesPC!: ParticleContainer;
  private overlayG!: Container; // halos, compass glyphs, lobe flashes, pulses
  private haloG!: Graphics;

  private circleTex!: Texture;
  private glowTex!: Texture;

  private records: VRec[] = [];
  private byAddr = new Map<string, number>();
  private edgesByDepth: Array<Array<[number, number]>> = [];
  private maxDepthPresent = BALL_DEPTH;

  private cam: Mobius = { a: 0, c: { re: 0, im: 0 } };
  private camIsIdentity = true;
  private camAnimating = false;

  private zoom = 1;
  private panX = 0;
  private panY = 0;
  private R0 = 400;
  private cx = 0;
  private cy = 0;

  private theme: Theme = DARK;
  private tweens: Tween[] = [];
  private introStart = -1;
  private introDone = false;
  private lodCutoff = Infinity;

  private trailEdges: Array<[string, string]> = []; // active (pending) meander
  private frozenTrailEdges: Array<[string, string]> = []; // synthesized geodesics
  private selectedAddr: string | null = null;
  private endpointAddr: string | null = null;
  private endpointCutoff = false; // draw the endpoint red (out of range) vs green

  onVertexClick: ((addr: string) => void) | null = null;
  onVertexDblClick: ((addr: string) => void) | null = null;

  // -- construction ----------------------------------------------------------

  static async create(host: HTMLElement): Promise<Renderer> {
    const r = new Renderer();
    await r.init(host);
    return r;
  }

  private async init(host: HTMLElement): Promise<void> {
    this.app = new Application();
    await this.app.init({
      antialias: true,
      resolution: window.devicePixelRatio || 1,
      autoDensity: true,
      background: DARK.bg,
      resizeTo: window,
    });
    host.appendChild(this.app.canvas);

    this.world = new Container();
    this.boundaryG = new Graphics();
    this.edgesG = new Graphics();
    this.trailG = new Graphics();
    this.liveTrailG = new Graphics();
    this.verticesPC = new ParticleContainer({
      dynamicProperties: { position: true, color: true, rotation: false, uvs: false, vertex: false },
    });
    this.overlayG = new Container();
    this.haloG = new Graphics();
    this.overlayG.addChild(this.haloG);

    this.world.addChild(this.boundaryG, this.edgesG, this.trailG, this.liveTrailG, this.verticesPC, this.overlayG);
    this.app.stage.addChild(this.world);

    this.makeTextures();
    this.buildBall();
    this.layoutViewport();
    this.refreshPositions();
    this.redrawEdges(this.maxDepthPresent);
    this.redrawBoundary();
    this.setThemeInternal(DARK);

    this.introStart = performance.now();
    this.app.ticker.add(() => this.tick());
    this.bindInput();
    window.addEventListener('resize', () => this.handleResize());
  }

  private makeTextures(): void {
    const g = new Graphics();
    g.circle(32, 32, 30).fill({ color: 0xffffff });
    this.circleTex = this.app.renderer.generateTexture({ target: g, resolution: 2 });
    g.destroy();

    const cv = document.createElement('canvas');
    cv.width = cv.height = 128;
    const ctx = cv.getContext('2d')!;
    const grad = ctx.createRadialGradient(64, 64, 0, 64, 64, 64);
    grad.addColorStop(0, 'rgba(255,255,255,1)');
    grad.addColorStop(0.35, 'rgba(255,255,255,0.55)');
    grad.addColorStop(1, 'rgba(255,255,255,0)');
    ctx.fillStyle = grad;
    ctx.fillRect(0, 0, 128, 128);
    this.glowTex = Texture.from(cv);
  }

  private buildBall(): void {
    this.edgesByDepth = [];
    for (let d = 0; d <= BALL_DEPTH; d++) this.edgesByDepth.push([]);
    for (const addr of generateBall(BALL_DEPTH)) this.addVertexRecord(addr);
  }

  private addVertexRecord(addr: string): number {
    const depth = addr.length;
    const z0 = addressToDisk(addr);
    const r = vertexRadius(depth);
    const s = r / 32;
    const particle = new Particle({
      texture: this.circleTex,
      anchorX: 0.5,
      anchorY: 0.5,
      scaleX: s,
      scaleY: s,
      tint: depth % 2 === 0 ? COLOR_PURE : COLOR_ALT,
      alpha: 0,
    });
    const rec: VRec = { addr, depth, z0, zx: z0.re, zy: z0.im, particle };
    const idx = this.records.length;
    this.records.push(rec);
    this.byAddr.set(addr, idx);
    this.verticesPC.addParticle(particle);
    if (depth > 0) {
      const p = parentOf(addr)!;
      const pi = this.byAddr.get(p);
      if (pi !== undefined) {
        while (this.edgesByDepth.length <= depth) this.edgesByDepth.push([]);
        this.edgesByDepth[depth].push([pi, idx]);
      }
      if (depth > this.maxDepthPresent) this.maxDepthPresent = depth;
    }
    return idx;
  }

  /** Add a vertex (and its ancestors) outside the prebuilt ball. */
  private ensureVertex(addr: string): number {
    const have = this.byAddr.get(addr);
    if (have !== undefined) return have;
    if (addr !== '') this.ensureVertex(addr.slice(0, -1));
    const idx = this.addVertexRecord(addr);
    const rec = this.records[idx];
    this.projectRecord(rec);
    rec.particle.alpha = this.targetAlpha(rec);
    return idx;
  }

  // -- viewport / camera ------------------------------------------------------

  private layoutViewport(): void {
    const w = this.app.screen.width;
    const h = this.app.screen.height;
    this.cx = w / 2;
    this.cy = h / 2;
    this.R0 = (Math.min(w, h) / 2) * 0.96;
    this.applyWorldTransform();
  }

  private applyWorldTransform(): void {
    this.world.position.set(this.cx + this.panX, this.cy + this.panY);
    this.world.scale.set(this.zoom);
  }

  private handleResize(): void {
    requestAnimationFrame(() => {
      this.layoutViewport();
      this.refreshPositions();
      this.redrawEdges(this.maxDepthPresent);
      this.redrawBoundary();
      this.redrawTrail();
    });
  }

  private projectRecord(rec: VRec): void {
    if (this.camIsIdentity) {
      rec.zx = rec.z0.re;
      rec.zy = rec.z0.im;
    } else {
      const z = mobiusApply(this.cam, rec.z0);
      rec.zx = z.re;
      rec.zy = z.im;
    }
    rec.particle.x = rec.zx * this.R0;
    rec.particle.y = rec.zy * this.R0;
  }

  private refreshPositions(): void {
    for (const rec of this.records) this.projectRecord(rec);
  }

  /** Animate the Möbius parameter c along the hyperbolic segment to addr. */
  focusOn(addr: string, durationMs = 700): Promise<void> {
    if (addr.length > MAX_DEPTH) return Promise.resolve(); // out of range: don't glitch
    const idx = this.ensureVertex(addr);
    const target = this.records[idx].z0;
    const path = poincareSegmentPoints(this.cam.c, target, 64);
    const pan0x = this.panX;
    const pan0y = this.panY;
    this.camAnimating = true;
    return this.animate(durationMs, (tRaw) => {
      const t = easeInOut(tRaw);
      const f = t * (path.length - 1);
      const i = Math.min(Math.floor(f), path.length - 2);
      const frac = f - i;
      const c = {
        re: path[i].re + (path[i + 1].re - path[i].re) * frac,
        im: path[i].im + (path[i + 1].im - path[i].im) * frac,
      };
      this.cam = { a: 0, c };
      this.camIsIdentity = Math.abs(c.re) < 1e-12 && Math.abs(c.im) < 1e-12;
      this.panX = pan0x * (1 - t);
      this.panY = pan0y * (1 - t);
      this.applyWorldTransform();
      this.refreshPositions();
      this.redrawEdges(Math.min(ANIM_EDGE_DETAIL, this.maxDepthPresent));
      this.redrawTrail();
    }).then(() => {
      this.cam = { a: 0, c: { re: target.re, im: target.im } };
      this.camIsIdentity = target.re === 0 && target.im === 0;
      this.camAnimating = false;
      this.refreshPositions();
      this.redrawEdges(this.maxDepthPresent);
      this.redrawTrail();
    });
  }

  /** Reset zoom and glide the camera to centre `addr`. */
  recenterOn(addr: string, durationMs = 700): Promise<void> {
    this.zoom = 1;
    this.applyWorldTransform();
    this.updateLod();
    return this.focusOn(addr, durationMs);
  }

  resetCamera(durationMs = 700): Promise<void> {
    return this.recenterOn('', durationMs);
  }

  // -- drawing ----------------------------------------------------------------

  private redrawBoundary(): void {
    this.boundaryG.clear();
    this.boundaryG
      .circle(0, 0, this.R0)
      .stroke({ width: 1.2, color: this.theme.boundary, alpha: this.theme.boundaryAlpha });
    // faint inner shell ring at depth 1 radius: ambient bipartite rhythm
    this.boundaryG
      .circle(0, 0, Math.tanh(Math.log(3) / 2) * this.R0)
      .stroke({ width: 0.6, color: this.theme.boundary, alpha: this.theme.boundaryAlpha * 0.3 });
  }

  private redrawEdges(maxDepth: number): void {
    const g = this.edgesG;
    g.clear();
    for (let d = 1; d < this.edgesByDepth.length && d <= maxDepth; d++) {
      const list = this.edgesByDepth[d];
      if (list.length === 0) continue;
      for (const [pi, ci] of list) {
        const a = this.records[pi];
        const b = this.records[ci];
        if (d <= ARC_DEPTH) {
          const pts = poincareSegmentPoints(
            { re: a.zx, im: a.zy },
            { re: b.zx, im: b.zy },
            ARC_POINTS,
          );
          g.moveTo(pts[0].re * this.R0, pts[0].im * this.R0);
          for (let i = 1; i < pts.length; i++) g.lineTo(pts[i].re * this.R0, pts[i].im * this.R0);
        } else {
          g.moveTo(a.zx * this.R0, a.zy * this.R0);
          g.lineTo(b.zx * this.R0, b.zy * this.R0);
        }
      }
      const width = Math.max(0.35, 2.4 * Math.pow(0.8, d));
      g.stroke({ width, color: this.theme.edge, alpha: this.theme.edgeAlpha, cap: 'round', join: 'round' });
    }
  }

  /** Polyline (px coords) along the hyperbolic edge a—b under the camera. */
  edgePolylinePx(a: string, b: string, n = ARC_POINTS): Array<{ x: number; y: number }> {
    const ra = this.records[this.ensureVertex(a)];
    const rb = this.records[this.ensureVertex(b)];
    const pts = poincareSegmentPoints({ re: ra.zx, im: ra.zy }, { re: rb.zx, im: rb.zy }, n);
    return pts.map((p) => ({ x: p.re * this.R0, y: p.im * this.R0 }));
  }

  private strokeTrailPass(
    g: Graphics,
    polys: Array<Array<{ x: number; y: number }>>,
    width: number,
    alpha: number,
    color: number,
  ): void {
    for (const pts of polys) {
      g.moveTo(pts[0].x, pts[0].y);
      for (let i = 1; i < pts.length; i++) g.lineTo(pts[i].x, pts[i].y);
    }
    g.stroke({ width, color, alpha, cap: 'round', join: 'round' });
  }

  private redrawTrail(): void {
    const g = this.trailG;
    g.clear();
    // frozen (already-synthesized geodesics) + active (current pending meander)
    const all = this.frozenTrailEdges.concat(this.trailEdges);
    if (all.length === 0) return;
    const polys = all.map(([a, b]) => this.edgePolylinePx(a, b));
    const c = this.theme.trail;
    // layered strokes: wide low-alpha under narrow high-alpha = cheap bloom
    this.strokeTrailPass(g, polys, 9, 0.14, c);
    this.strokeTrailPass(g, polys, 5, 0.4, c);
    this.strokeTrailPass(g, polys, 3, 0.95, c);
  }

  // -- StagePort --------------------------------------------------------------

  drawDeepPath(addrs: string[]): void {
    let added = false;
    for (const a of addrs) {
      if (a.length > MAX_DEPTH) continue; // beyond the renderable range
      if (!this.byAddr.has(a)) {
        this.ensureVertex(a);
        added = true;
      }
    }
    if (added) this.redrawEdges(this.maxDepthPresent);
  }

  private inRange(e: [string, string]): boolean {
    return e[0].length <= MAX_DEPTH && e[1].length <= MAX_DEPTH;
  }

  setTrailEdges(edges: Array<[string, string]>): void {
    this.trailEdges = edges.filter((e) => this.inRange(e));
    this.redrawTrail();
  }

  clearTrail(): void {
    this.trailEdges = [];
    this.liveTrailG.clear();
    this.redrawTrail();
  }

  /** Frozen geodesics = the parts already synthesized; they persist across
   *  appends and further synthesizes (never re-expanded). */
  addFrozenEdges(edges: Array<[string, string]>): void {
    const inRange = edges.filter((e) => this.inRange(e));
    for (const [a, b] of inRange) {
      this.ensureVertex(a);
      this.ensureVertex(b);
    }
    this.frozenTrailEdges.push(...inRange);
    this.redrawTrail();
  }

  clearFrozen(): void {
    this.frozenTrailEdges = [];
    this.redrawTrail();
  }

  appendTrailEdge(a: string, b: string): void {
    this.trailEdges.push([a, b]);
    this.liveTrailG.clear();
    this.redrawTrail();
  }

  /** S/R: 4-spoke compass glyph rotating 120° + a halo pulse at the vertex. */
  compassSpin(addr: string, ms: number): Promise<void> {
    const idx = this.ensureVertex(addr);
    const rec = this.records[idx];
    const glyph = new Graphics();
    const L = 16;
    for (let k = 0; k < 4; k++) {
      const ang = (k * Math.PI) / 2;
      glyph.moveTo(Math.cos(ang) * 5, Math.sin(ang) * 5);
      glyph.lineTo(Math.cos(ang) * L, Math.sin(ang) * L);
    }
    glyph.stroke({ width: 2, color: this.theme.pulse, alpha: 0.95, cap: 'round' });
    glyph.circle(0, 0, 5).stroke({ width: 1.5, color: this.theme.pulse, alpha: 0.8 });
    if (this.theme.additive) glyph.blendMode = 'add';
    const halo = new Graphics();
    if (this.theme.additive) halo.blendMode = 'add';
    this.overlayG.addChild(halo, glyph);
    return this.animate(ms, (tRaw) => {
      const t = easeInOut(tRaw);
      glyph.position.set(rec.zx * this.R0, rec.zy * this.R0);
      glyph.rotation = (t * 2 * Math.PI) / 3; // 120°: monomials permute branches
      glyph.alpha = tRaw < 0.15 ? tRaw / 0.15 : tRaw > 0.85 ? (1 - tRaw) / 0.15 : 1;
      halo.clear();
      halo
        .circle(rec.zx * this.R0, rec.zy * this.R0, 8 + 10 * t)
        .stroke({ width: 2, color: this.theme.trail, alpha: 0.5 * (1 - t) });
    }).then(() => {
      glyph.destroy();
      halo.destroy();
    });
  }

  /** H: a glowing pulse travels from → mid → to along the drawn edge arcs. */
  movePulse(from: string, mid: string, to: string, ms: number): Promise<void> {
    this.ensureVertex(from);
    this.ensureVertex(mid);
    this.ensureVertex(to);
    const glow = new Sprite({ texture: this.glowTex, anchor: 0.5 });
    glow.scale.set(0.55);
    glow.tint = this.theme.trail;
    glow.alpha = 0.85;
    const core = new Sprite({ texture: this.glowTex, anchor: 0.5 });
    core.scale.set(0.22);
    core.tint = this.theme.pulse;
    if (this.theme.additive) {
      glow.blendMode = 'add';
      core.blendMode = 'add';
    }
    this.overlayG.addChild(glow, core);
    let flourished = false;
    return this.animate(ms, (tRaw) => {
      const t = easeInOut(tRaw);
      // recompute each frame: the camera may be animating underneath us
      const p1 = this.edgePolylinePx(from, mid);
      const p2 = this.edgePolylinePx(mid, to);
      const pts = p1.concat(p2.slice(1));
      const f = t * (pts.length - 1);
      const i = Math.min(Math.floor(f), pts.length - 2);
      const frac = f - i;
      const x = pts[i].x + (pts[i + 1].x - pts[i].x) * frac;
      const y = pts[i].y + (pts[i + 1].y - pts[i].y) * frac;
      glow.position.set(x, y);
      core.position.set(x, y);
      // progressive reveal of the segment behind the pulse
      const lg = this.liveTrailG;
      lg.clear();
      if (i >= 1 || frac > 0) {
        lg.moveTo(pts[0].x, pts[0].y);
        for (let k = 1; k <= i; k++) lg.lineTo(pts[k].x, pts[k].y);
        lg.lineTo(x, y);
        lg.stroke({ width: 3, color: this.theme.trail, alpha: 0.9, cap: 'round', join: 'round' });
      }
      if (!flourished && tRaw >= 0.5) {
        flourished = true;
        this.ringBurst(mid, 26, 320); // 180° flourish marker at the midpoint
      }
    }).then(() => {
      glow.destroy();
      core.destroy();
      this.liveTrailG.clear();
    });
  }

  /** Expanding ring at a vertex (midpoint flourish / triumph accents). */
  private ringBurst(addr: string, maxR: number, ms: number, color?: number): void {
    const idx = this.byAddr.get(addr);
    if (idx === undefined) return;
    const rec = this.records[idx];
    const ring = new Graphics();
    if (this.theme.additive) ring.blendMode = 'add';
    this.overlayG.addChild(ring);
    void this.animate(ms, (t) => {
      ring.clear();
      ring
        .circle(rec.zx * this.R0, rec.zy * this.R0, 6 + (maxR - 6) * t)
        .stroke({ width: 2.5 * (1 - t) + 0.5, color: color ?? this.theme.pulse, alpha: 0.8 * (1 - t) });
    }).then(() => ring.destroy());
  }

  /** Camera follow: refocus if the vertex strays past ~70% of viewport radius. */
  maybeFollow(addr: string): void {
    if (addr.length > MAX_DEPTH) return; // never centre on an out-of-range vertex
    const idx = this.byAddr.get(addr);
    if (idx === undefined || this.camAnimating) return;
    const rec = this.records[idx];
    const sx = this.cx + this.panX + rec.zx * this.R0 * this.zoom;
    const sy = this.cy + this.panY + rec.zy * this.R0 * this.zoom;
    const d = Math.hypot(sx - this.cx, sy - this.cy);
    if (d > 0.7 * (this.R0 / 0.96)) void this.focusOn(addr, 600);
  }

  /** Straightening: flash the backtrack lobe bright, then fade it out. */
  flashLobe(lobe: [string, string, string], ms: number): Promise<void> {
    const [a, m, b] = lobe;
    const g = new Graphics();
    if (this.theme.additive) g.blendMode = 'add';
    this.overlayG.addChild(g);
    const draw = (alpha: number, width: number) => {
      g.clear();
      const p1 = this.edgePolylinePx(a, m);
      const p2 = this.edgePolylinePx(m, b);
      for (const pts of [p1, p2]) {
        g.moveTo(pts[0].x, pts[0].y);
        for (let i = 1; i < pts.length; i++) g.lineTo(pts[i].x, pts[i].y);
      }
      g.stroke({ width, color: this.theme.flash, alpha, cap: 'round', join: 'round' });
    };
    return this.animate(ms, (t) => {
      // quick flash up, then burn off
      const alpha = t < 0.3 ? t / 0.3 : 1 - (t - 0.3) / 0.7;
      draw(Math.max(alpha, 0) * 0.95, 6 - 3 * t);
    }).then(() => g.destroy());
  }

  /** The geodesic redrawn in a triumphant pulse: width 4, full gold. */
  triumph(geodesic: string[], ms: number): Promise<void> {
    const g = new Graphics();
    if (this.theme.additive) g.blendMode = 'add';
    this.overlayG.addChild(g);
    this.ringBurst('', 40, ms, this.theme.trail);
    const last = geodesic[geodesic.length - 1];
    if (last !== undefined && last !== '') this.ringBurst(last, 32, ms, this.theme.trail);
    return this.animate(ms, (tRaw) => {
      const t = easeInOut(tRaw);
      g.clear();
      if (geodesic.length < 2) return;
      const polys: Array<Array<{ x: number; y: number }>> = [];
      for (let i = 0; i + 1 < geodesic.length; i++) {
        polys.push(this.edgePolylinePx(geodesic[i], geodesic[i + 1]));
      }
      const swell = Math.sin(Math.PI * t);
      this.strokeTrailPass(g, polys, 12 + 10 * swell, 0.18 * (1 - 0.4 * t), this.theme.trail);
      this.strokeTrailPass(g, polys, 4 + 3 * swell, 1, this.theme.trail);
    }).then(() => g.destroy());
  }

  /**
   * Play a pulse that travels along `addrs` (the synthesized geodesic),
   * revealing the gold trail behind it and committing the clean path at the
   * end. The camera is static during this (synthesize awaits resetCamera
   * first), so screen coords are precomputed once.
   */
  tracePath(addrs: string[], msPerEdge = 320): Promise<void> {
    if (addrs.length < 2) {
      this.setTrailEdges([]);
      return Promise.resolve();
    }
    // disk-space (camera-independent) polyline through the geodesic vertices,
    // so the pulse position is well-defined even as the camera moves
    const verts = addrs.map((a) => this.records[this.ensureVertex(a)].z0);
    let diskPts: Complex[] = [];
    for (let i = 0; i + 1 < verts.length; i++) {
      const seg = poincareSegmentPoints(verts[i], verts[i + 1], ARC_POINTS);
      diskPts = i === 0 ? seg.slice() : diskPts.concat(seg.slice(1));
    }
    const edges: Array<[string, string]> = [];
    for (let i = 0; i + 1 < addrs.length; i++) edges.push([addrs[i], addrs[i + 1]]);
    this.setTrailEdges(edges); // base geodesic, redrawn each frame under the moving camera

    const glow = new Sprite({ texture: this.glowTex, anchor: 0.5 });
    glow.scale.set(0.7);
    glow.tint = this.theme.trail;
    glow.alpha = 0.95;
    const core = new Sprite({ texture: this.glowTex, anchor: 0.5 });
    core.scale.set(0.28);
    core.tint = this.theme.pulse;
    if (this.theme.additive) {
      glow.blendMode = 'add';
      core.blendMode = 'add';
    }
    this.overlayG.addChild(glow, core);

    const wp = (z: Complex): { x: number; y: number } => {
      const m = mobiusApply(this.cam, z);
      return { x: m.re * this.R0, y: m.im * this.R0 };
    };

    const dur = Math.max(msPerEdge * edges.length, 480);
    this.camAnimating = true; // suppress maybeFollow while we drive the camera
    return this.animate(dur, (t) => {
      const f = t * (diskPts.length - 1);
      const i = Math.min(Math.floor(f), diskPts.length - 2);
      const frac = f - i;
      const zp = {
        re: diskPts[i].re + (diskPts[i + 1].re - diskPts[i].re) * frac,
        im: diskPts[i].im + (diskPts[i + 1].im - diskPts[i].im) * frac,
      };
      // center the camera on the pulse → the camera follows it
      this.cam = { a: 0, c: zp };
      this.camIsIdentity = Math.abs(zp.re) < 1e-12 && Math.abs(zp.im) < 1e-12;
      this.refreshPositions();
      this.redrawEdges(Math.min(ANIM_EDGE_DETAIL, this.maxDepthPresent));
      this.redrawTrail();
      // bright reveal of the traversed portion, in current screen coords
      const head = wp(zp); // ≈ (0,0): the pulse sits at the centred camera focus
      const lg = this.liveTrailG;
      lg.clear();
      const drawTraversed = (width: number, alpha: number) => {
        const p0 = wp(diskPts[0]);
        lg.moveTo(p0.x, p0.y);
        for (let k = 1; k <= i; k++) {
          const p = wp(diskPts[k]);
          lg.lineTo(p.x, p.y);
        }
        lg.lineTo(head.x, head.y);
        lg.stroke({ width, color: this.theme.trail, alpha, cap: 'round', join: 'round' });
      };
      drawTraversed(11, 0.16);
      drawTraversed(4, 1);
      glow.position.set(head.x, head.y);
      core.position.set(head.x, head.y);
    }).then(() => {
      glow.destroy();
      core.destroy();
      this.liveTrailG.clear();
      // settle: camera centred on the green endpoint, where the pulse arrived
      const endZ = verts[verts.length - 1];
      this.cam = { a: 0, c: { re: endZ.re, im: endZ.im } };
      this.camIsIdentity = Math.abs(endZ.re) < 1e-12 && Math.abs(endZ.im) < 1e-12;
      this.camAnimating = false;
      this.refreshPositions();
      this.redrawEdges(this.maxDepthPresent);
      this.setTrailEdges(edges); // commit the full clean geodesic
      this.ringBurst(addrs[addrs.length - 1], 32, 700, this.theme.trail);
    });
  }

  // -- selection / theme ------------------------------------------------------

  selectVertex(addr: string | null): void {
    this.selectedAddr = addr;
  }

  /** Mark (or clear) the current walk endpoint with a pulsating node — green
   *  normally, red when `cutoff` (the path ran past the renderable range). */
  setEndpoint(addr: string | null, cutoff = false): void {
    if (addr !== null) this.ensureVertex(addr);
    this.endpointAddr = addr;
    this.endpointCutoff = cutoff;
  }

  get endpoint(): string | null {
    return this.endpointAddr;
  }

  setLightTheme(light: boolean): void {
    this.setThemeInternal(light ? LIGHT : DARK);
  }

  private setThemeInternal(t: Theme): void {
    this.theme = t;
    this.app.renderer.background.color = t.bg;
    this.trailG.blendMode = t.additive ? 'add' : 'normal';
    this.liveTrailG.blendMode = t.additive ? 'add' : 'normal';
    this.haloG.blendMode = t.additive ? 'add' : 'normal';
    this.redrawBoundary();
    this.redrawEdges(this.maxDepthPresent);
    this.redrawTrail();
  }

  // -- per-frame --------------------------------------------------------------

  private targetAlpha(rec: VRec): number {
    return rec.depth <= this.lodCutoff ? 1 : 0;
  }

  private updateLod(): void {
    // hide vertices whose projected radius would be < 0.8 px
    let cutoff = Infinity;
    for (let d = 0; d <= this.maxDepthPresent; d++) {
      if (vertexRadius(d) * this.zoom < 0.8) {
        cutoff = d - 1;
        break;
      }
    }
    if (cutoff !== this.lodCutoff) {
      this.lodCutoff = cutoff;
      if (this.introDone) {
        for (const rec of this.records) rec.particle.alpha = this.targetAlpha(rec);
      }
    }
  }

  private tick(): void {
    const now = performance.now();

    // tweens
    for (const tw of this.tweens) {
      const t = Math.min((now - tw.start) / tw.dur, 1);
      tw.update(t);
      if (t >= 1 && !tw.done) {
        tw.done = true;
        tw.resolve();
      }
    }
    this.tweens = this.tweens.filter((tw) => !tw.done);

    // staggered fade-in of the ball (~600 ms total, by depth)
    if (!this.introDone && this.introStart >= 0) {
      const el = now - this.introStart;
      let allDone = true;
      for (const rec of this.records) {
        const t = Math.min(Math.max((el - rec.depth * 55) / 200, 0), 1);
        rec.particle.alpha = t * this.targetAlpha(rec);
        if (t < 1) allDone = false;
      }
      if (allDone) this.introDone = true;
    }

    // origin halo: gentle idle breathing (gold) + selection ring
    const g = this.haloG;
    g.clear();
    const o = this.records[0];
    if (o) {
      const breathe = 1 + 0.18 * Math.sin((now / 2400) * 2 * Math.PI);
      const r = 11 * breathe;
      const x = o.zx * this.R0;
      const y = o.zy * this.R0;
      g.circle(x, y, r).stroke({ width: 2, color: COLOR_GOLD, alpha: 0.85 });
      g.circle(x, y, r + 4).stroke({ width: 5, color: COLOR_GOLD, alpha: 0.16 });
    }
    // current walk endpoint: pulsating green node (the target the geodesic reaches)
    if (this.endpointAddr !== null) {
      const i = this.byAddr.get(this.endpointAddr);
      if (i !== undefined) {
        const rec = this.records[i];
        const x = rec.zx * this.R0;
        const y = rec.zy * this.R0;
        const base = Math.max(vertexRadius(rec.depth), 5);
        const pulse = 1 + 0.28 * Math.sin((now / 620) * 2 * Math.PI);
        const c = this.endpointCutoff ? this.theme.cutoff : this.theme.endpoint;
        g.circle(x, y, base + 2).fill({ color: c, alpha: 0.95 });
        g.circle(x, y, (base + 6) * pulse).stroke({ width: 2.5, color: c, alpha: 0.9 });
        g.circle(x, y, (base + 12) * pulse).stroke({ width: 6, color: c, alpha: 0.18 });
      }
    }
    if (this.selectedAddr !== null) {
      const i = this.byAddr.get(this.selectedAddr);
      if (i !== undefined) {
        const rec = this.records[i];
        g.circle(rec.zx * this.R0, rec.zy * this.R0, vertexRadius(rec.depth) + 5).stroke({
          width: 1.5,
          color: this.theme.pulse,
          alpha: 0.9,
        });
      }
    }
  }

  private animate(dur: number, update: (t: number) => void): Promise<void> {
    return new Promise((resolve) => {
      this.tweens.push({ start: performance.now(), dur: Math.max(dur, 16), update, resolve, done: false });
    });
  }

  // -- input ------------------------------------------------------------------

  private hitTest(sx: number, sy: number): string | null {
    const wx = (sx - this.cx - this.panX) / this.zoom;
    const wy = (sy - this.cy - this.panY) / this.zoom;
    const thresh = HIT_RADIUS_PX / this.zoom;
    let best: string | null = null;
    let bestD = thresh;
    for (const rec of this.records) {
      if (rec.depth > this.lodCutoff) continue;
      const d = Math.hypot(rec.zx * this.R0 - wx, rec.zy * this.R0 - wy);
      if (d < bestD) {
        bestD = d;
        best = rec.addr;
      }
    }
    return best;
  }

  private bindInput(): void {
    const cv = this.app.canvas;
    cv.style.touchAction = 'none';

    cv.addEventListener('wheel', (e) => {
      e.preventDefault();
      const factor = Math.exp(-e.deltaY * 0.0014);
      const newZoom = Math.min(Math.max(this.zoom * factor, 0.15), 80);
      const sx = e.offsetX;
      const sy = e.offsetY;
      const wx = (sx - this.cx - this.panX) / this.zoom;
      const wy = (sy - this.cy - this.panY) / this.zoom;
      this.panX = sx - this.cx - wx * newZoom;
      this.panY = sy - this.cy - wy * newZoom;
      this.zoom = newZoom;
      this.applyWorldTransform();
      this.updateLod();
    }, { passive: false });

    let dragging = false;
    let moved = 0;
    let lastX = 0;
    let lastY = 0;
    cv.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      dragging = true;
      moved = 0;
      lastX = e.offsetX;
      lastY = e.offsetY;
      cv.setPointerCapture(e.pointerId);
    });
    cv.addEventListener('pointermove', (e) => {
      if (!dragging) return;
      const dx = e.offsetX - lastX;
      const dy = e.offsetY - lastY;
      moved += Math.abs(dx) + Math.abs(dy);
      lastX = e.offsetX;
      lastY = e.offsetY;
      this.panX += dx;
      this.panY += dy;
      this.applyWorldTransform();
    });
    cv.addEventListener('pointerup', (e) => {
      if (!dragging) return;
      dragging = false;
      cv.releasePointerCapture(e.pointerId);
      if (moved < 5) {
        const addr = this.hitTest(e.offsetX, e.offsetY);
        if (addr !== null && this.onVertexClick) this.onVertexClick(addr);
      }
    });
    cv.addEventListener('dblclick', (e) => {
      const addr = this.hitTest(e.offsetX, e.offsetY);
      if (addr !== null) {
        if (this.onVertexDblClick) this.onVertexDblClick(addr);
        else void this.focusOn(addr, 700);
      }
    });
  }
}
