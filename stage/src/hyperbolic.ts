/**
 * hyperbolic.ts — pure Poincaré-disk math. NO pixi imports (unit-testable in node).
 *
 * Everything here is *display* geometry. All exact arithmetic lives in the
 * Python engine (DATA_CONTRACT.md); addresses are the only shared language.
 */

// ---------------------------------------------------------------------------
// Complex arithmetic
// ---------------------------------------------------------------------------

export interface Complex {
  re: number;
  im: number;
}

export const C = (re: number, im = 0): Complex => ({ re, im });
export const ZERO: Complex = { re: 0, im: 0 };

export const cadd = (a: Complex, b: Complex): Complex => ({ re: a.re + b.re, im: a.im + b.im });
export const csub = (a: Complex, b: Complex): Complex => ({ re: a.re - b.re, im: a.im - b.im });
export const cneg = (a: Complex): Complex => ({ re: -a.re, im: -a.im });
export const conj = (a: Complex): Complex => ({ re: a.re, im: -a.im });
export const cmul = (a: Complex, b: Complex): Complex => ({
  re: a.re * b.re - a.im * b.im,
  im: a.re * b.im + a.im * b.re,
});
export const cdiv = (a: Complex, b: Complex): Complex => {
  const d = b.re * b.re + b.im * b.im;
  return { re: (a.re * b.re + a.im * b.im) / d, im: (a.im * b.re - a.re * b.im) / d };
};
export const cabs2 = (a: Complex): number => a.re * a.re + a.im * a.im;
export const cabs = (a: Complex): number => Math.hypot(a.re, a.im);
export const expi = (t: number): Complex => ({ re: Math.cos(t), im: Math.sin(t) });

// ---------------------------------------------------------------------------
// Möbius transforms of the unit disk:  M(z) = e^{i a} (z − c) / (1 − conj(c) z)
// ---------------------------------------------------------------------------

export interface Mobius {
  a: number; // rotation angle
  c: Complex; // the point sent to 0, |c| < 1
}

export const MOBIUS_ID: Mobius = { a: 0, c: { re: 0, im: 0 } };

export function mobiusApply(m: Mobius, z: Complex): Complex {
  const num = csub(z, m.c);
  const den = csub(C(1), cmul(conj(m.c), z));
  return cmul(expi(m.a), cdiv(num, den));
}

/** Inverse transform. Derivation: c' = M(0) = −e^{ia} c, a' = −a. */
export function mobiusInvert(m: Mobius): Mobius {
  return { a: -m.a, c: cneg(cmul(expi(m.a), m.c)) };
}

/**
 * Composition f∘g (apply g first, then f), via the 2×2 matrix representation
 * M = [[e^{ia}, −e^{ia}c], [−conj(c), 1]]  (a unit-modulus multiple of an
 * SU(1,1) matrix, so the product is again of the same shape up to scalar).
 */
export function mobiusCompose(f: Mobius, g: Mobius): Mobius {
  const ef = expi(f.a);
  const eg = expi(g.a);
  // f matrix entries
  const fA = ef, fB = cneg(cmul(ef, f.c)), fC = cneg(conj(f.c)), fD = C(1);
  // g matrix entries
  const gA = eg, gB = cneg(cmul(eg, g.c)), gC = cneg(conj(g.c)), gD = C(1);
  // product P = F · G
  const p = cadd(cmul(fA, gA), cmul(fB, gC));
  const q = cadd(cmul(fA, gB), cmul(fB, gD));
  const s = cadd(cmul(fC, gB), cmul(fD, gD));
  // (pz+q)/(rz+s) = e^{ia'}(z−c')/(1−conj(c')z) with c' = −q/p, e^{ia'} = p/s
  const c = cneg(cdiv(q, p));
  const ratio = cdiv(p, s);
  return { a: Math.atan2(ratio.im, ratio.re), c };
}

// ---------------------------------------------------------------------------
// Address → layout (the contract's pure layout function)
// ---------------------------------------------------------------------------

export const LAMBDA = Math.log(3);
const TAU = 2 * Math.PI;

export interface Layout {
  theta: number;
  rho: number;
}

/**
 * Root interval [0, 2π); child c of a vertex splits the parent interval into
 * k equal parts (k = 4 at root, else 3) and takes the c-th; theta = interval
 * midpoint; rho = depth · ln 3. Iterative over the digits.
 */
export function addressToLayout(addr: string): Layout {
  let lo = 0;
  let width = TAU;
  for (let i = 0; i < addr.length; i++) {
    const k = i === 0 ? 4 : 3;
    const d = addr.charCodeAt(i) - 48;
    if (d < 0 || d >= k) {
      throw new Error(`invalid address digit '${addr[i]}' at position ${i} of '${addr}'`);
    }
    width /= k;
    lo += d * width;
  }
  return { theta: lo + width / 2, rho: addr.length * LAMBDA };
}

/** The half-open angular interval [lo, hi) owned by an address. */
export function addressInterval(addr: string): { lo: number; hi: number } {
  let lo = 0;
  let width = TAU;
  for (let i = 0; i < addr.length; i++) {
    const k = i === 0 ? 4 : 3;
    width /= k;
    lo += (addr.charCodeAt(i) - 48) * width;
  }
  return { lo, hi: lo + width };
}

/** Poincaré coordinate z = tanh(rho/2) e^{i theta}. Origin ('' / rho 0) → 0. */
export function diskPoint(theta: number, rho: number): Complex {
  const r = Math.tanh(rho / 2);
  return { re: r * Math.cos(theta), im: r * Math.sin(theta) };
}

export function addressToDisk(addr: string): Complex {
  const { theta, rho } = addressToLayout(addr);
  return diskPoint(theta, rho);
}

// ---------------------------------------------------------------------------
// Geodesic arcs (circles orthogonal to the unit circle)
// ---------------------------------------------------------------------------

/**
 * n points (n ≥ 2, endpoints included) along the hyperbolic geodesic between
 * two disk points. The geodesic is the arc of the circle through z1, z2
 * orthogonal to the unit circle; when z1, z2, 0 are (near-)collinear the
 * geodesic is the straight chord (a diameter), handled separately.
 */
export function poincareSegmentPoints(z1: Complex, z2: Complex, n: number): Complex[] {
  if (n < 2) n = 2;
  const pts: Complex[] = new Array(n);
  const x1 = z1.re, y1 = z1.im, x2 = z2.re, y2 = z2.im;
  const cross = x1 * y2 - y1 * x2;
  const straight = () => {
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
      pts[i] = { re: x1 + (x2 - x1) * t, im: y1 + (y2 - y1) * t };
    }
    pts[0] = { re: x1, im: y1 };
    pts[n - 1] = { re: x2, im: y2 };
    return pts;
  };
  // coincident or (near-)diametral → straight chord
  const d2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
  if (d2 < 1e-24 || Math.abs(cross) < 1e-9) return straight();

  // Orthogonality |w|^2 = R^2 + 1 plus |z_i − w| = R gives the linear system
  //   2 x_i wx + 2 y_i wy = |z_i|^2 + 1   (i = 1, 2)
  const b1 = x1 * x1 + y1 * y1 + 1;
  const b2 = x2 * x2 + y2 * y2 + 1;
  const det = 4 * cross;
  const wx = (b1 * 2 * y2 - b2 * 2 * y1) / det;
  const wy = (2 * x1 * b2 - 2 * x2 * b1) / det;
  const R2 = wx * wx + wy * wy - 1;
  if (R2 <= 0) return straight(); // numerical guard
  const R = Math.sqrt(R2);
  const phi1 = Math.atan2(y1 - wy, x1 - wx);
  const phi2 = Math.atan2(y2 - wy, x2 - wx);
  let dphi = phi2 - phi1;
  while (dphi > Math.PI) dphi -= TAU;
  while (dphi < -Math.PI) dphi += TAU;
  for (let i = 0; i < n; i++) {
    const t = i / (n - 1);
    const phi = phi1 + dphi * t;
    pts[i] = { re: wx + R * Math.cos(phi), im: wy + R * Math.sin(phi) };
  }
  pts[0] = { re: x1, im: y1 };
  pts[n - 1] = { re: x2, im: y2 };
  return pts;
}

// ---------------------------------------------------------------------------
// Tiny assertion suite (used by test/math.test.mjs and runnable in-browser)
// ---------------------------------------------------------------------------

export interface MathCheckResult {
  passed: number;
  failures: string[];
}

export function runMathChecks(): MathCheckResult {
  let passed = 0;
  const failures: string[] = [];
  const ok = (cond: boolean, label: string) => {
    if (cond) passed++;
    else failures.push(label);
  };
  const near = (a: number, b: number, eps: number) => Math.abs(a - b) <= eps;

  // determinism
  const l1 = addressToLayout('2101');
  const l2 = addressToLayout('2101');
  ok(l1.theta === l2.theta && l1.rho === l2.rho, 'addressToLayout deterministic');

  // root children at 4 distinct equally spaced quarter midpoints
  const thetas = ['0', '1', '2', '3'].map((a) => addressToLayout(a).theta);
  ok(
    thetas.every((t, i) => near(t, Math.PI / 4 + (i * Math.PI) / 2, 1e-12)),
    'root children at quarter midpoints',
  );

  // child theta inside parent interval
  const iv = addressInterval('21');
  const th = addressToLayout('212').theta;
  ok(th >= iv.lo && th < iv.hi, 'nested interval containment');

  // |diskPoint| < 1 (deep valid address: first digit 0–3, later digits 0–2)
  ok(cabs(addressToDisk('3' + '012210120121012')) < 1, 'disk point inside unit disk');

  // Möbius round trip
  const m: Mobius = { a: 0.7, c: C(0.31, -0.42) };
  const id = mobiusCompose(m, mobiusInvert(m));
  ok(cabs(id.c) < 1e-12 && Math.abs(Math.sin(id.a)) < 1e-12, 'mobius compose/invert round trip');

  // compose = sequential application
  const m2: Mobius = { a: -1.2, c: C(-0.2, 0.55) };
  const z = C(0.27, 0.13);
  const seq = mobiusApply(m, mobiusApply(m2, z));
  const comp = mobiusApply(mobiusCompose(m, m2), z);
  ok(cabs(csub(seq, comp)) < 1e-12, 'mobius composition agrees with sequential apply');

  // geodesic endpoints
  const a = C(0.5, 0.1);
  const b = C(-0.2, 0.6);
  const seg = poincareSegmentPoints(a, b, 12);
  ok(
    cabs(csub(seg[0], a)) < 1e-9 && cabs(csub(seg[11], b)) < 1e-9,
    'geodesic segment endpoints',
  );

  // arc stays inside the disk
  ok(seg.every((p) => cabs(p) < 1.0 + 1e-12), 'geodesic arc inside closed disk');

  return { passed, failures };
}
