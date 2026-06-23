// Plain-node math tests for src/hyperbolic.ts (compiled to test/build by
// `npm run build:lib`). No framework; throws on first failure, exit code 1.

import {
  addressToLayout,
  addressInterval,
  addressToDisk,
  diskPoint,
  mobiusApply,
  mobiusCompose,
  mobiusInvert,
  poincareSegmentPoints,
  runMathChecks,
  cabs,
  csub,
  LAMBDA,
} from './build/hyperbolic.js';
import { generateBall, childCount, kindOf, parentOf } from './build/tree.js';

let n = 0;
function assert(cond, label) {
  n++;
  if (!cond) {
    console.error(`FAIL [${label}]`);
    process.exit(1);
  }
}
const near = (a, b, eps) => Math.abs(a - b) <= eps;

// --- the module's own assertion suite ---------------------------------------
const r = runMathChecks();
assert(r.failures.length === 0, `runMathChecks: ${r.failures.join('; ')}`);
console.log(`runMathChecks: ${r.passed} internal checks passed`);

// --- determinism -------------------------------------------------------------
{
  const a = addressToLayout('30121');
  const b = addressToLayout('30121');
  assert(a.theta === b.theta && a.rho === b.rho, 'addressToLayout deterministic');
  assert(near(a.rho, 5 * LAMBDA, 1e-15), 'rho = depth * ln 3');
}

// --- root children: 4 distinct equally spaced quarter midpoints ---------------
{
  const thetas = ['0', '1', '2', '3'].map((x) => addressToLayout(x).theta);
  const distinct = new Set(thetas.map((t) => t.toFixed(12)));
  assert(distinct.size === 4, 'root children thetas distinct');
  const sorted = [...thetas].sort((x, y) => x - y);
  for (let i = 0; i < 3; i++) {
    assert(near(sorted[i + 1] - sorted[i], Math.PI / 2, 1e-12), 'root children equally spaced');
  }
  // the canonical orientation: {pi/4, 3pi/4, 5pi/4, 7pi/4}
  for (let i = 0; i < 4; i++) {
    assert(near(sorted[i], Math.PI / 4 + (i * Math.PI) / 2, 1e-12), 'quarter midpoints');
  }
}

// --- sibling angle ordering ----------------------------------------------------
{
  const sibs = ['20', '21', '22'].map((x) => addressToLayout(x).theta);
  assert(sibs[0] < sibs[1] && sibs[1] < sibs[2], 'sibling angle ordering (digits ascending)');
  const root = ['0', '1', '2', '3'].map((x) => addressToLayout(x).theta);
  assert(
    root[0] < root[1] && root[1] < root[2] && root[2] < root[3],
    'root sibling angle ordering',
  );
}

// --- nested interval containment for 50 random addresses ----------------------
{
  let seed = 12345;
  const rand = () => {
    // LCG — deterministic across runs
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
  };
  for (let trial = 0; trial < 50; trial++) {
    const len = 1 + Math.floor(rand() * 11);
    let addr = String(Math.floor(rand() * 4));
    for (let i = 1; i < len; i++) addr += String(Math.floor(rand() * 3));
    const child = addr;
    const parent = child.slice(0, -1);
    const iv = addressInterval(parent);
    const th = addressToLayout(child).theta;
    assert(th >= iv.lo && th < iv.hi, `containment: theta(${child}) inside interval(${parent})`);
    // and the child's own interval nests inside the parent's
    const civ = addressInterval(child);
    assert(
      civ.lo >= iv.lo - 1e-15 && civ.hi <= iv.hi + 1e-15,
      `containment: interval(${child}) ⊂ interval(${parent})`,
    );
  }
}

// --- |diskPoint| < 1 always -----------------------------------------------------
{
  for (const addr of generateBall(6)) {
    assert(cabs(addressToDisk(addr)) < 1, `|diskPoint(${addr})| < 1`);
  }
  // very deep point still inside
  const deep = diskPoint(1.0, 30 * LAMBDA);
  assert(cabs(deep) < 1, '|diskPoint| < 1 at depth 30');
}

// --- Möbius compose/invert round trip to 1e-12 ----------------------------------
{
  let seed = 777;
  const rand = () => {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return seed / 2147483648;
  };
  for (let trial = 0; trial < 25; trial++) {
    const m1 = { a: rand() * 6 - 3, c: { re: rand() * 1.6 - 0.8, im: rand() * 1.6 - 0.8 } };
    const m2 = { a: rand() * 6 - 3, c: { re: rand() * 1.6 - 0.8, im: rand() * 1.6 - 0.8 } };
    if (cabs(m1.c) >= 0.95 || cabs(m2.c) >= 0.95) continue;
    // compose agrees with sequential application
    const z = { re: rand() * 1.2 - 0.6, im: rand() * 1.2 - 0.6 };
    const seq = mobiusApply(m1, mobiusApply(m2, z));
    const comp = mobiusApply(mobiusCompose(m1, m2), z);
    assert(cabs(csub(seq, comp)) < 1e-12, 'compose agrees with sequential apply');
    // m ∘ m^{-1} = id as a transform
    const id = mobiusCompose(m1, mobiusInvert(m1));
    const w = mobiusApply(id, z);
    assert(cabs(csub(w, z)) < 1e-12, 'compose(m, invert(m)) acts as identity');
    // and invert really inverts pointwise
    const back = mobiusApply(mobiusInvert(m1), mobiusApply(m1, z));
    assert(cabs(csub(back, z)) < 1e-12, 'invert round trip pointwise');
  }
}

// --- geodesic segment endpoints to 1e-9 -------------------------------------------
{
  const cases = [
    [{ re: 0.5, im: 0.1 }, { re: -0.2, im: 0.6 }], // generic arc
    [{ re: 0.3, im: 0.0 }, { re: -0.7, im: 0.0 }], // diametral (straight)
    [{ re: 0.0, im: 0.0 }, { re: 0.4, im: 0.4 }], // through the origin (straight)
    [{ re: 0.81, im: 0.02 }, { re: 0.8, im: 0.1 }], // short, near the rim
    [{ re: 0.25, im: 0.25 }, { re: 0.25, im: 0.25 }], // degenerate (same point)
  ];
  for (const [z1, z2] of cases) {
    const pts = poincareSegmentPoints(z1, z2, 12);
    assert(pts.length === 12, 'segment point count');
    assert(cabs(csub(pts[0], z1)) < 1e-9, 'segment start endpoint');
    assert(cabs(csub(pts[11], z2)) < 1e-9, 'segment end endpoint');
    for (const p of pts) assert(cabs(p) < 1 + 1e-12, 'segment stays in closed disk');
  }
}

// --- tree skeleton sanity ----------------------------------------------------------
{
  let count = 0;
  let depth8 = 0;
  for (const addr of generateBall(8)) {
    count++;
    if (addr.length === 8) depth8++;
  }
  assert(count === 13121, `|ball(8)| = 13121 (got ${count})`);
  assert(depth8 === 4 * 3 ** 7, 'sphere(8) = 4·3^7');
  assert(childCount('') === 4 && childCount('21') === 3, 'child counts');
  assert(kindOf('') === 'P' && kindOf('2') === 'A' && kindOf('21') === 'P', 'kind = depth parity');
  assert(parentOf('210') === '21' && parentOf('') === null, 'parentOf');
}

console.log(`math.test.mjs: all ${n} assertions passed`);
