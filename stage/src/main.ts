/**
 * main.ts — UI shell: top bar, examples, keyboard, captions, toasts, modal,
 * engine health. Wires renderer (pixi stage) + walk driver + scrubber +
 * inspector together. All math comes from the engine over /api (contract).
 */

import './style.css';
import { runMathChecks } from './hyperbolic.js';
import { Renderer } from './renderer.js';
import { Scrubber } from './scrubber.js';
import { Inspector } from './inspector.js';
import { WalkDriver, pathEdges, type Trajectory } from './walk.js';
import {
  apiWalk,
  apiRandomWord,
  apiHealth,
  backendName,
  onBackendStatus,
  warmup,
} from './transport.js';

// ---------------------------------------------------------------------------
// field notes (static for MVP; every on-screen claim carries its badge)
// ---------------------------------------------------------------------------

const NOTES = {
  idle: `Every circuit is a walk on the (4,4)-biregular Bruhat–Tits tree: d(e₀, U·e₀) = 2·sde(U). <span class="badge proved">proved</span>`,
  walk: `Each H crosses pure → alternating → pure. S and R never move the vertex — the compass spin is real: monomials permute the branch directions. The <b style="color:#2ee36e">green node</b> is where this circuit lands. <span class="badge proved">proved</span>`,
  append: `Appended from the <b style="color:#2ee36e">green node</b>: the walk continues along the same path. The endpoint is U·e₀ for the full accumulated word. <span class="badge proved">proved</span>`,
  oracle: `The <b style="color:#2ee36e">green node</b> marks U·e₀. Watch the 12 candidates at each step: exactly one descends. It always does — we checked every step. <span class="badge verified">machine-verified: 690/690</span>`,
  straighten: `Same endpoint, shortest route: the geodesic from e₀ to the <b style="color:#2ee36e">green node</b> is the optimal circuit. Circuit optimization is path straightening — and the tree counts H, the denominator gate (Clifford), not non-Clifford gates. <span class="badge verified">machine-verified</span>`,
  incremental: `Each appended segment is synthesized from where the last one ended — already-clean parts stay frozen. The committed circuit is the concatenation of per-segment-optimal geodesics. <span class="badge proved">proved</span>`,
  cutoff: `The tree is infinite — this circuit runs past the radius this flat view can draw exactly. The <b style="color:#ff5a5a">red node</b> marks where it leaves the rendered range.`,
} as const;

const CUTOFF_MSG = 'This path runs beyond the visualizer’s range — shown out to the red edge marker.';

// ---------------------------------------------------------------------------

function toast(msg: string): void {
  const holder = document.getElementById('toast-holder')!;
  const el = document.createElement('div');
  el.className = 'toast';
  el.textContent = msg;
  holder.appendChild(el);
  requestAnimationFrame(() => el.classList.add('show'));
  setTimeout(() => {
    el.classList.remove('show');
    setTimeout(() => el.remove(), 400);
  }, 4200);
}

function setCaption(html: string): void {
  document.getElementById('caption')!.innerHTML = html;
}

function showBanner(text: string): void {
  const b = document.getElementById('banner')!;
  b.innerHTML = `${text} <span class="badge verified">machine-verified</span>`;
  b.classList.add('show');
}

function hideBanner(): void {
  document.getElementById('banner')!.classList.remove('show');
}

async function postWalk(body: { word?: string; matrix?: string[][] }): Promise<Trajectory | null> {
  try {
    const data = (await apiWalk(body)) as Trajectory;
    if (!data || !data.ok) {
      toast(data?.error ?? 'rejected the input');
      return null;
    }
    return data;
  } catch {
    toast(
      backendName() === 'pyodide'
        ? 'compute failed — please retry'
        : 'engine unreachable — start it on 127.0.0.1:8137',
    );
    return null;
  }
}

// ---------------------------------------------------------------------------

async function boot(): Promise<void> {
  // display-math self-check (exact math lives engine-side; this is layout only)
  const checks = runMathChecks();
  if (checks.failures.length > 0) {
    console.warn('hyperbolic.ts self-checks FAILED:', checks.failures);
  } else {
    console.info(`hyperbolic.ts self-checks: ${checks.passed} passed`);
  }

  // In the deployed (Pyodide) build the verified Python boots in a Web Worker;
  // show a splash while it downloads. The tree renders immediately regardless
  // (its layout is procedural and needs no backend).
  const splash = document.getElementById('splash')!;
  const splashStatus = document.getElementById('splash-status')!;
  if (backendName() === 'pyodide') {
    splash.classList.remove('hidden');
    onBackendStatus((msg, ready, error) => {
      if (error) {
        splashStatus.textContent = `failed to load: ${error}`;
        splashStatus.classList.add('error');
        return;
      }
      splashStatus.textContent = ready ? 'ready' : msg;
      if (ready) {
        splash.classList.add('done');
        setTimeout(() => splash.classList.add('hidden'), 700);
      }
    });
    warmup();
  } else {
    splash.classList.add('hidden');
  }

  const renderer = await Renderer.create(document.getElementById('stage')!);
  const scrubber = new Scrubber(document.getElementById('scrubber')!);
  const inspector = new Inspector(document.getElementById('inspector')!);
  const driver = new WalkDriver(renderer, scrubber);

  const wordInput = document.getElementById('word') as HTMLInputElement;
  const speedInput = document.getElementById('speed') as HTMLInputElement;
  const examplesSel = document.getElementById('examples') as HTMLSelectElement;
  const statusEl = document.getElementById('engine-status')!;

  let light = false;
  // Incremental synthesis state. `anchorWord` is the already-synthesized prefix
  // (frozen on screen as clean geodesics); `pendingWord` is the appended,
  // not-yet-synthesized suffix (the active meander). committed = anchor+pending.
  // Play replaces everything from the origin; Append extends pending; Synthesize
  // freezes the pending meander into a geodesic and advances the anchor — so an
  // already-synthesized part never reverts or re-expands. matrixMode disables
  // append (a pasted matrix has no word suffix).
  let anchorWord = '';
  let pendingWord = '';
  let matrixMode = false;
  let matrixDone = false;
  let lockedH = 0; // H's in the frozen geodesics (= sum of segment sde)
  let segCount = 0; // number of synthesized segments
  const committedWord = (): string => anchorWord + pendingWord;
  const countH = (w: string): number => (w.match(/H/g) || []).length;

  // -- wiring ----------------------------------------------------------------

  renderer.onVertexClick = (addr) => {
    renderer.selectVertex(addr);
    void inspector.show(addr);
  };
  renderer.onVertexDblClick = (addr) => void renderer.focusOn(addr, 700);
  inspector.onClose = () => renderer.selectVertex(null);

  // when a walk finishes, mark the vertex it landed on with the green endpoint
  driver.onPlayDone = () => {
    setCaption(NOTES.oracle);
    const t = driver.trajectory;
    if (t && t.trail.length > 0) renderer.setEndpoint(t.trail[t.trail.length - 1]);
  };
  // a walk that runs past the renderable range stops at a red "out of range" node
  driver.onCutoff = (addr) => {
    renderer.setEndpoint(addr, true);
    setCaption(NOTES.cutoff);
    toast(CUTOFF_MSG);
  };

  driver.msPerStep = Number(speedInput.value);
  speedInput.addEventListener('input', () => {
    driver.msPerStep = Number(speedInput.value);
  });

  // live validation: letters HSR only, uppercased as you type
  wordInput.addEventListener('input', () => {
    const cleaned = wordInput.value.toUpperCase().replace(/[^HSR]/g, '');
    if (cleaned !== wordInput.value) {
      wordInput.classList.add('invalid');
      setTimeout(() => wordInput.classList.remove('invalid'), 350);
    }
    wordInput.value = cleaned;
  });

  // Play: wipe everything and start fresh from the origin with `word`.
  async function playWord(word: string): Promise<void> {
    if (!/^[HSR]+$/.test(word)) {
      toast('word must be non-empty letters H, S, R');
      return;
    }
    hideBanner();
    renderer.setEndpoint(null);
    renderer.clearFrozen();
    const traj = await postWalk({ word });
    if (!traj) return;
    anchorWord = '';
    pendingWord = word;
    matrixMode = false;
    matrixDone = false;
    lockedH = 0;
    segCount = 0;
    driver.load(traj); // clears the active trail
    scrubber.setLocked(0);
    scrubber.build(traj); // pending = the whole word
    setCaption(NOTES.walk);
    void renderer.resetCamera(450); // back to the origin, then walk out
    wordInput.value = '';
    void driver.play();
  }

  // Append: continue from the green endpoint (the current anchor/endpoint),
  // adding `seg` to the pending segment. Frozen (synthesized) parts are untouched.
  async function appendWord(seg: string): Promise<void> {
    if (matrixMode) {
      toast('append needs a word — press Play with H, S, R first');
      return;
    }
    if (!committedWord()) {
      await playWord(seg); // nothing to append to yet → behave like Play
      return;
    }
    if (!/^[HSR]+$/.test(seg)) {
      toast('append needs letters H, S, R');
      return;
    }
    hideBanner();
    renderer.setEndpoint(null);
    const oldLen = committedWord().length; // one walk step per letter
    const newWord = committedWord() + seg;
    const traj = await postWalk({ word: newWord });
    if (!traj) return; // state unchanged on failure
    pendingWord = pendingWord + seg;
    driver.load(traj, true); // keep frozen geodesics + the active meander on screen
    scrubber.build(traj, anchorWord.length); // pending-scoped tokens
    setCaption(NOTES.append);
    wordInput.value = '';
    // keepTrail: extend the active meander from the current endpoint; never
    // touch the frozen (already-synthesized) prefix
    void driver.play(oldLen, true);
  }

  async function playMatrix(rows: string[][]): Promise<void> {
    hideBanner();
    renderer.setEndpoint(null);
    renderer.clearFrozen();
    const traj = await postWalk({ matrix: rows });
    if (!traj) return;
    wordInput.value = '';
    anchorWord = '';
    pendingWord = '';
    matrixMode = true;
    matrixDone = false;
    lockedH = 0;
    segCount = 0;
    driver.load(traj);
    scrubber.setLocked(0);
    scrubber.build(traj);
    setCaption(NOTES.walk);
    void renderer.resetCamera(450);
    void driver.play();
  }

  // Synthesize: straighten ONLY the pending segment into its geodesic, freeze it,
  // and advance the anchor. Already-synthesized parts stay frozen (no re-expand,
  // no revert), so alternating append/synthesize builds the circuit incrementally.
  async function synthesize(): Promise<void> {
    if (!driver.trajectory) {
      const word = committedWord() || wordInput.value.trim();
      if (!word) {
        toast('type a word (or pick an example) first');
        return;
      }
      const traj = await postWalk({ word });
      if (!traj) return;
      anchorWord = '';
      pendingWord = word;
      matrixMode = false;
      matrixDone = false;
      lockedH = 0;
      segCount = 0;
      driver.load(traj);
      scrubber.setLocked(0);
      scrubber.build(traj);
    }
    if (matrixMode && matrixDone) return; // a matrix has no further segments
    if (!matrixMode && !pendingWord) return; // nothing new to synthesize

    const t = driver.trajectory!;
    const fromIndex = matrixMode ? 0 : 2 * countH(anchorWord);
    const anchorAddr = t.trail[fromIndex] ?? '';
    const G = t.trail[t.trail.length - 1];
    renderer.setEndpoint(G);
    hideBanner();
    setCaption(segCount > 0 ? NOTES.incremental : NOTES.straighten);

    await renderer.recenterOn(anchorAddr, 700); // centre where this segment starts
    const { geodesic: segGeo, cut } = await driver.straighten(fromIndex);

    if (cut) {
      // the path runs past the renderable range — show it out to a red marker and
      // leave the incremental state unchanged (this segment isn't fully synthesized)
      renderer.setTrailEdges(segGeo.length >= 2 ? pathEdges(segGeo) : []);
      renderer.setEndpoint(segGeo[segGeo.length - 1] ?? anchorAddr, true);
      setCaption(NOTES.cutoff);
      toast(CUTOFF_MSG);
      return;
    }

    if (segGeo.length >= 2) renderer.addFrozenEdges(pathEdges(segGeo));
    renderer.setTrailEdges([]); // the segment is now frozen; clear the active meander

    const segSde = Math.max(0, (segGeo.length - 1) / 2);
    lockedH += segSde;
    segCount += 1;
    scrubber.setLocked(lockedH);
    scrubber.clearPending();

    if (matrixMode) {
      matrixDone = true;
    } else {
      anchorWord = committedWord();
      pendingWord = '';
    }

    if (segCount === 1) {
      showBanner(`H-count = distance/2 = sde = ${lockedH}`);
    } else {
      showBanner(`circuit H-count = ${lockedH} (built in ${segCount} optimal pieces)`);
    }
    setCaption(segCount > 1 ? NOTES.incremental : NOTES.straighten);
  }

  async function random20(): Promise<void> {
    try {
      const seed = Math.floor(Math.random() * 1e6);
      const data = (await apiRandomWord(20, seed)) as { word?: string };
      if (!data.word) throw new Error('bad response');
      wordInput.value = data.word;
      await playWord(data.word);
    } catch {
      toast(
        backendName() === 'pyodide'
          ? 'still loading — try again in a moment'
          : 'engine unreachable — start it on 127.0.0.1:8137',
      );
    }
  }

  // A random Clifford word uses H and S only. ⟨H,S⟩ is the finite Clifford group
  // (order 1296), so the orbit is bounded — these never run out of range. The
  // word is plain input (no arithmetic), so it's generated client-side.
  async function randomClifford(): Promise<void> {
    let w = '';
    for (let i = 0; i < 20; i++) w += Math.random() < 0.5 ? 'H' : 'S';
    wordInput.value = w;
    await playWord(w);
  }

  function clearAll(): void {
    driver.cancel();
    renderer.clearTrail();
    renderer.clearFrozen();
    renderer.setEndpoint(null);
    renderer.selectVertex(null);
    scrubber.clear();
    hideBanner();
    void renderer.resetCamera(600); // glide the view back to the origin
    anchorWord = '';
    pendingWord = '';
    matrixMode = false;
    matrixDone = false;
    lockedH = 0;
    segCount = 0;
    setCaption(NOTES.idle);
  }

  document.getElementById('btn-play')!.addEventListener('click', () => {
    void playWord(wordInput.value.trim());
  });
  document.getElementById('btn-append')!.addEventListener('click', () => {
    void appendWord(wordInput.value.trim());
  });
  document.getElementById('btn-synth')!.addEventListener('click', () => void synthesize());
  document.getElementById('btn-random')!.addEventListener('click', () => void random20());
  document.getElementById('btn-clifford')!.addEventListener('click', () => void randomClifford());
  document.getElementById('btn-clear')!.addEventListener('click', clearAll);

  examplesSel.addEventListener('change', () => {
    const v = examplesSel.value;
    examplesSel.selectedIndex = 0;
    if (!v) return;
    if (v === '__random20__') {
      void random20();
    } else if (v === '__clifford__') {
      void randomClifford();
    } else {
      wordInput.value = v;
      void playWord(v);
    }
  });

  // -- matrix modal ------------------------------------------------------------

  const modal = document.getElementById('modal')!;
  const matrixText = document.getElementById('matrix-text') as HTMLTextAreaElement;
  const openModal = () => {
    modal.classList.remove('hidden');
    matrixText.focus();
  };
  const closeModal = () => modal.classList.add('hidden');
  document.getElementById('btn-matrix')!.addEventListener('click', openModal);
  document.getElementById('modal-cancel')!.addEventListener('click', closeModal);
  modal.addEventListener('click', (e) => {
    if (e.target === modal) closeModal();
  });
  document.getElementById('modal-go')!.addEventListener('click', () => {
    const lines = matrixText.value
      .split('\n')
      .map((l) => l.trim())
      .filter((l) => l.length > 0);
    const rows = lines.map((l) => l.split(/[,\s]+/).filter((x) => x.length > 0));
    if (rows.length !== 3 || rows.some((r) => r.length !== 3)) {
      toast('need exactly 3 lines of 3 entries each');
      return;
    }
    closeModal();
    void playMatrix(rows);
  });

  // -- theme / keyboard ----------------------------------------------------------

  function toggleTheme(): void {
    light = !light;
    document.body.classList.toggle('light', light);
    renderer.setLightTheme(light);
  }

  window.addEventListener('keydown', (e) => {
    const target = e.target as HTMLElement;
    const typing =
      target.tagName === 'INPUT' || target.tagName === 'TEXTAREA' || target.tagName === 'SELECT';
    if (e.key === 'Escape') {
      closeModal();
      inspector.hide();
      return;
    }
    if (typing) {
      if (e.key === 'Enter' && target === wordInput) {
        // Enter = Play (fresh); Shift+Enter = Append (continue from green node)
        if (e.shiftKey) void appendWord(wordInput.value.trim());
        else void playWord(wordInput.value.trim());
      }
      return;
    }
    switch (e.key) {
      case ' ':
        e.preventDefault();
        driver.togglePlay();
        break;
      case 'a':
      case 'A':
        void appendWord(wordInput.value.trim());
        break;
      case 'g':
      case 'G':
        void synthesize();
        break;
      case 'l':
      case 'L':
        toggleTheme();
        break;
      case 'r':
      case 'R':
        void random20();
        break;
      case 'c':
      case 'C':
        void randomClifford();
        break;
      case '0':
        void renderer.resetCamera(700);
        break;
    }
  });

  // -- engine health -----------------------------------------------------------

  async function pollHealth(): Promise<void> {
    const pyo = backendName() === 'pyodide';
    try {
      const data = (await apiHealth()) as { ok?: boolean };
      const ok = data.ok === true;
      statusEl.textContent = ok
        ? pyo
          ? 'exact arithmetic: in-browser'
          : 'engine: connected'
        : 'backend: error';
      statusEl.className = ok ? 'connected' : 'offline';
    } catch {
      statusEl.textContent = pyo ? 'loading…' : 'engine: offline';
      statusEl.className = 'offline';
    }
  }
  void pollHealth();
  setInterval(() => void pollHealth(), 5000);

  setCaption(NOTES.idle);
}

void boot();
