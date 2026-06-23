/**
 * walk.ts — animation driver for /api/walk trajectories, plus the PURE
 * reduction-bookkeeping logic (computeTokenFates), kept dependency-free so it
 * compiles standalone for the node contract test.
 *
 * The engine emits the full straightening schedule (DATA_CONTRACT.md); this
 * module only *replays* it. No arithmetic is ever re-derived here.
 */

// ---------------------------------------------------------------------------
// Contract types (subset the stage consumes)
// ---------------------------------------------------------------------------

export interface WalkStep {
  letter: string;
  type: 'move' | 'fix';
  from?: string;
  mid?: string;
  to?: string;
  at?: string;
}

export interface ReductionEvent {
  index: number; // remove trail[i], trail[i+1] from the *evolving* trail copy
  removed: string[];
}

export interface NormalForm {
  blocks: { monomial: string; h: boolean }[];
  tailMonomial: string;
}

export interface MatrixEntryDisplay {
  re: number;
  im: number;
  str: string;
  vpi: number | null;
}

export interface MatrixDisplay {
  entries: MatrixEntryDisplay[][];
  sde: number;
}

export interface Trajectory {
  ok: boolean;
  word: string | null;
  sde: number;
  steps: WalkStep[];
  trail: string[];
  edgeOwner: number[];
  geodesic: string[];
  reduction: ReductionEvent[];
  normalForm: NormalForm;
  matrix: MatrixDisplay;
  error: string | null;
}

// ---------------------------------------------------------------------------
// Pure token-lifetime bookkeeping (unit-tested against the engine's JSON)
// ---------------------------------------------------------------------------

export type TokenFate = 'alive' | 'half' | 'dead';

export interface FateEvent {
  index: number; // index into the working trail at the time of the event
  lobe: [string, string, string]; // work[i-1], work[i], work[i+1] (a backtrack lobe)
  edgeOwners: [number, number]; // H ordinals owning the two removed edges
  died: number[]; // H ordinals whose last edge was removed by this event
  halved: number[]; // H ordinals reduced from 2 edges to 1 by this event
  survivingEdges: number; // trail edges remaining after this event
}

export interface TokenFates {
  events: FateEvent[];
  fates: TokenFate[]; // final per-H fate, indexed by H ordinal
  survivingFull: number; // count of H tokens with both edges intact (= sde)
  survivingEdges: number; // = 2 · sde
  hCount: number;
  finalTrail: string[]; // must equal the engine's geodesic
}

/**
 * Replay the engine's reduction schedule on working copies of trail/edgeOwner
 * and track per-H token fates. Event {index: i} removes trail vertices i, i+1,
 * i.e. working edges i−1 and i (the lobe around the backtrack).
 */
export function computeTokenFates(
  trail: string[],
  edgeOwner: number[],
  reduction: ReductionEvent[],
): TokenFates {
  if (edgeOwner.length !== Math.max(0, trail.length - 1)) {
    throw new Error(
      `contract violation: ${trail.length} trail vertices but ${edgeOwner.length} edge owners`,
    );
  }
  const workTrail = trail.slice();
  const workOwner = edgeOwner.slice();
  const hCount = edgeOwner.length === 0 ? 0 : Math.max(...edgeOwner) + 1;
  const remaining = new Array<number>(hCount).fill(0);
  for (const h of edgeOwner) remaining[h]++;

  const events: FateEvent[] = [];
  let survivingEdges = workOwner.length;

  for (const ev of reduction) {
    const i = ev.index;
    if (i < 1 || i + 1 >= workTrail.length) {
      throw new Error(`contract violation: reduction index ${i} out of range`);
    }
    if (workTrail[i + 1] !== workTrail[i - 1]) {
      throw new Error(
        `contract violation: event at ${i} is not a backtrack ` +
          `(${workTrail[i - 1]} vs ${workTrail[i + 1]})`,
      );
    }
    if (ev.removed && (ev.removed[0] !== workTrail[i] || ev.removed[1] !== workTrail[i + 1])) {
      throw new Error(`contract violation: removed vertices mismatch at index ${i}`);
    }
    const lobe: [string, string, string] = [workTrail[i - 1], workTrail[i], workTrail[i + 1]];
    const owners: [number, number] = [workOwner[i - 1], workOwner[i]];
    const died: number[] = [];
    const halved: number[] = [];
    if (owners[0] === owners[1]) {
      remaining[owners[0]] -= 2;
      if (remaining[owners[0]] === 0) died.push(owners[0]);
    } else {
      for (const h of owners) {
        remaining[h] -= 1;
        if (remaining[h] === 1) halved.push(h);
        else if (remaining[h] === 0) died.push(h);
      }
    }
    workTrail.splice(i, 2);
    workOwner.splice(i - 1, 2);
    survivingEdges -= 2;
    events.push({ index: i, lobe, edgeOwners: owners, died, halved, survivingEdges });
  }

  const fates: TokenFate[] = remaining.map((r) =>
    r >= 2 ? 'alive' : r === 1 ? 'half' : 'dead',
  );
  return {
    events,
    fates,
    survivingFull: fates.filter((f) => f === 'alive').length,
    survivingEdges,
    hCount,
    finalTrail: workTrail,
  };
}

/** Consecutive [a,b] pairs of a vertex path. */
export function pathEdges(path: string[]): Array<[string, string]> {
  const out: Array<[string, string]> = [];
  for (let i = 0; i + 1 < path.length; i++) out.push([path[i], path[i + 1]]);
  return out;
}

/**
 * Reduce a tree walk to its geodesic by removing adjacent backtracks — the same
 * algorithm the engine uses (vizgen/walks.py `_reduction_schedule`). Returns the
 * reduced path and the event schedule. Used to straighten a *segment* of the
 * walk (the path between two pure vertices) independently of the whole trail.
 */
export function reducePath(trail: string[]): { geodesic: string[]; events: ReductionEvent[] } {
  const work = trail.slice();
  const events: ReductionEvent[] = [];
  let i = 1;
  while (i < work.length) {
    if (i + 1 < work.length && work[i + 1] === work[i - 1]) {
      events.push({ index: i, removed: [work[i], work[i + 1]] });
      work.splice(i, 2);
      i = Math.max(1, i - 1);
    } else {
      i++;
    }
  }
  return { geodesic: work, events };
}

// ---------------------------------------------------------------------------
// Animation driver (renderer/scrubber accessed through structural ports so
// this file stays pixi-free and standalone-compilable)
// ---------------------------------------------------------------------------

export interface StagePort {
  drawDeepPath(addrs: string[]): void;
  setTrailEdges(edges: Array<[string, string]>): void;
  clearTrail(): void;
  appendTrailEdge(a: string, b: string): void;
  compassSpin(addr: string, ms: number): Promise<void>;
  movePulse(from: string, mid: string, to: string, ms: number): Promise<void>;
  maybeFollow(addr: string): void;
  flashLobe(lobe: [string, string, string], ms: number): Promise<void>;
  triumph(geodesic: string[], ms: number): Promise<void>;
  tracePath(addrs: string[], msPerEdge: number): Promise<void>;
}

export interface ScrubPort {
  setActive(step: number): void;
  onTokenDead(hOrdinal: number): void;
  onTokenHalf(hOrdinal: number): void;
  setMeter(currentH: number, optimal: number | null): void;
  resetTokens(): void;
}

const sleep = (ms: number) => new Promise<void>((r) => setTimeout(r, ms));

export class WalkDriver {
  msPerStep = 550;

  private traj: Trajectory | null = null;
  private gen = 0; // generation counter: bumping cancels in-flight playback
  private paused = false;
  private playing = false;
  private resumeAt = 0;

  onStraightenDone: ((sde: number) => void) | null = null;
  onPlayDone: (() => void) | null = null;

  constructor(
    private stage: StagePort,
    private scrub: ScrubPort,
  ) {}

  get trajectory(): Trajectory | null {
    return this.traj;
  }

  get isPlaying(): boolean {
    return this.playing;
  }

  get isPaused(): boolean {
    return this.paused;
  }

  /**
   * Load a trajectory; cancels any running animation. Clears the trail unless
   * `keepTrail` is set (used by append, where the existing prefix stays on
   * screen and the new trajectory extends it seamlessly).
   */
  load(traj: Trajectory, keepTrail = false): void {
    this.gen++;
    this.playing = false;
    this.paused = false;
    this.resumeAt = 0;
    this.traj = traj;
    this.stage.drawDeepPath(traj.trail);
    if (!keepTrail) this.stage.clearTrail();
    this.scrub.resetTokens();
    this.scrub.setMeter(this.hCount(), null);
  }

  hCount(): number {
    return this.traj ? this.traj.steps.filter((s) => s.type === 'move').length : 0;
  }

  /** Trail edges contributed by steps 0..stepIndex inclusive. */
  private edgesUpToStep(stepIndex: number): Array<[string, string]> {
    const out: Array<[string, string]> = [];
    if (!this.traj) return out;
    for (let i = 0; i <= stepIndex && i < this.traj.steps.length; i++) {
      const s = this.traj.steps[i];
      if (s.type === 'move') {
        out.push([s.from!, s.mid!]);
        out.push([s.mid!, s.to!]);
      }
    }
    return out;
  }

  private async pauseGate(myGen: number): Promise<void> {
    while (this.paused && myGen === this.gen) await sleep(60);
  }

  /**
   * Animate steps from `fromStep`. With `keepTrail`, the currently displayed
   * path is left intact and only new edges are appended (used by append, so
   * continuing after a Synthesize extends the *geodesic* on screen rather than
   * snapping back to the meander prefix).
   */
  async play(fromStep = 0, keepTrail = false): Promise<void> {
    const t = this.traj;
    if (!t) return;
    const myGen = ++this.gen;
    this.playing = true;
    this.paused = false;
    if (!keepTrail) this.stage.setTrailEdges(this.edgesUpToStep(fromStep - 1));
    for (let i = fromStep; i < t.steps.length; i++) {
      await this.pauseGate(myGen);
      if (myGen !== this.gen) return;
      this.scrub.setActive(i);
      this.resumeAt = i + 1;
      const s = t.steps[i];
      if (s.type === 'fix') {
        await this.stage.compassSpin(s.at!, this.msPerStep);
      } else {
        await this.stage.movePulse(s.from!, s.mid!, s.to!, this.msPerStep);
        if (myGen !== this.gen) return;
        this.stage.appendTrailEdge(s.from!, s.mid!);
        this.stage.appendTrailEdge(s.mid!, s.to!);
        this.stage.maybeFollow(s.to!);
      }
      if (myGen !== this.gen) return;
    }
    this.playing = false;
    this.resumeAt = t.steps.length;
    this.scrub.setActive(-1);
    if (this.onPlayDone) this.onPlayDone();
  }

  /** Space bar: pause if playing, resume (or restart) otherwise. */
  togglePlay(): void {
    if (!this.traj) return;
    if (this.playing) {
      this.paused = !this.paused;
    } else {
      const from = this.resumeAt >= this.traj.steps.length ? 0 : this.resumeAt;
      void this.play(from);
    }
  }

  /** Scrubber click: jump instantly to the state *after* step i. */
  jumpTo(stepIndex: number): void {
    const t = this.traj;
    if (!t) return;
    this.gen++; // cancel any animation
    this.playing = false;
    this.paused = false;
    this.resumeAt = stepIndex + 1;
    this.stage.setTrailEdges(this.edgesUpToStep(stepIndex));
    this.scrub.setActive(stepIndex);
    const s = t.steps[stepIndex];
    if (s && s.type === 'move') this.stage.maybeFollow(s.to!);
  }

  /**
   * The signature moment: straighten the PENDING segment `t.trail[fromIndex:]`
   * into its geodesic and return that geodesic (addresses). Only the pending
   * meander is retracted/redrawn — the frozen prefix (already-synthesized
   * geodesics, kept as renderer frozen edges by the caller) is never touched,
   * so synthesizing after an append does not re-expand earlier work.
   * fromIndex = 0 is the whole walk (first synthesize / global case).
   */
  async straighten(fromIndex = 0): Promise<string[]> {
    const t = this.traj;
    if (!t) return [];
    const myGen = ++this.gen;
    this.playing = false;
    this.paused = false;

    const seg = t.trail.slice(fromIndex);
    if (seg.length < 2) return seg; // nothing pending (segment is sde 0)

    const anchorH = fromIndex / 2; // H's frozen before this segment (trail: 2 verts per H)
    const segOwners = t.edgeOwner.slice(fromIndex).map((o) => o - anchorH); // local ordinals
    const { geodesic: segGeo, events } = reducePath(seg);
    const fates = computeTokenFates(seg, segOwners, events);
    const segSde = (segGeo.length - 1) / 2;

    this.stage.drawDeepPath(seg);
    this.stage.setTrailEdges(pathEdges(seg)); // active = the pending meander
    this.scrub.setActive(-1);
    this.scrub.setMeter(fates.hCount, segSde);
    await sleep(300);
    if (myGen !== this.gen) return segGeo;

    const work = seg.slice();
    const lobeMs = Math.min(Math.max(this.msPerStep * 0.8, 180), 480);
    for (const ev of fates.events) {
      await this.pauseGate(myGen);
      if (myGen !== this.gen) return segGeo;
      await this.stage.flashLobe(ev.lobe, lobeMs);
      if (myGen !== this.gen) return segGeo;
      work.splice(ev.index, 2);
      this.stage.setTrailEdges(pathEdges(work));
      for (const h of ev.halved) this.scrub.onTokenHalf(h);
      for (const h of ev.died) this.scrub.onTokenDead(h);
      this.scrub.setMeter(ev.survivingEdges / 2, segSde);
    }
    // pulse following the synthesized segment, anchor → green endpoint (camera follows)
    await this.stage.tracePath(segGeo, 280);
    return segGeo;
  }

  cancel(): void {
    this.gen++;
    this.playing = false;
    this.paused = false;
  }
}
