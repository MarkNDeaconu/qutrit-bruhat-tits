/**
 * scrubber.ts — bottom timeline (plain DOM): one token per letter of the word.
 * H = large gold square, S = azure dot, R = vermilion dot. Implements the
 * ScrubPort consumed by walk.ts; annihilation visuals are CSS transitions.
 */

import type { ScrubPort, Trajectory } from './walk.js';

export class Scrubber implements ScrubPort {
  private tokensEl: HTMLElement;
  private meterEl: HTMLElement;
  private tokens: HTMLElement[] = [];
  private hToToken: number[] = []; // H ordinal -> token index
  private activeIdx = -1;

  onScrub: ((step: number) => void) | null = null;

  private root: HTMLElement;

  constructor(root: HTMLElement) {
    this.root = root;
    root.classList.add('empty');
    root.innerHTML = `
      <div class="scrub-tokens" id="scrub-tokens"></div>
      <div class="scrub-meter" id="scrub-meter"></div>`;
    this.tokensEl = root.querySelector('#scrub-tokens')!;
    this.meterEl = root.querySelector('#scrub-meter')!;
  }

  /** One token per letter (from traj.steps, which covers matrix input too). */
  build(traj: Trajectory): void {
    this.root.classList.remove('empty');
    this.tokensEl.innerHTML = '';
    this.tokens = [];
    this.hToToken = [];
    this.activeIdx = -1;
    traj.steps.forEach((s, i) => {
      const el = document.createElement('button');
      const letter = s.letter.toUpperCase();
      el.className = `tok tok-${letter.toLowerCase()}`;
      el.title = `${letter} — step ${i + 1}`;
      el.textContent = letter === 'H' ? 'H' : '';
      el.addEventListener('click', () => {
        if (this.onScrub) this.onScrub(i);
      });
      this.tokensEl.appendChild(el);
      this.tokens.push(el);
      if (s.type === 'move') this.hToToken.push(i);
    });
    this.setMeter(this.hToToken.length, null);
  }

  clear(): void {
    this.root.classList.add('empty');
    this.tokensEl.innerHTML = '';
    this.tokens = [];
    this.hToToken = [];
    this.activeIdx = -1;
    this.meterEl.textContent = '';
  }

  setActive(step: number): void {
    if (this.activeIdx >= 0 && this.tokens[this.activeIdx]) {
      this.tokens[this.activeIdx].classList.remove('active');
    }
    this.activeIdx = step;
    const el = step >= 0 ? this.tokens[step] : undefined;
    if (el) {
      el.classList.add('active');
      el.scrollIntoView({ inline: 'nearest', block: 'nearest', behavior: 'smooth' });
    }
  }

  /** Both edges of this H were removed: the token annihilates. */
  onTokenDead(hOrdinal: number): void {
    const el = this.tokens[this.hToToken[hOrdinal]];
    if (el) {
      el.classList.remove('half');
      el.classList.add('dead');
    }
  }

  /** One of its two edges was removed: the token merges (half-alpha). */
  onTokenHalf(hOrdinal: number): void {
    const el = this.tokens[this.hToToken[hOrdinal]];
    if (el && !el.classList.contains('dead')) el.classList.add('half');
  }

  resetTokens(): void {
    for (const el of this.tokens) el.classList.remove('dead', 'half', 'active');
    this.activeIdx = -1;
  }

  /** Optimality meter: 'H letters: w → optimal: n', counting down live. */
  setMeter(currentH: number, optimal: number | null): void {
    if (optimal === null) {
      this.meterEl.innerHTML = `H letters: <b>${currentH}</b>`;
    } else {
      this.meterEl.innerHTML = `H letters: <b>${currentH}</b> &nbsp;→&nbsp; optimal: <b>${optimal}</b>`;
    }
  }
}
