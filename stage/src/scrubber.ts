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
  private hToToken: number[] = []; // (pending) H ordinal -> token index
  private activeIdx = -1;
  private locked = 0; // H's already frozen in synthesized segments

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

  /** One token per letter of the PENDING segment (traj.steps from `fromStep`). */
  build(traj: Trajectory, fromStep = 0): void {
    this.root.classList.remove('empty');
    this.tokensEl.innerHTML = '';
    this.tokens = [];
    this.hToToken = [];
    this.activeIdx = -1;
    traj.steps.slice(fromStep).forEach((s, i) => {
      const el = document.createElement('button');
      const letter = s.letter.toUpperCase();
      el.className = `tok tok-${letter.toLowerCase()}`;
      el.title = letter;
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

  /** H's locked into already-synthesized segments (shown alongside pending). */
  setLocked(n: number): void {
    this.locked = n;
    this.setMeter(this.hToToken.length, null);
  }

  /** Drop the pending tokens once the segment has been frozen. */
  clearPending(): void {
    this.tokensEl.innerHTML = '';
    this.tokens = [];
    this.hToToken = [];
    this.activeIdx = -1;
    this.root.classList.toggle('empty', this.locked === 0);
    this.setMeter(0, null);
  }

  clear(): void {
    this.root.classList.add('empty');
    this.tokensEl.innerHTML = '';
    this.tokens = [];
    this.hToToken = [];
    this.activeIdx = -1;
    this.locked = 0;
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

  /** Meter: 'locked M H · pending w → optimal n', counting down live. */
  setMeter(currentH: number, optimal: number | null): void {
    const parts: string[] = [];
    if (this.locked > 0) parts.push(`locked <b>${this.locked}</b> H`);
    if (currentH > 0 || optimal !== null) {
      parts.push(
        optimal === null
          ? `pending <b>${currentH}</b> H`
          : `pending <b>${currentH}</b> &nbsp;→&nbsp; optimal <b>${optimal}</b>`,
      );
    }
    this.meterEl.innerHTML = parts.join(' &nbsp;·&nbsp; ') || '&nbsp;';
  }
}
