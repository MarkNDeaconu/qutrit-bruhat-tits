/**
 * inspector.ts — right-side vertex inspector (plain DOM).
 * Shows the vertex card (HNF basis + Gram matrices with v_π badges, residue
 * matrix mod χ as 0/1/2 chips, kind badge, and the representative unitary for
 * pure vertices). Data comes from the transport (HTTP engine or Pyodide
 * worker — same contract). Lazy, cached per address.
 */

import { apiVertex } from './transport.js';

interface EntryDisplay {
  re: number;
  im: number;
  str: string;
  vpi: number | null;
}

interface MatDisplay {
  entries: EntryDisplay[][];
  sde?: number;
}

interface VertexCard {
  addr: string;
  kind: 'P' | 'A';
  depth: number;
  basis: MatDisplay;
  gram: MatDisplay;
  residue: number[][];
  unitary: MatDisplay | null;
  [k: string]: unknown;
}

const esc = (s: string): string =>
  s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

function vpiBadge(v: number | null): string {
  if (v === null) return `<span class="vpi vpi-none">·</span>`;
  const cls = v < 0 ? 'vpi-neg' : v > 0 ? 'vpi-pos' : 'vpi-zero';
  const sign = v > 0 ? '+' : '';
  return `<span class="vpi ${cls}">v<sub>π</sub>=${sign}${v}</span>`;
}

function matTable(m: MatDisplay): string {
  const rows = m.entries
    .map(
      (row) =>
        `<tr>${row
          .map(
            (c) => `<td><span class="cell-str">${esc(c.str)}</span>${vpiBadge(c.vpi)}</td>`,
          )
          .join('')}</tr>`,
    )
    .join('');
  return `<table class="mat"><tbody>${rows}</tbody></table>`;
}

function residueChips(r: number[][]): string {
  const rows = r
    .map(
      (row) =>
        `<tr>${row.map((v) => `<td><span class="chip chip-${v}">${v}</span></td>`).join('')}</tr>`,
    )
    .join('');
  return `<table class="residue"><tbody>${rows}</tbody></table>`;
}

export class Inspector {
  private el: HTMLElement;
  private body: HTMLElement;
  private cache = new Map<string, VertexCard>();
  private current: string | null = null;

  onClose: (() => void) | null = null;

  constructor(root: HTMLElement) {
    this.el = root;
    root.innerHTML = `
      <div class="insp-head">
        <span class="insp-title" id="insp-title"></span>
        <button class="insp-close" id="insp-close" title="close (Esc)">×</button>
      </div>
      <div class="insp-body" id="insp-body"></div>`;
    this.body = root.querySelector('#insp-body')!;
    root.querySelector('#insp-close')!.addEventListener('click', () => this.hide());
  }

  get isOpen(): boolean {
    return this.el.classList.contains('open');
  }

  hide(): void {
    this.el.classList.remove('open');
    this.current = null;
    if (this.onClose) this.onClose();
  }

  async show(addr: string): Promise<void> {
    this.current = addr;
    this.el.classList.add('open');
    const title = this.el.querySelector('#insp-title')!;
    title.textContent = addr === '' ? 'origin e₀' : `vertex ${addr}`;
    const cached = this.cache.get(addr);
    if (cached) {
      this.render(cached);
      return;
    }
    this.body.innerHTML = `<p class="insp-loading">computing…</p>`;
    try {
      const { status, data } = await apiVertex(addr);
      // worker returns {error}, the HTTP engine returns {detail}; accept either
      if (status !== 200) throw new Error(data?.error ?? data?.detail ?? `status ${status}`);
      const card = data as VertexCard;
      this.cache.set(addr, card);
      if (this.current === addr) this.render(card);
    } catch (err) {
      if (this.current === addr) {
        this.body.innerHTML = `<p class="insp-error">unavailable — ${esc(String(err))}</p>`;
      }
    }
  }

  private render(card: VertexCard): void {
    const kindName = card.kind === 'P' ? 'pure' : 'alternating';
    const parts: string[] = [];
    parts.push(
      `<div class="insp-badges">
         <span class="badge kind-${card.kind}">${kindName}</span>
         <span class="badge depth">depth ${card.depth}</span>
       </div>`,
    );
    parts.push(`<h4>canonical HNF basis</h4>${matTable(card.basis)}`);
    parts.push(`<h4>Gram matrix</h4>${matTable(card.gram)}`);
    parts.push(
      `<h4>residue mod χ <span class="muted">(${
        card.kind === 'P' ? 'Gram' : 'χ·Gram'
      } over 𝔽₃)</span></h4>${residueChips(card.residue)}`,
    );
    if (card.unitary) {
      parts.push(
        `<h4>representative unitary <span class="muted">U·e₀ = v</span></h4>${matTable(card.unitary)}`,
      );
    }
    parts.push(
      `<p class="insp-note">exact data from the frozen verified core — every cell is exact arithmetic <span class="badge verified">machine-verified</span></p>`,
    );
    this.body.innerHTML = parts.join('');
  }
}
