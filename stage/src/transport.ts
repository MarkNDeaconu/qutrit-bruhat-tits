/**
 * transport.ts — backend abstraction. Two modes, chosen at build time:
 *   • 'http'    — talk to the local FastAPI engine over /api (dev, `make demo`)
 *   • 'pyodide' — run the verified Python in a Web Worker (deployed static site)
 * Both answer the identical DATA_CONTRACT.md, so the rest of the stage is
 * backend-agnostic.
 */

// Injected by vite (see vite.config.ts `define`). Defaults to http.
declare const __BACKEND__: string;
const BACKEND: 'http' | 'pyodide' = __BACKEND__ === 'pyodide' ? 'pyodide' : 'http';

export function backendName(): 'http' | 'pyodide' {
  return BACKEND;
}

export interface ApiResult {
  status: number;
  data: any;
}

type StatusCb = (msg: string, ready: boolean, error?: string) => void;
let statusCb: StatusCb | null = null;
/** Subscribe to Pyodide boot progress (no-op in http mode). */
export function onBackendStatus(cb: StatusCb): void {
  statusCb = cb;
}

// ---------------------------------------------------------------------------
// Pyodide worker client
// ---------------------------------------------------------------------------

let worker: Worker | null = null;
let nextId = 1;
const pending = new Map<number, (r: ApiResult) => void>();

function ensureWorker(): Worker {
  if (worker) return worker;
  // py-worker.js lives next to index.html (public/); resolve via the document
  // base so it works at the site root OR a project subpath (GitHub Pages).
  worker = new Worker(new URL('py-worker.js', document.baseURI));
  worker.onmessage = (e: MessageEvent) => {
    const m = e.data;
    if (m.type === 'progress') return void statusCb?.(m.msg, false);
    if (m.type === 'ready') return void statusCb?.('ready', true);
    if (m.type === 'error') return void statusCb?.('failed', false, m.msg);
    const resolve = pending.get(m.id);
    if (resolve) {
      pending.delete(m.id);
      resolve({ status: m.status, data: m.data });
    }
  };
  worker.onerror = (e) => {
    e.preventDefault?.();
    statusCb?.('failed', false, e?.message || 'worker error (see console)');
  };
  return worker;
}

/** Start Pyodide booting now (so the splash can show progress). */
export function warmup(): void {
  if (BACKEND === 'pyodide') ensureWorker();
}

function requestPyodide(op: string, payload: unknown): Promise<ApiResult> {
  const w = ensureWorker();
  return new Promise((resolve) => {
    const id = nextId++;
    pending.set(id, resolve);
    w.postMessage({ id, op, payload });
  });
}

// ---------------------------------------------------------------------------
// HTTP client (local engine)
// ---------------------------------------------------------------------------

async function requestHttp(op: string, payload: any): Promise<ApiResult> {
  let url = '';
  let init: RequestInit | undefined;
  const post = (b: unknown): RequestInit => ({
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(b),
  });
  switch (op) {
    case 'health':
      url = '/api/health';
      break;
    case 'walk':
      url = '/api/walk';
      init = post(payload);
      break;
    case 'synthesis':
      url = '/api/synthesis';
      init = post(payload);
      break;
    case 'vertex':
      url = payload.addr ? `/api/vertex/${encodeURIComponent(payload.addr)}` : '/api/vertex';
      break;
    case 'random_word':
      url = `/api/random_word?length=${payload.length}&seed=${payload.seed}`;
      break;
    case 'equal':
      url = '/api/equal';
      init = post(payload);
      break;
  }
  const resp = await fetch(url, init);
  let data: any = {};
  try {
    data = await resp.json();
  } catch {
    data = {};
  }
  return { status: resp.status, data };
}

// ---------------------------------------------------------------------------

function request(op: string, payload: any = {}): Promise<ApiResult> {
  return BACKEND === 'pyodide' ? requestPyodide(op, payload) : requestHttp(op, payload);
}

// -- convenience wrappers used by the stage --------------------------------

export async function apiWalk(body: { word?: string; matrix?: string[][] }): Promise<any> {
  return (await request('walk', body)).data;
}
export async function apiRandomWord(length: number, seed: number): Promise<any> {
  return (await request('random_word', { length, seed })).data;
}
export async function apiHealth(): Promise<any> {
  return (await request('health')).data;
}
export async function apiVertex(addr: string): Promise<ApiResult> {
  return request('vertex', { addr });
}
