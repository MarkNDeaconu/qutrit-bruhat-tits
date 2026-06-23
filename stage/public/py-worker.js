/*
 * py-worker.js — runs the verified arithmetic in the browser.
 *
 * Loads Pyodide (CPython → WASM), writes btlib.py + vizgen/*.py into its
 * virtual filesystem verbatim, and dispatches requests through the SAME
 * vizgen.webapi.handle that the local FastAPI engine uses. The frontend talks
 * to this worker exactly as it would the HTTP engine (DATA_CONTRACT.md).
 *
 * Classic worker (importScripts) so Pyodide loads straight from the CDN with no
 * bundler involvement. The .py files are served as static assets next to this
 * script (copy-py.mjs puts them in public/py/).
 */

const PYODIDE_VERSION = '0.26.4';
const PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;

const PY_FILES = [
  'btlib.py',
  'vizgen/__init__.py',
  'vizgen/treegen.py',
  'vizgen/serialize.py',
  'vizgen/walks.py',
  'vizgen/webapi.py',
];

let pyHandle = null; // PyProxy of the dispatch function

async function boot() {
  postMessage({ type: 'progress', msg: 'downloading Pyodide…' });
  importScripts(`${PYODIDE_URL}pyodide.js`);
  const pyodide = await loadPyodide({ indexURL: PYODIDE_URL });

  postMessage({ type: 'progress', msg: 'loading verified arithmetic…' });
  pyodide.FS.mkdirTree('/py/vizgen');
  for (const f of PY_FILES) {
    // resolve relative to THIS worker's URL → works at site root or a subpath
    const url = new URL('py/' + f, self.location.href).href;
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`failed to fetch ${f} (${resp.status})`);
    pyodide.FS.writeFile('/py/' + f, await resp.text());
  }

  pyodide.runPython('import sys; sys.path.insert(0, "/py")');
  pyHandle = pyodide.runPython(`
import json
from vizgen import webapi

def _handle(op, payload_json):
    payload = json.loads(payload_json) if payload_json else {}
    status, data = webapi.handle(op, payload)
    return json.dumps({"status": status, "data": data})

_handle
`);

  // warm the lazy caches (coset reps, address tree) so the first real call is fast
  pyHandle('health', '{}');
  pyHandle('walk', JSON.stringify({ word: 'H' }));

  postMessage({ type: 'ready' });
}

const bootPromise = boot().catch((err) => {
  postMessage({ type: 'error', msg: String(err && err.message ? err.message : err) });
  throw err;
});

self.onmessage = async (e) => {
  const { id, op, payload } = e.data || {};
  if (id === undefined) return;
  try {
    await bootPromise;
    const out = pyHandle(op, JSON.stringify(payload || {}));
    const parsed = JSON.parse(out);
    postMessage({ id, status: parsed.status, data: parsed.data });
  } catch (err) {
    postMessage({
      id,
      status: 500,
      data: { ok: false, error: String(err && err.message ? err.message : err) },
    });
  }
};
