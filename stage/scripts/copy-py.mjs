// Copy the verified Python (btlib.py + vizgen/*.py) into stage/public/py/ so
// the Pyodide worker can load them at runtime. Run before `vite build` for the
// deployable build. Paths resolve relative to THIS file, so cwd doesn't matter.
//
// This is a verbatim copy at build time — the deployed site runs the exact same
// frozen, verified arithmetic as the local engine. No fork, ever.

import { mkdirSync, copyFileSync, rmSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url)); // stage/scripts
const repo = resolve(here, '..', '..'); // project root
const stage = resolve(here, '..'); // stage/
const dest = resolve(stage, 'public', 'py');

const FILES = [
  'btlib.py',
  'vizgen/__init__.py',
  'vizgen/treegen.py',
  'vizgen/serialize.py',
  'vizgen/walks.py',
  'vizgen/webapi.py',
];

rmSync(dest, { recursive: true, force: true });
mkdirSync(resolve(dest, 'vizgen'), { recursive: true });

for (const f of FILES) {
  copyFileSync(resolve(repo, f), resolve(dest, f));
}

console.log(`copy-py: ${FILES.length} files → ${dest}`);
