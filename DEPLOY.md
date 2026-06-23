# Deploying "The Tree Behind the Gates"

The deployable build is a **fully static site** that runs the verified Python
(`btlib.py` + `vizgen/`) in the browser via Pyodide (WASM CPython) in a Web
Worker. No backend, no server, no maintenance. It computes **byte-identical**
results to the local engine — proven by `npm run test:web` (16 fixtures, WASM
CPython vs native CPython).

Your local laptop/HDMI demo is unchanged: `make demo` still runs the fast
CPython FastAPI engine. Deployment is a separate build target.

```
make web          # build the deployable static site → stage/dist
make preview-web  # build it, open the browser, serve locally (Ctrl-C to stop)
make verify-web   # headless end-to-end check: builds, serves, drives it in
                  # real Chromium, fails on any page/console error
```

`make preview-web` builds, **opens your browser automatically**, and serves the
deployable build. Wait ~2–4 s for the splash ("downloading Pyodide…" → "ready"),
then it works exactly as the local demo — it's the same artifact you deploy.

`make verify-web` is the automated proof it's flawless (used to validate this
build: Pyodide boots, a walk plays, Synthesize lands the geodesic, zero errors).
One-time browser install for it: `cd stage && npx playwright install chromium`.

## Deploy to Vercel (recommended)

`vercel.json` (repo root) is already configured: it builds `stage` with
`build:web` and serves `stage/dist`.

**Dashboard:** push this repo to GitHub, then on vercel.com → *Add New Project*
→ import the repo → Deploy. Vercel reads `vercel.json`; no other settings
needed. (Leave the Root Directory as the repo root.)

**CLI:** from the repo root,
```
npx vercel          # preview deploy (first run links/creates the project)
npx vercel --prod   # production deploy
```

That's it — you get a permanent `*.vercel.app` URL (add a custom domain in the
dashboard if you like).

## Deploy to GitHub Pages

Pages serves under a project subpath (`https://<user>.github.io/<repo>/`), so
build with a matching base, then publish `stage/dist`:

```
cd stage
VITE_BASE=/<repo>/ npm run build:web      # e.g. VITE_BASE=/qutrits_v2/
npx gh-pages -d dist                       # publish to the gh-pages branch
```

(Then enable Pages → branch `gh-pages` in the repo settings.) A `<user>.github.io`
*user* repo serves at the root, so there you can use the default
`make web` with no `VITE_BASE`.

Or automate it with a workflow (`.github/workflows/pages.yml`):
```yaml
name: deploy-pages
on: { push: { branches: [main] } }
permissions: { contents: read, pages: write, id-token: write }
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with: { node-version: 20 }
      - run: npm --prefix stage install
      - run: VITE_BASE=/${{ github.event.repository.name }}/ npm --prefix stage run build:web
      - uses: actions/upload-pages-artifact@v3
        with: { path: stage/dist }
  deploy:
    needs: build
    runs-on: ubuntu-latest
    environment: { name: github-pages }
    steps: [{ uses: actions/deploy-pages@v4 }]
```

## Notes

- **First load needs internet + a CDN trust.** The deployed worker
  `importScripts` Pyodide (~6.5 MB, browser-cached) from the jsDelivr CDN at
  runtime. A classic worker cannot use Subresource Integrity, so the deploy
  trusts jsDelivr (a jsDelivr compromise could run code in the worker). jsDelivr
  serves immutable versioned mirrors, so availability is not the concern —
  integrity is. **For a hardened/offline/auditable deploy, vendor Pyodide** (this
  is also the supply-chain fix): `npm i pyodide`, copy these files from
  `node_modules/pyodide/` into `stage/public/pyodide/` — `pyodide.js`,
  `pyodide.asm.js`, `pyodide.asm.wasm`, `python_stdlib.zip`,
  `pyodide-lock.json` — and point `PYODIDE_URL` in `stage/public/py-worker.js`
  at `./pyodide/`. (Your local `make demo` is already fully offline.)
- **Don't set `NODE_ENV=production` in Vercel project settings.** It makes `npm`
  skip devDependencies (tsc/vite/pyodide) and the build fails. `vercel.json`
  already passes `--include=dev` to guard against this.
- **GitHub Pages project sites: use the trailing slash.** The worker resolves
  its assets via `document.baseURI`, so visit `…github.io/<repo>/` (with the
  slash; Pages redirects to add it). Never add a `<base>` tag.
- **Pyodide version** is pinned to 0.26.4 in `stage/public/py-worker.js` (CDN)
  and matched by the `pyodide` devDependency used by `test:web`. Bump both
  together.
- **Math stays unforked.** `stage/scripts/copy-py.mjs` copies the verified
  `.py` files verbatim into the build at deploy time; the browser runs the same
  bytes as the local engine. `make web` always re-copies, so the deployed site
  can never drift from the verified source.
- **No special headers** are required (Pyodide 0.26 runs without
  cross-origin isolation), so it works on plain static hosts including Pages.
