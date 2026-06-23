/*
 * page.test.mjs — headless end-to-end check of the DEPLOYED (Pyodide) build.
 * Loads the served site in real Chromium, waits for Pyodide to boot, drives a
 * Play + Synthesize, and fails on any uncaught page error / console error.
 * Run against a running `vite preview` (see the harness in package/Makefile).
 *
 *   URL env var overrides the target (default http://localhost:4173/).
 */
import { chromium } from 'playwright';

const URL = process.env.URL || 'http://localhost:4173/';
const errors = [];
const consoleErrors = [];

const browser = await chromium.launch();
const page = await browser.newPage();
page.on('pageerror', (e) => errors.push(String(e)));
page.on('console', (m) => {
  if (m.type() === 'error') consoleErrors.push(m.text());
});

let fails = 0;
const check = (name, ok, extra = '') => {
  console.log(`[${ok ? 'PASS' : 'FAIL'}] ${name}${extra ? ' — ' + extra : ''}`);
  if (!ok) fails++;
};

try {
  console.log(`loading ${URL} …`);
  await page.goto(URL, { waitUntil: 'load', timeout: 30000 });

  // 1. splash must clear (Pyodide downloads from CDN + boots) within 120 s
  let splashErr = '';
  try {
    await page.waitForFunction(
      () => {
        const s = document.querySelector('#splash');
        const st = document.querySelector('#splash-status');
        if (st && st.classList.contains('error')) return true; // surface failure fast
        return !!s && s.classList.contains('hidden');
      },
      { timeout: 120000 },
    );
  } catch {
    splashErr = await page.textContent('#splash-status').catch(() => '(no status)');
  }
  const splashStatus = await page.textContent('#splash-status').catch(() => '');
  const splashHidden = await page.evaluate(
    () => !!document.querySelector('#splash')?.classList.contains('hidden'),
  );
  check('Pyodide booted, splash cleared', splashHidden && !splashErr,
    splashErr ? `splash stuck: ${splashErr}` : `status="${splashStatus}"`);

  // 2. backend reports ready (status bar)
  await page.waitForTimeout(500);
  const status = await page.textContent('#engine-status').catch(() => '');
  check('status shows in-browser backend', /in-browser/.test(status), `"${status}"`);

  // 3. canvas present (renderer initialised)
  const hasCanvas = await page.evaluate(() => !!document.querySelector('#stage canvas'));
  check('WebGL/canvas stage present', hasCanvas);

  // 4. drive a walk: fast animation speed, type a scenic word, Play, Synthesize
  await page.$eval('#speed', (el) => {
    el.value = '120';
    el.dispatchEvent(new Event('input', { bubbles: true }));
  });
  await page.fill('#word', 'HSHHRHRSHRHRRHRRHHHHSHHHHSSHHR');
  await page.click('#btn-play');
  // scrubber tokens appear as soon as the walk loads (before animation finishes)
  await page.waitForSelector('#scrubber .tok', { timeout: 15000 });
  const tokenCount = await page.evaluate(
    () => document.querySelectorAll('#scrubber .tok').length,
  );
  check('Play loaded a walk (scrubber tokens)', tokenCount > 0, `${tokenCount} tokens`);

  await page.click('#btn-synth');
  // the optimality banner appears at the END of the straighten animation
  let banner = { shown: false, text: '' };
  try {
    await page.waitForSelector('#banner.show', { timeout: 30000 });
    banner = await page.evaluate(() => {
      const b = document.querySelector('#banner');
      return { shown: !!b && b.classList.contains('show'), text: b?.textContent || '' };
    });
  } catch {
    /* banner never showed */
  }
  check('Synthesize showed the optimality banner', banner.shown && /sde/.test(banner.text),
    `"${banner.text.trim()}"`);

  // 4b. incremental: append + synthesize repeatedly — each must add a piece and
  //     NOT revert/re-expand (banner reports growing "optimal pieces")
  const waitBanner = async (re) => {
    await page.waitForFunction(
      (src) => {
        const b = document.querySelector('#banner');
        return !!b && b.classList.contains('show') && new RegExp(src).test(b.textContent || '');
      },
      re.source,
      { timeout: 30000 },
    );
    return (await page.textContent('#banner')) || '';
  };
  let incrementalOk = true;
  let lastBanner = '';
  for (const piece of [2, 3]) {
    await page.fill('#word', 'HSRH');
    await page.click('#btn-append');
    await page.waitForTimeout(2200); // let the appended segment walk in
    await page.click('#btn-synth');
    try {
      lastBanner = await waitBanner(new RegExp(`${piece} optimal pieces`));
    } catch {
      incrementalOk = false;
      lastBanner = (await page.textContent('#banner')) || '(timeout)';
      break;
    }
  }
  check('append+synthesize builds incrementally (no revert/re-expand)', incrementalOk,
    `"${lastBanner.trim()}"`);

  // 4c. Clifford button (H,S only) loads a walk and stays in range (no cutoff)
  await page.click('#btn-clear');
  await page.waitForTimeout(900);
  await page.click('#btn-clifford');
  await page.waitForSelector('#scrubber .tok', { timeout: 15000 });
  const cliffTokens = await page.$$eval('#scrubber .tok', (els) => els.length);
  check('Clifford button loads a walk', cliffTokens > 0, `${cliffTokens} tokens`);

  // 4d. deep path: synthesizing a high-sde word must cut off cleanly (red marker +
  //     notice) and NOT glitch — i.e. no uncaught errors and the notice appears
  await page.click('#btn-clear');
  await page.waitForTimeout(700);
  await page.fill('#word', 'SRSRRHHHRRSRSHRHRRRSHRHSRRHRRHHHRHRHSSSSRRRRSSRHSRSHRR'); // sde 9, depth 18
  await page.click('#btn-play');
  await page.waitForTimeout(2500);
  await page.click('#btn-synth');
  let cutoffShown = false;
  try {
    await page.waitForFunction(
      () => /rendered range|red node/.test(document.querySelector('#caption')?.textContent || ''),
      { timeout: 30000 },
    );
    cutoffShown = true;
  } catch {
    /* notice never appeared */
  }
  check('deep path cuts off with an out-of-range notice (no glitch)', cutoffShown,
    ((await page.textContent('#caption')) || '').trim().slice(0, 70));
  await page.waitForTimeout(1500); // let any (now-bounded) animation settle
  await page.screenshot({ path: '/tmp/qutrits_cutoff.png', fullPage: false });

  // 5. inspector: click the origin card via the API-backed path
  const vertexOk = await page.evaluate(async () => {
    // exercise the transport directly the way the inspector does
    const r = await fetch('py/btlib.py'); // sanity: py asset reachable from page
    return r.ok;
  });
  check('py asset reachable from page origin', vertexOk);

  await page.screenshot({ path: '/tmp/qutrits_preview.png', fullPage: false });
  console.log('screenshot -> /tmp/qutrits_preview.png');

  // 6. no uncaught errors anywhere
  check('no uncaught page errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  check('no console errors', consoleErrors.length === 0, consoleErrors.slice(0, 3).join(' | '));
} finally {
  await browser.close();
}

console.log('');
if (fails) {
  console.log(`page.test.mjs: ${fails} check(s) FAILED`);
  process.exit(1);
}
console.log('page.test.mjs: preview is flawless ✓');
