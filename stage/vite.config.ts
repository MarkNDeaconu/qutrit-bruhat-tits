import { defineConfig } from 'vite';

// Two build modes:
//   • default              → backend 'http'    (local engine on :8137, `make demo`)
//   • VITE_BACKEND=pyodide → backend 'pyodide' (verified Python in a worker; deploy)
// base './' keeps assets relative so the deployed build works both at a site
// root (Vercel) and under a project subpath (GitHub Pages). Override with
// VITE_BASE for an absolute base if desired.
export default defineConfig(() => {
  const backend = process.env.VITE_BACKEND === 'pyodide' ? 'pyodide' : 'http';
  return {
    base: process.env.VITE_BASE || './',
    define: {
      __BACKEND__: JSON.stringify(backend),
    },
    server: {
      proxy: {
        '/api': 'http://127.0.0.1:8137',
      },
    },
    build: {
      target: 'es2022',
    },
  };
});
