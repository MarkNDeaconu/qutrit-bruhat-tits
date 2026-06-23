# "The Tree Behind the Gates" — local engine + stage
# All Python runs from the project .venv; nothing touches the system.

ROOT    := /Users/markdeaconu/projects/qutrits_v2
PY      := $(ROOT)/.venv/bin/python
UVICORN := $(ROOT)/.venv/bin/uvicorn
HOST    := 127.0.0.1
PORT    := 8137
URL     := http://$(HOST):$(PORT)

.PHONY: engine stage-dev stage demo test web preview-web verify-web

# Run the FastAPI engine in the foreground (Ctrl-C stops it).
engine:
	$(UVICORN) engine.app:app --host $(HOST) --port $(PORT)

# Build the DEPLOYABLE static site: the verified Python runs in-browser via
# Pyodide (no backend). Output in stage/dist — upload anywhere static.
web:
	cd stage && npm install && npm run build:web

# Build the deployable site and serve it locally to test before deploying
# (this is the Pyodide build, NOT the local engine). Frees a stale server on
# 4173, builds, opens the browser, then serves in the foreground (Ctrl-C stops).
preview-web:
	cd stage && npm run build:web
	@lsof -ti tcp:4173 | xargs kill 2>/dev/null || true
	@( sleep 1.5; open http://localhost:4173 ) &
	cd stage && npx vite preview --port 4173 --strictPort

# Headless end-to-end verification of the deployable build: builds it, serves it,
# loads it in real Chromium, drives Play+Synthesize, fails on any page/console
# error. One-time browser install: cd stage && npx playwright install chromium
verify-web:
	cd stage && npm run build:web
	@lsof -ti tcp:4173 | xargs kill 2>/dev/null || true
	@cd stage && ( npx vite preview --port 4173 --strictPort >/tmp/qutrits_vp.log 2>&1 & ) ; \
	 sleep 3 ; \
	 node test/page.test.mjs ; S=$$? ; \
	 lsof -ti tcp:4173 | xargs kill 2>/dev/null || true ; \
	 exit $$S

# Vite dev server for the stage (hot reload; talks to the engine over CORS).
stage-dev:
	cd stage && npm run dev

# Install deps and build the stage into stage/dist (served by the engine at /).
stage:
	cd stage && npm install && npm run build

# One-command demo: build the stage if needed, start the engine, open the
# browser. The engine runs in the FOREGROUND, so Ctrl-C in this terminal
# stops the demo cleanly; the browser-opener is a one-shot background
# subshell that exits on its own after launching the default browser.
demo:
	@if [ ! -d stage/dist ]; then \
		if [ -d stage ]; then $(MAKE) stage; \
		else echo "note: no stage/ yet — the engine will serve a JSON hint at /"; fi; \
	fi
	@( sleep 1; open $(URL) ) &
	$(UVICORN) engine.app:app --host $(HOST) --port $(PORT)

# Verification + smoke test. Runs verify/v09_vizgen.py if present, then curls
# /api/health (against a running engine if one is on $(PORT), otherwise
# against a temporary engine started and killed just for the smoke).
test:
	@if [ -f verify/v09_vizgen.py ]; then \
		$(PY) verify/v09_vizgen.py; \
	else \
		echo "skip: verify/v09_vizgen.py not present"; \
	fi
	@if lsof -ti :$(PORT) >/dev/null 2>&1; then \
		echo "smoke: using running engine on :$(PORT)"; \
		curl -sf $(URL)/api/health && echo; \
	else \
		echo "smoke: starting temporary engine on :$(PORT)"; \
		$(UVICORN) engine.app:app --host $(HOST) --port $(PORT) & \
		ENGINE_PID=$$!; \
		sleep 2; \
		curl -sf $(URL)/api/health; STATUS=$$?; echo; \
		kill $$ENGINE_PID 2>/dev/null; \
		exit $$STATUS; \
	fi
