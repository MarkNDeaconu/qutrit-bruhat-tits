"""Local FastAPI engine for "The Tree Behind the Gates".

Implements exactly the REST endpoints of DATA_CONTRACT.md. All exact
arithmetic lives in btlib (frozen) + vizgen; this module only validates
inputs, dispatches, and shapes errors:

- /api/walk and /api/synthesis return HTTP 200 with {ok: false, error: str}
  on any compute-layer failure (the stage displays the message);
- everything else returns HTTP 400 on bad input.

Endpoints are plain `def` handlers: FastAPI runs them on its threadpool, so
the event loop is never blocked by the (sub-second) exact computations.
"""

from __future__ import annotations
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from fastapi import Body, FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from vizgen import webapi  # the dispatch shared with the Pyodide worker

STAGE_DIST = os.path.join(_ROOT, "stage", "dist")
_HAS_STAGE = os.path.isdir(STAGE_DIST)

app = FastAPI(title="qutrits engine", docs_url=None, redoc_url=None)

# The Vite dev server runs on another localhost port: allow all localhost
# origins (any port, http or https).
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=r"^https?://(localhost|127\.0\.0\.1)(:\d+)?$",
    allow_methods=["*"],
    allow_headers=["*"],
)


# --- endpoints ------------------------------------------------------------------
# Thin HTTP shells over webapi.handle (the same dispatch the Pyodide worker
# uses). walk/synthesis surface compute errors as 200 + {ok:false}; the rest
# raise HTTP 400 on bad input.

@app.get("/api/health")
def health():
    return webapi.handle("health", {})[1]


@app.post("/api/walk")
def walk(payload: dict = Body(...)):
    return webapi.handle("walk", payload)[1]


@app.post("/api/synthesis")
def synthesis(payload: dict = Body(...)):
    return webapi.handle("synthesis", payload)[1]


@app.get("/api/vertex/{addr}")
def vertex(addr: str):
    status, data = webapi.handle("vertex", {"addr": addr})
    if status != 200:
        raise HTTPException(status_code=status, detail=data.get("error"))
    return data


@app.get("/api/vertex")
def vertex_origin():
    """The origin's inspector card (the empty address cannot be a path param)."""
    return webapi.handle("vertex", {"addr": ""})[1]


@app.get("/api/random_word")
def random_word(length: int = Query(40), seed: int = Query(0)):
    return webapi.handle("random_word", {"length": length, "seed": seed})[1]


@app.post("/api/equal")
def equal(payload: dict = Body(...)):
    status, data = webapi.handle("equal", payload)
    if status != 200:
        raise HTTPException(status_code=status, detail=data.get("error"))
    return data


# --- static stage ---------------------------------------------------------------
# Mounted last so the /api/* routes above take precedence.

if _HAS_STAGE:
    app.mount("/", StaticFiles(directory=STAGE_DIST, html=True), name="stage")
else:
    @app.get("/")
    def stage_missing():
        return {
            "ok": False,
            "hint": "stage/dist not found — build the stage with `make stage`, "
                    "or use the API directly under /api/*",
        }
