"""vizgen/webapi.py — the single request dispatch shared by BOTH backends:

  * the local FastAPI engine (engine/app.py), and
  * the in-browser Pyodide worker (stage/public/py-worker.js).

One implementation means the deployed static site and the local demo run
identical code paths over the same frozen, verified arithmetic (btlib + vizgen).
Pure stdlib so it runs unchanged inside Pyodide.

`handle(op, payload) -> (status, data)`:
  status is an HTTP-style int (200 ok, 400 bad input); data is a JSON-safe dict.
Mirrors DATA_CONTRACT.md exactly. walk/synthesis return 200 with
{ok: False, error} on compute failure (the stage displays it); vertex/equal/
random_word return 400 on bad input.
"""

from __future__ import annotations
import os
import re
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import btlib
from vizgen import walks
from vizgen.serialize import parse_matrix

_ADDR_RE = re.compile(r"[0-3][0-2]{0,23}")  # first digit 0-3, rest 0-2, len<=24

# Live sanity check of the frozen core, computed once at import:
# the radius-3 ball must have sphere sizes [1, 4, 12, 36].
_levels, _parent, _cycle_edges = btlib.bfs_tree(3)
BALL_CHECK = [len(level) for level in _levels]


def _extract_input(payload: dict):
    if not isinstance(payload, dict):
        raise ValueError("body must be a JSON object with 'word' or 'matrix'")
    word = payload.get("word")
    matrix = payload.get("matrix")
    if (word is None) == (matrix is None):
        raise ValueError("provide exactly one of 'word' or 'matrix'")
    if word is not None and not isinstance(word, str):
        raise ValueError("'word' must be a string over H, S, R")
    return word, matrix


def _build_unitary(word, matrix):
    if word is not None:
        return walks.build_unitary(word.strip().upper())
    return parse_matrix(matrix)


def handle(op: str, payload: dict | None = None):
    """Dispatch one request. Returns (status:int, data:dict)."""
    payload = payload or {}

    if op == "health":
        return 200, {"ok": True, "btlib": "frozen", "ball_check": BALL_CHECK}

    if op == "walk":
        try:
            word, matrix = _extract_input(payload)
            return 200, walks.trajectory(word=word, matrix=matrix)
        except Exception as e:  # contract: stage displays, never crashes
            return 200, {"ok": False, "error": str(e)}

    if op == "synthesis":
        try:
            word, matrix = _extract_input(payload)
            traj = walks.trajectory(word=word, matrix=matrix)
            # trajectory() already validated the input; rebuild U for the oracle
            traj["oracle"] = walks.synthesis_oracle(_build_unitary(word, matrix))
            return 200, traj
        except Exception as e:
            return 200, {"ok": False, "error": str(e)}

    if op == "vertex":
        addr = payload.get("addr", "")
        if addr != "" and not _ADDR_RE.fullmatch(addr):
            return 400, {"error": "bad address: digits only, first in 0-3, "
                                  "rest in 0-2, length <= 24"}
        try:
            return 200, walks.vertex_card(addr)
        except Exception as e:
            return 400, {"error": str(e)}

    if op == "random_word":
        try:
            length = max(1, min(200, int(payload.get("length", 40))))
            return 200, {"word": walks.random_word(length, int(payload.get("seed", 0)))}
        except Exception as e:
            return 400, {"error": str(e)}

    if op == "equal":
        try:
            a, b = payload.get("a"), payload.get("b")
            if a is None or b is None:
                raise ValueError("body must contain 'a' and 'b'")
            if isinstance(a, str):
                a = {"word": a}
            if isinstance(b, str):
                b = {"word": b}
            return 200, walks.equality(a, b)
        except Exception as e:
            return 400, {"error": str(e)}

    return 400, {"error": f"unknown op {op!r}"}
