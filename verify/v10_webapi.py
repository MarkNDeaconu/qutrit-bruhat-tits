"""v10_webapi.py — the shared dispatch (vizgen.webapi) must be faithful to the
underlying verified functions (vizgen.walks) AND must be the single source of
truth for both backends (engine + Pyodide worker). Deterministic; exits nonzero
on any failure.
"""

import json
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from vizgen import webapi, walks

FAILS = []
def check(name, cond, extra=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")
    if not cond:
        FAILS.append(name)


# 1. every response is JSON-serializable (must cross the worker boundary)
for op, payload in [
    ("health", {}),
    ("walk", {"word": "HSHHHSH"}),
    ("synthesis", {"word": "HSRHSR"}),
    ("vertex", {"addr": "21"}),
    ("vertex", {"addr": ""}),
    ("random_word", {"length": 40, "seed": 7}),
    ("equal", {"a": "HSH", "b": "HSH"}),
]:
    status, data = webapi.handle(op, payload)
    try:
        json.loads(json.dumps(data))
        ok = True
    except Exception:
        ok = False
    check(f"JSON-serializable: {op}", ok and isinstance(status, int))

# 2. walk == walks.trajectory exactly (the dispatch must not alter the math)
import random
rng = random.Random(2026)
mismatch = 0
for _ in range(40):
    w = "".join(rng.choice("HSR") for _ in range(rng.randint(1, 40)))
    s, d = webapi.handle("walk", {"word": w})
    direct = walks.trajectory(word=w)
    if s != 200 or json.dumps(d, sort_keys=True) != json.dumps(direct, sort_keys=True):
        mismatch += 1
check("walk dispatch == walks.trajectory (40 random words)", mismatch == 0,
      f"{mismatch} mismatches")

# 3. health ball_check is the verified [1,4,12,36]
s, d = webapi.handle("health", {})
check("health ball_check == [1,4,12,36]", s == 200 and d.get("ball_check") == [1, 4, 12, 36],
      str(d.get("ball_check")))

# 4. error shaping matches the contract
s, d = webapi.handle("walk", {"word": "XYZ"})
check("walk bad letter -> 200 + ok:false", s == 200 and d.get("ok") is False and "error" in d)
s, d = webapi.handle("walk", {})
check("walk no input -> 200 + ok:false", s == 200 and d.get("ok") is False)
s, d = webapi.handle("walk", {"word": "H", "matrix": [["1"]]})
check("walk both inputs -> 200 + ok:false", s == 200 and d.get("ok") is False)
s, d = webapi.handle("vertex", {"addr": "9"})
check("vertex bad addr -> 400", s == 400 and "error" in d)
s, d = webapi.handle("vertex", {"addr": "4"})
check("vertex first-digit 4 -> 400", s == 400)
s, d = webapi.handle("unknown_op", {})
check("unknown op -> 400", s == 400)

# 5. vertex ok cases
s, d = webapi.handle("vertex", {"addr": ""})
check("origin vertex card ok", s == 200 and d.get("kind") == "P" and d.get("depth") == 0)
s, d = webapi.handle("vertex", {"addr": "21"})
check("vertex 21 card ok (pure, 4 f3 lines)", s == 200 and d.get("kind") == "P"
      and len(d.get("f3Data", [])) == 4)
s, d = webapi.handle("vertex", {"addr": "2"})
check("vertex 2 card ok (alternating, 4 planes)", s == 200 and d.get("kind") == "A"
      and len(d.get("f3Data", [])) == 4)

# 6. matrix input path (the H gate, round-tripped through serialize)
from vizgen.serialize import mat_display
from btlib import H_GATE
rows = [[e["str"] for e in r] for r in mat_display(H_GATE)["entries"]]
s, d = webapi.handle("walk", {"matrix": rows})
check("matrix walk (H) -> sde 1", s == 200 and d.get("sde") == 1 and d.get("ok") is not False)

# 7. synthesis adds an oracle with the verified 1-2-9 split at every step
s, d = webapi.handle("synthesis", {"word": walks.random_word(24, 3)})
ok = s == 200 and "oracle" in d
for step in d.get("oracle", []):
    ch = [c["lChange"] for c in step["candidates"]]
    if not (ch.count(-2) == 1 and ch.count(0) == 2 and ch.count(2) == 9):
        ok = False
check("synthesis oracle 1-2-9 split at every step", ok, f"{len(d.get('oracle', []))} steps")

# 8. random_word clamps and is deterministic
s1, d1 = webapi.handle("random_word", {"length": 40, "seed": 7})
s2, d2 = webapi.handle("random_word", {"length": 40, "seed": 7})
check("random_word deterministic", d1 == d2 and len(d1["word"]) == 40)
_, dlo = webapi.handle("random_word", {"length": 0, "seed": 1})
_, dhi = webapi.handle("random_word", {"length": 9999, "seed": 1})
check("random_word clamps to [1,200]", len(dlo["word"]) == 1 and len(dhi["word"]) == 200)

print()
if FAILS:
    print(f"FAILURES ({len(FAILS)}): {FAILS}")
    sys.exit(1)
print("ALL WEBAPI CHECKS PASSED")
