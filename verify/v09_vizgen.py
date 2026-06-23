#!/usr/bin/env python
"""
v09_vizgen.py -- Verification that the NEW vizgen layer (treegen / walks /
serialize) is faithful to the frozen, verified btlib core. The data contract
(DATA_CONTRACT.md) is binding: addresses are digit strings (first digit 0-3,
later digits 0-2), kind = parity of depth, child order = canonical HNF key.

Checks (one PASS/FAIL line per claim, nonzero exit on any FAIL):

  1. Address determinism: the addr -> vertex-key map of ball(4) is built in
     TWO SEPARATE SUBPROCESSES (with different PYTHONHASHSEED) and must be
     byte-identical, and identical to the in-process map; children_of
     ordering is strictly increasing in canonical HNF key.
  2. Address <-> vertex faithfulness on ALL 161 addresses of depth <= 4
     (procedural enumeration, 4*3^(d-1) per depth): kind == parity of depth;
     even depth => btlib d_tilde(I, basis) == depth; odd depth => neighbor of
     its parent; parent/child relations agree with btlib neighbors();
     addresses of depth d biject onto the btlib BFS sphere(d), no cycles.
  3. Trajectory faithfulness vs btlib on 40 seeded random {H,S,R} words
     (lengths 5..60): sde == btlib sde(U); len(trail) == 2*#H + 1; every
     consecutive trail pair is parent/child; geodesic ==
     treegen.geodesic_addresses(U) with len == 2*sde + 1; replaying the
     reduction events on the trail yields the geodesic EXACTLY (each event a
     genuine backtrack, 'removed' fields exact); edgeOwner = each H ordinal
     exactly twice, consecutively; normal-form block count == sde.
  4. Geodesic vs btlib synthesis on 20 seeded words: the pure vertices of the
     geodesic (even positions) match the prefix products of the btlib
     synthesis_steps word, as Vertex keys, in order.
  5. Oracle vs v08 semantics on 15 seeded words: every step shows the exact
     1-2-9 split, the chosen candidate is the unique lChange == -2, #steps ==
     sde; chosen rep indices AND full per-step lChange vectors are recomputed
     independently with raw btlib calls (s0_coset_reps loop) and must match
     exactly; 'at' addresses run down the geodesic pures from U e0 to S_0.
  6. step_address candidate forms: 200 seeded random (vertex, H-step) pairs
     generated from random words; the incremental step_address result equals
     the slow address_of_pure(U) recomputation; all three candidate forms
     (grandparent / sibling / grandchild) are exercised.
  7. vertex_card sanity on all 53 addresses of depth <= 3: card builds,
     f3Data has exactly 4 entries, branches/structure match the contract;
     pure cards: parse the unitary back from the 'str' fields via
     serialize.parse_matrix and check in_Gamma and Vertex('P', parsed) ==
     vertex_at(addr); parse_fw/fw_str round-trip on 300 random Fw elements
     (denominators in {1,3,9,27,5,15}).
  8. Performance guard: trajectory(random word of length 60) under 5 s,
     measured in a FRESH subprocess (cold caches); sde cross-checked in-process.

Deterministic: all RNGs seeded; subprocesses get fixed (different) hash seeds.
"""

import json
import os
import random
import subprocess
import sys
import tempfile
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from btlib import (
    Zw, Fw, Mat, Vertex, ORIGIN,
    H_GATE, S_GATE, R_GATE,
    ell, sde, d_tilde, in_Gamma, neighbors, bfs_tree,
    s0_coset_reps, synthesis_steps,
)
from vizgen import treegen, walks, serialize

SEED = 20260610
PYTHON = sys.executable           # the .venv interpreter running this script
FAILURES = []

GATES = {"H": H_GATE, "S": S_GATE, "R": R_GATE}


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)
    return ok


def addresses_up_to(maxd):
    """All addresses of depth <= maxd (procedural: 4*3^(d-1) per depth d>=1)."""
    addrs, frontier = [""], [""]
    for _d in range(1, maxd + 1):
        frontier = [a + str(c) for a in frontier
                    for c in range(4 if a == "" else 3)]
        addrs += frontier
    return addrs


def word_unitary(w):
    """Product of btlib gate matrices (independent of walks.build_unitary)."""
    U = Mat.identity()
    for ch in w:
        U = U * GATES[ch]
    return U


# ===========================================================================
# Section 1: address determinism across separate subprocesses
# ===========================================================================
print("=" * 76)
print("Section 1: address determinism (two subprocesses) + child ordering")
print("=" * 76)

SNIPPET = r'''
import json, sys
sys.path.insert(0, sys.argv[2])
from vizgen import treegen
addrs, frontier = [""], [""]
for _d in range(1, 5):
    frontier = [a + str(c) for a in frontier
                for c in range(4 if a == "" else 3)]
    addrs += frontier
amap = {a: treegen.vertex_at(a).key() for a in addrs}
with open(sys.argv[1], "w") as f:
    json.dump(amap, f, sort_keys=True, separators=(",", ":"))
'''

blobs = []
sub_ok = True
t0 = time.perf_counter()
for run, hashseed in ((1, "0"), (2, "424242")):
    fd, path = tempfile.mkstemp(suffix=f".v09run{run}.json")
    os.close(fd)
    try:
        env = dict(os.environ, PYTHONHASHSEED=hashseed)
        r = subprocess.run([PYTHON, "-c", SNIPPET, path, ROOT],
                           env=env, capture_output=True, text=True)
        if r.returncode != 0:
            sub_ok = False
            print(f"    subprocess {run} failed:\n{r.stderr}")
            blobs.append(None)
        else:
            with open(path, "rb") as f:
                blobs.append(f.read())
    finally:
        os.unlink(path)
t1 = time.perf_counter()

maps = [json.loads(b) if b else None for b in blobs]
check("1a ball(4) address map identical across two subprocesses "
      "(different PYTHONHASHSEED), byte-for-byte",
      sub_ok and blobs[0] is not None and blobs[0] == blobs[1]
      and maps[0] is not None and len(maps[0]) == 161,
      f"entries={len(maps[0]) if maps[0] else '?'}, "
      f"bytes={len(blobs[0]) if blobs[0] else '?'}  ({t1-t0:.1f}s)")

ALL4 = addresses_up_to(4)
amap_local = {a: treegen.vertex_at(a).key() for a in ALL4}
amap_local_json = json.loads(json.dumps(amap_local))   # tuples -> lists
check("1b in-process ball(4) address map identical to the subprocess maps",
      maps[0] is not None and amap_local_json == maps[0],
      f"{len(amap_local)} addresses (1+4+12+36+108)")

n_ok, n_tot = 0, 0
for a in addresses_up_to(3):
    keys = [w.key() for w in treegen.children_of(a)]
    n_tot += 1
    n_ok += (len(keys) == (4 if a == "" else 3)
             and all(keys[i] < keys[i + 1] for i in range(len(keys) - 1)))
check("1c children_of ordering strictly increasing in canonical HNF key "
      "(all 53 vertices of depth <= 3)",
      n_ok == n_tot, f"{n_ok}/{n_tot}")


# ===========================================================================
# Section 2: address <-> vertex faithfulness on ball(4)
# ===========================================================================
print()
print("=" * 76)
print("Section 2: address <-> vertex faithfulness (all 161 addresses, depth <= 4)")
print("=" * 76)

depth_counts = Counter(len(a) for a in ALL4)
check("2a procedural enumeration sizes: depth d has 4*3^(d-1) addresses",
      all(depth_counts[d] == (1 if d == 0 else 4 * 3 ** (d - 1))
          for d in range(5)),
      f"counts={dict(sorted(depth_counts.items()))}")

n_kind = sum(treegen.vertex_at(a).kind == ("P" if len(a) % 2 == 0 else "A")
             and treegen.kind_of(a) == treegen.vertex_at(a).kind
             for a in ALL4)
check("2b kind == parity of depth (P even / A odd) for all 161 addresses",
      n_kind == len(ALL4), f"{n_kind}/{len(ALL4)}")

I3 = Mat.identity()
n_even, t_even = 0, 0
n_odd, t_odd = 0, 0
for a in ALL4:
    v = treegen.vertex_at(a)
    if len(a) % 2 == 0:
        t_even += 1
        n_even += (d_tilde(I3, v.g) == len(a))
    else:
        t_odd += 1
        p = treegen.vertex_at(a[:-1])
        n_odd += (p.key() in {w.key() for w in neighbors(v)})
check("2c even depth: btlib d_tilde(I, basis) == depth (tree distance to e0)",
      n_even == t_even, f"{n_even}/{t_even} pure addresses")
check("2d odd depth: vertex is a btlib neighbor of its parent (symmetric edge)",
      n_odd == t_odd, f"{n_odd}/{t_odd} alternating addresses")

n_ok, n_tot = 0, 0
for a in addresses_up_to(3):
    v = treegen.vertex_at(a)
    nbr_keys = {w.key() for w in neighbors(v)}
    child_keys = {w.key() for w in treegen.children_of(a)}
    if a:
        pk = treegen.vertex_at(a[:-1]).key()
        ok = (len(nbr_keys) == 4 and pk not in child_keys
              and nbr_keys == child_keys | {pk})
    else:
        ok = (len(nbr_keys) == 4 and nbr_keys == child_keys)
    n_tot += 1
    n_ok += ok
check("2e neighbors(v) == {parent} U children exactly (4 distinct), depth <= 3",
      n_ok == n_tot, f"{n_ok}/{n_tot}")

t0 = time.perf_counter()
levels, _parent, cycles = bfs_tree(4)
t1 = time.perf_counter()
by_depth = {}
for a in ALL4:
    by_depth.setdefault(len(a), []).append(treegen.vertex_at(a).key())
bij_ok = (cycles == [])
for d in range(5):
    ks = by_depth[d]
    bij_ok = bij_ok and (len(set(ks)) == len(ks)
                         and set(ks) == {v.key() for v in levels[d]})
check("2f addresses of depth d biject onto btlib BFS sphere(d), d <= 4; "
      "no cycle edges",
      bij_ok,
      f"sphere sizes {[len(l) for l in levels]}  ({t1-t0:.1f}s)")


# ===========================================================================
# Section 3: trajectory faithfulness vs btlib (40 seeded words, len 5..60)
# ===========================================================================
print()
print("=" * 76)
print("Section 3: trajectory() vs btlib on 40 seeded random words (len 5..60)")
print("=" * 76)

rng = random.Random(SEED + 3)
words3 = []
for _ in range(40):
    L = rng.randint(5, 60)
    words3.append("".join(rng.choice("HSR") for _ in range(L)))


def replay_reduction(trail, events):
    """Independent replay: each event must remove a genuine backtrack pair
    work[i], work[i+1] with work[i+1] == work[i-1], matching 'removed'."""
    work = list(trail)
    for e in events:
        i = e["index"]
        if not (1 <= i <= len(work) - 2):
            return None
        if work[i + 1] != work[i - 1]:
            return None
        if e.get("removed") != [work[i], work[i + 1]]:
            return None
        del work[i:i + 2]
    return work


def is_parent_child(x, y):
    return ((len(y) == len(x) + 1 and y[:-1] == x)
            or (len(x) == len(y) + 1 and x[:-1] == y))


n = Counter()
n_err = 0
sdes3 = []
t0 = time.perf_counter()
for w in words3:
    U = word_unitary(w)
    s = sde(U)            # btlib, directly
    sdes3.append(s)
    nH = w.count("H")
    try:
        traj = walks.trajectory(word=w)
    except Exception as exc:                                  # noqa: BLE001
        n_err += 1
        print(f"    trajectory({w!r}) raised: {exc!r}")
        continue
    trail, geo = traj["trail"], traj["geodesic"]
    n["sde"] += (traj["sde"] == s)
    n["trail_len"] += (len(trail) == 2 * nH + 1)
    n["parent_child"] += all(is_parent_child(trail[i], trail[i + 1])
                             for i in range(len(trail) - 1))
    n["geodesic"] += (geo == treegen.geodesic_addresses(U)
                      and len(geo) == 2 * s + 1)
    n["reduce"] += (replay_reduction(trail, traj["reduction"]) == geo)
    n["owner"] += (traj["edgeOwner"] == [i // 2 for i in range(2 * nH)])
    n["blocks"] += (len(traj["normalForm"]["blocks"]) == s)
    n["letters"] += ("".join(st["letter"] for st in traj["steps"]) == w)
t1 = time.perf_counter()

N3 = len(words3)
check("3a no exceptions; sde == btlib sde(U) for all 40 words",
      n_err == 0 and n["sde"] == N3,
      f"{n['sde']}/{N3}; sde range [{min(sdes3)},{max(sdes3)}]  ({t1-t0:.1f}s)")
check("3b len(trail) == 2*(#H letters) + 1", n["trail_len"] == N3,
      f"{n['trail_len']}/{N3}")
check("3c every consecutive trail pair is parent/child (one-digit extension)",
      n["parent_child"] == N3, f"{n['parent_child']}/{N3}")
check("3d geodesic == treegen.geodesic_addresses(U), len == 2*sde + 1",
      n["geodesic"] == N3, f"{n['geodesic']}/{N3}")
check("3e replaying reduction events on the trail yields the geodesic EXACTLY "
      "(every event a genuine backtrack, 'removed' fields exact)",
      n["reduce"] == N3, f"{n['reduce']}/{N3}")
check("3f edgeOwner: each H ordinal exactly twice, consecutively",
      n["owner"] == N3, f"{n['owner']}/{N3}")
check("3g normal-form block count == sde", n["blocks"] == N3,
      f"{n['blocks']}/{N3}")
check("3h step letters reproduce the input word in order",
      n["letters"] == N3, f"{n['letters']}/{N3}")


# ===========================================================================
# Section 4: geodesic pures == prefix products of btlib synthesis word
# ===========================================================================
print()
print("=" * 76)
print("Section 4: geodesic vs btlib synthesis_steps (20 seeded words)")
print("=" * 76)

rng = random.Random(SEED + 4)
n_ok, n_tot = 0, 0
sdes4 = []
for _ in range(20):
    L = rng.randint(5, 50)
    w = "".join(rng.choice("HSR") for _ in range(L))
    U = word_unitary(w)
    s = sde(U)
    sdes4.append(s)
    geo = treegen.geodesic_addresses(U)
    pures = geo[::2]
    word_m, _M = synthesis_steps(U)
    prefix_keys = [Vertex("P", Mat.identity()).key()]
    acc = Mat.identity()
    for m in word_m:
        acc = acc * (m * H_GATE)
        prefix_keys.append(Vertex("P", acc).key())
    addr_keys = [treegen.vertex_at(a).key() for a in pures]
    n_tot += 1
    n_ok += (len(word_m) == s and len(geo) == 2 * s + 1
             and addr_keys == prefix_keys
             and addr_keys[-1] == Vertex("P", U).key())
check("4a geodesic pure vertices (even positions) == prefix products of the "
      "btlib synthesis word, as Vertex keys, in order",
      n_ok == n_tot,
      f"{n_ok}/{n_tot}; sde range [{min(sdes4)},{max(sdes4)}]")


# ===========================================================================
# Section 5: synthesis_oracle vs raw btlib descent (v08 semantics)
# ===========================================================================
print()
print("=" * 76)
print("Section 5: synthesis_oracle vs raw btlib recomputation (15 seeded words)")
print("=" * 76)

REPS = s0_coset_reps()
MHinv = [(m * H_GATE).inv() for m in REPS]

rng = random.Random(SEED + 5)
n = Counter()
n_tot, steps_total = 0, 0
for _ in range(15):
    L = rng.randint(5, 45)
    w = "".join(rng.choice("HSR") for _ in range(L))
    U = word_unitary(w)
    s = sde(U)
    oracle = walks.synthesis_oracle(U)
    n_tot += 1
    n["count"] += (len(oracle) == s)

    # independent recomputation with raw btlib calls
    cur = U
    my_dls, my_chosen, unique = [], [], True
    while ell(cur) > 0:
        l0 = ell(cur)
        dls = [ell(MHinv[i] * cur) - l0 for i in range(12)]
        idxs = [i for i, d in enumerate(dls) if d == -2]
        if len(idxs) != 1:
            unique = False
            break
        my_dls.append(dls)
        my_chosen.append(idxs[0])
        cur = MHinv[idxs[0]] * cur
    steps_total += len(my_dls)

    split_ok = rep_ok = chosen_ok = match_ok = True
    or_chosen = []
    for st in oracle:
        cands = st["candidates"]
        rep_ok &= ([c["rep"] for c in cands] == list(range(12)))
        split_ok &= (Counter(c["lChange"] for c in cands)
                     == Counter({-2: 1, 0: 2, 2: 9}))
        ch = [c["rep"] for c in cands if c["chosen"]]
        chosen_ok &= (len(ch) == 1
                      and cands[ch[0]]["lChange"] == -2)
        or_chosen += ch
    match_ok = (unique and len(oracle) == len(my_dls)
                and or_chosen == my_chosen
                and all([c["lChange"] for c in oracle[k]["candidates"]]
                        == my_dls[k] for k in range(len(oracle))))
    n["split"] += split_ok
    n["reps"] += rep_ok
    n["chosen"] += chosen_ok
    n["match"] += match_ok

    pures = treegen.geodesic_addresses(U)[::2]
    n["at"] += ([st["at"] for st in oracle]
                == [pures[len(pures) - 1 - k] for k in range(len(oracle))])

check("5a oracle has exactly sde steps for all 15 words",
      n["count"] == n_tot, f"{n['count']}/{n_tot}; {steps_total} steps total")
check("5b every step lists reps 0..11 in order with the exact 1-2-9 split "
      "{-2:1, 0:2, +2:9}",
      n["split"] == n_tot and n["reps"] == n_tot,
      f"split {n['split']}/{n_tot}, rep order {n['reps']}/{n_tot}")
check("5c the chosen candidate is the unique lChange == -2 at every step",
      n["chosen"] == n_tot, f"{n['chosen']}/{n_tot}")
check("5d chosen rep indices AND full lChange vectors match an independent "
      "raw-btlib descent (s0_coset_reps loop) exactly",
      n["match"] == n_tot, f"{n['match']}/{n_tot}")
check("5e 'at' addresses run down the geodesic pures from U e0 toward e0",
      n["at"] == n_tot, f"{n['at']}/{n_tot}")


# ===========================================================================
# Section 6: step_address (incremental) == address_of_pure (slow recompute)
# ===========================================================================
print()
print("=" * 76)
print("Section 6: step_address candidate forms, 200 seeded (vertex, H-step) pairs")
print("=" * 76)

rng = random.Random(SEED + 6)
n_pairs, n_ok = 0, 0
forms = Counter()
t0 = time.perf_counter()
while n_pairs < 200:
    L = rng.randint(8, 30)
    w = "".join(rng.choice("HSR") for _ in range(L))
    acc = Mat.identity()
    cur_addr = ""
    for ch in w:
        acc = acc * GATES[ch]
        if ch != "H":
            continue
        got = treegen.step_address(cur_addr, Vertex("P", acc))
        slow = treegen.address_of_pure(acc)
        if len(got) == len(cur_addr) - 2:
            forms["grandparent"] += 1
        elif len(got) == len(cur_addr):
            forms["sibling"] += 1
        else:
            forms["grandchild"] += 1
        n_pairs += 1
        n_ok += (got == slow)
        cur_addr = got
        if n_pairs >= 200:
            break
t1 = time.perf_counter()
check("6a incremental step_address == slow address_of_pure recomputation "
      "on 200 H-steps",
      n_ok == n_pairs == 200, f"{n_ok}/{n_pairs}  ({t1-t0:.1f}s)")
check("6b all three candidate forms exercised (grandparent/sibling/grandchild)",
      set(forms) == {"grandparent", "sibling", "grandchild"},
      f"forms={dict(forms)}")


# ===========================================================================
# Section 7: vertex_card sanity + serializer round trips
# ===========================================================================
print()
print("=" * 76)
print("Section 7: vertex_card (depth <= 3) + parse_fw/fw_str round trip")
print("=" * 76)

ALL3 = addresses_up_to(3)
n = Counter()
n_err = 0
n_pure = 0
for a in ALL3:
    try:
        card = walks.vertex_card(a)
    except Exception as exc:                                  # noqa: BLE001
        n_err += 1
        print(f"    vertex_card({a!r}) raised: {exc!r}")
        continue
    v = treegen.vertex_at(a)
    n["build"] += 1
    n["f3"] += (len(card["f3Data"]) == 4)
    res = card["residue"]
    struct_ok = (card["addr"] == a and card["kind"] == v.kind
                 and card["depth"] == len(a)
                 and len(res) == 3 and all(len(r) == 3 for r in res)
                 and all(x in (0, 1, 2) for r in res for x in r)
                 and card["f3Kind"] == ("line" if v.kind == "P" else "plane"))
    br = card["branches"]
    if a:
        struct_ok &= (len(br) == 4 and br[0]["child"] is None
                      and br[0]["toward"] == a[:-1]
                      and [b["child"] for b in br[1:]] == [0, 1, 2]
                      and [b["toward"] for b in br[1:]]
                      == [a + str(i) for i in range(3)])
    else:
        struct_ok &= (len(br) == 4
                      and [b["child"] for b in br] == [0, 1, 2, 3]
                      and [b["toward"] for b in br]
                      == [str(i) for i in range(4)])
    n["struct"] += struct_ok
    if v.kind == "P":
        n_pure += 1
        u = card["unitary"]
        ok = u is not None
        if ok:
            strs = [[e["str"] for e in row] for row in u["entries"]]
            parsed = serialize.parse_matrix(strs)
            ok = (in_Gamma(parsed)
                  and Vertex("P", parsed).key() == v.key())
        n["unitary"] += ok

check("7a vertex_card builds without error on all 53 addresses of depth <= 3",
      n_err == 0 and n["build"] == len(ALL3), f"{n['build']}/{len(ALL3)}")
check("7b f3Data has exactly 4 entries (isotropic lines / self-dual planes)",
      n["f3"] == len(ALL3), f"{n['f3']}/{len(ALL3)}")
check("7c card structure matches the contract (addr/kind/depth/residue/branches)",
      n["struct"] == len(ALL3), f"{n['struct']}/{len(ALL3)}")
check("7d pure cards: unitary parses back via parse_matrix('str' fields), "
      "in_Gamma, and Vertex('P', parsed) == vertex_at(addr)",
      n["unitary"] == n_pure, f"{n['unitary']}/{n_pure} pure cards (depth 0, 2)")

rng = random.Random(SEED + 7)
n_ok, n_tot = 0, 0
for _ in range(300):
    x = Fw(Zw(rng.randint(-99, 99), rng.randint(-99, 99)),
           rng.choice([1, 3, 9, 27, 5, 15]))
    n_tot += 1
    n_ok += (serialize.parse_fw(serialize.fw_str(x)) == x)
check("7e parse_fw(fw_str(x)) == x on 300 random Fw "
      "(denominators in {1,3,9,27,5,15})",
      n_ok == n_tot, f"{n_ok}/{n_tot}")


# ===========================================================================
# Section 8: performance guard
# ===========================================================================
print()
print("=" * 76)
print("Section 8: performance guard")
print("=" * 76)

PERF_SNIPPET = r'''
import json, sys, time
sys.path.insert(0, sys.argv[2])
from vizgen import walks
w = sys.argv[1]
t0 = time.perf_counter()
traj = walks.trajectory(word=w)
dt = time.perf_counter() - t0
print(json.dumps({"dt": dt, "ok": traj["ok"], "sde": traj["sde"]}))
'''

w60 = walks.random_word(60, SEED + 8)
r = subprocess.run([PYTHON, "-c", PERF_SNIPPET, w60, ROOT],
                   env=dict(os.environ, PYTHONHASHSEED="0"),
                   capture_output=True, text=True)
perf = json.loads(r.stdout) if r.returncode == 0 else None
if perf is None:
    print(f"    perf subprocess failed:\n{r.stderr}")
check("8a trajectory(random word of length 60) under 5 s in a FRESH "
      "subprocess (cold caches); sde matches in-process btlib",
      perf is not None and perf["dt"] < 5.0 and perf["ok"]
      and perf["sde"] == sde(word_unitary(w60)),
      (f"{perf['dt']*1000:.0f} ms cold, sde={perf['sde']}, "
       f"#H={w60.count('H')}") if perf else "subprocess failure")


# ===========================================================================
print()
print("=" * 76)
if FAILURES:
    print(f"RESULT: {len(FAILURES)} FAILURE(S): {FAILURES}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
sys.exit(0)
