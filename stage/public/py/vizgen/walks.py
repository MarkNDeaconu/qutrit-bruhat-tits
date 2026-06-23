"""Trajectories, straightening schedules, descent oracle, circuit tools.

The straightening model: the walk is recorded at TREE level (every vertex,
midpoints included), where path reduction = adjacent-backtrack removal is
confluent and terminates at the unique geodesic (treeness, verified). Each H
letter contributes two tree edges; S and R contribute none (they fix the
vertex). The engine emits the full reduction schedule so the stage only
animates — it never reasons about the math.
"""

from __future__ import annotations
import random as _random
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from btlib import (Mat, Vertex, H_GATE, S_GATE, R_GATE, ell, sde, d_tilde,
                   in_Gamma, s0_coset_reps, monomial_matrices, neighbors,
                   gram, hnf_local, residue_matrix, isotropic_lines,
                   selfdual_planes, CHI, synthesis_steps)
from . import treegen
from .serialize import mat_display, monomial_str, parse_matrix

GATES = {"H": H_GATE, "S": S_GATE, "R": R_GATE}

_s0_reps = None
def _reps():
    global _s0_reps
    if _s0_reps is None:
        _s0_reps = s0_coset_reps()
    return _s0_reps


def build_unitary(word: str) -> Mat:
    U = Mat.identity()
    for ch in word.upper():
        if ch not in GATES:
            raise ValueError(f"invalid gate letter {ch!r} (use H, S, R)")
        U = U * GATES[ch]
    return U


def random_word(length: int, seed: int) -> str:
    rng = _random.Random(seed)
    return "".join(rng.choice("HSR") for _ in range(length))


def trajectory(word: str | None = None, matrix=None,
               max_len: int = 200, max_sde: int = 40) -> dict:
    """The /api/walk payload. Exactly one of word/matrix."""
    if word is not None:
        word = word.strip().upper()
        if len(word) > max_len:
            raise ValueError(f"word too long (max {max_len})")
        U = build_unitary(word)
    else:
        U = parse_matrix(matrix)
        if not in_Gamma(U):
            raise ValueError("matrix is not in U_3(Z[omega,1/3]) "
                             "(unitarity or ring membership fails)")
        word = ""
    if sde(U) > max_sde:
        raise ValueError(f"sde {sde(U)} exceeds cap {max_sde}")

    # --- letter-by-letter walk (tree-level trail) ----------------------------
    steps = []
    tree_trail = [""]            # tree-level vertex addresses (midpoints incl.)
    edge_owner = []              # per tree edge: ordinal of the H that made it
    cur_addr = ""
    acc = Mat.identity()
    h_ord = 0
    for ch in word:
        g = GATES[ch]
        acc = acc * g
        if ch in ("S", "R"):
            steps.append({"letter": ch, "type": "fix", "at": cur_addr})
            continue
        new_v = Vertex("P", acc)
        mid_addr, new_addr = _step_with_mid(cur_addr, new_v)
        steps.append({"letter": "H", "type": "move",
                      "from": cur_addr, "mid": mid_addr, "to": new_addr})
        tree_trail += [mid_addr, new_addr]
        edge_owner += [h_ord, h_ord]
        h_ord += 1
        cur_addr = new_addr

    if matrix is not None:
        # no letter walk: trail is just the geodesic itself (played as a walk)
        geo = treegen.geodesic_addresses(U)
        tree_trail = list(geo)
        edge_owner = [i // 2 for i in range(len(geo) - 1)]
        steps = [{"letter": "H", "type": "move",
                  "from": geo[2 * k], "mid": geo[2 * k + 1], "to": geo[2 * k + 2]}
                 for k in range((len(geo) - 1) // 2)]

    # --- geodesic + straightening schedule -----------------------------------
    geodesic = treegen.geodesic_addresses(U)
    reduction = _reduction_schedule(tree_trail)
    reduced = _apply_reduction(tree_trail, reduction)
    assert reduced == geodesic, "reduction did not reach the geodesic"

    blocks, tail = synthesis_steps(U)
    return {
        "ok": True,
        "word": word or None,
        "sde": int(sde(U)),
        "steps": steps,
        "trail": tree_trail,
        "edgeOwner": edge_owner,
        "geodesic": geodesic,
        "reduction": reduction,
        "normalForm": {
            "blocks": [{"monomial": monomial_str(m), "h": True} for m in blocks],
            "tailMonomial": monomial_str(tail),
        },
        "matrix": mat_display(U),
        "error": None,
    }


def _step_with_mid(cur_addr: str, new_vertex: Vertex) -> tuple[str, str]:
    """(midpoint address, new pure address) for a distance-2 H step."""
    new_addr = treegen.step_address(cur_addr, new_vertex)
    if new_addr == cur_addr[:-2]:
        return cur_addr[:-1], new_addr
    if len(new_addr) == len(cur_addr):
        return cur_addr[:-1], new_addr            # sibling via parent midpoint
    return new_addr[:-1], new_addr                # grandchild via child midpoint


def _reduction_schedule(trail: list[str]) -> list[dict]:
    """Adjacent-backtrack removal events, innermost-first, replayable on the
    evolving trail: event {"index": i} removes trail[i], trail[i+1]."""
    work = list(trail)
    events = []
    i = 1
    while i < len(work) - 0:
        if i + 1 < len(work) and work[i + 1] == work[i - 1]:
            events.append({"index": i, "removed": [work[i], work[i + 1]]})
            del work[i:i + 2]
            i = max(1, i - 1)
        else:
            i += 1
    return events


def _apply_reduction(trail: list[str], events: list[dict]) -> list[str]:
    work = list(trail)
    for e in events:
        del work[e["index"]:e["index"] + 2]
    return work


def synthesis_oracle(U: Mat) -> list[dict]:
    """Per descent step: all 12 candidates with exact l-change, plus address."""
    geodesic = treegen.geodesic_addresses(U)
    pures = geodesic[::2]              # origin ... target
    out = []
    cur = U
    k = len(pures) - 1
    while ell(cur) > 0:
        l0 = ell(cur)
        cands, chosen_done = [], False
        for i, m in enumerate(_reps()):
            cand = (m * H_GATE).inv() * cur
            dl = ell(cand) - l0
            chosen = (dl == -2 and not chosen_done)
            cands.append({"rep": i, "lChange": dl, "chosen": chosen})
            if chosen:
                nxt, chosen_done = cand, True
        assert chosen_done
        out.append({"at": pures[k], "candidates": cands})
        cur, k = nxt, k - 1
    return out


def equality(a: dict, b: dict) -> dict:
    """a, b: {"word": ...} or {"matrix": ...}."""
    def to_U(x):
        return build_unitary(x["word"]) if "word" in x and x["word"] is not None \
            else parse_matrix(x["matrix"])
    Ua, Ub = to_U(a), to_U(b)
    return {
        "equal": Ua == Ub,
        "dTilde": int(d_tilde(Ua, Ub)),
        "vertexA": treegen.address_of_pure(Ua),
        "vertexB": treegen.address_of_pure(Ub),
    }


# --- inspector ----------------------------------------------------------------

def descend_vertex(v: Vertex) -> Mat:
    """A unitary U in Gamma with U*e0 = v (pure vertices only)."""
    assert v.kind == "P"
    I = Mat.identity()
    word = []
    cur = v.g
    while d_tilde(I, cur) > 0:
        d0 = d_tilde(I, cur)
        for m in _reps():
            cand = (m * H_GATE).inv() * cur
            if d_tilde(I, cand) == d0 - 2:
                word.append(m)
                cur = cand
                break
        else:
            raise RuntimeError("vertex descent failed")
    U = I
    for m in word:
        U = U * (m * H_GATE)
    return U


def vertex_card(addr: str) -> dict:
    v = treegen.vertex_at(addr)
    g = hnf_local(v.g)
    G = gram(v.g)
    if v.kind == "P":
        res = residue_matrix(G)
        f3 = [list(u) for u in isotropic_lines(res)]
        f3_kind = "line"
    else:
        res = residue_matrix(G * CHI)
        f3 = [[list(p[0]), list(p[1])] for p in selfdual_planes(res)]
        f3_kind = "plane"
    kids = treegen.children_of(addr)
    branches = []
    if addr:
        branches.append({"child": None, "toward": addr[:-1]})
    for i, _w in enumerate(kids):
        branches.append({"child": i, "toward": addr + str(i)})
    card = {
        "addr": addr, "kind": v.kind, "depth": len(addr),
        "basis": mat_display(g),
        "gram": mat_display(G),
        "residue": res, "f3Kind": f3_kind, "f3Data": f3,
        "branches": branches,
        "unitary": None,
    }
    if v.kind == "P" and len(addr) <= 12:
        card["unitary"] = mat_display(descend_vertex(v))
    return card
