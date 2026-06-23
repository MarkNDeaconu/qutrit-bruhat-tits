"""Canonical addressing of tree vertices.

Every vertex of the Bruhat-Tits tree gets a digit-string address:
  ""        the origin e0 (pure)
  first digit in {0,1,2,3}   (the 4 neighbors of the origin)
  later digits in {0,1,2}    (the 3 non-parent children of any other vertex)

Child order is canonical: btlib neighbors minus the parent, sorted by the
canonical HNF vertex key. This makes addresses deterministic across runs.
"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from btlib import (Mat, Vertex, ORIGIN, neighbors, d_tilde, ell, sde,
                   synthesis_steps, H_GATE, monomial_matrices)

_vertex_cache: dict[str, Vertex] = {"": ORIGIN}
_children_cache: dict[str, list[Vertex]] = {}


def children_of(addr: str) -> list[Vertex]:
    """Sorted non-parent neighbors of the vertex at `addr` (4 at root, else 3)."""
    if addr in _children_cache:
        return _children_cache[addr]
    v = vertex_at(addr)
    parent_key = vertex_at(addr[:-1]).key() if addr else None
    kids = sorted((w for w in neighbors(v) if w.key() != parent_key),
                  key=lambda w: w.key())
    expected = 4 if addr == "" else 3
    assert len(kids) == expected, (addr, len(kids))
    _children_cache[addr] = kids
    return kids


def vertex_at(addr: str) -> Vertex:
    """The Vertex object at a given address (cached, built incrementally)."""
    if addr in _vertex_cache:
        return _vertex_cache[addr]
    parent = addr[:-1]
    kids = children_of(parent)
    v = kids[int(addr[-1])]
    _vertex_cache[addr] = v
    return v


def kind_of(addr: str) -> str:
    return "P" if len(addr) % 2 == 0 else "A"


def geodesic_addresses(U: Mat) -> list[str]:
    """Addresses of ALL vertices (pure and alternating, alternating order
    P,A,P,A,...,P) on the geodesic from e0 to U*e0.

    Walks down from the origin: at each vertex pick the child whose subtree
    contains the target, which is the child at distance (remaining-1) from it.
    Uses exact d_tilde for pure vertices; alternating vertices are the unique
    midpoints, identified as the common neighbor key.
    """
    target = Vertex("P", U)
    total = d_tilde(Mat.identity(), U)
    path = [""]
    if total == 0:
        return path
    # Pure vertices along the geodesic via synthesis (descent from U side):
    # synthesis_steps gives U = (m1 H)(m2 H)...(mn H) M; prefix products are
    # the geodesic pure vertices from the origin outward.
    word, _M = synthesis_steps(U)
    pures = [Vertex("P", Mat.identity())]
    acc = Mat.identity()
    for m in word:
        acc = acc * (m * H_GATE)
        pures.append(Vertex("P", acc))
    assert pures[-1] == target
    # Convert to addresses by walking down, inserting alternating midpoints.
    addr = ""
    for k in range(1, len(pures)):
        nxt = pures[k]
        # midpoint: the neighbor of current vertex adjacent to nxt
        kids_mid = children_of(addr)
        mid_idx = None
        nxt_nbr_keys = {w.key() for w in neighbors(nxt)}
        for i, a in enumerate(kids_mid):
            if a.key() in nxt_nbr_keys:
                mid_idx = i
                break
        assert mid_idx is not None, "geodesic midpoint not among children"
        addr_mid = addr + str(mid_idx)
        path.append(addr_mid)
        kids_p = children_of(addr_mid)
        p_idx = None
        for i, p in enumerate(kids_p):
            if p == nxt:
                p_idx = i
                break
        assert p_idx is not None, "geodesic pure vertex not among children"
        addr = addr_mid + str(p_idx)
        path.append(addr)
        _vertex_cache.setdefault(addr, nxt)
    return path


def address_of_pure(U: Mat) -> str:
    """Address of the pure vertex U*e0."""
    return geodesic_addresses(U)[-1]


def step_address(cur_addr: str, new_vertex: Vertex) -> str:
    """Address of a pure vertex at distance 2 from the pure vertex at cur_addr.

    The 12 candidates have address forms:
      cur[:-2]                      (grandparent)
      cur[:-1] + c, c != last digit (siblings through the parent midpoint)
      cur + c + c'                  (grandchildren)
    Matches by exact vertex key. O(1) local neighbor computations.
    """
    nk = new_vertex.key()
    if len(cur_addr) >= 2 and vertex_at(cur_addr[:-2]).key() == nk:
        return cur_addr[:-2]
    if cur_addr:
        last = cur_addr[-1]
        for i, w in enumerate(children_of(cur_addr[:-1])):
            if str(i) != last and w.key() == nk:
                return cur_addr[:-1] + str(i)
    for c, _mid in enumerate(children_of(cur_addr)):
        mid_addr = cur_addr + str(c)
        for c2, w in enumerate(children_of(mid_addr)):
            if w.key() == nk:
                a = mid_addr + str(c2)
                _vertex_cache.setdefault(a, w)
                return a
    raise ValueError("new vertex is not at distance 2 from cur_addr")
