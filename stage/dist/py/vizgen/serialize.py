"""Exact btlib objects -> display-ready JSON (floats/strings only)."""

from __future__ import annotations
import math, sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from btlib import Mat, Fw, gram, sde, hnf_local, residue_matrix, CHI

_SQ3_2 = math.sqrt(3) / 2.0


def fw_floats(x: Fw):
    """(re, im) floats of (a + b*omega)/d with omega = (-1 + i sqrt3)/2."""
    a, b, d = x.num.a, x.num.b, x.den
    return ((a - b / 2.0) / d, (b * _SQ3_2) / d)


def fw_str(x: Fw) -> str:
    a, b, d = x.num.a, x.num.b, x.den
    if a == 0 and b == 0:
        return "0"
    if b == 0:
        core = f"{a}"
    elif a == 0:
        core = f"{b}w" if b not in (1, -1) else ("w" if b == 1 else "-w")
    else:
        bw = f"{'+' if b > 0 else '-'}{abs(b) if abs(b) != 1 else ''}w"
        core = f"{a}{bw}"
    if d == 1:
        return core
    paren = core if (b == 0 or a == 0) and "+" not in core[1:] else f"({core})"
    return f"{paren}/{d}"


def fw_display(x: Fw) -> dict:
    re, im = fw_floats(x)
    v = x.v_pi()
    return {"re": round(re, 9), "im": round(im, 9), "str": fw_str(x),
            "vpi": None if v == math.inf else int(v)}


def mat_display(M: Mat) -> dict:
    return {"entries": [[fw_display(M[i, j]) for j in range(M.n)]
                        for i in range(M.n)],
            "sde": int(sde(M))}


def monomial_str(M: Mat) -> str:
    """Human-readable form of a monomial matrix: phases + permutation."""
    unit_names = {(1, 0): "1", (0, 1): "w", (-1, -1): "w2",
                  (-1, 0): "-1", (0, -1): "-w", (1, 1): "-w2"}
    perm, phases = [None] * 3, [None] * 3
    for i in range(3):
        for j in range(3):
            x = M[i, j]
            if not x.is_zero():
                perm[i] = j
                phases[i] = unit_names[(x.num.a, x.num.b)]
    if perm == [0, 1, 2] and phases == ["1", "1", "1"]:
        return "I"
    p = ""
    if perm != [0, 1, 2]:
        p = "·P(" + "".join(str(j + 1) for j in perm) + ")"
    return f"diag({','.join(phases)}){p}"


# --- exact-string parser (matrix input from the stage) -----------------------

def parse_fw(s: str) -> Fw:
    """Parse entries like '0', '-1', 'w', '2w', '(1+2w)/3', '1-w', '(-1+w)/9'.

    Grammar: optional parenthesized a+bw core, optional /d with d a positive
    integer (any value; engine checks 3-power separately via membership test).
    """
    import re as _re
    s = s.strip().replace(" ", "").replace("ω", "w").replace("omega", "w")
    if s in ("0", "+0", "-0"):
        return Fw(0)
    m = _re.fullmatch(r"\(?([^()/]+)\)?(?:/(\d+))?", s)
    if not m:
        raise ValueError(f"cannot parse entry {s!r}")
    core, den = m.group(1), int(m.group(2) or 1)
    a = b = 0
    for term in _re.findall(r"[+-]?[^+-]+", core):
        t = term
        sign = -1 if t.startswith("-") else 1
        t = t.lstrip("+-")
        if t.endswith("w2"):
            c = int(t[:-2]) if t[:-2] else 1
            a -= sign * c          # w^2 = -1 - w
            b -= sign * c
        elif t.endswith("w"):
            c = int(t[:-1]) if t[:-1] else 1
            b += sign * c
        else:
            a += sign * int(t)
    from btlib import Zw
    return Fw(Zw(a, b), den)


def parse_matrix(rows: list[list[str]]) -> Mat:
    assert len(rows) == 3 and all(len(r) == 3 for r in rows), "need 3x3"
    return Mat([[parse_fw(x) for x in row] for row in rows])
