"""Quick smoke test of btlib.py before the full verification suite."""
import sys
sys.path.insert(0, "/Users/markdeaconu/projects/qutrits_v2")
from btlib import *

# --- ring basics ---
assert (CHI * CHI.conj()) == Fw(3), "chi * conj(chi) = 3"
assert CHI.conj() == -(OMEGA * OMEGA) * CHI, "conj(chi) = -omega^2 chi"
assert Fw(3).v_pi() == 2 and CHI.v_pi() == 1
assert OMEGA.res3() == 1, "omega == 1 mod chi"
assert I_SQRT3 * I_SQRT3.conj() == Fw(1, 3), "|i/sqrt3|^2 = 1/3"
assert I_SQRT3.v_pi() == -1

# --- gates ---
assert H_GATE.is_unitary(), "H unitary"
assert S_GATE.is_unitary() and R_GATE.is_unitary()
assert ell(H_GATE) == 2, "l(H) = 2"
assert in_Gamma(H_GATE) and in_Gamma(S_GATE) and in_Gamma(R_GATE)
HH = H_GATE * H_GATE
P23 = Mat([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
assert HH == -P23, "H^2 = -P23"

# --- monomials ---
mons = monomial_matrices()
assert len(mons) == 1296
assert all(m.is_unitary() for m in mons[:20])

# --- lattices / duality ---
assert is_self_dual(Mat.identity())
assert is_self_dual(H_GATE)
g = Mat([[CHI, 1, 1], [0, -1, 1], [0, 0, 1]])  # arbitrary
dd = dual_basis(dual_basis(g))
assert lattice_eq(g, dd), "double dual"

# --- vertex machinery ---
assert ORIGIN.kind == "P"
nbrs = neighbors_of_pure(ORIGIN)
assert len(nbrs) == 4, f"origin degree {len(nbrs)}"
back = neighbors_of_alternating(nbrs[0])
assert len(back) == 4
assert sum(1 for w in back if w == ORIGIN) == 1, "origin among neighbors-of-neighbor"

# --- canonical keys: H*O^3 as a vertex, distance 2 ---
vH = Vertex("P", H_GATE)
assert vH != ORIGIN
assert d_tilde(Mat.identity(), H_GATE) == 2

# --- BFS radius 3 ---
levels, parent, cycles = bfs_tree(3)
sizes = [len(l) for l in levels]
assert sizes == [1, 4, 12, 36], sizes
assert cycles == [], f"cycle edges found: {cycles[:3]}"

# --- synthesis on a random word ---
import random
random.seed(7)
U = Mat.identity()
gates = [H_GATE, S_GATE, R_GATE]
for _ in range(25):
    U = U * random.choice(gates)
word, M = synthesis_steps(U)
assert reconstruct(word, M) == U, "synthesis round-trip"
print(f"synthesis: sde={sde(U)}, word length {len(word)} H-blocks, OK")

print("ALL SMOKE TESTS PASSED")
