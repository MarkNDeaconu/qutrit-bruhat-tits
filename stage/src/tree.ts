/**
 * tree.ts — procedural skeleton of the (4,4)-biregular Bruhat–Tits tree.
 *
 * Addresses are the contract's digit strings: '' = origin; first digit 0–3,
 * later digits 0–2. The skeleton needs NO server data: the entire combinatorial
 * tree is a pure function of the addressing scheme.
 */

/** 4 children at the root, 3 everywhere else (1 parent + 3 children). */
export function childCount(addr: string): number {
  return addr === '' ? 4 : 3;
}

export function parentOf(addr: string): string | null {
  return addr === '' ? null : addr.slice(0, -1);
}

export function depthOf(addr: string): number {
  return addr.length;
}

/** Bipartite class: even depth = pure ('P'), odd depth = alternating ('A'). */
export function kindOf(addr: string): 'P' | 'A' {
  return addr.length % 2 === 0 ? 'P' : 'A';
}

/**
 * BFS enumeration of the ball of radius maxDepth: yields '' first, then each
 * depth in order. |ball(8)| = 1 + 4·(3^8 − 1)/2 = 13121.
 */
export function* generateBall(maxDepth: number): Generator<string, void, void> {
  yield '';
  let frontier: string[] = [''];
  for (let d = 1; d <= maxDepth; d++) {
    const next: string[] = [];
    for (const v of frontier) {
      const k = childCount(v);
      for (let c = 0; c < k; c++) {
        const a = v + String(c);
        yield a;
        next.push(a);
      }
    }
    frontier = next;
  }
}
