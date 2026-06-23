/**
 * tree.ts — procedural skeleton of the (4,4)-biregular Bruhat–Tits tree.
 *
 * Addresses are the contract's digit strings: '' = origin; first digit 0–3,
 * later digits 0–2. The skeleton needs NO server data: the entire combinatorial
 * tree is a pure function of the addressing scheme.
 */
/** 4 children at the root, 3 everywhere else (1 parent + 3 children). */
export function childCount(addr) {
    return addr === '' ? 4 : 3;
}
export function parentOf(addr) {
    return addr === '' ? null : addr.slice(0, -1);
}
export function depthOf(addr) {
    return addr.length;
}
/** Bipartite class: even depth = pure ('P'), odd depth = alternating ('A'). */
export function kindOf(addr) {
    return addr.length % 2 === 0 ? 'P' : 'A';
}
/**
 * BFS enumeration of the ball of radius maxDepth: yields '' first, then each
 * depth in order. |ball(8)| = 1 + 4·(3^8 − 1)/2 = 13121.
 */
export function* generateBall(maxDepth) {
    yield '';
    let frontier = [''];
    for (let d = 1; d <= maxDepth; d++) {
        const next = [];
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
