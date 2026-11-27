"""Check the bound and the rounding against exact answers on small dictionaries.

For each random dictionary, the shortest complete sentence is found exactly by
Dijkstra over (state, set of words covered). Then:

    LP optimum  <=  LP with parity cuts  <=  exact shortest  <=  sentence built

and the sentence built must be a palindrome holding every includable word.
Every parity cut found must also hold for the whole-number flow of the
sentence built.
"""
import heapq
import math
import random
import sys

import numpy as np

import cuts as C
import lp as L
import shortest as S
from palgraph import Graph


def exact(g, inc):
    bit = {w: 1 << i for i, w in enumerate(sorted(inc))}
    full = (1 << len(bit)) - 1
    live = g.coreachable()
    dist = {(0, 0): 0}
    pq = [(0, 0, 0)]
    while pq:
        d, u, mask = heapq.heappop(pq)
        if d > dist[(u, mask)]:
            continue
        if mask == full and g.accept[u]:
            return d
        for a in g.out(u):
            v = g.dst[a]
            if not live[v]:
                continue
            key = (v, mask | bit[g.lab[a]])
            nd = d + len(g.words[g.lab[a]])
            if nd < dist.get(key, math.inf):
                dist[key] = nd
                heapq.heappush(pq, (nd, v, key[1]))
    raise AssertionError("no complete sentence")


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    done = tight = gap_lp = joined = ncuts = closed = searched = better = 0
    while done < 150:
        alpha = "ab" if rng.random() < 0.4 else "abc"
        d = set()
        for _ in range(rng.randint(3, 9)):
            d.add("".join(rng.choice(alpha) for _ in range(rng.randint(1, 5))))
        g = Graph(d)
        inc = g.includable()
        if len(inc) < 2:
            continue
        best = exact(g, inc)
        cols = L.build(g)
        model, _, _ = L.to_highs(g, cols)
        m = L.Model(model)
        r = m.solve("simplex")
        assert r["status"] == "Optimal"
        assert r["objective"] <= best + 1e-6, (d, r["objective"], best)
        x = S.whole_flow(g, cols, S.dive(m, cols, r["x"], 0.25, lambda msg: None, "simplex"), "simplex")
        S.check_flow(g, cols, x)
        before = int((x * cols["cost"]).sum())
        x = S.connect(g, cols, x, lambda msg: None)
        S.check_flow(g, cols, x)
        keys = S.sentence(g, cols, S.euler(cols, x))
        s = "".join(keys)
        assert s == s[::-1], (d, keys)
        assert {g.wid[k] for k in keys} == inc, (d, keys)
        assert len(s) == int((x * cols["cost"]).sum()) >= best, (d, keys, best)
        assert np.array_equal(S.flow_of_sentence(g, cols, keys), x), (d, keys)
        # the other way to settle the split words: branch and bound on them
        pins, fixed = S.settle(g, cols, r["x"], [], 30, lambda msg: None, show=False)
        if fixed is not None:
            x2 = S.whole_flow(g, cols, pins, "simplex", fixed=fixed)
            S.check_flow(g, cols, x2)
            x2 = S.connect(g, cols, x2, lambda msg: None)
            k2 = S.sentence(g, cols, S.euler(cols, x2))
            s2 = "".join(k2)
            assert s2 == s2[::-1] and {g.wid[k] for k in k2} == inc, (d, k2)
            assert len(s2) >= best, (d, k2, best)
            searched += 1
            better += len(s2) < len(s)
        # parity cuts: still a lower bound, and valid for the flow just built
        fresh, _, _ = L.to_highs(g, cols)
        rc, found = C.cutting_planes(g, cols, L.Model(fresh), 20, lambda msg: None, exact=5)
        assert r["objective"] - 1e-6 <= rc["objective"] <= best + 1e-6, (d, rc["objective"], best)
        index = C.Index(cols, g.n)
        for Sset, W in found:
            idx, rhs = index.row(Sset, W)
            assert len(W) % 2 == 1 and x[idx].sum() >= rhs, (d, Sset, W)
        ncuts += len(found)
        closed += math.ceil(rc["objective"] - 1e-6) == best > math.ceil(r["objective"] - 1e-6)
        done += 1
        tight += len(s) == best
        gap_lp += math.ceil(r["objective"] - 1e-6) < best
        joined += len(s) > before
    print(f"{done} dictionaries: LP <= exact <= built, every sentence checked; "
          f"built == exact in {tight}, LP bound below exact in {gap_lp}, "
          f"stray pieces joined in {joined}; {ncuts} parity cuts, all valid, "
          f"closing the bound to exact in {closed}; branch and bound on the split "
          f"words gave a valid sentence in {searched}, shorter than the dive in {better}")


if __name__ == "__main__":
    main()
