"""Check the bound and the rounding against exact answers on small dictionaries.

For each random dictionary, the shortest complete sentence is found exactly by
Dijkstra over (state, set of words covered). Then:

    LP optimum  <=  exact shortest  <=  the sentence shortest.py builds

and the sentence built must be a palindrome holding every includable word.
"""
import heapq
import math
import random
import sys

import numpy as np

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
    done = tight = gap_lp = joined = 0
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
        x = S.round_lp(m, cols, r["x"], 0.25, lambda msg: None, "simplex")
        S.check_flow(g, cols, x)
        before = int((x * cols["cost"]).sum())
        x = S.connect(g, cols, x, lambda msg: None)
        S.check_flow(g, cols, x)
        keys = S.sentence(g, cols, S.euler(cols, x))
        s = "".join(keys)
        assert s == s[::-1], (d, keys)
        assert {g.wid[k] for k in keys} == inc, (d, keys)
        assert len(s) == int((x * cols["cost"]).sum()) >= best, (d, keys, best)
        done += 1
        tight += len(s) == best
        gap_lp += math.ceil(r["objective"] - 1e-6) < best
        joined += len(s) > before
    print(f"{done} dictionaries: LP <= exact <= built, every sentence checked; "
          f"built == exact in {tight}, LP bound below exact in {gap_lp}, "
          f"stray pieces joined in {joined}")


if __name__ == "__main__":
    main()
