"""Check the exact separation against brute force on small dictionaries.

For each LP solution (before and after a few rounds of cuts) with few enough
states carrying flow, try every set of those states with cuts.Round.evaluate
and compare the most violated cut with the first one exact.Separator returns.
"""
import itertools
import random
import sys

import cuts as C
import exact as X
import lp as L
from palgraph import Graph


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    solutions = violated = 0
    while solutions < 250:
        alpha = "ab" if rng.random() < 0.4 else "abc"
        d = set()
        for _ in range(rng.randint(3, 9)):
            d.add("".join(rng.choice(alpha) for _ in range(rng.randint(1, 5))))
        g = Graph(d)
        if len(g.includable()) < 2:
            continue
        cols = L.build(g)
        index = C.Index(cols, g.n)
        model, _, _ = L.to_highs(g, cols)
        m = L.Model(model)
        for _ in range(4):
            r = m.solve("simplex")
            rnd = C.Round(index, r["x"])
            states = sorted(rnd.inc)
            if len(states) > 14:
                break
            best = 0.0
            for k in range(1, len(states) + 1):
                for S in itertools.combinations(states, k):
                    best = max(best, rnd.evaluate(set(S))[0])
            found, proved = X.Separator(rnd).cuts(max_cuts=1, min_violation=0.001)
            got = found[0][0] if found else 0.0
            if best > 0.001 + 1e-9:
                assert abs(got - best) < 1e-6, (d, got, best)
            else:
                assert not found and proved, (d, got, best)
            solutions += 1
            violated += bool(found)
            if not found:
                break
            m.add_rows([index.row(S, W) for _, S, W in found])
    print(f"{solutions} LP solutions: the exact separation agrees with brute force over "
          f"every set of states; {violated} had a violated cut")


if __name__ == "__main__":
    main()
