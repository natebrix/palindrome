"""A short complete palindromic sentence, with a lower bound to hold it against.

    python3 shortest.py [--dict data/npdict.txt] [--out out/short_palindrome.txt]
                        [--share 0.25]

1. Solve the linear program of lp.py: a lower bound on every complete sentence.
2. Round it by diving. For each word pick one arc to carry it (the arc with
   the most flow) and require a unit of flow there. With the arcs chosen, what
   is left is a minimum-cost circulation with lower bounds, whose vertices are
   whole numbers. Words are settled in rounds, surest first, and the LP is
   solved again between rounds.
3. If the flow does not hang together, join each stray piece to the rest by
   its cheapest closed walk through the level state E.
4. Read the sentence off an Eulerian circuit (Hierholzer) and check it by
   reversing the string.

Needs numpy and highspy (see requirements.txt).
"""
import argparse
import heapq
import math
import os
import time

import numpy as np

import lp as L
from palgraph import Graph, letters, load


def best_arcs(cols, x, free):
    """For each free word, its arc with the most flow: [(flow, word, column)],
    most flow first."""
    word = cols["word"]
    order = np.lexsort((-x, word))          # by word, most flow first
    first = np.ones(len(order), dtype=bool)
    first[1:] = word[order][1:] != word[order][:-1]
    best = [(float(x[j]), int(word[j]), int(j)) for j in order[first]
            if word[j] >= 0 and int(word[j]) in free]
    best.sort(key=lambda t: (-t[0], t[1]))
    return best


def round_lp(model, cols, x, share, log, solver="ipm"):
    """Dive. Each round requires a unit of flow on one arc for every free word
    that already has a unit there (which costs nothing) and for the surest
    `share` of the rest, then solves again. With every
    word settled the problem is a circulation with lower bounds, and its
    vertices are whole numbers."""
    free = set(int(w) for w in cols["includable"])
    rounds = 0
    while free:
        best = best_arcs(cols, x, free)
        sure = [t for t in best if t[0] >= 1 - 1e-7]
        rest = best[len(sure):]
        take = sure + rest[: max(10, math.ceil(share * len(rest)))]
        model.require(j for _, _, j in take)
        free.difference_update(w for _, w, _ in take)
        r = model.solve(solver)
        assert r["status"] == "Optimal", r["status"]
        x = r["x"]
        rounds += 1
        log(f"  round {rounds}: {len(sure):,} words at a full unit, {len(take) - len(sure):,} "
            f"more required, {len(free):,} free; objective {r['objective']:,.2f}")
    xi = np.round(x)
    assert np.abs(x - xi).max() < 1e-6, "the rounded problem did not give whole numbers"
    return xi.astype(np.int64)


def check_flow(g, cols, x):
    """Whole-number flows must conserve flow, cover every word and close at
    most one walk."""
    n = g.n
    net = np.bincount(cols["head"], weights=x, minlength=n) - np.bincount(cols["tail"], weights=x, minlength=n)
    assert not net.any(), "flow is not conserved"
    arcs = cols["word"] >= 0
    cov = np.bincount(cols["word"][arcs], weights=x[arcs], minlength=len(g.words))
    assert (cov[cols["includable"]] >= 1).all(), "a word is not covered"
    assert x[~arcs].sum() <= 1, "more than one walk"


def components(cols, x):
    """Label the states touched by flow: 0 for the piece holding E."""
    used = np.nonzero(x > 0)[0]
    adj = {}
    for j in used:
        u, v = int(cols["tail"][j]), int(cols["head"][j])
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    comp = {}
    k = 0
    for s in [0] + list(adj):
        if s in comp:
            continue
        comp[s] = k
        stack = [s]
        while stack:
            u = stack.pop()
            for v in adj.get(u, ()):
                if v not in comp:
                    comp[v] = k
                    stack.append(v)
        k += 1
    return comp, k


def connect(g, cols, x, log):
    """Join stray pieces to E's piece with closed walks through E. Returns the
    flows with the extra walks added."""
    comp, k = components(cols, x)
    if k == 1:
        return x
    col_of_arc = {int(a): j for j, a in enumerate(cols["arc"]) if a >= 0}
    wlen = [len(w) for w in g.words]
    INF = float("inf")

    def dijkstra(forward):
        dist = [INF] * g.n
        via = [-1] * g.n
        dist[0] = 0
        pq = [(0, 0)]
        while pq:
            d, u = heapq.heappop(pq)
            if d > dist[u]:
                continue
            for a in (g.out(u) if forward else g.inn(u)):
                if a not in col_of_arc:
                    continue
                v = g.dst[a] if forward else g.src[a]
                nd = d + wlen[g.lab[a]]
                if nd < dist[v]:
                    dist[v] = nd
                    via[v] = a
                    heapq.heappush(pq, (nd, v))
        return dist, via

    df, pf = dijkstra(True)
    db, pb = dijkstra(False)
    x = x.copy()
    added = 0
    for c in range(1, k):
        u = min((s for s in comp if comp[s] == c), key=lambda s: df[s] + db[s])
        added += df[u] + db[u]
        v = u
        while v != 0:
            x[col_of_arc[pf[v]]] += 1
            v = g.src[pf[v]]
        v = u
        while v != 0:
            x[col_of_arc[pb[v]]] += 1
            v = g.dst[pb[v]]
    log(f"  joined {k - 1:,} stray pieces at a cost of {added:,} letters")
    return x


def euler(cols, x):
    """Hierholzer: the columns of an Eulerian circuit from E, in order. If a
    return arc is in use the circuit is turned so that it comes last."""
    tail, head = cols["tail"], cols["head"]
    out = {}
    for j in np.nonzero(x > 0)[0]:
        out.setdefault(int(tail[j]), []).append([int(j), int(x[j])])
    stack = [(0, -1)]
    circuit = []
    while stack:
        u, via = stack[-1]
        lst = out.get(u)
        if lst:
            e = lst[-1]
            e[1] -= 1
            if e[1] == 0:
                lst.pop()
            stack.append((int(head[e[0]]), e[0]))
        else:
            stack.pop()
            if via >= 0:
                circuit.append(via)
    circuit.reverse()
    assert len(circuit) == int(x.sum()), "the flow does not hang together"
    ret = [i for i, j in enumerate(circuit) if cols["arc"][j] < 0]
    if ret:
        i = ret[0]
        circuit = circuit[i + 1:] + circuit[:i + 1]
        assert tail[circuit[0]] == 0
        circuit.pop()
    return circuit


def sentence(g, cols, circuit):
    left, right = [], []
    for j in circuit:
        a = int(cols["arc"][j])
        (right if g.side(a) == "R" else left).append(g.words[g.lab[a]])
    return left + right[::-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dict", default="data/npdict.txt")
    ap.add_argument("--out", default="out/short_palindrome.txt")
    ap.add_argument("--resolve", default="ipm", choices=["ipm", "simplex"],
                    help="solver for the rounds after the first LP")
    ap.add_argument("--share", type=float, default=0.25,
                    help="share of the unsettled words fixed in each round")
    args = ap.parse_args()
    t0 = time.time()

    def log(msg):
        print(f"{msg}   [{time.time() - t0:.0f}s]", flush=True)

    D = load(args.dict)
    g = Graph(D.keys())
    cols = L.build(g)
    model, nstate, nword = L.to_highs(g, cols)
    log(f"LP: {model.num_col_:,} columns, {model.num_row_:,} rows")
    m = L.Model(model)
    r = m.solve("ipm")
    assert r["status"] == "Optimal", r["status"]
    bound = math.ceil(r["objective"] - 1e-6)
    log(f"LP optimum {r['objective']:,.2f}: no complete sentence has fewer than {bound:,} letters")

    x = round_lp(m, cols, r["x"], args.share, log, args.resolve)
    check_flow(g, cols, x)
    x = connect(g, cols, x, log)
    check_flow(g, cols, x)
    keys = sentence(g, cols, euler(cols, x))

    s = "".join(keys)
    assert s == s[::-1], "not a palindrome"
    assert {g.wid[k] for k in keys} == g.includable(), "a phrase is missing"
    assert len(s) == int((x * cols["cost"]).sum())
    text = ", ".join(D[k] for k in keys)
    t = letters(text)
    assert t == t[::-1] == s
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(text)
    log(f"sentence: {len(keys):,} phrases, {len(set(keys)):,} distinct, {len(s):,} letters")
    print(f"lower bound {bound:,}, sentence {len(s):,}: gap {len(s) - bound:,} letters "
          f"({100 * (len(s) - bound) / bound:.3f}%)")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
