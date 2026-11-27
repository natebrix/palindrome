"""Parity cuts for the linear program of lp.py.

Take a set S of states and an odd set W of words. Add up the covering rows of
the words in W and the flow balance of the states in S:

    x(arcs labelled by W) + x(arcs into S) - x(arcs out of S) >= |W|

Halve it and round up, coefficients and right-hand side (flows are whole and
not negative):

    x(arcs labelled by W that do not leave S)
      + x(arcs into S not labelled by W)        >= (|W| + 1) / 2

The cut bites when each word of W is covered exactly once, by flow that
crosses the boundary of S, and nothing else crosses it: then the flow across
the boundary is |W|, odd, yet as much must enter S as leaves it.

Separation is a heuristic in two parts. A cheap exact evaluation gives, for any
S, the best W and the violation. Candidate sets S come from the structure of
the LP solution: contract every arc that carries flow for a word covered more
than once (it must not cross), then ask that all the arcs of a split word cross
together. That is a linear system over GF(2); its null space, and the single
contracted states, are the candidates.
"""
import math
from collections import defaultdict

import numpy as np

EPS = 1e-6


class Index:
    """Lookups over the columns of the LP that do not change between rounds."""

    def __init__(self, cols, nnode):
        self.cols = cols
        self.nnode = nnode
        tail, head, word = cols["tail"], cols["head"], cols["word"]
        ncol = len(tail)
        self.by_word = _group(word, ncol)      # word -> columns (word -1: return arcs)
        self.by_head = _group(head, ncol)
        self.by_tail = _group(tail, ncol)

    def row(self, S, W):
        """Columns of the cut for states S and words W (all with coefficient 1)
        and its right-hand side. Of the two equivalent orientations, the one
        with fewer columns."""
        tail, head, word = self.cols["tail"], self.cols["head"], self.cols["word"]
        S = set(S)
        W = set(W)
        inS = np.zeros(self.nnode, dtype=bool)
        inS[list(S)] = True
        label = np.concatenate([self.by_word.get(w, _EMPTY) for w in W]) if W else _EMPTY
        best = None
        for into in (True, False):
            # words of W on arcs that do not leave S (or, mirrored, do not enter)
            away = (inS[tail[label]] & ~inS[head[label]]) if into else (~inS[tail[label]] & inS[head[label]])
            keep = label[~away]
            side = self.by_head if into else self.by_tail
            other = self.by_tail if into else self.by_head
            cross = np.concatenate([side.get(u, _EMPTY) for u in S])
            far = tail[cross] if into else head[cross]
            cross = cross[~inS[far]]
            cross = cross[~np.isin(word[cross], list(W))]
            idx = np.unique(np.concatenate([keep, cross]))
            if best is None or len(idx) < len(best):
                best = idx
        return best, (len(W) + 1) // 2


_EMPTY = np.zeros(0, dtype=np.int64)


def _group(key, n):
    order = np.argsort(key, kind="stable")
    k = key[order]
    cut = np.nonzero(np.diff(k))[0] + 1
    return {int(key[s[0]]): s for s in np.split(order, cut) if len(s)}


class Round:
    """One LP solution, prepared for separation."""

    def __init__(self, index, x):
        self.index = index
        cols = index.cols
        self.x = x
        self.tail, self.head, self.word = cols["tail"], cols["head"], cols["word"]
        self.pos = np.nonzero(x > EPS)[0]
        nword = int(self.word.max()) + 1
        arcs = self.pos[self.word[self.pos] >= 0]
        self.flow = np.bincount(self.word[arcs], weights=x[arcs], minlength=nword)
        inc = defaultdict(list)                 # state -> columns with flow at it
        for j in self.pos:
            inc[int(self.tail[j])].append(int(j))
            inc[int(self.head[j])].append(int(j))
        self.inc = inc

    def evaluate(self, S):
        """Best odd word set for the states S: (violation, words). The cut is
        violated when the violation is positive; at most 0.5."""
        S = S if isinstance(S, (set, frozenset)) else set(S)
        x, tail, head, word = self.x, self.tail, self.head, self.word
        cross = defaultdict(float)
        seen = set()
        for u in S:
            for j in self.inc.get(u, ()):
                if j in seen:
                    continue
                seen.add(j)
                if (int(tail[j]) in S) != (int(head[j]) in S):
                    cross[int(word[j])] += x[j]
        cost = cross.pop(-1, 0.0)               # return arcs belong to no word
        W = []
        flip = (math.inf, None, None)
        for w, c in cross.items():
            c0 = c                              # leave w out: its crossing flow is slack
            c1 = (self.flow[w] - 1) + (self.flow[w] - c)   # take w: surplus and flow not crossing
            if c1 < c0:
                W.append(w)
                cost += c1
            else:
                cost += c0
            d = abs(c1 - c0)
            if d < flip[0]:
                flip = (d, w, c1 < c0)
        if len(W) % 2 == 0:
            if flip[1] is None:
                return -1.0, []
            cost += flip[0]
            if flip[2]:
                W.remove(flip[1])
            else:
                W.append(flip[1])
        return (1 - cost) / 2, W

    def candidates(self):
        """Sets of states worth evaluating."""
        x, tail, head, word = self.x, self.tail, self.head, self.word
        parent = {}

        def find(u):
            r = u
            while parent.get(r, r) != r:
                r = parent[r]
            while parent.get(u, u) != u:
                parent[u], u = r, parent[u]
            return r

        tight = self.flow < 1 + EPS
        arcs_of = defaultdict(list)
        for j in self.pos:
            w = int(word[j])
            if w < 0 or not tight[w]:
                a, b = find(int(tail[j])), find(int(head[j]))
                if a != b:
                    parent[a] = b
            else:
                arcs_of[w].append(int(j))
        members = defaultdict(list)
        for u in self.inc:
            members[find(u)].append(u)
        giant = max(members, key=lambda r: len(members[r]))
        out = [members[r] for r in members if r != giant]

        # split words: all their arcs cross together. Rows over GF(2), one
        # variable per contracted state, the giant fixed outside S.
        rows = []
        for w, js in arcs_of.items():
            if len(js) < 2:
                continue
            first = {find(int(tail[js[0]])), find(int(head[js[0]]))} if find(int(tail[js[0]])) != find(int(head[js[0]])) else set()
            for j in js[1:]:
                a, b = find(int(tail[j])), find(int(head[j]))
                r = set(first) ^ ({a, b} if a != b else set())
                r.discard(giant)
                if r:
                    rows.append(r)
        # components of the system
        link = {}

        def root(u):
            while link.get(u, u) != u:
                link[u] = link.get(link[u], link[u])
                u = link[u]
            return u

        for r in rows:
            it = iter(r)
            a = root(next(it))
            for b in it:
                b = root(b)
                if a != b:
                    link[b] = a
        groups = defaultdict(list)
        for r in rows:
            groups[root(next(iter(r)))].append(r)
        for rs in groups.values():
            var = sorted(set().union(*rs))
            pos = {v: i for i, v in enumerate(var)}
            for vec in _null_space([sum(1 << pos[v] for v in r) for r in rs], len(var)):
                S = [u for i, v in enumerate(var) if vec >> i & 1 for u in members[v]]
                if S:
                    out.append(S)
        return out

    def separate(self, limit=None, min_violation=0.05):
        """Violated cuts, most violated first: [(violation, S, W)]."""
        found = {}
        for S in self.candidates():
            v, W = self.evaluate(S)
            if v > min_violation:
                key = (frozenset(S), frozenset(W))
                found[key] = v
        cuts = sorted(((v, sorted(S), sorted(W)) for (S, W), v in found.items()),
                      key=lambda t: (-t[0], len(t[1])))
        return cuts[:limit] if limit else cuts


def _null_space(rows, n):
    """Basis of the null space over GF(2) of rows given as bit masks."""
    pivots = {}                                  # pivot bit -> reduced row
    for r in rows:
        for p, q in pivots.items():
            if r >> p & 1:
                r ^= q
        if r:
            p = r.bit_length() - 1
            for k in list(pivots):
                if pivots[k] >> p & 1:
                    pivots[k] ^= r
            pivots[p] = r
    basis = []
    for f in range(n):
        if f in pivots:
            continue
        v = 1 << f
        for p, q in pivots.items():
            if q >> f & 1:
                v |= 1 << p
        basis.append(v)
    return basis


def cutting_planes(g, cols, model, rounds, log, known=(), keep=None, exact=0.0):
    """Solve, separate, add, repeat. Returns the last LP result and every cut
    added, as (states, words). `known` cuts are added before the first solve.
    keep(result, cuts) is called after every solve, with the cuts that LP had.
    exact: seconds a round to spend on the exact separation of exact.py, after
    the heuristic one (0: heuristic only)."""
    index = Index(cols, g.n)
    cuts = [(list(S), list(W)) for S, W in known]
    model.add_rows([index.row(S, W) for S, W in cuts])
    have = {(frozenset(S), frozenset(W)) for S, W in cuts}
    r = None
    for k in range(rounds + 1):
        r = model.solve("ipm")
        assert r["status"] == "Optimal", r["status"]
        x = r["x"]
        frac = int((np.abs(x - np.round(x)) > EPS).sum())
        if keep:
            keep(r, cuts)
        note = ""
        found = []
        if k < rounds:
            rnd = Round(index, x)
            found = rnd.separate()
            if exact:
                from exact import Separator
                more, proved = Separator(rnd).cuts(budget=exact)
                found += more
                note = f" ({len(more):,} by exact separation" + (", which found no more)" if proved else ")")
        new, seen = [], set()
        for _, S, W in found:
            key = (frozenset(S), frozenset(W))
            if key not in have and key not in seen:
                seen.add(key)
                new.append((S, W))
        log(f"  round {k}: objective {r['objective']:,.3f}, {len(cuts):,} cuts, "
            f"{frac:,} fractional columns, {len(new):,} new cuts{note}")
        if not new:
            break
        model.add_rows([index.row(S, W) for S, W in new])
        cuts += new
        have.update((frozenset(S), frozenset(W)) for S, W in new)
    return r, cuts


def main():
    import argparse
    import os
    import pickle
    import time

    import lp as L
    from palgraph import Graph, load

    ap = argparse.ArgumentParser(description="Tighten the LP bound with parity cuts.")
    ap.add_argument("--dict", default="data/npdict.txt")
    ap.add_argument("--rounds", type=int, default=30)
    ap.add_argument("--cuts", default="out/cuts.pkl", help="where to keep the cuts found")
    ap.add_argument("--resume", action="store_true", help="start from the cuts already kept")
    ap.add_argument("--exact", type=float, default=0.0, metavar="SECONDS",
                    help="also separate exactly by integer programming, this long a round")
    args = ap.parse_args()
    t0 = time.time()

    def log(msg):
        print(f"{msg}   [{time.time() - t0:.0f}s]", flush=True)

    D = load(args.dict)
    g = Graph(D.keys())
    cols = L.build(g)
    lp, _, _ = L.to_highs(g, cols)
    known = pickle.load(open(args.cuts, "rb")) if args.resume else []

    def keep(r, found):
        with open(args.cuts + ".tmp", "wb") as f:
            pickle.dump(found, f)
        os.replace(args.cuts + ".tmp", args.cuts)
        np.savez_compressed(args.cuts + ".x.npz", x=r["x"], objective=r["objective"])

    os.makedirs(os.path.dirname(args.cuts) or ".", exist_ok=True)
    r, found = cutting_planes(g, cols, L.Model(lp), args.rounds, log, known, keep, args.exact)
    print(f"lower bound with {len(found):,} cuts: {math.ceil(r['objective'] - 1e-6):,} letters")


if __name__ == "__main__":
    main()
