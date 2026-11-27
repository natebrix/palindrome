"""Exact separation of the parity cuts of cuts.py, by integer programming.

A parity cut is a choice of states S and words W with |W| odd. Against an LP
solution x its slack is

    sum over w in W of (flow of w - 1)
      + sum over arcs a of x[a] * [ (a's word is in W) xor (a crosses S) ]

and it is violated, by (1 - slack) / 2, when the slack is under one. Finding
the least slack is an integer program in 0-1 variables y (state in S) and z
(word in W). Three reductions, none of which loses a violated cut, make it
small:

- A word covered by a single arc at flow exactly one must be in W exactly when
  that arc crosses S, or the slack is already one. Its z is eliminated; it
  only feeds the parity of |W|.
- An arc at flow one or more whose word is covered at least twice must not
  cross. Its two states are merged.
- A state touched only by arcs of the first kind changes neither slack nor
  parity. It is left out of S.

What remains has a few thousand variables. It is solved repeatedly, each time
with the states of the cuts already found kept out of S, until no cut is
violated by more than the threshold. When the first solve finds none, no
parity cut of this family is violated at all.
"""
import time
from collections import defaultdict

import highspy
import numpy as np

TOL = 1e-6


class Separator:
    def __init__(self, rnd):
        """rnd: a cuts.Round (one LP solution, prepared)."""
        self.rnd = rnd
        x, tail, head, word = rnd.x, rnd.tail, rnd.head, rnd.word
        arcs_of = defaultdict(list)
        for j in rnd.pos:
            arcs_of[int(word[j])].append(int(j))

        def kind(w):
            js = arcs_of[w]
            if w < 0:
                return "heavy"                      # return arcs: no covering row
            s = rnd.flow[w] - 1
            if len(js) == 1 and abs(x[js[0]] - 1) < TOL:
                return "unit"
            return "heavy" if s >= 1 - TOL else "soft"

        kinds = {w: kind(w) for w in arcs_of}

        # merge the states of arcs that must not cross
        parent = {}

        def find(u):
            r = u
            while parent.get(r, r) != r:
                r = parent[r]
            while parent.get(u, u) != u:
                parent[u], u = r, parent[u]
            return r

        for w, js in arcs_of.items():
            if kinds[w] == "heavy":
                for j in js:
                    if x[j] >= 1 - TOL:
                        a, b = find(int(tail[j])), find(int(head[j]))
                        if a != b:
                            parent[a] = b
        members = defaultdict(list)
        for u in rnd.inc:
            members[find(u)].append(u)
        giant = max(members, key=lambda r: len(members[r]))

        # variables
        self.cost, self.lower, self.upper, self.integer = [], [], [], []
        self.rows = []                               # (indices, values, lower, upper)

        def var(cost=0.0, upper=1.0, integer=True):
            self.cost.append(cost)
            self.lower.append(0.0)
            self.upper.append(upper)
            self.integer.append(integer)
            return len(self.cost) - 1

        yvar, zvar = {}, {}

        def y(p):
            if p not in yvar:
                yvar[p] = var()
            return yvar[p]

        for w, js in arcs_of.items():
            if kinds[w] == "unit":
                continue
            soft = kinds[w] == "soft"
            if soft:
                zvar[w] = var(cost=max(rnd.flow[w] - 1, 0.0))
            for j in js:
                p, q = find(int(tail[j])), find(int(head[j]))
                terms = [] if p == q else [y(r) for r in (p, q) if r != giant]
                if soft:
                    terms.append(zvar[w])
                if not terms:
                    continue
                if x[j] >= 1 - TOL:                  # must come out even
                    if len(terms) == 1:
                        self.upper[terms[0]] = 0.0
                    elif len(terms) == 2:
                        self.rows.append((terms, [1.0, -1.0], 0.0, 0.0))
                    else:
                        m = var()
                        self.rows.append((terms + [m], [1.0, 1.0, 1.0, -2.0], 0.0, 0.0))
                elif len(terms) == 1:
                    self.cost[terms[0]] += x[j]
                else:                                # e >= xor of the terms
                    e = var(cost=float(x[j]), integer=False)
                    if len(terms) == 2:
                        a, b = terms
                        self.rows.append(([e, a, b], [1.0, -1.0, 1.0], 0.0, highspy.kHighsInf))
                        self.rows.append(([e, a, b], [1.0, 1.0, -1.0], 0.0, highspy.kHighsInf))
                    else:
                        a, b, c = terms
                        for sa, sb, sc, lo in ((-1, 1, 1, 0), (1, -1, 1, 0), (1, 1, -1, 0), (-1, -1, -1, -2)):
                            self.rows.append(([e, a, b, c], [1.0, float(sa), float(sb), float(sc)],
                                              float(lo), highspy.kHighsInf))

        # parity of |W|: the soft words taken, and the unit words whose arc crosses
        odd = defaultdict(int)
        for w, js in arcs_of.items():
            if kinds[w] == "unit":
                p, q = find(int(tail[js[0]])), find(int(head[js[0]]))
                if p != q:
                    odd[p] ^= 1
                    odd[q] ^= 1
        # a state outside the model changes nothing: its unit arcs come in pairs
        self.stray = sum(1 for p, b in odd.items() if b and p != giant and p not in yvar)
        k = var(upper=highspy.kHighsInf)
        idx = [v for p, v in yvar.items() if odd.get(p)] + list(zvar.values()) + [k]
        self.rows.append((idx, [1.0] * (len(idx) - 1) + [-2.0], 1.0, 1.0))
        self.cutoff_row = len(self.rows)
        nz = [i for i, c in enumerate(self.cost) if c]
        self.rows.append((nz, [self.cost[i] for i in nz], -highspy.kHighsInf, 1.0))

        self.yvar, self.zvar, self.members = yvar, zvar, members
        self.h = None

    def _load(self, min_violation, time_limit):
        lp = highspy.HighsLp()
        n = len(self.cost)
        lp.num_col_, lp.num_row_ = n, len(self.rows)
        lp.col_cost_ = np.array(self.cost)
        lp.col_lower_ = np.array(self.lower)
        lp.col_upper_ = np.array(self.upper)
        lp.integrality_ = [highspy.HighsVarType.kInteger if i else highspy.HighsVarType.kContinuous
                           for i in self.integer]
        lo = [r[2] for r in self.rows]
        hi = [r[3] for r in self.rows]
        hi[self.cutoff_row] = 1.0 - 2 * min_violation
        lp.row_lower_, lp.row_upper_ = np.array(lo), np.array(hi)
        lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
        lp.a_matrix_.start_ = np.cumsum([0] + [len(r[0]) for r in self.rows]).astype(np.int32)
        lp.a_matrix_.index_ = np.concatenate([r[0] for r in self.rows]).astype(np.int32)
        lp.a_matrix_.value_ = np.concatenate([r[1] for r in self.rows]).astype(np.float64)
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("time_limit", float(time_limit))
        h.setOptionValue("mip_rel_gap", 0.0)
        h.setOptionValue("mip_abs_gap", 1e-7)
        h.passModel(lp)
        self.h = h

    def cuts(self, max_cuts=300, budget=120.0, min_violation=0.02, solve_limit=20.0):
        """Violated cuts, each the most violated among those whose states avoid
        the cuts before it: [(violation, states, words)]. Also returns whether
        the search ended by proof (no further cut above the threshold) rather
        than by the limits."""
        self._load(min_violation, solve_limit)
        h = self.h
        t0 = time.time()
        found = []
        proved = False
        while len(found) < max_cuts and time.time() - t0 < budget:
            h.run()
            status = h.modelStatusToString(h.getModelStatus())
            if status == "Infeasible":
                proved = True
                break
            if status != "Optimal":
                break                                # a time limit: stop, proving nothing
            v = np.array(h.getSolution().col_value)
            chosen = [p for p, i in self.yvar.items() if v[i] > 0.5]
            S = [u for p in chosen for u in self.members[p]]
            violation, W = self.rnd.evaluate(S)
            slack = h.getInfo().objective_function_value
            assert abs((1 - slack) / 2 - violation) < 1e-5, (slack, violation)
            found.append((violation, sorted(S), sorted(W)))
            for p in chosen:                         # next: a cut on other states
                h.changeColBounds(self.yvar[p], 0.0, 0.0)
        return found, proved
