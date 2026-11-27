"""Solve an exported model with SCIP.

    python3 scip_run.py out/palindrome.mps [--start SENTENCE.txt] [--time-limit 1800]
                        [--separate] [--out out/scip_palindrome.txt]

Reads the MPS file written by export.py.

With --separate the parity cuts are also separated inside the search, by the
same code as cuts.py (heuristic) and, at the root, exact.py. SCIP's own
parameters are set with --set, e.g. --set separating/maxroundsroot=5.

Needs pyscipopt as well as numpy and highspy.
"""
import argparse
import math
import time

import numpy as np
from pyscipopt import SCIP_RESULT, Model, Sepa

import cuts as C
import lp as L
import shortest as S
from palgraph import Graph, letters, load


def flow(model, variables, sol=None):
    """The value of every column, in an LP solution (sol None) or a stored one."""
    return np.fromiter((model.getSolVal(sol, v) for v in variables), dtype=np.float64,
                       count=len(variables))


class Parity(Sepa):
    """The parity cuts of cuts.py as a SCIP separator."""

    def __init__(self, index, variables, exact_seconds, max_depth):
        self.index = index
        self.vars = variables
        self.exact_seconds = exact_seconds
        self.max_depth = max_depth
        self.have = set()
        self.calls = self.added = 0
        self.seconds = 0.0

    def sepaexeclp(self):
        m = self.model
        depth = m.getDepth()
        if depth > self.max_depth:
            return {"result": SCIP_RESULT.DIDNOTRUN}
        t0 = time.time()
        self.calls += 1
        rnd = C.Round(self.index, flow(m, self.vars))
        found = rnd.separate()
        if self.exact_seconds and depth == 0:
            from exact import Separator
            more, _ = Separator(rnd).cuts(budget=self.exact_seconds)
            found += more
        result = SCIP_RESULT.DIDNOTFIND
        for _, Sset, W in found:
            key = (frozenset(Sset), frozenset(W))
            if key in self.have:
                continue
            self.have.add(key)
            idx, rhs = self.index.row(Sset, W)
            row = m.createEmptyRowSepa(self, f"parity{len(self.have)}", lhs=float(rhs), rhs=None,
                                       local=False)
            m.cacheRowExtensions(row)
            for j in idx:
                m.addVarToRow(row, self.vars[int(j)], 1.0)
            m.flushRowExtensions(row)
            if m.isCutEfficacious(row):
                infeasible = m.addCut(row, forcecut=True)
                self.added += 1
                result = SCIP_RESULT.CUTOFF if infeasible else SCIP_RESULT.SEPARATED
            m.releaseRow(row)
        self.seconds += time.time() - t0
        return {"result": result}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mps")
    ap.add_argument("--dict", default="data/npdict.txt")
    ap.add_argument("--start", default=None)
    ap.add_argument("--time-limit", type=float, default=1800)
    ap.add_argument("--separate", action="store_true")
    ap.add_argument("--exact", type=float, default=30.0,
                    help="seconds of exact separation per round at the root (with --separate)")
    ap.add_argument("--depth", type=int, default=0,
                    help="deepest node at which parity cuts are separated (with --separate)")
    ap.add_argument("--set", action="append", default=[], metavar="PARAM=VALUE",
                    help="a SCIP parameter, e.g. separating/maxroundsroot=5")
    ap.add_argument("--out", default="out/scip_palindrome.txt")
    args = ap.parse_args()
    t0 = time.time()

    D = load(args.dict)
    g = Graph(D.keys())
    cols = L.build(g)
    ncol = len(cols["arc"])
    m = Model()
    m.readProblem(args.mps)
    by_name = {v.name: v for v in m.getVars()}
    variables = [by_name[f"c{j}"] for j in range(ncol)]
    print(f"read {args.mps}: {m.getNVars():,} variables, {m.getNConss():,} rows "
          f"[{time.time() - t0:.0f}s]", flush=True)
    m.setParam("limits/time", args.time_limit)
    for kv in args.set:
        k, v = kv.split("=", 1)
        cur = m.getParam(k)
        m.setParam(k, v.lower() in ("1", "true") if isinstance(cur, bool) else type(cur)(v))
    sepa = None
    if args.separate:
        sepa = Parity(C.Index(cols, g.n), variables, args.exact, args.depth)
        m.includeSepa(sepa, "parity", "parity cuts over covering rows and flow balance",
                      priority=1000000, freq=1)
    if args.start:
        keys = [letters(p) for p in open(args.start, encoding="utf-8").read().split(", ")]
        x0 = S.flow_of_sentence(g, cols, keys)
        sol = m.createSol()
        for j in np.nonzero(x0)[0]:
            m.setSolVal(sol, variables[int(j)], float(x0[j]))
        ok = m.addSol(sol)
        print(f"start {args.start}: {int((x0 * cols['cost']).sum()):,} letters, accepted: {ok}",
              flush=True)
    m.optimize()
    bound = m.getDualbound()
    print(f"status: {m.getStatus()}   lower bound: {math.ceil(bound - 1e-6):,} ({bound:,.3f})   "
          f"best: {m.getPrimalbound():,.0f}   nodes: {m.getNNodes():,}   "
          f"[{time.time() - t0:.0f}s]")
    if sepa:
        print(f"parity separator: {sepa.calls} calls, {sepa.added:,} cuts, {sepa.seconds:.0f}s")
    if m.getNSols():
        x = np.round(flow(m, variables, m.getBestSol())).astype(np.int64)
        try:
            S.check_flow(g, cols, x)
        except AssertionError as e:            # a pinned export has continuous columns
            print(f"best solution is not a whole-number flow ({e}); no sentence written")
            return
        x = S.connect(g, cols, x, print)
        keys = S.sentence(g, cols, S.euler(cols, x))
        s = "".join(keys)
        assert s == s[::-1] and {g.wid[k] for k in keys} == g.includable()
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(", ".join(D[k] for k in keys))
        print(f"sentence: {len(keys):,} phrases, {len(s):,} letters; wrote {args.out}")


if __name__ == "__main__":
    main()
