"""Write the tightened model as an MPS file, for another solver.

    python3 export.py [--cuts out/cuts.pkl] [--pin] [--out out/palindrome.mps]

The model of lp.py with whole-number flows and the parity cuts kept by cuts.py
as extra rows. Columns are named c0, c1, ... in the order of lp.build().

With --pin the LP is solved first and every word it gives a full unit on one
arc is pinned there, as in `shortest.py --search`: whole numbers are then asked
only of the arcs of the other words. That model is a restriction, so its
optimum is a sentence length, not a lower bound.
"""
import argparse
import os
import pickle

import highspy
import numpy as np

import cuts as C
import lp as L
import shortest as S
from palgraph import Graph, load


def build(g, cols, rows, pin=False, log=print):
    """A HiGHS instance holding the model, and the mask of whole-number columns."""
    ncol = len(cols["arc"])
    lower = np.zeros(ncol)
    whole = np.ones(ncol, dtype=bool)
    if pin:
        relaxed, _, _ = L.to_highs(g, cols)
        m = L.Model(relaxed)
        m.add_rows(rows)
        r = m.solve("ipm")
        best = S.best_arcs(cols, r["x"], set(int(w) for w in cols["includable"]))
        pinned = [j for f, _, j in best if f >= 1 - 1e-7]
        split = np.array(sorted(w for f, w, _ in best if f < 1 - 1e-7))
        lower[pinned] = 1.0
        whole = np.isin(cols["word"], split)
        log(f"LP {r['objective']:,.2f}; {len(pinned):,} words pinned, {len(split):,} left, "
            f"{int(whole.sum()):,} whole-number arcs")
    model, _, _ = L.to_highs(g, cols, lower=lower)
    model.integrality_ = [highspy.HighsVarType.kInteger if b else highspy.HighsVarType.kContinuous
                          for b in whole]
    m = L.Model(model)
    m.add_rows(rows)
    return m.h, whole


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dict", default="data/npdict.txt")
    ap.add_argument("--cuts", default="out/cuts.pkl")
    ap.add_argument("--pin", action="store_true")
    ap.add_argument("--out", default="out/palindrome.mps")
    args = ap.parse_args()
    D = load(args.dict)
    g = Graph(D.keys())
    cols = L.build(g)
    rows = []
    if os.path.exists(args.cuts):
        with open(args.cuts, "rb") as f:
            index = C.Index(cols, g.n)
            rows = [index.row(Sset, W) for Sset, W in pickle.load(f)]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    h, _ = build(g, cols, rows, args.pin)
    h.writeModel(args.out)
    print(f"wrote {args.out}: {h.getNumCol():,} columns, {h.getNumRow():,} rows, "
          f"{len(rows):,} of them parity cuts")


if __name__ == "__main__":
    main()
