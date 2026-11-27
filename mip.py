"""Branch and bound on the cut-tightened model, to close the last of the gap.

    python3 mip.py --cuts out/cuts.pkl [--start out/short_palindrome.txt]
                   [--time-limit 1800] [--out out/mip_palindrome.txt]

The model of lp.py with whole-number flows and the parity cuts of cuts.py as
extra rows, solved by HiGHS with interior point at the root. A sentence given
with --start is handed over as the first incumbent. The model does not ask the
flow to hang together, so its optimum is a lower bound; the flow it returns is
joined up if need be and read off as a sentence, which is an upper bound.

The lower bound printed is the solver's: it rests on HiGHS's own cutting
planes as well as the parity cuts, and is not checked independently here.
"""
import argparse
import math
import os
import pickle
import time

import highspy
import numpy as np

import cuts as C
import lp as L
import shortest as S
from palgraph import Graph, letters, load


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dict", default="data/npdict.txt")
    ap.add_argument("--cuts", default="out/cuts.pkl")
    ap.add_argument("--start", default=None, help="a sentence to start from")
    ap.add_argument("--out", default="out/mip_palindrome.txt")
    ap.add_argument("--time-limit", type=float, default=1800)
    ap.add_argument("--threads", type=int, default=0)
    args = ap.parse_args()
    t0 = time.time()

    D = load(args.dict)
    g = Graph(D.keys())
    cols = L.build(g)
    model, _, _ = L.to_highs(g, cols)
    model.integrality_ = [highspy.HighsVarType.kInteger] * model.num_col_
    m = L.Model(model, log=True, time_limit=args.time_limit)
    h = m.h
    if os.path.exists(args.cuts):
        with open(args.cuts, "rb") as f:
            known = pickle.load(f)
        index = C.Index(cols, g.n)
        m.add_rows([index.row(Sset, W) for Sset, W in known])
        print(f"added {len(known):,} parity cuts", flush=True)
    h.setOptionValue("mip_lp_solver", "ipm")
    h.setOptionValue("mip_rel_gap", 0.0)
    h.setOptionValue("mip_abs_gap", 0.999)
    if args.threads:
        h.setOptionValue("threads", args.threads)
    if args.start:
        text = open(args.start, encoding="utf-8").read()
        keys = [letters(p) for p in text.split(", ")]
        x0 = S.flow_of_sentence(g, cols, keys)
        S.check_flow(g, cols, x0)
        sol = highspy.HighsSolution()
        sol.col_value = x0.astype(np.float64)
        h.setSolution(sol)
        print(f"starting from {args.start}: {int((x0 * cols['cost']).sum()):,} letters", flush=True)
    h.run()
    info = h.getInfo()
    status = h.modelStatusToString(h.getModelStatus())
    bound = math.ceil(info.mip_dual_bound - 1e-6)
    print(f"status: {status}   lower bound: {bound:,}   best flow: "
          f"{info.objective_function_value:,.0f}   ({time.time() - t0:.0f}s)")
    x = np.round(np.array(h.getSolution().col_value)).astype(np.int64)
    S.check_flow(g, cols, x)
    # trust the flow returned, not the value reported for it
    assert int((x * cols["cost"]).sum()) == round(info.objective_function_value), \
        "the solver's reported value does not match the flow it returned"
    x = S.connect(g, cols, x, print)
    keys = S.sentence(g, cols, S.euler(cols, x))
    s = "".join(keys)
    assert s == s[::-1] and {g.wid[k] for k in keys} == g.includable()
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(", ".join(D[k] for k in keys))
    print(f"sentence: {len(keys):,} phrases, {len(s):,} letters; lower bound {bound:,}; "
          f"gap {len(s) - bound:,}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
