"""A lower bound on the shortest complete palindromic sentence.

    python3 lp.py [--solver simplex|ipm|pdlp] [--dict data/npdict.txt]
                  [--out out/lp_solution.npz] [--time-limit SECONDS]

The linear program, on the states of the leftover graph that can still finish:

    minimize    sum over arcs a of letters(a) * x[a]
    subject to  flow in = flow out            at every state
                sum of x over arcs labelled w >= 1   for every includable word w
                sum of x over return arcs <= 1
                x >= 0

A return arc runs from an accepting state back to the level state E at no
cost: it closes the one walk that a sentence is. A complete sentence is a
feasible integer point with cost equal to its letters, so the optimum is a
lower bound on every complete sentence. Two things are relaxed: the flows need
not be whole numbers, and the flow need not hang together as one walk.

Needs numpy and highspy (see requirements.txt).
"""
import argparse
import math
import time

import highspy
import numpy as np

from palgraph import Graph, load


def build(g):
    """Columns of the LP. Returns a dict of arrays:
    arc[j]   graph arc of column j, or -1 for a return arc
    tail[j], head[j]   states (graph node ids)
    word[j]  word id covered by column j, or -1
    cost[j]  letters
    and the list of live states and includable word ids."""
    co = g.coreachable()
    src, dst, lab = g.src, g.dst, g.lab
    live = np.frombuffer(co, dtype=np.uint8).astype(bool)
    src = np.frombuffer(src, dtype=np.int32)
    dst = np.frombuffer(dst, dtype=np.int32)
    lab = np.frombuffer(lab, dtype=np.int32)
    keep = np.nonzero(live[dst])[0]
    wlen = np.array([len(w) for w in g.words], dtype=np.int64)
    acc = np.nonzero(np.array(g.accept, dtype=bool) & live)[0]
    acc = acc[acc != 0]
    nret = len(acc)
    return {
        "arc": np.concatenate([keep, -np.ones(nret, dtype=np.int64)]),
        "tail": np.concatenate([src[keep], acc]),
        "head": np.concatenate([dst[keep], np.zeros(nret, dtype=np.int64)]),
        "word": np.concatenate([lab[keep], -np.ones(nret, dtype=np.int64)]),
        "cost": np.concatenate([wlen[lab[keep]], np.zeros(nret, dtype=np.int64)]),
        "states": np.nonzero(live)[0],
        "includable": np.unique(lab[keep]),
    }


def to_highs(g, cols, lower=None):
    """The HighsLp for the columns from build(). lower: optional per-column
    lower bounds (used to fix chosen arcs when rounding)."""
    ncol = len(cols["arc"])
    states, inc = cols["states"], cols["includable"]
    nstate, nword = len(states), len(inc)
    row_of_state = np.full(g.n, -1, dtype=np.int64)
    row_of_state[states] = np.arange(nstate)
    row_of_word = np.full(len(g.words), -1, dtype=np.int64)
    row_of_word[inc] = nstate + np.arange(nword)
    ret_row = nstate + nword
    nrow = ret_row + 1

    is_ret = cols["arc"] < 0
    # three entries per column: -1 at the tail, +1 at the head, +1 in its
    # cover row (an arc) or in the return row (a return arc)
    third = np.where(is_ret, ret_row, row_of_word[np.maximum(cols["word"], 0)])
    index = np.stack([row_of_state[cols["tail"]], row_of_state[cols["head"]], third], axis=1)
    value = np.tile(np.array([-1.0, 1.0, 1.0]), (ncol, 1))
    assert (index >= 0).all()

    inf = highspy.kHighsInf
    lp = highspy.HighsLp()
    lp.num_col_ = ncol
    lp.num_row_ = nrow
    lp.col_cost_ = cols["cost"].astype(np.float64)
    lp.col_lower_ = np.zeros(ncol) if lower is None else np.asarray(lower, dtype=np.float64)
    lp.col_upper_ = np.full(ncol, inf)
    row_lower = np.zeros(nrow)
    row_upper = np.zeros(nrow)
    row_lower[nstate:ret_row] = 1.0
    row_upper[nstate:ret_row] = inf
    row_lower[ret_row] = -inf
    row_upper[ret_row] = 1.0
    lp.row_lower_ = row_lower
    lp.row_upper_ = row_upper
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = np.arange(0, 3 * ncol + 1, 3, dtype=np.int32)
    lp.a_matrix_.index_ = index.reshape(-1).astype(np.int32)
    lp.a_matrix_.value_ = value.reshape(-1)
    return lp, nstate, nword


class Model:
    """The LP held in HiGHS, so that it can be solved again after a change of
    bounds from the basis it stopped at."""

    def __init__(self, lp, log=False, time_limit=None):
        self.h = highspy.Highs()
        self.h.setOptionValue("output_flag", bool(log))
        if time_limit:
            self.h.setOptionValue("time_limit", float(time_limit))
        self.h.passModel(lp)

    def solve(self, solver="ipm"):
        """solver: 'ipm' (interior point, then crossover to a vertex),
        'simplex' (starts from the current basis when there is one) or 'pdlp'."""
        self.h.setOptionValue("solver", solver)
        self.h.run()
        info = self.h.getInfo()
        sol = self.h.getSolution()
        return {
            "status": self.h.modelStatusToString(self.h.getModelStatus()),
            "objective": info.objective_function_value,
            "x": np.array(sol.col_value),
            "row_dual": np.array(sol.row_dual),
            "col_dual": np.array(sol.col_dual),
        }

    def add_rows(self, rows):
        """Add cuts: rows is a list of (columns, right-hand side), each cut
        the sum of its columns >= its right-hand side."""
        if not rows:
            return
        inf = highspy.kHighsInf
        starts = np.cumsum([0] + [len(c) for c, _ in rows[:-1]]).astype(np.int32)
        index = np.concatenate([c for c, _ in rows]).astype(np.int32)
        lower = np.array([float(r) for _, r in rows])
        self.h.addRows(len(rows), lower, np.full(len(rows), inf), len(index),
                       starts, index, np.ones(len(index)))

    def require(self, columns):
        """Require at least one unit of flow on each of these columns."""
        inf = highspy.kHighsInf
        for j in columns:
            self.h.changeColBounds(int(j), 1.0, inf)


def solve(lp, solver="simplex", time_limit=None, log=True):
    return Model(lp, log=log, time_limit=time_limit).solve(solver)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dict", default="data/npdict.txt")
    ap.add_argument("--solver", default="simplex", choices=["simplex", "ipm", "pdlp"])
    ap.add_argument("--out", default="out/lp_solution.npz")
    ap.add_argument("--time-limit", type=float, default=None)
    args = ap.parse_args()

    t0 = time.time()
    D = load(args.dict)
    g = Graph(D.keys())
    cols = build(g)
    lp, nstate, nword = to_highs(g, cols)
    nret = int((cols["arc"] < 0).sum())
    print(f"LP: {lp.num_col_:,} columns ({lp.num_col_ - nret:,} arcs, {nret:,} return arcs), "
          f"{lp.num_row_:,} rows ({nstate:,} states, {nword:,} cover rows, 1 return row)",
          flush=True)
    r = solve(lp, args.solver, args.time_limit)
    x = r["x"]
    print(f"status: {r['status']}   objective: {r['objective']:,.4f}   ({time.time() - t0:.0f}s)")
    if r["status"] == "Optimal":
        print(f"lower bound on any complete sentence: {math.ceil(r['objective'] - 1e-6):,} letters")
    frac = np.abs(x - np.round(x)) > 1e-6
    print(f"columns in use: {(x > 1e-9).sum():,}   fractional: {frac.sum():,}")
    import os
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    np.savez_compressed(args.out, x=x, objective=r["objective"], status=r["status"],
                        row_dual=r["row_dual"], col_dual=r["col_dual"],
                        **{k: v for k, v in cols.items()})
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
