"""Write the tightened model as AMPL files, for solvers reached through AMPL.

    python3 ampl_export.py [--cuts out/cuts.pkl] [--start SENTENCE.txt] [--dir out/ampl]

Writes palindrome.mod, palindrome.dat and palindrome.run. The model is stated
once and the data is the arc list and the parity cuts as their sets of states
and words, so the files stay small: about 9 MB for the full dictionary, where
the same model as an MPS file is 83 MB. (NEOS takes no file over 16 MB.)

States are renumbered 0, 1, ... with 0 the level state, and words 1, 2, ...,
both in the order of lp.build(). The run file prints the flow as lines
"x tail word head value" and "r state value", which `solfile.py read` turns
back into a sentence. With --start, that sentence's flow is given as the
starting point.
"""
import argparse
import os
import pickle

import numpy as np

import lp as L
import shortest as S
from palgraph import Graph, letters, load

MODEL = """\
# The shortest complete palindromic sentence: a walk in the leftover graph that
# uses every word at least once. See README.md.
param nstate integer > 0;
param nword integer > 0;
set STATES := 0 .. nstate - 1;          # 0 is the level state
set WORDS := 1 .. nword;
set ARCS within {STATES, WORDS, STATES};  # (tail, word, head)
set RET within STATES;                    # accepting states: a free arc back to 0
param len {WORDS} > 0;

param ncut integer >= 0 default 0;        # parity cuts
set CS {1 .. ncut} within STATES;
set CW {1 .. ncut} within WORDS;
param rhs {1 .. ncut};

var x {ARCS} integer >= 0;
var r {RET} integer >= 0;

minimize Letters: sum {(t,w,h) in ARCS} len[w] * x[t,w,h];

subject to Balance {s in STATES}:
    sum {(t,w,s) in ARCS} x[t,w,s] + (if s = 0 then sum {t in RET} r[t])
  = sum {(s,w,h) in ARCS} x[s,w,h] + (if s in RET then r[s]);

subject to Cover {w in WORDS}: sum {(t,w,h) in ARCS} x[t,w,h] >= 1;

subject to OneWalk: sum {t in RET} r[t] <= 1;

subject to Parity {c in 1 .. ncut}:
    sum {w in CW[c], (t,w,h) in ARCS: not (t in CS[c] and h not in CS[c])} x[t,w,h]
  + sum {h in CS[c], (t,w,h) in ARCS: t not in CS[c] and w not in CW[c]} x[t,w,h]
  + (if 0 in CS[c] then sum {t in RET: t not in CS[c]} r[t])
  >= rhs[c];
"""

RUN = """\
option gurobi_options '{options}';
# make the starting point complete: every flow not given in the data is zero
let {{(t,w,h) in ARCS}} x[t,w,h] := x[t,w,h];
let {{t in RET}} r[t] := r[t];
solve;
display Letters;
printf {{(t,w,h) in ARCS: x[t,w,h] > 0.5}} "x %d %d %d %d\\n", t, w, h, x[t,w,h];
printf {{t in RET: r[t] > 0.5}} "r %d %d\\n", t, r[t];
display Letters.bestbound;
"""

OPTIONS = ("outlev=1 method=2 threads=4 nodefilestart=0.5 mipfocus=3 mipgap=0 "
           "mipgapabs=0.99 timelim=27000 bestbound=1")


class Numbering:
    """AMPL's numbers for the states and words of lp.build()."""

    def __init__(self, g, cols):
        self.states = cols["states"]
        assert self.states[0] == 0
        self.state = np.full(g.n, -1, dtype=np.int64)
        self.state[self.states] = np.arange(len(self.states))
        self.words = cols["includable"]
        self.word = np.full(len(g.words), -1, dtype=np.int64)
        self.word[self.words] = 1 + np.arange(len(self.words))
        arc = cols["arc"] >= 0
        self.t = self.state[cols["tail"]]
        self.h = self.state[cols["head"]]
        self.w = np.where(arc, self.word[np.maximum(cols["word"], 0)], 0)
        self.is_arc = arc

    def column(self):
        """{(tail, word): column} for arcs and {state: column} for return arcs."""
        arcs, rets = {}, {}
        for j in range(len(self.t)):
            if self.is_arc[j]:
                arcs[(int(self.t[j]), int(self.w[j]))] = j
            else:
                rets[int(self.t[j])] = j
        return arcs, rets


def wrap(items, per=12):
    items = [str(i) for i in items]
    return "\n".join(" ".join(items[k:k + per]) for k in range(0, len(items), per))


def write(g, cols, cuts, dirname, x0=None, options=OPTIONS):
    num = Numbering(g, cols)
    os.makedirs(dirname, exist_ok=True)
    with open(os.path.join(dirname, "palindrome.mod"), "w") as f:
        f.write(MODEL)
    with open(os.path.join(dirname, "palindrome.run"), "w") as f:
        f.write(RUN.format(options=options))
    arcs = np.nonzero(num.is_arc)[0]
    with open(os.path.join(dirname, "palindrome.dat"), "w") as f:
        f.write(f"data;\nparam nstate := {len(num.states)};\nparam nword := {len(num.words)};\n")
        f.write("param len :=\n")
        f.write(wrap((f"{i + 1} {len(g.words[w])}" for i, w in enumerate(num.words)), 8))
        f.write(";\nset RET :=\n" + wrap(num.t[~num.is_arc]) + ";\n")
        f.write("set ARCS :=\n")
        f.write("\n".join(f"{num.t[j]} {num.w[j]} {num.h[j]}" for j in arcs))
        f.write(";\n")
        f.write(f"param ncut := {len(cuts)};\n")
        for c, (Sset, W) in enumerate(cuts, 1):
            f.write(f"set CS[{c}] :=\n" + wrap(num.state[list(Sset)]) + ";\n")
            f.write(f"set CW[{c}] :=\n" + wrap(num.word[list(W)]) + ";\n")
        if cuts:
            f.write("param rhs :=\n")
            f.write(wrap((f"{c} {(len(W) + 1) // 2}" for c, (_, W) in enumerate(cuts, 1)), 8) + ";\n")
        if x0 is not None:
            used = arcs[x0[arcs] > 0]
            f.write("var x :=\n")
            f.write("\n".join(f"{num.t[j]} {num.w[j]} {num.h[j]} {x0[j]}" for j in used) + ";\n")
            ret = np.nonzero(~num.is_arc & (x0 > 0))[0]
            if len(ret):
                f.write("var r :=\n" + wrap(f"{num.t[j]} {x0[j]}" for j in ret) + ";\n")
    return num


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dict", default="data/npdict.txt")
    ap.add_argument("--cuts", default="out/cuts.pkl")
    ap.add_argument("--start", default=None)
    ap.add_argument("--dir", default="out/ampl")
    ap.add_argument("--options", default=OPTIONS, help="the gurobi_options string")
    args = ap.parse_args()
    D = load(args.dict)
    g = Graph(D.keys())
    cols = L.build(g)
    cuts = []
    if os.path.exists(args.cuts):
        with open(args.cuts, "rb") as f:
            cuts = pickle.load(f)
    x0 = None
    if args.start:
        keys = [letters(p) for p in open(args.start, encoding="utf-8").read().split(", ")]
        x0 = S.flow_of_sentence(g, cols, keys)
        S.check_flow(g, cols, x0)
    write(g, cols, cuts, args.dir, x0, args.options)
    for name in ("palindrome.mod", "palindrome.dat", "palindrome.run"):
        p = os.path.join(args.dir, name)
        print(f"wrote {p}: {os.path.getsize(p):,} bytes")


if __name__ == "__main__":
    main()
