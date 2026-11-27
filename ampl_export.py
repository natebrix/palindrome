"""Write the tightened model as AMPL files, for solvers reached through AMPL.

    python3 ampl_export.py [--cuts out/cuts.pkl] [--start SENTENCE.txt] [--dir out/ampl]

Writes palindrome.mod, palindrome.dat and palindrome.run. The model is stated
once and the data is the arc list and the parity cuts as their sets of states
and words, so the files stay small: about 9 MB for the full dictionary, where
the same model as an MPS file is 83 MB. (NEOS takes no file over 16 MB.)

States are renumbered 0, 1, ... with 0 the level state, and words 1, 2, ...,
both in the order of lp.build(). With --start, that sentence's flow is given
as the starting point.

The run file solves the root node again and again, each time from the best
flow so far, and after each solve prints the flows that changed as lines
"x tail word head value" and "r state value", then a line "flow N complete".
`solfile.py read` turns that output back into a sentence. The lower bound is
the best bound in the solver's log.
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
# Rounds of the root node alone, each from the best flow so far and with its
# own random seed. The search past the root does not fit in the memory NEOS
# allows (3 GB, of which AMPL itself holds 0.5 GB); the solver's heuristics at
# the root do, and they are what shortens the flow. Each round prints what it
# changed, so a job that dies has already printed the rounds before.
param x1 {{ARCS}} default 0;
param r1 {{RET}} default 0;
# make the starting point complete: every flow not given in the data is zero
let {{(t,w,h) in ARCS}} x[t,w,h] := x[t,w,h];
let {{t in RET}} r[t] := r[t];
for {{k in 1 .. {rounds}}} {{
    if _ampl_elapsed_time > {total} - {round_seconds} - 150 then break;
    option gurobi_options ('{options} timelim={round_seconds} seed=' & k);
    solve;
    display Letters;
    printf {{(t,w,h) in ARCS: round(x[t,w,h]) <> x1[t,w,h]}} "x %d %d %d %d\\n", t, w, h, round(x[t,w,h]);
    printf {{t in RET: round(r[t]) <> r1[t]}} "r %d %d\\n", t, round(r[t]);
    printf "flow %d complete\\n", k;
    let {{(t,w,h) in ARCS}} x1[t,w,h] := round(x[t,w,h]);
    let {{t in RET}} r1[t] := round(r[t]);
}}
"""

# Gurobi through AMPL: the settings of a NEOS run that got through the root
# inside the memory allowed and found shorter flows there. ("bestbound=1" is
# one of them: NEOS's driver reads it as cvt:pre:boundsbest, which passes the
# bound of 1 on the return arcs to the solver.) That run was killed 1,300
# seconds in, just after its root, hence the time limit of a round.
OPTIONS = ("outlev=1 method=2 mipfocus=3 mipgap=0 mipgapabs=0.99 bestbound=1 "
           "threads=4 nodefilestart=0.5 nodelim=1")
ROUND = 1150             # seconds a round
TOTAL = 27000            # seconds for the whole job; NEOS allows 28,800


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


def write(g, cols, cuts, dirname, x0=None, options=OPTIONS, round_seconds=ROUND, total=TOTAL):
    num = Numbering(g, cols)
    os.makedirs(dirname, exist_ok=True)
    with open(os.path.join(dirname, "palindrome.mod"), "w") as f:
        f.write(MODEL)
    with open(os.path.join(dirname, "palindrome.run"), "w") as f:
        f.write(RUN.format(options=options, round_seconds=round_seconds, total=total,
                           rounds=total // round_seconds + 1))
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
    ap.add_argument("--options", default=OPTIONS, help="the gurobi_options of a round")
    ap.add_argument("--round", type=int, default=ROUND, help="seconds a round")
    ap.add_argument("--total", type=int, default=TOTAL, help="seconds for the whole job")
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
    write(g, cols, cuts, args.dir, x0, args.options, args.round, args.total)
    for name in ("palindrome.mod", "palindrome.dat", "palindrome.run"):
        p = os.path.join(args.dir, name)
        print(f"wrote {p}: {os.path.getsize(p):,} bytes")


if __name__ == "__main__":
    main()
