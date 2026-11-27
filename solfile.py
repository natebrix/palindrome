"""Move flows between sentences and the solution files other solvers use.

    python3 solfile.py write SENTENCE.txt START.mst     a MIP start from a sentence
    python3 solfile.py read SOLUTION.sol SENTENCE.txt   a sentence from a solution

Both file kinds are lines of "name value", with # comments, as Gurobi writes
and reads them; columns are named c0, c1, ... as in export.py. `write` lists
only the arcs in use. `read` checks the flow, joins it up if it is in pieces
(which adds letters), and checks the sentence by reversing it.

`read` also takes the output of the run file written by ampl_export.py: lines
"x tail word head value" and "r state value" anywhere in the text, such as a
NEOS result saved to a file.
"""
import sys

import numpy as np

import lp as L
import shortest as S
from palgraph import Graph, letters, load


def main():
    mode, src, dst = sys.argv[1:4]
    dict_path = sys.argv[4] if len(sys.argv) > 4 else "data/npdict.txt"
    D = load(dict_path)
    g = Graph(D.keys())
    cols = L.build(g)
    if mode == "write":
        keys = [letters(p) for p in open(src, encoding="utf-8").read().split(", ")]
        x = S.flow_of_sentence(g, cols, keys)
        S.check_flow(g, cols, x)
        with open(dst, "w") as f:
            f.write(f"# MIP start: {int((x * cols['cost']).sum())} letters\n")
            for j in np.nonzero(x)[0]:
                f.write(f"c{j} {x[j]}\n")
        print(f"wrote {dst}: {int((x > 0).sum()):,} arcs in use, "
              f"{int((x * cols['cost']).sum()):,} letters")
    elif mode == "read":
        x = np.zeros(len(cols["arc"]))
        arcs = rets = None
        for line in open(src, errors="replace"):
            parts = line.split()
            if len(parts) == 2 and parts[0][0] == "c" and parts[0][1:].isdigit():
                x[int(parts[0][1:])] = float(parts[1])
            elif parts and parts[0] in ("x", "r") and len(parts) == (5 if parts[0] == "x" else 3):
                if arcs is None:
                    from ampl_export import Numbering
                    arcs, rets = Numbering(g, cols).column()
                try:
                    v = [int(p) for p in parts[1:]]
                except ValueError:
                    continue
                x[arcs[(v[0], v[1])] if parts[0] == "x" else rets[v[0]]] = v[-1]
        xi = np.round(x).astype(np.int64)
        assert np.abs(x - xi).max() < 1e-6, "the solution is not in whole numbers"
        S.check_flow(g, cols, xi)
        print(f"flow: {int((xi * cols['cost']).sum()):,} letters")
        xi = S.connect(g, cols, xi, print)
        keys = S.sentence(g, cols, S.euler(cols, xi))
        s = "".join(keys)
        assert s == s[::-1] and {g.wid[k] for k in keys} == g.includable()
        with open(dst, "w", encoding="utf-8") as f:
            f.write(", ".join(D[k] for k in keys))
        print(f"sentence: {len(keys):,} phrases, {len(s):,} letters; wrote {dst}")
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
