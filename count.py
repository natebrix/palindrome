"""How many distinct phrases can one palindromic sentence hold?

    python3 count.py [data/npdict.txt]

Prints the size of the leftover-state graph, the number of phrases that appear
in at least one palindromic sentence, and the check that all of them lie on a
closed walk through the level state E, so that a single sentence holds them all.
"""
import sys
import time
from collections import Counter

from palgraph import Graph, load

path = sys.argv[1] if len(sys.argv) > 1 else "data/npdict.txt"
t0 = time.time()
D = load(path)
print(f"phrases (distinct keys): {len(D):,}   letters: {sum(map(len, D)):,}")
g = Graph(D.keys())
print(f"graph: {g.n:,} states, {g.m:,} arcs")

co = g.coreachable()
print(f"states that can still finish a palindrome: {sum(co):,}")
U = g.includable()
print(f"phrases in at least one palindromic sentence: {len(U):,}")
print(f"phrases in none: {len(D) - len(U):,}")

# The check. If T is a palindromic sentence then so is T T, level in the
# middle, so every state that can finish should share E's component.
comp, ncomp = g.scc(co)
c0 = comp[0]
size = Counter(c for c in comp if c >= 0)
W0 = {g.lab[a] for a in range(g.m) if comp[g.src[a]] == c0 and comp[g.dst[a]] == c0}
print(f"strongly connected components among those states: {ncomp:,}")
print(f"states in E's component: {size[c0]:,}")
print(f"phrases on a closed walk through E: {len(W0):,}")
assert W0 == U, "some includable phrase is not on a closed walk through E"
print(f"maximum distinct phrases in one palindromic sentence: {len(U):,}")
print(f"letters if each were used once (a floor on any complete sentence): "
      f"{sum(len(g.words[i]) for i in U):,}")
print(f"({time.time() - t0:.0f}s)")
