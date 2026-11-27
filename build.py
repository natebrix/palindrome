"""Build one palindromic sentence that holds every includable phrase.

    python3 build.py [data/npdict.txt] [out/complete_palindrome.txt]

Greedy, and not the shortest: for each phrase not yet covered, costliest first,
splice in its cheapest closed walk through the level state E (cheapest in
letters). The closed walks nest: left halves in order, right halves in reverse.
The result is checked by reversing the string.
"""
import heapq
import os
import sys
import time

from palgraph import Graph, letters, load

dict_path = sys.argv[1] if len(sys.argv) > 1 else "data/npdict.txt"
out_path = sys.argv[2] if len(sys.argv) > 2 else "out/complete_palindrome.txt"

t0 = time.time()
D = load(dict_path)
g = Graph(D.keys())
src, dst, lab, words = g.src, g.dst, g.lab, g.words
wlen = [len(w) for w in words]
INF = float("inf")


def dijkstra(forward):
    """Fewest letters from E to each state (forward) or from each state back
    to E, with the arc that achieves it."""
    dist = [INF] * g.n
    via = [-1] * g.n
    dist[0] = 0
    pq = [(0, 0)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist[u]:
            continue
        for a in (g.out(u) if forward else g.inn(u)):
            v = dst[a] if forward else src[a]
            nd = d + wlen[lab[a]]
            if nd < dist[v]:
                dist[v] = nd
                via[v] = a
                heapq.heappush(pq, (nd, v))
    return dist, via


df, pf = dijkstra(True)
db, pb = dijkstra(False)

best = {}  # word id -> (letters, arc) of its cheapest closed walk through E
for a in range(g.m):
    c = df[src[a]] + wlen[lab[a]] + db[dst[a]]
    if c < INF and c < best.get(lab[a], (INF,))[0]:
        best[lab[a]] = (c, a)


def closed_walk(a):
    arcs = []
    u = src[a]
    while u != 0:
        arcs.append(pf[u])
        u = src[pf[u]]
    arcs.reverse()
    arcs.append(a)
    v = dst[a]
    while v != 0:
        arcs.append(pb[v])
        v = dst[pb[v]]
    return arcs


covered = set()
shells = []
for w in sorted(best, key=lambda w: -best[w][0]):
    if w in covered:
        continue
    arcs = closed_walk(best[w][1])
    left = [lab[a] for a in arcs if g.side(a) == "L"]
    right = [lab[a] for a in arcs if g.side(a) == "R"][::-1]
    shells.append((left, right))
    covered.update(left)
    covered.update(right)

left = [i for L, _ in shells for i in L]
right = [i for _, R in reversed(shells) for i in R]
keys = [words[i] for i in left + right]
s = "".join(keys)
assert s == s[::-1], "not a palindrome"
assert set(left + right) == g.includable(), "a phrase is missing"

text = ", ".join(D[k] for k in keys)
# read the text back the way a reader would: strip to letters and reverse
t = letters(text)
assert t == t[::-1] == s

os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
with open(out_path, "w", encoding="utf-8") as f:
    f.write(text)
print(f"{len(keys):,} phrases, {len(set(keys)):,} distinct, {len(s):,} letters, "
      f"{len(shells):,} nested closed walks")
print(f"floor (each includable phrase once): {sum(wlen[i] for i in best):,} letters")
print(f"wrote {out_path} ({time.time() - t0:.0f}s)")
