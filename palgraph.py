"""Leftover-state graph of palindromic sentences over a phrase dictionary.

A sentence is a sequence of dictionary keys (letters only, lowercase) whose
concatenation is a palindrome. Build it from the outside in, adding the next
word to the side that is behind (the left on a tie). The state is the unmatched
leftover:
  E        both sides level
  (L, s)   the left is ahead by s, a suffix of the last left word
  (R, p)   the right is ahead by p, a prefix of the last right word
Each arc places one word. A state accepts when its leftover is a palindrome.
Every sentence is exactly one walk from E to an accepting state, and every such
walk is a sentence.
"""
import re
from array import array
from collections import defaultdict, deque

_sub = re.compile(r"[\W]+").sub


def letters(phrase):
    return _sub("", phrase).lower()


def load(path):
    d = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.strip()
            k = letters(p)
            if k:
                d[k] = p
    return d


class Graph:
    def __init__(self, keys):
        self.words = sorted(keys)
        self.wid = {w: i for i, w in enumerate(self.words)}
        wid = self.wid
        by_suffix = defaultdict(list)  # proper suffix -> words ending in it
        by_prefix = defaultdict(list)  # proper prefix -> words starting with it
        for w, i in wid.items():
            for k in range(1, len(w)):
                by_suffix[w[k:]].append(i)
                by_prefix[w[:k]].append(i)
        # node 0 is E
        self.kind = [0]
        self.text = [""]
        ids = ({}, {"": 0}, {"": 0})  # index 1: L strings, 2: R strings
        ids[1][""] = 0
        src, dst, lab = array("i"), array("i"), array("i")
        queue = deque([0])
        words = self.words

        def node(kind, s):
            if not s:
                return 0
            m = ids[kind]
            n = m.get(s)
            if n is None:
                n = len(self.kind)
                m[s] = n
                self.kind.append(kind)
                self.text.append(s)
                queue.append(n)
            return n

        while queue:
            u = queue.popleft()
            kind, s = self.kind[u], self.text[u]
            if kind == 0:
                for w, i in wid.items():
                    src.append(u); dst.append(node(1, w)); lab.append(i)
            elif kind == 1:
                # a word goes on the right; its reverse must line up with s
                n = len(s)
                for k in range(1, n + 1):
                    i = wid.get(s[:k][::-1])
                    if i is not None:
                        src.append(u); dst.append(node(1, s[k:])); lab.append(i)
                for i in by_suffix.get(s[::-1], ()):
                    w = words[i]
                    src.append(u); dst.append(node(2, w[: len(w) - n])); lab.append(i)
            else:
                # a word goes on the left
                n = len(s)
                for k in range(1, n + 1):
                    i = wid.get(s[n - k:][::-1])
                    if i is not None:
                        src.append(u); dst.append(node(2, s[: n - k])); lab.append(i)
                for i in by_prefix.get(s[::-1], ()):
                    w = words[i]
                    src.append(u); dst.append(node(1, w[n:])); lab.append(i)
        self.src, self.dst, self.lab = src, dst, lab
        self.n = len(self.kind)
        self.m = len(src)
        self.accept = [t == t[::-1] for t in self.text]
        # CSR, forward and backward (arc indices)
        self.fptr, self.farc = _csr(self.n, src)
        self.bptr, self.barc = _csr(self.n, dst)

    def out(self, u):
        return self.farc[self.fptr[u]: self.fptr[u + 1]]

    def inn(self, v):
        return self.barc[self.bptr[v]: self.bptr[v + 1]]

    def side(self, a):
        """'L' if the arc's word goes on the left, 'R' if on the right."""
        return "R" if self.kind[self.src[a]] == 1 else "L"

    def includable(self):
        """Ids of the words that appear in at least one palindromic sentence:
        the labels of arcs whose head can still finish. (Every state here is
        reachable from E, since the graph is built by search from E.)"""
        co = self.coreachable()
        dst, lab = self.dst, self.lab
        return {lab[a] for a in range(self.m) if co[dst[a]]}

    def coreachable(self):
        """Nodes from which an accepting state can be reached."""
        seen = bytearray(self.n)
        stack = [v for v in range(self.n) if self.accept[v]]
        for v in stack:
            seen[v] = 1
        src = self.src
        while stack:
            v = stack.pop()
            for a in self.inn(v):
                u = src[a]
                if not seen[u]:
                    seen[u] = 1
                    stack.append(u)
        return seen

    def scc(self, alive):
        """Strongly connected components of the subgraph on alive nodes
        (iterative Tarjan). Returns comp[] (-1 for dead nodes) and the count;
        components are numbered in reverse topological order."""
        n = self.n
        index = [-1] * n
        low = [0] * n
        comp = [-1] * n
        onstack = bytearray(n)
        st = []
        counter = 0
        ncomp = 0
        fptr, farc, dst = self.fptr, self.farc, self.dst
        for root in range(n):
            if not alive[root] or index[root] != -1:
                continue
            work = [(root, fptr[root])]
            index[root] = low[root] = counter; counter += 1
            st.append(root); onstack[root] = 1
            while work:
                u, pos = work[-1]
                end = fptr[u + 1]
                advanced = False
                while pos < end:
                    v = dst[farc[pos]]
                    pos += 1
                    if not alive[v]:
                        continue
                    if index[v] == -1:
                        work[-1] = (u, pos)
                        index[v] = low[v] = counter; counter += 1
                        st.append(v); onstack[v] = 1
                        work.append((v, fptr[v]))
                        advanced = True
                        break
                    if onstack[v] and index[v] < low[u]:
                        low[u] = index[v]
                if advanced:
                    continue
                work.pop()
                if work:
                    p = work[-1][0]
                    if low[u] < low[p]:
                        low[p] = low[u]
                if low[u] == index[u]:
                    while True:
                        x = st.pop(); onstack[x] = 0
                        comp[x] = ncomp
                        if x == u:
                            break
                    ncomp += 1
        return comp, ncomp


def _csr(n, key):
    ptr = [0] * (n + 1)
    for k in key:
        ptr[k + 1] += 1
    for i in range(n):
        ptr[i + 1] += ptr[i]
    fill = ptr[:]
    arcs = [0] * len(key)
    for a, k in enumerate(key):
        arcs[fill[k]] = a
        fill[k] += 1
    return ptr, arcs
