"""Check the graph against brute force on small random dictionaries.

For each dictionary: enumerate every word sequence of at most K tokens whose
concatenation is a palindrome, and compare with the walks of the graph.
  soundness:    every walk E -> accepting of at most K arcs is a palindrome
  completeness: every palindromic sequence is such a walk (same multiset count)
"""
import itertools
import random
import sys

from palgraph import Graph


def walks(g, K):
    """All (left tokens, right tokens) of walks from E to an accepting state
    with 1..K arcs."""
    out = []

    def rec(u, left, right):
        if (left or right) and g.accept[u]:
            out.append(tuple(left) + tuple(reversed(right)))
        if len(left) + len(right) == K:
            return
        for a in g.out(u):
            w = g.words[g.lab[a]]
            if g.side(a) == "L":
                left.append(w); rec(g.dst[a], left, right); left.pop()
            else:
                right.append(w); rec(g.dst[a], left, right); right.pop()

    rec(0, [], [])
    return out


def main():
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    trials = 0
    total = 0
    for _ in range(400):
        alpha = "ab" if rng.random() < 0.5 else "abc"
        nw = rng.randint(2, 6)
        d = set()
        while len(d) < nw:
            d.add("".join(rng.choice(alpha) for _ in range(rng.randint(1, 5))))
        K = 5 if nw <= 4 else 4
        brute = sorted(
            seq
            for k in range(1, K + 1)
            for seq in itertools.product(sorted(d), repeat=k)
            if "".join(seq) == "".join(seq)[::-1]
        )
        g = Graph(d)
        got = sorted(walks(g, K))
        assert got == brute, (d, K, set(got) ^ set(brute))
        trials += 1
        total += len(brute)
    print(f"{trials} dictionaries, {total} palindromic sentences: graph walks == brute force")


if __name__ == "__main__":
    main()
