# palindrome

How many distinct phrases can one palindromic sentence hold, and how short can
such a sentence be?

The setting is Peter Norvig's [palindrome dictionary](https://norvig.com/palindrome.html)
of noun phrases (`npdict.txt`) and Roland Meertens's 2018 post
[World's longest palindrome?](https://meertens.dev/blog/worlds-longest-palindrome/),
which builds a sentence of 64,810 distinct phrases and says it is safe to assume
that no more can be included. This repository settles that count exactly.

## Results

A sentence is a sequence of dictionary phrases, repeats allowed, whose letters
read the same in both directions. A phrase's key is its text lowercased with
non-word characters removed, as in Norvig's and Meertens's code (so digits stay).

| | Phrases | Distinct | Letters |
|---|---|---|---|
| Meertens's published sentence | 954,401 | 64,810 | 4,836,010 |
| Sentence built by `build.py` | 309,748 | 64,811 | 1,406,982 |
| Floor: each includable phrase once | 64,811 | 64,811 | 488,137 |

- The dictionary has 125,512 distinct keys. Exactly 64,811 of them appear in at
  least one palindromic sentence, and a single sentence can hold all of them.
  The other 60,701 can appear in none.
- Meertens's sentence lacks one includable phrase, "M2M".
- The shortest sentence that holds all 64,811 has between 488,137 and 1,406,982
  letters. That gap is open.

Letters here are key characters. The post's figure of 5,813,576 also counts
commas and other punctuation.

## How it works

Build a sentence from the outside in, always adding the next word to the side
that is behind, and to the left on a tie. The state is the unmatched leftover:

- `E`: both sides are level.
- `(L, s)`: the left is ahead by `s`, a suffix of the last word placed there.
- `(R, p)`: the right is ahead by `p`, a prefix of the last word placed there.

An arc places one word. A state accepts when its leftover is itself a
palindrome. Every palindromic sentence is exactly one walk from `E` to an
accepting state, and every such walk is a sentence. For `npdict.txt` the graph
has 716,216 states and 1,263,316 arcs.

A word is includable when it labels an arc whose head can still reach an
accepting state. One sentence holds them all because, if `T` is a palindromic
sentence, then `T T` is one too and is level in the middle: every includable
word lies on a closed walk through `E`, and closed walks through `E` nest.
`count.py` checks this on the data rather than relying on the argument: the
135,455 states that can finish form a single strongly connected component.

## Running it

Python 3 with the standard library only. Fetch the dictionary first; it is
Norvig's file and is not kept in this repository.

```bash
mkdir -p data && curl -L -o data/npdict.txt https://norvig.com/npdict.txt
```

```bash
python3 count.py
```

```bash
python3 build.py
```

```bash
python3 check_sentence.py out/complete_palindrome.txt
```

```bash
python3 test_bruteforce.py
```

| File | Purpose |
|---|---|
| `palgraph.py` | The leftover-state graph: construction, backward reachability, strongly connected components |
| `count.py` | The exact count, about 5 seconds |
| `build.py` | A greedy sentence holding every includable phrase, checked by reversing the string, about 8 seconds |
| `check_sentence.py` | Checks any sentence in the `", "`-joined format and lists the includable phrases it lacks |
| `test_bruteforce.py` | Compares the graph's walks with brute-force enumeration on small random dictionaries |

To check Meertens's sentence, pass `longest_palindrome.txt` from
[rmeertens/palindromes](https://github.com/rmeertens/palindromes) to
`check_sentence.py`.

## What is not done

`build.py` is a heuristic: for each phrase not yet covered, costliest first, it
splices in that phrase's cheapest closed walk through `E`. Nothing here bounds
the shortest complete sentence from below beyond the trivial floor.

The natural model is a directed rural postman problem on this graph: integer
arc flows, flow conservation, a cost of one per letter placed, and for each
includable word a row requiring at least one unit over the arcs it labels.
Unlike the portmantout, where each required word is a single arc and the linear
program is a pure network problem, those covering rows break total
unimodularity, so this one is a genuine integer program.
