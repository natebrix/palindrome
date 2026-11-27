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
| Greedy sentence, `build.py` | 309,748 | 64,811 | 1,406,982 |
| Short sentence, `shortest.py` | 153,781 | 64,811 | 815,407 |
| Lower bound from the linear program, `lp.py` | | | 813,106 |
| Floor: each includable phrase once | 64,811 | 64,811 | 488,137 |

- The dictionary has 125,512 distinct keys. Exactly 64,811 of them appear in at
  least one palindromic sentence, and a single sentence can hold all of them.
  The other 60,701 can appear in none.
- Meertens's sentence lacks one includable phrase, "M2M".
- The shortest sentence that holds all 64,811 has between 813,106 and 815,407
  letters, a gap of 2,301 letters or 0.28%.

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

## The shortest complete sentence

A complete sentence is a walk in the graph that uses every includable word at
least once, and its length is the letters of the words on its arcs. That is a
directed rural postman problem with one twist: a word labels many arcs, one for
each state it can be played from, and any one of them will do.

`lp.py` solves the linear relaxation on the 135,455 states that can finish:

- a flow on each of 426,676 arcs, at a cost of the word's letters;
- flow in equals flow out at every state;
- for each of the 64,811 includable words, at least one unit of flow over the
  arcs it labels;
- 903 free return arcs from the accepting states to `E`, carrying at most one
  unit between them, which close the single walk.

Its optimum is 813,105.25, so no complete sentence has fewer than 813,106
letters. Two things are relaxed, whole-number flows and the flow hanging
together as one walk, and only the first one bites: the optimal flow is already
in one piece, with 2,020 of its 104,420 arcs in use at a fractional value,
nearly all at one half. 64,180 words have a full unit on a single arc and 631
are split.

`shortest.py` turns that into a sentence by diving. Each round it requires a
unit of flow on one arc for every word that already has one, and for the surest
quarter of the words still split, and solves again. Once every word has its
arc, the covering rows are redundant and what is left is a circulation with
lower bounds, whose vertices are whole numbers. The flow is then read off as an
Eulerian circuit by Hierholzer's algorithm and checked by reversing the string.

| Rounding | Letters | Gap | Time |
|---|---|---|---|
| All split words at once, `--share 1` | 815,933 | 2,827 | 2 min |
| A quarter per round, the default | 815,451 | 2,345 | 13 min |
| A tenth per round, `--share 0.1 --resolve simplex` | 815,407 | 2,301 | 27 min |

Times are on four cores with HiGHS 1.16, and the last was measured with other
solves running beside it. Each interior point solve takes about a minute; dual
simplex from a cold start is far slower on this model.

`test_lp.py` checks all of this on small random dictionaries, where the
shortest complete sentence can be found exactly: the LP optimum is never above
it and the sentence built is never below it. On those small cases the flow
often does fall into pieces, and the join step is exercised.

## Running it

Python 3. The count, the greedy sentence and the checker use the standard
library only; `lp.py`, `shortest.py` and `test_lp.py` need `numpy` and
`highspy`. Fetch the dictionary first; it is Norvig's file and is not kept in
this repository.

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

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

```bash
.venv/bin/python lp.py --solver ipm
```

```bash
.venv/bin/python shortest.py
```

```bash
.venv/bin/python test_lp.py
```

| File | Purpose |
|---|---|
| `palgraph.py` | The leftover-state graph: construction, backward reachability, strongly connected components |
| `count.py` | The exact count, about 5 seconds |
| `build.py` | A greedy sentence holding every includable phrase, checked by reversing the string, about 8 seconds |
| `check_sentence.py` | Checks any sentence in the `", "`-joined format and lists the includable phrases it lacks |
| `lp.py` | The linear program and its lower bound, about a minute with `--solver ipm` |
| `shortest.py` | LP, rounding, joining, Hierholzer: a short complete sentence and its gap to the bound |
| `test_bruteforce.py` | Compares the graph's walks with brute-force enumeration on small random dictionaries |
| `test_lp.py` | Compares the LP bound and the sentence built with exact answers on small random dictionaries |

To check Meertens's sentence, pass `longest_palindrome.txt` from
[rmeertens/palindromes](https://github.com/rmeertens/palindromes) to
`check_sentence.py`.

## What is not done

The 2,301 letters between the bound and the best sentence are open. Closing
them needs either a stronger bound or a better sentence:

- The integer program itself, the model of `lp.py` with whole-number flows. HiGHS
  did not finish its root relaxation in 15 minutes, because its branch and
  bound solves that relaxation by simplex.
- Cuts. The fractional flows are nearly all one half, which suggests a parity
  argument of the kind that gives blossom inequalities for matching. No such
  cut has been tried.
- A search over the 631 split words alone, with the other 64,180 held at the
  arcs the LP gives them. That would be a heuristic too: nothing says the
  optimum uses those arcs.
