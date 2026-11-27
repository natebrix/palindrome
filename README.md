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
| Shortest sentence found, `shortest.py --cuts` | 153,658 | 64,811 | 814,515 |
| Lower bound, branch and bound at its root, `mip.py` | | | 814,399 |
| Lower bound, linear program with parity cuts, `cuts.py` | | | 814,371 |
| Lower bound, linear program alone, `lp.py` | | | 813,106 |
| Floor: each includable phrase once | 64,811 | 64,811 | 488,137 |

- The dictionary has 125,512 distinct keys. Exactly 64,811 of them appear in at
  least one palindromic sentence, and a single sentence can hold all of them.
  The other 60,701 can appear in none.
- Meertens's sentence lacks one includable phrase, "M2M".
- The shortest sentence that holds all 64,811 has between 814,399 and 814,515
  letters, a gap of 116 letters or 0.014%.

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

### Parity cuts

The split words are split for a reason of parity. Take the palindromic word
"ababa". Placed on the left from `E` it leaves `(L, ababa)`, and placed on the
right from there it leads back to `E`. The LP sends half a unit around that
loop: the word is covered once, for 5 letters. A whole-number flow has to pay
for both copies, the sentence "ababa, ababa", or reach the word another way.
Longer loops through several words do the same.

The cut that removes this takes a set `S` of states and an odd set `W` of
words. Add the covering rows of `W` to the flow balance of `S`, halve, and
round up:

    x(arcs labelled by W that do not leave S)
      + x(arcs into S not labelled by W)        >= (|W| + 1) / 2

It is violated when every word of `W` is covered exactly once by flow crossing
the boundary of `S` and nothing else crosses: the crossing flow is then `|W|`,
an odd number, although as much must enter `S` as leave it. For "ababa", `S` is
the one state `(L, ababa)` and `W` is the word alone, and the cut asks for a
whole unit where the LP has a half.

`cuts.py` finds such cuts in the LP solution. For any `S` the best `W` and the
violation come from one pass over the arcs at `S`. Candidates for `S` come from
contracting every arc that carries flow for a word covered more than once, then
asking that all the arcs of each split word cross together, which is a linear
system over GF(2): its null space and the single contracted states are tried.
This is a heuristic, not an exact separation.

| | Bound | Cuts | Fractional arcs |
|---|---|---|---|
| LP alone | 813,105.25 | 0 | 2,020 |
| After 1 round | 814,185.75 | 468 | 1,705 |
| After 3 rounds | 814,317.75 | 730 | 1,291 |
| After 8 rounds | 814,350.83 | 1,076 | 1,304 |
| After 16 rounds | 814,370.74 | 1,460 | 1,277 |

The first round does most of the work, and by round 16 each round adds about
45 cuts for half a letter. Many of the fractions that remain are thirds and
stranger values, which cuts of this one kind reach less and less.

The cuts help the sentence more than the bound. Rounding from the tightened LP
hardly moves the objective until the last hundred or so words:

| Rounding from | Letters |
|---|---|
| LP alone, a tenth per round | 815,407 |
| 888 cuts (5 rounds), a quarter per round | 814,607 |
| 1,460 cuts (16 rounds), a quarter per round | 814,521 |
| 1,460 cuts (16 rounds), a tenth per round, `--share 0.1` | 814,515 |

`mip.py` gives the tightened model to HiGHS as an integer program, with a
sentence as its starting point and interior point at the root. Run with the
1,416 cuts of round 15 and a 40-minute limit, it never left the root node, but
its own cutting planes raised the bound to 814,398.99, that is 814,399. That figure rests on the solver's cuts and is not
checked independently here; 814,371 rests only on the cuts above.

### Checks

`test_lp.py` checks all of this on small random dictionaries, where the
shortest complete sentence can be found exactly: the LP optimum, with and
without parity cuts, is never above it, and the sentence built is never below
it. Every parity cut found there holds for the whole-number flow of the
sentence built. On those small cases the flow often does fall into pieces, and
the join step is exercised.

## Running it

Python 3. The count, the greedy sentence and the checker use the standard
library only; `lp.py`, `cuts.py`, `shortest.py`, `mip.py` and `test_lp.py` need
`numpy` and `highspy`. Fetch the dictionary first; it is Norvig's file and is not kept in
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
.venv/bin/python cuts.py --rounds 16
```

```bash
.venv/bin/python shortest.py --cuts out/cuts.pkl --share 0.1
```

```bash
.venv/bin/python mip.py --cuts out/cuts.pkl --start out/short_palindrome.txt --time-limit 2400
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
| `cuts.py` | Parity cuts: separation, and the loop that tightens the LP, about a minute a round |
| `shortest.py` | LP (with cuts if given), rounding, joining, Hierholzer: a short complete sentence and its gap to the bound |
| `mip.py` | The tightened model as an integer program in HiGHS, started from a sentence |
| `test_bruteforce.py` | Compares the graph's walks with brute-force enumeration on small random dictionaries |
| `test_lp.py` | Compares the LP bound, the parity cuts and the sentence built with exact answers on small random dictionaries |

To check Meertens's sentence, pass `longest_palindrome.txt` from
[rmeertens/palindromes](https://github.com/rmeertens/palindromes) to
`check_sentence.py`.

## What is not done

The last 116 letters are open.

- The parity cuts here are one family, found by a heuristic. An exact
  separation, or cuts that combine the parity cuts themselves, might lift the
  bound past the plateau at 814,371.
- Branch and bound has not branched. HiGHS spent its 40 minutes on cutting
  planes at the root; closing the gap that way needs far more time or a faster
  way to solve each node.
- Rounding is still a dive without backtracking. The last hundred words it
  settles account for most of what the sentence costs above the bound.
