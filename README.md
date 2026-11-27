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
| Shortest sentence found, Gurobi's root heuristics from the one below | 153,710 | 64,811 | 814,419 |
| Shortest sentence found here, `shortest.py --cuts --search` | 153,656 | 64,811 | 814,455 |
| Lower bound, branch and bound at its root, `mip.py` | | | 814,402 |
| Lower bound, linear program with parity cuts, `cuts.py` | | | 814,391 |
| Lower bound, linear program alone, `lp.py` | | | 813,106 |
| Floor: each includable phrase once | 64,811 | 64,811 | 488,137 |

- The dictionary has 125,512 distinct keys. Exactly 64,811 of them appear in at
  least one palindromic sentence, and a single sentence can hold all of them.
  The other 60,701 can appear in none.
- Meertens's sentence lacks one includable phrase, "M2M".
- The shortest sentence that holds all 64,811 has between 814,402 and 814,419
  letters, a gap of 17 letters or 0.002%.

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
This is a heuristic.

| | Bound | Cuts | Fractional arcs |
|---|---|---|---|
| LP alone | 813,105.25 | 0 | 2,020 |
| After 1 round | 814,185.75 | 468 | 1,705 |
| After 3 rounds | 814,317.75 | 730 | 1,291 |
| After 8 rounds | 814,350.83 | 1,076 | 1,304 |
| After 16 rounds | 814,370.74 | 1,460 | 1,277 |

The first round does most of the work, and by round 16 each round adds about
45 cuts for half a letter.

`exact.py` separates the same cuts exactly. The least slack over all `S` and
odd `W` is an integer program in one 0-1 variable per state and per word, and
three reductions that lose no violated cut make it small: a word covered by a
single arc at flow one is in `W` exactly when that arc crosses, so it only
feeds the parity; an arc at flow one or more of a word covered at least twice
must not cross, so its states merge; and a state touched only by arcs of the
first kind changes nothing. For the LP after 16 rounds that leaves 570 states
and 342 words. It is solved again and again, each time keeping the states of
the cuts already found out of `S`, and it finds what the heuristic misses:
cuts violated by less than a half, on larger sets of states.

| Resumed with `--exact 90` | Bound | Cuts |
|---|---|---|
| The 16 rounds above | 814,370.74 | 1,460 |
| After 1 more round | 814,385.74 | 1,524 |
| After 3 | 814,389.57 | 1,624 |
| After 5 | 814,389.90 | 1,727 |
| After 9 | 814,390.74 | 1,891 |

One round of exact separation gains 15 letters where the heuristic was gaining
a half, and then this stalls too, near 814,391. What is still fractional then
is no longer a set of small loops: it is one cluster of some 460 words, most of
them short bridging words used many times over, with flows at thirds and odder
values. Cuts of this one family do not reach it.

The cuts help the sentence more than the bound. Rounding from the tightened LP
hardly moves the objective until the last hundred or so words:

| Rounding from | Letters |
|---|---|
| LP alone, a tenth per round | 815,407 |
| 888 cuts (5 rounds), a quarter per round | 814,607 |
| 1,460 cuts (16 rounds), a quarter per round | 814,521 |
| 1,460 cuts (16 rounds), a tenth per round, `--share 0.1` | 814,515 |
| 1,727 cuts, a quarter per round | 814,503 |
| 1,891 cuts, branch and bound on the split words, `--search` | 814,459 |
| 1,460 cuts, branch and bound on the split words, `--search` | 814,455 |

The last two rows replace the dive. The LP with 1,891 cuts gives 64,575 words
a full unit on one arc. Pinning them there costs nothing, and only the arcs of
the other 236 words, 1,506 of them, then need whole-number flows: with those
fixed the rest is a circulation. HiGHS solves that to optimality in about 20
minutes and 382 nodes: 814,459 is the best sentence that keeps those 64,575
words where that LP puts them. The LP with 1,460 cuts pins fewer words, 64,479,
and leaves 332 to the search, which is harder and reaches further: 814,455 at
its 25-minute limit, with 814,454 not ruled out for that pinning. A shorter
sentence than these must move some of the pinned words.

`mip.py` gives the whole tightened model to HiGHS as an integer program, with
a sentence as its starting point and interior point at the root. It does not
get past the root node, but its own cutting planes raise the bound: to 814,399
with the 1,416 cuts of round 15 in 40 minutes, and to 814,402 with
the 1,891 cuts in 30. That figure rests on
the solver's cuts and is not checked independently here; 814,391 rests only on
the parity cuts.

### Another solver

`export.py` writes the tightened model as an MPS file, for any solver that
reads one, and `ampl_export.py` writes it as AMPL files, which are small
enough for the NEOS server's 16 MB limit on a file where the MPS file is not.
`scip_run.py` solves the MPS file with SCIP 10, where the parity cuts can
also be separated inside the search by a plugin, which HiGHS does not allow.

SCIP's search is at least as good as HiGHS's here and its first LP is far
worse. On the pinned model with 1,891 cuts, the one HiGHS solves in 21 minutes
and 382 nodes, SCIP proves the same 814,459 in 25 minutes and 167 nodes, with a
root bound of 814,455.4 against 814,452.8. But 550 seconds of that go to the
first LP, which HiGHS solves by interior point in 36: SCIP's LP solver,
SoPlex, has only the simplex method, and this LP is badly degenerate, with
274,000 of its 427,579 arcs at reduced cost zero in the optimum. After that
first solve SCIP's cut rounds take about a second each.

On the full model that first solve is the whole story: SCIP never got through
it. Three runs of 200 to 210 minutes each, with SoPlex's default, steepest edge
and devex pricing, all ended at their time limits still inside the first LP,
after 750,057, 747,484 and 967,698 simplex iterations on the presolved model
(337,000 columns, 112,000 rows). None produced a bound. A solver for this model
needs an interior point method for its first LP.

Gurobi 13, on the NEOS server through the AMPL files, has one. Its barrier
method solved the first LP of the full model in 194 seconds (814,393.07, with
the 1,891 cuts and its own presolve), its cuts took the root bound to 814,401.5,
the same 814,402 as HiGHS, and its heuristics found flows of 814,449, 814,443,
814,437, 814,423 and 814,419 letters at the root. Then the job was killed, nine
nodes into the search, for passing NEOS's 3 GB of memory, and a killed job
prints no solution.

A second job stopped at the root, as soon as the flow of 814,419 was found, and
printed it. It is one walk as it stands, so it is a sentence of 814,419 letters,
36 shorter than the best found here, and `check_sentence.py` passes it. The
same job then tried the search on one thread under a memory limit of the
solver's own, 1.9 GB, and reached that limit in the root's third round of cuts:
branch and bound on the whole model does not fit on NEOS. What does fit is the
root, where the shorter flows came from, so the run file now repeats the root
alone, each round from the best flow so far with a new random seed, and prints
what changed after each.

The parity plugin did not pay on the pinned model: 199 cuts in 31 calls, a root
bound of 814,456.2, and no node left after the same 25 minutes. Rewriting the
model as a deviation from the best sentence, so that the simplex method starts
feasible and 64 letters from the LP optimum, did not help either: 134,000
primal iterations in nine minutes without progress.

### Checks

`test_lp.py` checks all of this on small random dictionaries, where the
shortest complete sentence can be found exactly: the LP optimum, with and
without parity cuts, is never above it, and the sentence built is never below
it. Every parity cut found there holds for the whole-number flow of the
sentence built. On those small cases the flow often does fall into pieces, and
the join step is exercised.

`test_exact.py` checks the exact separation against brute force: for small LP
solutions it tries every set of states and compares the most violated cut with
the one the integer program returns.

## Running it

Python 3. The count, the greedy sentence and the checker use the standard
library only; the rest need `numpy` and `highspy`, and `scip_run.py` needs
`pyscipopt` as well. Fetch the dictionary first; it is Norvig's file and is not kept in
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
.venv/bin/python shortest.py --cuts out/cuts.pkl --search 1500
```

```bash
.venv/bin/python cuts.py --resume --rounds 9 --exact 90
```

```bash
.venv/bin/python mip.py --cuts out/cuts.pkl --start out/short_palindrome.txt --time-limit 1800
```

```bash
.venv/bin/python export.py --cuts out/cuts.pkl --out out/palindrome.mps
```

```bash
.venv/bin/python scip_run.py out/palindrome.mps --start out/short_palindrome.txt --time-limit 12600
```

```bash
.venv/bin/python test_lp.py
```

```bash
.venv/bin/python test_exact.py
```

The exact separation works to time limits, so a rerun will not reproduce the
cuts, or the last digits, exactly.

| File | Purpose |
|---|---|
| `palgraph.py` | The leftover-state graph: construction, backward reachability, strongly connected components |
| `count.py` | The exact count, about 5 seconds |
| `build.py` | A greedy sentence holding every includable phrase, checked by reversing the string, about 8 seconds |
| `check_sentence.py` | Checks any sentence in the `", "`-joined format and lists the includable phrases it lacks |
| `lp.py` | The linear program and its lower bound, about a minute with `--solver ipm` |
| `cuts.py` | Parity cuts: heuristic separation, and the loop that tightens the LP, about a minute a round |
| `exact.py` | Exact separation of the parity cuts by integer programming (`cuts.py --exact`) |
| `shortest.py` | LP (with cuts if given), rounding or branch and bound on the split words, joining, Hierholzer: a short complete sentence and its gap to the bound |
| `mip.py` | The tightened model as an integer program in HiGHS, started from a sentence |
| `export.py` | The tightened model, or its pinned restriction, as an MPS file |
| `scip_run.py` | Solves an exported model with SCIP, with the parity cuts as a separator plugin if asked |
| `ampl_export.py` | The tightened model as AMPL model, data and run files: 12 MB where the MPS file is 83 MB |
| `solfile.py` | Writes a sentence as a MIP start file, and turns a solver's solution file, or the AMPL run's output, back into a checked sentence |
| `test_bruteforce.py` | Compares the graph's walks with brute-force enumeration on small random dictionaries |
| `test_lp.py` | Compares the LP bound, the parity cuts and the sentence built with exact answers on small random dictionaries |
| `test_exact.py` | Compares the exact separation with brute force over every set of states |

To check Meertens's sentence, pass `longest_palindrome.txt` from
[rmeertens/palindromes](https://github.com/rmeertens/palindromes) to
`check_sentence.py`.

## What is not done

The last 17 letters are open.

- The parity cuts are exhausted for practical purposes: separated exactly,
  they stall at 814,391. A stronger bound needs another family, or cuts that
  combine the parity cuts with each other.
- Branch and bound on the whole model has not branched. Asking for whole
  numbers only on the arcs the LP leaves fractional keeps a valid bound and is
  much faster a round, but the fractions move elsewhere: it reached 814,394 in
  ten minutes, less than the all-integer root.
- The best sentence found here, 814,455, is at most a letter from optimal for
  its pinning, and Gurobi's 814,419 shows how much the pinning costs. Freeing
  more words is the way down, but the obvious wider neighbourhood, words whose arc touches a
  state of a split word, is half the dictionary.
