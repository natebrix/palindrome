"""Build the explainer page from the repository's data.

    npm install --prefix docs/explainer      once, for Temml (LaTeX to MathML)
    python3 docs/explainer/build.py [--fragment PATH]

Writes docs/explainer/explainer.html, a page that stands alone. With
--fragment it also writes the same page without its document skeleton, the
form a claude.ai artifact is published from.

Every example and count on the page is computed here from data/npdict.txt and
the sentence in sentences/, and the asserts fail if the text's claims about
them stop being true. The bounds in the chart are the README's.
"""
import argparse
import collections
import html
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
from palgraph import Graph, letters, load  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--dict", default=os.path.join(ROOT, "data", "npdict.txt"))
ap.add_argument("--sentence", default=os.path.join(ROOT, "sentences", "shortest_814403.txt"))
ap.add_argument("--fragment", default=None, metavar="PATH")
args = ap.parse_args()

D = load(args.dict)
g = Graph(D.keys())
inc = g.includable()
live = g.coreachable()
esc = html.escape


def n(v):
    return f"{v:,}"


# ---------------------------------------------------------------- formulas
TEX = {
    "obj": (r"\sum_{a \in A} \ell_{w(a)}\, x_a", True),
    "bal": (r"""\begin{aligned}
 & \sum_{a \in \delta^-(v)} x_a \;+\; [v = E] \sum_{f \in F} r_f \\
 & \qquad =\; \sum_{a \in \delta^+(v)} x_a \;+\; [v \in F]\, r_v
\end{aligned}""", True),
    "cov": (r"\sum_{a \,:\, w(a) = w} x_a \;\ge\; 1", True),
    "one": (r"\sum_{f \in F} r_f \;\le\; 1", True),
    "int": (r"x_a,\ r_f \in \{0, 1, 2, \dots\}", True),
    "half": (r"""\tfrac12 \sum_{w \in W} \; \sum_{a \,:\, w(a) = w} x_a
\;+\; \tfrac12 \Bigl( x\bigl(\delta^-(S)\bigr) - x\bigl(\delta^+(S)\bigr) \Bigr)
\;\ge\; \frac{|W|}{2}""", True),
    "cut": (r"""x\bigl(A_W \setminus \delta^+(S)\bigr) \;+\; x\bigl(\delta^-(S) \setminus A_W\bigr)
\;\ge\; \frac{|W| + 1}{2}""", True),
    "xa": (r"x_a", False), "rf": (r"r_f", False), "lw": (r"\ell_w", False),
    "wa": (r"w(a)", False), "S": (r"S", False), "W": (r"W", False), "AW": (r"A_W", False),
    "din": (r"\delta^-(S)", False), "dout": (r"\delta^+(S)", False),
    "TT": (r"T\,T", False), "T": (r"T", False), "E": (r"E", False),
}
M = json.loads(subprocess.run(["node", os.path.join(HERE, "tex.js")], cwd=HERE, input=json.dumps(TEX),
                              capture_output=True, text=True, check=True).stdout)


# ---------------------------------------------------------------- the Panama walk
def state_name(u):
    return "E" if u == 0 else f"({'L' if g.kind[u] == 1 else 'R'}, {g.text[u]})"


def step(u, key, side):
    for a in g.out(u):
        if g.words[g.lab[a]] == key and g.side(a) == side:
            return int(g.dst[a])
    raise AssertionError((state_name(u), key, side))


walk = [("aman", "L"), ("panama", "R"), ("aplan", "L"), ("acanal", "R")]
u = 0
states = [0]
for key, side in walk:
    u = step(u, key, side)
    states.append(u)
assert [state_name(s) for s in states] == ["E", "(L, aman)", "(R, pa)", "(L, lan)", "(R, aca)"]
assert g.accept[states[-1]] and not any(g.accept[s] for s in states[1:-1])
ways = [(len(list(g.out(s))), sum(1 for a in g.out(s) if live[g.dst[a]])) for s in states]
fan = sorted((D[g.words[g.lab[a]]], state_name(int(g.dst[a]))) for a in g.out(states[1]))
assert fan == [("Ama", "(L, n)"), ("Manama", "(R, ma)"), ("Nama", "E"), ("Panama", "(R, pa)")], fan
aba = step(0, "ababa", "L")
assert state_name(aba) == "(L, ababa)" and step(aba, "ababa", "R") == 0
assert "afghanistan" in D and g.wid["afghanistan"] not in inc
assert "m2m" in D and g.wid["m2m"] in inc


def cells(left, right, n_left_over, side):
    """Letter cells of a sentence in progress; the leftover is on `side`."""
    def cell(ch, cls):
        return f'<span class="c {cls}">{ch}</span>'
    if side == "L":
        lc = [cell(ch, "m") for ch in left[:len(left) - n_left_over]] + \
             [cell(ch, "o") for ch in left[len(left) - n_left_over:]]
        rc = [cell(ch, "m") for ch in right]
    else:
        lc = [cell(ch, "m") for ch in left]
        rc = [cell(ch, "o") for ch in right[:n_left_over]] + [cell(ch, "m") for ch in right[n_left_over:]]
    return (f'<span class="side">{"".join(lc)}</span><span class="gapdots" aria-hidden="true">···</span>'
            f'<span class="side">{"".join(rc)}</span>')


build_rows = [
    ("Start level.", "", "", 0, "L", "E"),
    ("<b>a man</b> goes on the left.", "aman", "", 4, "L", "(L, aman)"),
    ("<b>Panama</b> goes on the right. Its last four letters, read backwards, match "
     "<i>aman</i>; its first two are now unmatched.", "aman", "panama", 2, "R", "(R, pa)"),
    ("<b>a plan</b> goes on the left. Its <i>ap</i> answers the <i>pa</i>.", "amanaplan", "panama", 3, "L",
     "(L, lan)"),
    ("<b>a canal</b> goes on the right. The leftover <i>aca</i> reads the same both ways, "
     "so the sentence can stop here.", "amanaplan", "acanalpanama", 3, "R", "(R, aca)"),
]
build_html = []
for i, (text, left, right, k, side, st) in enumerate(build_rows):
    total, good = ways[i]
    last = i == len(build_rows) - 1
    note = (f"{n(total)} phrases fit next, {n(good)} of them without hitting a dead end"
            if not last else "accepting state")
    build_html.append(f"""
    <li class="bstep">
      <p class="btext">{text}</p>
      <div class="brow" aria-label="{esc(left)} ... {esc(right)}">{cells(left, right, k, side)}</div>
      <p class="bstate"><code class="state{' ok' if last else ''}">{st}</code><span>{note}</span></p>
    </li>""")

# ---------------------------------------------------------------- the sentence
text = open(args.sentence, encoding="utf-8").read().split(", ")
keys = [letters(p) for p in text]
s = "".join(keys)
assert s == s[::-1] and len(s) == 814403 and {g.wid[k] for k in keys} == inc
count = collections.Counter(keys)
once = sum(1 for v in count.values() if v == 1)
reused = len(count) - once
placements = len(keys)
c = len(s) // 2
strip = "".join(
    f'<span class="c {"mid" if i == c else "hl"}">{s[i]}</span>' for i in range(c - 12, c + 13))
top = [(D[k], v) for k, v in count.most_common(8)]
once_letters = sum(len(k) for k in count)
assert once_letters == 488137

# ---------------------------------------------------------------- the funnel
LB, UB = 814402, 814403
rows = [
    ("ub", "Meertens&rsquo;s 2018 sentence", 4836010),
    ("lb", "Each phrase once, nothing else", 488137),
    ("ub", "Greedy: one small loop per phrase", 1406982),
    ("lb", "Linear program", 813106),
    ("ub", "Rounding the linear program", 815407),
    ("lb", "Parity cuts, found by a heuristic", 814371),
    ("ub", "Rounding after the cuts", 814503),
    ("lb", "Parity cuts, separated exactly", 814391),
    ("ub", "Pin most words, search the rest", 814455),
    ("lb", "Solver&rsquo;s own cuts at the root", 814402),
    ("ub", "Gurobi&rsquo;s root heuristics, first job", 814419),
    ("ub", "Gurobi&rsquo;s root heuristics, second job", 814403),
]
DEC = 7.0


def pos(kind, v):
    d = (UB - v) if kind == "lb" else (v - LB)
    off = 3 + math.log10(d) / DEC * 45
    return 50 - off if kind == "lb" else 50 + off, d


funnel = []
for kind, label, v in rows:
    p, d = pos(kind, v)
    lo, hi = (p, 50 - 3) if kind == "lb" else (50 + 3, p)
    word = "below the best sentence" if kind == "lb" else "above the bound"
    funnel.append(f"""
    <div class="frow {kind}" role="row">
      <div class="flabel" role="cell"><span class="fname">{label}</span>
        <span class="fval">{n(v)}</span></div>
      <div class="ftrack" role="cell" title="{n(v)} letters: {n(d)} {word}">
        <span class="fstem" style="left:{lo:.2f}%;width:{hi - lo:.2f}%"></span>
        <span class="fdot" style="left:{p:.2f}%"></span>
      </div>
    </div>""")
ticks = []
for e in range(0, 8, 2):
    for kind in ("lb", "ub"):
        off = 3 + e / DEC * 45
        p = 50 - off if kind == "lb" else 50 + off
        lab = {0: "1", 2: "100", 4: "10k", 6: "1M"}[e]
        ticks.append(f'<span class="ftick" style="left:{p:.2f}%">{lab}</span>')

# ---------------------------------------------------------------- the page
share_inc = 100 * len(inc) / len(D)
page = f"""<title>Shortest Complete Palindrome</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Hanken+Grotesk:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;700&display=swap">
<style>
/* Layout: one reading column; figures set in the mono face, letters in cells. */
:root {{
  --bg: #fdfefe; --bg-subtle: #f5f1ea; --surface: #ffffff;
  --text: #2c2825; --text-2: #6b635a; --text-3: #9c948a;
  --accent: #b7714a; --accent-soft: #f0ddd1; --link: #9a5c38;
  --border: #e2ddd5; --border-strong: #c8c1b7;
  --lb: #257de0; --ub: #f47230; --ok: #236146;
  --font: "Hanken Grotesk", "HK Grotesk", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "JetBrains Mono", ui-monospace, "SFMono-Regular", Menlo, Consolas, monospace;
  --measure: 40rem;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg: #1b1917; --bg-subtle: #24211e; --surface: #211e1b;
    --text: #ece7e0; --text-2: #b3aa9f; --text-3: #857c72;
    --accent: #d69068; --accent-soft: #4a3327; --link: #e0a27d;
    --border: #37322d; --border-strong: #4d463f;
    --lb: #338ced; --ub: #e76624; --ok: #6fbf96; color-scheme: dark;
  }}
}}
:root[data-theme="dark"] {{
  --bg: #1b1917; --bg-subtle: #24211e; --surface: #211e1b;
  --text: #ece7e0; --text-2: #b3aa9f; --text-3: #857c72;
  --accent: #d69068; --accent-soft: #4a3327; --link: #e0a27d;
  --border: #37322d; --border-strong: #4d463f;
  --lb: #338ced; --ub: #e76624; --ok: #6fbf96; color-scheme: dark;
}}
body {{
  background: var(--bg); color: var(--text); font-family: var(--font);
  font-size: 1.0625rem; line-height: 1.68;
  padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px));
  padding-block: 2.5rem 4rem;
}}
main {{ max-width: var(--measure); margin-inline: auto; display: flex; flex-direction: column; gap: 1.15rem; min-width: 0; }}
main > * {{ margin: 0; min-width: 0; }}
h1 {{ font-size: clamp(1.9rem, 5.2vw, 2.7rem); line-height: 1.12; font-weight: 700; letter-spacing: -0.015em; text-wrap: balance; }}
h2 {{ font-size: 1.42rem; line-height: 1.25; font-weight: 700; margin-top: 1.9rem; text-wrap: balance; }}
h3 {{ font-size: 1.05rem; font-weight: 700; margin-top: 0.6rem; }}
p, li {{ text-wrap: pretty; }}
a {{ color: var(--link); text-underline-offset: 2px; }}
a:focus-visible, summary:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
code, .mono {{ font-family: var(--mono); font-size: 0.9em; }}
.kicker {{ font-size: 0.8rem; letter-spacing: 0.09em; text-transform: uppercase; color: var(--text-2); font-weight: 600; }}
.lede {{ font-size: 1.2rem; line-height: 1.55; color: var(--text); }}
.num {{ font-variant-numeric: tabular-nums; }}
ul.plain {{ padding-left: 1.2rem; display: flex; flex-direction: column; gap: 0.35rem; }}

/* letter cells */
.c {{ display: inline-block; width: 1.32ch; text-align: center; font-family: var(--mono); }}
.c.m {{ color: var(--text-3); }}
.c.o {{ color: var(--text); background: var(--accent-soft); font-weight: 700; }}
.c.hl {{ color: var(--text); }}
.c.mid {{ color: var(--bg); background: var(--accent); font-weight: 700; border-radius: 2px; }}

.hero-strip {{ border-block: 1px solid var(--border); padding-block: 1rem; display: flex; flex-direction: column; gap: 0.35rem; align-items: center; }}
.hero-strip .letters {{ font-size: clamp(0.82rem, 3.3vw, 1.25rem); white-space: nowrap; }}
.hero-strip .words {{ font-size: 0.86rem; color: var(--text-2); text-align: center; }}

figure {{ display: flex; flex-direction: column; gap: 0.7rem; padding: 1.1rem 0; }}
figcaption {{ font-size: 0.88rem; line-height: 1.5; color: var(--text-2); }}
figcaption b {{ color: var(--text); font-weight: 600; }}
.fig-title {{ font-size: 0.8rem; letter-spacing: 0.08em; text-transform: uppercase; color: var(--text-2); font-weight: 600; }}

/* building a palindrome */
ol.build {{ list-style: none; padding: 0; margin: 0; display: flex; flex-direction: column; gap: 0; counter-reset: b -1; }}
.bstep {{ counter-increment: b; display: grid; grid-template-columns: 1.6rem minmax(0, 1fr); column-gap: 0.6rem; row-gap: 0.3rem; padding: 0.75rem 0; border-top: 1px solid var(--border); }}
.bstep::before {{ content: counter(b); font-family: var(--mono); font-size: 0.85rem; color: var(--text-3); grid-row: 1 / span 3; padding-top: 0.1rem; }}
.bstep > * {{ grid-column: 2; margin: 0; min-width: 0; }}
.btext {{ font-size: 0.95rem; line-height: 1.45; }}
.brow {{ display: flex; justify-content: space-between; align-items: baseline; font-size: clamp(0.78rem, 3.4vw, 1.05rem); background: var(--bg-subtle); padding: 0.45rem 0.6rem; border-radius: 4px; min-height: 1.6em; white-space: nowrap; }}
.gapdots {{ color: var(--text-3); letter-spacing: 0.2em; }}
.bstate {{ display: flex; flex-wrap: wrap; gap: 0.2rem 0.7rem; align-items: baseline; font-size: 0.84rem; color: var(--text-2); }}
code.state {{ color: var(--text); border: 1px solid var(--border-strong); border-radius: 999px; padding: 0.02rem 0.55rem; white-space: nowrap; }}
code.state.ok {{ border-color: var(--ok); color: var(--ok); border-width: 2px; }}

svg.dia {{ max-width: 100%; height: auto; color: var(--text); font-family: var(--mono); }}
.dia .node {{ fill: var(--surface); stroke: currentColor; stroke-width: 1.2; }}
.dia .dim {{ fill: var(--text-2); }}
.dia .flow {{ fill: var(--accent); font-weight: 700; }}

/* proportion bars */
.split {{ display: flex; height: 2.1rem; gap: 2px; }}
.split > span {{ display: flex; align-items: center; padding-inline: 0.6rem; font-size: 0.82rem; font-weight: 600; white-space: nowrap; overflow: hidden; min-width: 0; }}
.split .yes {{ background: var(--accent); color: var(--bg); border-radius: 4px 0 0 4px; }}
.split .no {{ background: var(--bg-subtle); color: var(--text-2); border: 1px solid var(--border); border-radius: 0 4px 4px 0; }}
.split-legend {{ display: flex; justify-content: space-between; gap: 1rem; font-size: 0.85rem; color: var(--text-2); flex-wrap: wrap; }}
.split-legend b {{ color: var(--text); }}

/* math */
.math {{ overflow-x: auto; padding: 0.9rem 1rem; background: var(--bg-subtle); border-radius: 6px; }}
math {{ font-size: 1.02em; }}
math[display="block"] {{ font-size: 1em; margin: 0; }}
.model {{ background: var(--bg-subtle); border-radius: 6px; padding: 0.5rem 1rem; display: flex; flex-direction: column; }}
.mrow {{ display: grid; grid-template-columns: 5.4rem minmax(0, 1fr) minmax(0, 11rem); column-gap: 0.8rem; align-items: center; padding: 0.55rem 0; }}
.mrow + .mrow {{ border-top: 1px solid var(--border); }}
.mlab {{ font-size: 0.85rem; color: var(--text-2); }}
.mf {{ overflow-x: auto; min-width: 0; }}
.mf math[display="block"] {{ text-align: left; }}
.mnote {{ font-size: 0.82rem; color: var(--text-2); line-height: 1.35; }}
@media (max-width: 34rem) {{
  .mrow {{ grid-template-columns: minmax(0, 1fr); row-gap: 0.25rem; }}
  .mlab:empty {{ display: none; }}
}}

/* funnel */
.funnel {{ display: flex; flex-direction: column; }}
.flegend {{ display: flex; flex-wrap: wrap; gap: 0.3rem 1.2rem; font-size: 0.85rem; color: var(--text-2); }}
.flegend i {{ display: inline-block; width: 0.7rem; height: 0.7rem; border-radius: 50%; margin-right: 0.4rem; vertical-align: -0.03rem; }}
.frow, .faxis {{ display: grid; grid-template-columns: minmax(0, 15.5rem) minmax(0, 1fr); column-gap: 0.9rem; align-items: center; }}
.frow {{ padding: 0.32rem 0; border-top: 1px solid var(--border); }}
.flabel {{ display: flex; justify-content: space-between; gap: 0.6rem; align-items: baseline; font-size: 0.86rem; line-height: 1.3; }}
.fname {{ color: var(--text); min-width: 0; }}
.fval {{ font-family: var(--mono); font-size: 0.8rem; font-variant-numeric: tabular-nums; color: var(--text-2); white-space: nowrap; }}
.ftrack {{ position: relative; height: 1.3rem; background:
  linear-gradient(var(--border-strong), var(--border-strong)) 50% 0 / 1px 100% no-repeat; }}
.fstem {{ position: absolute; top: calc(50% - 1px); height: 2px; }}
.fdot {{ position: absolute; top: 50%; width: 10px; height: 10px; border-radius: 50%; transform: translate(-50%, -50%); box-shadow: 0 0 0 2px var(--bg); }}
.lb .fstem, .lb .fdot, .flegend .lb {{ background: var(--lb); }}
.ub .fstem, .ub .fdot, .flegend .ub {{ background: var(--ub); }}
.faxis {{ border-top: 1px solid var(--border-strong); padding-top: 0.25rem; }}
.faxis .ftrack {{ background: none; height: 2.3rem; }}
.ftick {{ position: absolute; top: 0; transform: translateX(-50%); font-family: var(--mono); font-size: 0.7rem; color: var(--text-2); }}
.fside {{ position: absolute; top: 1.15rem; font-size: 0.72rem; color: var(--text-2); white-space: nowrap; }}
@media (max-width: 34rem) {{
  .frow, .faxis {{ grid-template-columns: minmax(0, 1fr); row-gap: 0.15rem; }}
  .faxis .flabel {{ display: none; }}
}}

/* usage */
.stack {{ display: grid; grid-template-columns: minmax(0, 8.5rem) minmax(0, 1fr); gap: 0.5rem 0.8rem; align-items: center; font-size: 0.86rem; }}
.stack .bar {{ display: flex; height: 1.7rem; gap: 2px; }}
.stack .bar span {{ display: flex; align-items: center; padding-inline: 0.5rem; font-size: 0.78rem; font-weight: 600; white-space: nowrap; overflow: hidden; min-width: 0; }}
.stack .a {{ background: var(--bg-subtle); border: 1px solid var(--border); color: var(--text-2); border-radius: 4px 0 0 4px; }}
.stack .b {{ background: var(--accent); color: var(--bg); border-radius: 0 4px 4px 0; }}
.chips {{ display: flex; flex-wrap: wrap; gap: 0.4rem; padding: 0; margin: 0; list-style: none; }}
.chips li {{ font-family: var(--mono); font-size: 0.82rem; border: 1px solid var(--border-strong); border-radius: 4px; padding: 0.1rem 0.5rem; }}
.chips li span {{ color: var(--text-2); margin-left: 0.4rem; }}
.ends {{ font-family: var(--mono); font-size: 0.84rem; line-height: 1.6; background: var(--bg-subtle); border-radius: 6px; padding: 0.8rem 1rem; overflow-wrap: anywhere; }}
.ends .el {{ color: var(--text-3); }}

table.t {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
.t th, .t td {{ text-align: left; padding: 0.4rem 0.6rem 0.4rem 0; border-top: 1px solid var(--border); vertical-align: top; }}
.t th {{ font-weight: 600; color: var(--text-2); font-size: 0.78rem; letter-spacing: 0.06em; text-transform: uppercase; border-top: 0; }}
.t td.r, .t th.r {{ text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; padding-right: 0; }}
.tablewrap {{ overflow-x: auto; }}
.note {{ border-left: 3px solid var(--accent); padding: 0.2rem 0 0.2rem 1rem; color: var(--text); font-size: 0.97rem; }}
.sources {{ font-size: 0.92rem; }}
footer {{ margin-top: 2rem; padding-top: 1rem; border-top: 1px solid var(--border); font-size: 0.85rem; color: var(--text-2); }}
</style>

<main>
<p class="kicker">Palindromes &middot; integer programming</p>
<h1>The shortest palindrome that says everything it can</h1>
<p class="lede">Take Peter Norvig&rsquo;s dictionary of noun phrases. Exactly <span class="num">{n(len(inc))}</span> of its phrases can appear in a palindromic sentence at all. The shortest sentence that contains every one of them has <span class="num">{n(LB)}</span> or <span class="num">{n(UB)}</span> letters, and a sentence of {n(UB)} is in hand.</p>

<div class="hero-strip" role="img" aria-label="The 25 letters at the centre of the sentence: comanecimmer, e, remmicenamoc">
  <div class="letters">{strip}</div>
  <div class="words">the middle of the sentence: &hellip; Comaneci, Mme, <b>reremmice</b>, Nam, OCTU &hellip;</div>
</div>

<h2>The puzzle</h2>
<p>&ldquo;A man, a plan, a canal: Panama&rdquo; reads the same in both directions once you drop the spaces and punctuation. In 2002 <a href="https://norvig.com/palindrome.html">Peter Norvig</a> asked a program to stretch it: keep the opening and the ending, and insert as many different noun phrases as possible in between, each used once. His dictionary, <code>npdict.txt</code>, has 126,342 entries drawn from the Moby word list. His first palindrome ran to 15,319 words, and the 2016 version reaches 21,012 words and 16,111 phrases.</p>
<p>In 2018 <a href="https://meertens.dev/blog/worlds-longest-palindrome/">Roland Meertens</a> changed the rules. Let phrases repeat, and ask how many <em>different</em> phrases one sentence can hold. He called such a sentence complete. He built a small palindromic sub-sentence around each dictionary phrase by depth-first search and joined them, reaching 64,810 distinct phrases in 954,401 placements, and wrote that it is safe to assume no more can be included.</p>
<p>That leaves two questions with exact answers:</p>
<ul class="plain">
  <li><b>How many phrases can be included?</b> Exactly {n(len(inc))} of the {n(len(D))} distinct ones. Meertens&rsquo;s sentence has all but one of them, &ldquo;M2M&rdquo;.</li>
  <li><b>How short can a complete sentence be?</b> Meertens&rsquo;s has 4,836,010 letters. The answer is {n(LB)} or {n(UB)}, about a sixth of that.</li>
</ul>
<p>A phrase counts by its letters and digits only, lowercased, as in Norvig&rsquo;s and Meertens&rsquo;s code: &ldquo;a man&rdquo; is <code>aman</code>. Lengths below are in those characters.</p>

<h2>Build it from the outside in</h2>
<p>A palindrome is easiest to grow from both ends at once. Put a phrase at the left end. Its letters now need a mirror image at the right end, so put a phrase there that supplies it. That phrase usually overshoots and leaves letters of its own unmatched, which the left side must answer next. Always add to the side that is behind.</p>
<figure>
  <div class="fig-title">Four phrases, five states</div>
  <ol class="build">{"".join(build_html)}</ol>
  <figcaption><b>Highlighted letters are the leftover:</b> the part of one side that the other side has not yet mirrored. Grey letters are already matched. Everything that matters about the sentence so far is which side is ahead and by which letters. The counts are from the full dictionary.</figcaption>
</figure>
<p>The leftover is the whole state of the construction. There are three kinds:</p>
<ul class="plain">
  <li><code>E</code>: the two sides are level.</li>
  <li><code>(L, s)</code>: the left is ahead by <code>s</code>, the end of the last phrase placed there.</li>
  <li><code>(R, p)</code>: the right is ahead by <code>p</code>, the start of the last phrase placed there.</li>
</ul>
<p>A sentence may stop whenever the leftover is itself a palindrome, because the leftover is then the middle of the sentence. Such a state is called accepting.</p>

<h2>A graph of leftovers</h2>
<p>Make every possible leftover a node and draw an arc for every way of placing one phrase. Each arc carries the phrase it places.</p>
<figure>
  <div class="fig-title">Every arc out of one state</div>
  <svg class="dia" viewBox="0 0 380 212" width="460" role="img" aria-label="From state (L, aman) four phrases can be placed on the right: Ama leads to (L, n), Nama leads to E, Manama leads to (R, ma), Panama leads to (R, pa).">
    <defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>
    <rect class="node" x="8" y="88" width="104" height="34" rx="17"/>
    <text x="60" y="110" text-anchor="middle" font-size="13" fill="currentColor">(L, aman)</text>
    <g fill="none" stroke="currentColor" stroke-width="1.2" marker-end="url(#ah)">
      <path d="M112,97 C170,80 190,32 268,26"/>
      <path d="M112,102 C180,98 200,82 268,80"/>
      <path d="M112,108 C180,112 200,128 268,131"/>
      <path d="M112,113 C170,130 190,180 268,185"/>
    </g>
    <g font-size="12" fill="currentColor">
      <text x="176" y="40" text-anchor="middle">Ama</text>
      <text x="196" y="82" text-anchor="middle">Nama</text>
      <text x="196" y="138" text-anchor="middle">Manama</text>
      <text x="176" y="178" text-anchor="middle">Panama</text>
    </g>
    <rect class="node" x="272" y="10" width="96" height="32" rx="16"/>
    <rect class="node" x="272" y="64" width="96" height="32" rx="16"/>
    <rect class="node" x="272" y="116" width="96" height="32" rx="16"/>
    <rect class="node" x="272" y="170" width="96" height="32" rx="16"/>
    <g font-size="13" fill="currentColor" text-anchor="middle">
      <text x="320" y="31">(L, n)</text><text x="320" y="85">E</text>
      <text x="320" y="137">(R, ma)</text><text x="320" y="191">(R, pa)</text>
    </g>
  </svg>
  <figcaption><b>After &ldquo;a man&rdquo; on the left, only four phrases in the dictionary can go on the right.</b> Each either ends in <i>nama</i> or is the mirror image of the first letters of <i>aman</i>. &ldquo;Nama&rdquo; matches exactly and returns to the level state.</figcaption>
</figure>
<p>Every palindromic sentence is exactly one walk in this graph from <code>E</code> to an accepting state, and every such walk is a sentence. For Norvig&rsquo;s dictionary the graph has {n(int(g.n))} states and {n(len(g.src))} arcs.</p>
<p>The first question is now a question about reachability. A phrase can be included when it labels an arc that some walk from <code>E</code> to an accepting state uses. Of the states, {n(sum(live))} can still reach an accepting one, and the arcs among them carry {n(len(inc))} different phrases.</p>
<figure>
  <div class="split" role="img" aria-label="{n(len(inc))} phrases can be included, {n(len(D) - len(inc))} cannot.">
    <span class="yes" style="flex:{len(inc)}">{n(len(inc))} can appear</span>
    <span class="no" style="flex:{len(D) - len(inc)}">{n(len(D) - len(inc))} never can</span>
  </div>
  <figcaption><b>{share_inc:.1f}% of the dictionary&rsquo;s {n(len(D))} distinct phrases can appear in some palindromic sentence.</b> The rest, &ldquo;Afghanistan&rdquo; among them, contain a run of letters that no sequence of dictionary phrases spells backwards.</figcaption>
</figure>
<p>One sentence can hold all of them at once. If {M["T"]} is a palindromic sentence then so is {M["TT"]}, the sentence said twice, and {M["TT"]} is level in the middle. So every includable phrase lies on a closed walk that leaves <code>E</code> and returns to it, and closed walks through <code>E</code> can be nested one inside another. The computation confirms it directly: the {n(sum(live))} useful states form a single strongly connected component.</p>

<h2>How short? A postman problem</h2>
<p>A complete sentence is a walk that uses every includable phrase at least once, and its length is the sum of the letters on its arcs. Finding the shortest one is a directed rural postman problem: cover a required set of streets with the shortest closed route. Norvig&rsquo;s recent <a href="https://github.com/norvig/pytudes/blob/main/ipynb/Portmantout-optimal.ipynb">Portmantout notebook</a> solves a problem of the same family, the shortest string that contains every word of a list, and there the linear program lands on a whole-number answer by itself. Here it does not, for one reason: a phrase labels many arcs, one for each state it can be played from, and the walk needs only one of them.</p>
<p>Count how often the walk uses each arc. Let {M["xa"]} be that count for arc <i>a</i>, let {M["wa"]} be the phrase on the arc and {M["lw"]} its length. A walk that ends in an accepting state <i>f</i> is closed up by a free return arc to <code>E</code>, counted by {M["rf"]}; <i>F</i> is the set of accepting states.</p>
<div class="model" role="group" aria-label="The integer program">
  <div class="mrow"><span class="mlab">minimize</span><div class="mf">{M["obj"]}</div><span class="mnote">the letters of the sentence</span></div>
  <div class="mrow"><span class="mlab">subject to</span><div class="mf">{M["bal"]}</div><span class="mnote">for every state <i>v</i>: in equals out</span></div>
  <div class="mrow"><span class="mlab"></span><div class="mf">{M["cov"]}</div><span class="mnote">for every includable phrase <i>w</i></span></div>
  <div class="mrow"><span class="mlab"></span><div class="mf">{M["one"]}</div><span class="mnote">one sentence</span></div>
  <div class="mrow"><span class="mlab"></span><div class="mf">{M["int"]}</div><span class="mnote">whole numbers</span></div>
</div>
<p>The first family of constraints says the walk leaves every state as often as it enters it. The second says every phrase is used. The third allows one return to <code>E</code>, so the flow describes a single sentence. One requirement is left out: that the arcs in use hang together as one walk. It turns out not to matter here, since every flow found was already in one piece, and a flow in pieces can be joined through <code>E</code> at a small cost.</p>
<p>For this dictionary the model has 135,455 states, 426,676 arcs, 64,811 phrases and 903 return arcs. A whole-number flow is turned back into a sentence by Hierholzer&rsquo;s algorithm for Eulerian circuits, and the sentence is then checked from scratch: reverse the string, look up every phrase, count the distinct ones.</p>

<h2>Half a palindrome, and the cut that forbids it</h2>
<p>Drop the whole-number requirement and the model is a linear program. Its optimum is 813,105.25, so no complete sentence has fewer than 813,106 letters. The best sentence built by rounding that solution has 815,407. The gap comes from flows of one half.</p>
<figure>
  <div class="fig-title">How the linear program cheats</div>
  <svg class="dia" viewBox="0 0 380 150" width="460" role="img" aria-label="Two states, E and (L, ababa), joined by two arcs that both place the phrase ababa. The linear program sends half a unit around the loop.">
    <defs><marker id="ah2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="currentColor"/></marker></defs>
    <rect class="node" x="14" y="58" width="64" height="34" rx="17"/>
    <text x="46" y="80" text-anchor="middle" font-size="13" fill="currentColor">E</text>
    <rect class="node" x="256" y="58" width="110" height="34" rx="17"/>
    <text x="311" y="80" text-anchor="middle" font-size="13" fill="currentColor">(L, ababa)</text>
    <g fill="none" stroke="currentColor" stroke-width="1.2" marker-end="url(#ah2)">
      <path d="M70,58 C120,18 220,18 268,56"/>
      <path d="M268,94 C220,132 120,132 70,94"/>
    </g>
    <g font-size="12" fill="currentColor" text-anchor="middle">
      <text x="169" y="20">ababa, on the left</text>
      <text x="169" y="146">ababa, on the right</text>
    </g>
    <g font-size="14" text-anchor="middle" class="flow">
      <text x="169" y="50">1/2</text><text x="169" y="116">1/2</text>
    </g>
  </svg>
  <figcaption><b>Half a unit around a loop of two arcs covers the phrase once, for 5 letters.</b> A real sentence has to go around the loop a whole number of times, which costs 10, or reach &ldquo;ababa&rdquo; some other way. In the optimal solution 631 phrases are split like this, mostly around longer loops through several phrases.</figcaption>
</figure>
<p>The remedy is a family of inequalities that every whole-number flow satisfies and the half-unit loops violate. Take any set {M["S"]} of states and any set {M["W"]} of phrases with an odd number of members. Write {M["AW"]} for the arcs labelled by a phrase in {M["W"]}, and {M["din"]} and {M["dout"]} for the arcs entering and leaving {M["S"]}. Add the covering constraints of {M["W"]} to the balance of flow across the boundary of {M["S"]}, and halve:</p>
<div class="math" role="group" aria-label="Half the sum of the covering constraints and the flow balance">{M["half"]}</div>
<p>Arcs of {M["W"]} that enter {M["S"]} now have coefficient 1, those that leave have 0, and the other arcs of {M["W"]} and the other entering arcs have one half. Round the halves up, drop the negative terms, and the left side is a whole number for any whole-number flow. The right side, {M["W"]} being odd, is not, so it rounds up too:</p>
<div class="math" role="group" aria-label="The parity cut">{M["cut"]}</div>
<p>In words: if every phrase of an odd set were covered exactly once by flow crossing the boundary of {M["S"]}, and nothing else crossed, an odd amount would cross, although as much must enter as leave. For &ldquo;ababa&rdquo;, {M["S"]} is the one state <code>(L, ababa)</code> and {M["W"]} is the phrase alone, and the cut demands a whole unit where the linear program had a half.</p>
<p>Finding violated cuts is its own search problem. A heuristic based on linear algebra over the two-element field found 1,460 of them and raised the bound to 814,371. Separating the same family exactly, as a small integer program, added 431 more and reached 814,391. There the family runs out: what stays fractional is one tangle of about 460 short bridging phrases, with flows at thirds and odder values.</p>

<h2>Closing the gap</h2>
<p>Two numbers move toward each other. A lower bound says no complete sentence can be shorter. An upper bound is the length of a sentence actually built and checked. The answer lies between the best of each.</p>
<figure>
  <div class="fig-title">Bounds on the shortest complete sentence, in order of discovery</div>
  <div class="flegend"><span><i class="lb"></i>lower bound: nothing shorter exists</span><span><i class="ub"></i>upper bound: a checked sentence</span></div>
  <div class="funnel" role="table" aria-label="Lower and upper bounds">
    {"".join(funnel)}
    <div class="faxis"><div class="flabel"></div>
      <div class="ftrack">{"".join(ticks)}
        <span class="fside" style="left:0">&larr; letters below {n(UB)}</span>
        <span class="fside" style="right:0">letters above {n(LB)} &rarr;</span>
      </div></div>
  </div>
  <figcaption><b>The two sides now stand one letter apart.</b> Distances are on a logarithmic scale, measured from the best bound on the other side: lower bounds to the left, sentences to the right. Each row is a method, with the bound it gave.</figcaption>
</figure>
<p>The upper side moved in three ways. Rounding fixes one arc per phrase a few hundred phrases at a time and re-solves; once every phrase has its arc, what remains is an ordinary network flow, whose solutions are whole numbers. Starting that from the cut-tightened program instead of the plain one took a sentence from 815,407 letters to 814,503. Pinning the 64,479 phrases that the program already places on a single arc and running branch and bound on the other 332 gave 814,455.</p>
<p>The last 52 letters came from a commercial solver. Gurobi, run on the public <a href="https://neos-server.org">NEOS server</a>, improved the sentence at the root of its search to 814,419 and then to 814,403 by its own heuristics. It could not search further in the 3&nbsp;GB of memory NEOS allows. Its cutting planes also lift the lower bound from 814,391 to 814,401.5, which rounds up to 814,402, the same figure the open-source solver HiGHS reaches.</p>
<p class="note">How firm are the two ends? The sentence of {n(UB)} letters is checked independently of any solver. The lower bounds of 813,106 and 814,391 rest on a linear program plus cuts with a short proof. The final step to {n(LB)} rests on cuts generated inside two solvers that agree, and has not been verified separately.</p>

<h2>The sentence</h2>
<p>The best sentence has {n(placements)} phrase placements for its {n(len(count))} distinct phrases. Said once each, the phrases would take {n(once_letters)} letters. The other {n(UB - once_letters)} letters are glue: short phrases repeated to make the mirror images line up.</p>
<figure>
  <div class="fig-title">Most phrases appear once; a few do the joining</div>
  <div class="stack">
    <span>Distinct phrases</span>
    <div class="bar" role="img" aria-label="{n(once)} phrases used once, {n(reused)} used more than once">
      <span class="a" style="flex:{once}">{n(once)} used once</span><span class="b" style="flex:{reused}" title="{n(reused)} used more than once"></span></div>
    <span>Placements</span>
    <div class="bar" role="img" aria-label="{n(once)} placements of phrases used once, {n(placements - once)} placements of reused phrases">
      <span class="a" style="flex:{once}">{n(once)}</span><span class="b" style="flex:{placements - once}">{n(placements - once)} by {n(reused)} reused phrases</span></div>
  </div>
  <ul class="chips" aria-label="The most reused phrases">
    {"".join(f"<li>{esc(p)}<span>&times;{n(v)}</span></li>" for p, v in top)}
  </ul>
  <figcaption><b>{100 * once / len(count):.0f}% of the phrases are used exactly once, and {n(reused)} reused phrases supply {100 * (placements - once) / placements:.0f}% of the placements.</b> The chips are the eight most reused phrases with their counts.</figcaption>
</figure>
<p>It begins and ends like this:</p>
<div class="ends">{esc(", ".join(text[:9]))}, <span class="el">&hellip; {n(placements - 18)} more &hellip;</span> {esc(", ".join(text[-9:]))}</div>
<p>and its single middle letter is the second <i>e</i> of &ldquo;reremmice&rdquo;, a kind of bat.</p>

<h2>What is settled and what is open</h2>
<div class="tablewrap"><table class="t">
  <thead><tr><th>Question</th><th>Answer</th><th>Status</th></tr></thead>
  <tbody>
    <tr><td>Phrases that can appear in a palindromic sentence</td><td class="num">{n(len(inc))} of {n(len(D))}</td><td>Exact</td></tr>
    <tr><td>Phrases one sentence can hold</td><td class="num">{n(len(inc))}</td><td>Exact, with a sentence as witness</td></tr>
    <tr><td>Letters in the shortest complete sentence</td><td class="num">{n(LB)} or {n(UB)}</td><td>One letter open</td></tr>
  </tbody>
</table></div>
<p>Deciding the last letter means either finding a complete sentence of {n(LB)} letters or proving none exists. Both need a solver that can search the full model past its first node, which takes more memory than the free server offers and more speed than the open-source solvers have shown on this problem.</p>

<h2>Sources</h2>
<ul class="plain sources">
  <li>Peter Norvig, <a href="https://norvig.com/palindrome.html">World&rsquo;s Longest Palindrome Sentence</a>, 2002 to 2016, and the dictionary <code>npdict.txt</code>.</li>
  <li>Roland Meertens, <a href="https://meertens.dev/blog/worlds-longest-palindrome/">World&rsquo;s longest palindrome?</a>, 2 August 2018, and <a href="https://github.com/rmeertens/palindromes">his code and sentence</a>.</li>
  <li>The <a href="https://github.com/norvig/pytudes/blob/main/ipynb/Portmantout-optimal.ipynb">Portmantout-optimal notebook</a> in Norvig&rsquo;s pytudes, for the postman formulation and its whole-number linear program.</li>
  <li>Code, checks and the sentence: <a href="https://github.com/natebrix/palindrome">github.com/natebrix/palindrome</a>.</li>
</ul>
<footer>Solvers: HiGHS 1.16, SCIP 10 and Gurobi 13 (through NEOS). Figures as of 28 November 2026.</footer>
</main>
"""
SKELETON = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>body {{ margin: 0 }} img {{ max-width: 100% }}</style>
</head>
<body>
{page}</body>
</html>
"""
out = os.path.join(HERE, "explainer.html")
open(out, "w", encoding="utf-8").write(SKELETON.format(page=page))
print(f"wrote {out}")
if args.fragment:
    open(args.fragment, "w", encoding="utf-8").write(page)
    print(f"wrote {args.fragment}")
