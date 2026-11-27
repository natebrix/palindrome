"""Check a published palindromic sentence against the dictionary.

    python3 check_sentence.py SENTENCE.txt [data/npdict.txt]

The sentence is phrases joined by ", " (the format of Roland Meertens's
longest_palindrome.txt). Reports whether it is a palindrome, how many distinct
phrases it holds, and which includable phrases it lacks.
"""
import sys

from palgraph import Graph, letters, load

sentence_path = sys.argv[1]
dict_path = sys.argv[2] if len(sys.argv) > 2 else "data/npdict.txt"

D = load(dict_path)
g = Graph(D.keys())
U = {g.words[i] for i in g.includable()}

text = open(sentence_path, encoding="utf-8").read().strip()
s = letters(text)
print(f"letters: {len(s):,}   palindrome: {s == s[::-1]}")

# Split into phrases. A few dictionary phrases contain ", " themselves, so
# rejoin neighbours when a piece is not a phrase on its own.
pieces = [letters(p) for p in text.split(", ")]
keys = []
i = 0
unknown = 0
while i < len(pieces):
    k = pieces[i]
    j = i
    while k not in D and j + 1 < len(pieces) and j - i < 3:
        j += 1
        k += pieces[j]
    if k in D:
        keys.append(k)
        i = j + 1
    else:
        unknown += 1
        keys.append(pieces[i])
        i += 1
print(f"phrases: {len(keys):,}   not in the dictionary: {unknown:,}")
have = set(keys) & set(D)
print(f"distinct phrases: {len(have):,} of {len(U):,} includable")
extra = sorted(have - U)
if extra:
    print("in the sentence but not includable (a bug in the graph if any):", extra[:20])
missing = sorted(U - have)
print(f"includable phrases it lacks: {len(missing):,}")
for k in missing[:50]:
    print(f"  {D[k]}  ({k})")
