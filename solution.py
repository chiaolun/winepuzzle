#!/usr/bin/env python

from itertools import combinations
import numpy as np

nmice = 10
nbottles = 1000
nlabels = 2**nmice
npoisoned = 2
labels = range(nlabels)

# The resulting losses for every possible outcome would be
def worst_loss(alloc):
    nonzeros = [l for l in labels if alloc[l] > 0]
    loss_sets = [set() for _ in range(nlabels)]
    for N in range(1, npoisoned + 1):
        for label_set in combinations(nonzeros, N):
            scenario = 0
            for label in label_set:
                scenario |= label
            for label in label_set:
                loss_sets[scenario].add(label)
    losses = [sum([alloc[i] for i in loss_set]) for loss_set in loss_sets]
    return max([(x, i) for i, x in enumerate(losses)])

alloc = """
Kw0NBQ0EBQUMBQUFBQUGAAwFBQUFBQUABQUFAAUAAAAMBgUFBQYFAAUEBQAFAAAABQUFAAUAAAAF
AAAAAAAAAA0FBQUFBQUABQYFAAUAAAAFBQUABQAAAAUAAAAAAAAABQUFAAUAAAAFAAAAAAAAAAUA
AAAAAAAAAAAAAAAAAAAMBQUFBQUFAAUFBAAFAAAABAUGAAUAAAAFAAAAAAAAAAUFBQAFAAAABgAA
AAAAAAAFAAAAAAAAAAAAAAAAAAAABQUGAAUAAAAFAAAAAAAAAAUAAAAAAAAAAAAAAAAAAAAFAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA0FBQUFBQUABQUFAAUAAAAGBQUABQAAAAUAAAAA
AAAABQUFAAUAAAAGAAAAAAAAAAUAAAAAAAAAAAAAAAAAAAAFBQUABQAAAAUAAAAAAAAABQAAAAAA
AAAAAAAAAAAAAAUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABQUFAAUAAAAFAAAAAAAA
AAUAAAAAAAAAAAAAAAAAAAAFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAM
BQUFBQUFAAUFBQAFAAAABQUFAAUAAAAGAAAAAAAAAAUFBQAFAAAABQAAAAAAAAAFAAAAAAAAAAAA
AAAAAAAABQUFAAUAAAAFAAAAAAAAAAUAAAAAAAAAAAAAAAAAAAAFAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAYFBQAFAAAABQAAAAAAAAAFAAAAAAAAAAAAAAAAAAAABQAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABQUFAAYAAAAEAAAAAAAAAAUAAAAAAAAAAAAAAAAA
AAAFAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAFAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==
"""

alloc = map(ord,alloc.decode("base64"))
alloc = np.array(alloc)

assert sum(alloc) == 1000
assert min(alloc) >= 0
assert len(alloc) == 1024
loss, pattern = worst_loss(alloc)
print "Worst loss is", loss, "when pattern is", "{:010b}".format(73)
# Worst loss is 102 when pattern is 0001001001
