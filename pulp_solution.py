#!/usr/bin/env python

from pulp.constants import (
    LpMinimize,
)
from pulp import (
    LpProblem,
    LpVariable,
    LpInteger,
    lpSum,
    allcombinations,
)

nmice = 10
nbottles = 1000
npoisoned = 2

# Hamming weight - counts number of set bits in an integer
def popcount(x):
    x -= (x >> 1) & 0x5555555555555555
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
    return ((x * 0x0101010101010101) & 0xffffffffffffffff ) >> 56

labels = range(2**nmice)
# Only use labels with a Hamming weight <= 3
nonzeros = [i for i in labels if popcount(i) in [0,1,2,3]]

# Initialize problem
prob = LpProblem("drunk_mice", LpMinimize)

# This is the answer, the allocation of bottles to each label:
allocs = LpVariable.dicts(
    "alloc",
    nonzeros,
    0, nbottles,
    LpInteger,
)

# They must sum to nbottles:
prob += lpSum(allocs) == nbottles

# This is the number we are minimizing
worst_loss = LpVariable(
    "worst_loss",
    npoisoned,
    nbottles,
)

prob += worst_loss, "obj"

# This is the loss for each label
scenario_loss = LpVariable.dicts(
    "scenario_loss",
    labels,
    0, nbottles,
)

# The worst loss is worse than each scenario loss
for i in labels:
    prob += worst_loss >= scenario_loss[i]

# An scenario is generated by the overlapping of underlying
# labels. For each such set, if every underlying label has a non-zero
# allocation, then the allocations to all the labels must be added to
# the scenario loss

scenario2labels = {i : set() for i in labels}
for label_set in allcombinations(nonzeros, npoisoned):
    scenario = 0
    for label in label_set:
        scenario |= label
    for label in label_set:
        scenario2labels[scenario].add(label)

for scenario in labels:
    prob += scenario_loss[scenario] == lpSum(allocs[label0] for label0 in scenario2labels[scenario])

prob.writeLP("drunk_mice.lp")
prob.solve()

alloc = [0] * len(labels)
for l in nonzeros:
    alloc[l] = int(allocs[l].value())
print "".join([chr(i) for i in alloc]).encode("base64")

"""
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
