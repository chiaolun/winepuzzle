#!/usr/bin/env python

from itertools import combinations

nmice = 10
nbottles = 1000
npoisoned = 2
nlabels = 2**nmice

tothrow = [set() for i in range(nlabels)]
for n in range(npoisoned):
    for labelset in combinations(range(nlabels), n+1):
        k = 0
        for label in labelset:
            k |= label
        tothrow[k].add(labelset)

# The resulting losses for every possible outcome would be

def worst_loss(alloc):
    worst = 0, 0
    for j, labelsets in enumerate(tothrow):
        j_loss = 0
        thrown = set()
        for labelset in labelsets:
            if any(alloc[i] == 0 for i in labelset):
                continue
            for i in labelset:
                if i not in thrown:
                    j_loss += alloc[i]
                    thrown.add(i)
        if j_loss > worst[0]:
            worst = j_loss, j
    return worst

alloc = map(ord, 'x\x9c\xe3\xe5\xe5e\xe0e``\xe0E\xa3a\x00\x9d\x8f\x0e\xd0\xf5Q\xaa\x9f\x87\x87\x87\x81\x07I\x1e\x9dOm\x003\x1f\x9d\xa6\x97\xfd\x03\rF\xba\xffG:\x00\x00\x8f~\x03\xe9'.decode('zip'))

assert sum(alloc) == 1000
assert min(alloc) >= 0
assert len(alloc) == 1024
loss, pattern = worst_loss(alloc)
print "Worst loss is", loss, "when pattern is", "{:010b}".format(73)
# Worst loss is 104 when pattern is 0001001001
