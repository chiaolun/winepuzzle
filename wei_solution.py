#!/usr/bin/env python

from itertools import combinations
import numpy as np

nmice = 10
nbottles = 1000
nlabels = 2**nmice
labels = range(nlabels)

tothrow = np.zeros((nlabels,nlabels,nlabels),dtype=bool)
for i, j in combinations(labels, 2):
    k = i | j
    tothrow[k,i,j] = True
    tothrow[k,j,i] = True
for i in labels:
    tothrow[i,i,i] = True

# The resulting losses for every possible outcome would be
def worst_loss(alloc):
    nonzero = alloc > 0
    nonzero = np.outer(nonzero, nonzero)
    losses = (tothrow * nonzero).any(axis=2).dot(alloc)
    kmax = losses.argmax()
    return losses[kmax], kmax

alloc = map(ord, 'x\x9c\xe3\xe5\xe5e\xe0e``\xe0E\xa3a\x00\x9d\x8f\x0e\xd0\xf5Q\xaa\x9f\x87\x87\x87\x81\x07I\x1e\x9dOm\x003\x1f\x9d\xa6\x97\xfd\x03\rF\xba\xffG:\x00\x00\x8f~\x03\xe9'.decode('zip'))
alloc = np.array(alloc)

assert sum(alloc) == 1000
assert min(alloc) >= 0
assert len(alloc) == 1024
loss, pattern = worst_loss(alloc)
print "Worst loss is", loss, "when pattern is", "{:010b}".format(pattern)
# Worst loss is 104 when pattern is 0001001001
