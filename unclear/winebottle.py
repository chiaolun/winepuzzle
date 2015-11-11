#!/usr/bin/env python

# This is the solution to the poisoned wine problem for 2 poisoned
# bottles of wine and 1000 testers

from scipy.misc import comb
from numpy import array, dot, outer, eye
from itertools import groupby, product, chain, combinations
from math import floor, ceil
from cvxopt.solvers import lp
from cvxopt import matrix

import cvxopt.solvers
cvxopt.solvers.options['show_progress'] = False

Nmice = 10
Nbottles = 1000
Norm  = array([int(comb(Nmice,i,1)) for i in range(Nmice+1)])

def flatten(listOfLists):
    "Flatten one level of nesting"
    return chain.from_iterable(listOfLists)

def opt_allocation(labels,verbose=0):
    "label given in format of a list of numbers e.g., [1,4,5]"
    Nlabels = len(labels)
    cases = []
    # constructs all possible numbers of dead mice
    pairs = sorted(product(labels,repeat=2),key=sum)
    cases_detail = {}
    for g, l in groupby(pairs,sum):
        case = [0.] * Nlabels
        for x in set(flatten(l)):
             case[labels.index(x)] = float(comb(g,x,1))
        cases.append(tuple(case))
    for g, l in groupby(pairs,sum):
        cases_detail[g] = set(flatten(l))
    # normalization condition
    norm    = array([float(comb(Nmice,x,1)) for x in labels])
    cases   = array(cases)
    # optimizes category weights
    score, sol = opt_weights(norm,cases,labels)
    if verbose:
        print_solution(labels, sol, score, cases_detail)
    return score

def partition(n,k,aux=[]):
    if k == 1:
        return aux + [n]
    else:
        parts = [partition(n-i, k-1, aux + [i]) for i in range(n+1)]
        if k == 2:
            return parts
        else:
            return flatten(parts)

def opt_weights(norm,cases,labels):
    "Optimizes distribution of wine bottles to minimax cases given norm constraint"
    Nlabels = len(labels)
    # First does float version as a start
    sol = lp(-matrix(norm),
              matrix(array(list(cases)+list(-eye(Nlabels)))),
              matrix([1.]*len(cases)+[0.]*Nlabels))['x']
    sol = array(sol)
    sol /= dot(norm,sol)
    # Determine the floored allocation according to the float version
    alloc = array([int(i*j*Nbottles) for i,j in zip(norm,sol)])
    excess = Nbottles - sum(alloc)
    combs = [alloc + array(x) for x in partition(excess, Nlabels)]
    def score(comb):
        assert sum(comb) == Nbottles
        baserem = [divmod(combi,int(n)) for combi, n in zip(comb,norm)]
        return max([sum([(base*int(case0) + min(rem,int(case0)))
                         for (base, rem), case0 in zip(baserem,case)])
                    for case in cases])
    return min([(score(comb), list(comb)) for comb in combs])

def print_solution(labels, sol, score, cases_detail):
    solaux = {}
    print "Group the %i bottles as follows:" % Nbottles
    for l, s in zip(labels, sol):
        Ngroups = comb(Nmice,l,1)
        n1, rem = divmod(s,Ngroups)
        print "%4i %iC%i labels get %3i groups with %4i bottles each, and %3i groups with %4i bottles each"\
            % (Ngroups,Nmice,l,Ngroups-rem,n1,rem,n1+1)
        solaux[l] = rem, n1
    worst = 0
    for g, l in cases_detail.items():
        print "If %i mice die," % g
        totloss = 0
        for l0, ll in zip(labels,l):
            Ndiscard = comb(g,ll,1)
            print "  discard %4i %iC%i labels" % (Ndiscard, Nmice, ll),
            rem, n1 = solaux[ll]
            Nbdiscard = min(rem, Ndiscard)*(n1+1) + max(Ndiscard-rem,0)*n1
            totloss += Nbdiscard
            print "which at worst have %4i bottles" % Nbdiscard
        print "  totalling to %i bottles of wine lost" % totloss
        worst = max(worst,totloss)
    print "The worst possible loss is %4i bottles." % worst

#label_sets = flatten([combinations(range(Nmice),i) for i in range(1,Nmice+1)])
label_sets = [range(i) for i in range(1,Nmice+1)]

best_sol = Nbottles, None
for l in label_sets:
    waste = opt_allocation(l)
    if waste < best_sol[0]:
        best_sol = waste, l
#        print l, waste

opt_allocation(best_sol[1],verbose=1)

# Group the 1000 bottles as follows:
#    1 10C0 labels get   1 groups with   42 bottles each, and   0 groups with   43 bottles each
#   10 10C1 labels get   4 groups with   12 bottles each, and   6 groups with   13 bottles each
#   45 10C2 labels get  44 groups with    5 bottles each, and   1 groups with    6 bottles each
#  120 10C3 labels get 114 groups with    5 bottles each, and   6 groups with    6 bottles each
# If 0 mice die,
#   discard    1 10C0 labels which at worst have   42 bottles
#   totalling to 42 bottles of wine lost
# If 1 mice die,
#   discard    1 10C0 labels which at worst have   42 bottles
#   discard    1 10C1 labels which at worst have   13 bottles
#   totalling to 55 bottles of wine lost
# If 2 mice die,
#   discard    1 10C0 labels which at worst have   42 bottles
#   discard    2 10C1 labels which at worst have   26 bottles
#   discard    1 10C2 labels which at worst have    6 bottles
#   totalling to 74 bottles of wine lost
# If 3 mice die,
#   discard    1 10C0 labels which at worst have   42 bottles
#   discard    3 10C1 labels which at worst have   39 bottles
#   discard    3 10C2 labels which at worst have   16 bottles
#   discard    1 10C3 labels which at worst have    6 bottles
#   totalling to 103 bottles of wine lost
# If 4 mice die,
#   discard    4 10C1 labels which at worst have   52 bottles
#   discard    6 10C2 labels which at worst have   31 bottles
#   discard    4 10C3 labels which at worst have   24 bottles
#   totalling to 107 bottles of wine lost
# If 5 mice die,
#   discard   10 10C2 labels which at worst have   51 bottles
#   discard   10 10C3 labels which at worst have   56 bottles
#   totalling to 107 bottles of wine lost
# If 6 mice die,
#   discard   20 10C3 labels which at worst have  106 bottles
#   totalling to 106 bottles of wine lost
# The worst possible loss is  107 bottles.
