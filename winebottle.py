#!/usr/bin/env python

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

