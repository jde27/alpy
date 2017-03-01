#!/usr/bin/python

import finite_fields as ff
import rationals as QQ
import a_infinity as ainf
from itertools import permutations,combinations_with_replacement

K=ff.FF(2)
A=ainf.BP(2,3,K,2,1)
B=ainf.BP(2,3,K,3,1)
LA=A.total_yoneda()
LB=B.total_yoneda()

def word_twist(w,P):
    for n in w:
        P=P.twist(n)
    return P.width()[1]

results_A={(): A.total_yoneda()}
results_B={(): B.total_yoneda()}
for n in range(1,3):
    for results in [results_A,results_B]:
        new_results={}
        for w in results:
            new_results.update({w+(X,): results[w].twist(X)
                                for X in results[()].cat.objects})
        results.update(new_results)

for m in results_A:
    print('word',m,'garside',results_A[m].width()[1]-2,'new',results_B[m].width()[1])



    
'''words= [[1],[2],[1,1],[1,2],[2,1],[2,2],
        [1,1,1],[1,1,2],[1,2,1],[2,1,1],[1,2,2],[2,1,2],[2,2,1],[2,2,2],
        [1,1,1,1],[1,1,1,2],[1,1,2,1],[1,2,1,1],[2,1,1,1],
        [1,1,2,2],[1,2,1,2],[2,1,1,2],[1,2,2,1],[2,1,2,1],
        [2,2,1,1],[1,2,2,2],[2,1,2,2],[2,2,1,2],[2,2,2,1],
        [2,2,2,2]]'''
'''
items=[1,2]
words=[]
for n in range(1,8):
    for c in combinations_with_replacement(items,n):
        for d in permutations(c):
            if d not in words:
                words.append(d)

for w in words:
    print('word ',w,' garside ',word_twist(w,MA,NA)-2,
          ' other ',word_twist(w,MB,NB))

'''
