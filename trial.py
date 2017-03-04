#!/usr/bin/python

import finite_fields as ff
import rationals as QQ
import a_infinity as ainf
from itertools import permutations,combinations_with_replacement

K=ff.FF(2)
N=3
depth=4
milnor_number=3
A=ainf.BP(milnor_number+1,2,K,2,1)
B=ainf.BP(milnor_number+1,2,K,N,1)

def word_twist(w,P):
    for n in w:
        P=P.twist(n)
    return P

def depth_twist():
    results_A={(): A.total_yoneda()}
    results_B={(): B.total_yoneda()}
    for n in range(1,depth+1):
        for results in [results_A,results_B]:
            new_results={}
            for w in results:
                for X in results[()].cat.objects:
                    new_results.update({w+(X,): results[w].twist(X)})
                    #new_results[w+(X,)].verify()
            results.update(new_results)

    garsides={w: results_A[w].width()[1]-2 for w in results_A}
    new_answers={w: results_B[w].width()[1] for w in results_B}
    for w in garsides:
        if new_answers[w] not in range((N-1)*garsides[w]+N,(N-1)*garsides[w]+2*N-1):
            print('Exception:',w,'garside',garsides[w],'new',new_answers[w])
            print(results_A[w].total())
            print(results_B[w].total())
        else:
            print('Fine:',w,'garside',garsides[w],'new',new_answers[w])
            #print(results_A[w].total())
            #print(results_B[w].total())

depth_twist()
#word=[1,2,3,1,3,2]
#print(word_twist(word,A.total_yoneda())-2)
#print(word_twist(word,B.total_yoneda()))
