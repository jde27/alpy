#!/usr/bin/python

import numpy as np
import rings as ri
import rationals as QQ
import finite_fields as FF
import graded_linear_algebra as gla
import a8_categories as cat

'''print("Testing A8Category")'''
K=FF.FF(2)
A=cat.DynkinGraph.BP(3,2,K,2)
M=A.yoneda(1)
N=A.yoneda(2)
def disp(M):
    print(M.cpx(1).cohomology())
    print(M.cpx(2).cohomology())
E=N.twist(1).twist(1).twist(1).twist(1)
print(E.mod(1))
print(E.mod(2))
print(E.operations)
disp(E)
#disp(N)
#disp(N.twist(1).twist(1))


'''print(A.objects)
for X in A.objects:
    for Y in A.objects:
        print(A.hom(X,Y))
for X in A.objects:
    for Y in A.objects:
        for Z in A.objects:
            print(A.mu((X,Y,Z)))'''
'''V=gla.sph(K,2)
W=gla.pt(K,0)
X=gla.pt(K,2)
F=gla.sph_op(V)
G1=gla.simple(W.otimes(V),W,0)
G2=gla.simple(X.otimes(W),V,2)
G3=gla.simple(V.otimes(X),X,2)
G4=gla.simple(V.otimes(W),W,0)
G5=gla.simple(W.otimes(X),V,2)
G6=gla.simple(X.otimes(V),X,2)
objects={1,2}
morphisms={(1,1):V,(2,2):V,(1,2):W,(2,1):X}
operations={(1,1,1):F,(2,2,2):F,
            (1,1,2):G1,(1,2,1):G2,(2,1,1):G3,
            (1,2,2):G4,(2,1,2):G5,(2,2,1):G6}
A=cat.A8Category(K,objects,morphisms,operations)'''
'''print(A.base)
print(A.objects)
print(A.hom(1,1))
print(A.mu((1,1,1)).gr_map(0)[0,0])
print("Testing A8Modules")'''
'''print(M.mod(1))'''
### Need to test ltimes
#M4=M2.twist(1)
'''for n in M.modules:
    print(n,M.modules[n])'''
'''for n in M1.modules:
    print(n,M1.modules[n])'''
#print(M1.operations)
#print(M1.mu((2,)))
#print("The twist of 1 around 1")
#M1=M.twist(1)
'''print(M1.cpx(1).cohomology())
print(M1.cpx(2).cohomology())'''

'''
print("The sphere 1")
disp(M)
print("The twist of 1 around 2")
disp(M.twist(2))
print("The twist of 1 around 2 around 1")
disp(M.twist(2).twist(1))
print("The twist of 1 around 2 around 1 around 2")
disp(M.twist(2).twist(1).twist(2))
print("The twist of 1 around 2 around 1 around 2 around 1")
disp(M.twist(2).twist(1).twist(2).twist(1))
print("The twist of 1 around 2 around 1 around 2 around 1 around 2")
disp(M.twist(2).twist(1).twist(2).twist(1).twist(2))
print("The twist of 1 around 2 around 1 around 2 around 1 around 2 around 1")
disp(M.twist(2).twist(1).twist(2).twist(1).twist(2).twist(1))'''





#for n in M2.modules:
#    print(n,M2.modules[n])
#M3=M1.twist(2)
#for n in M3.modules:
#    print(n,M3.modules[n])
'''for n in M4.modules:
    print(n,M4.modules[n])'''

    
    
