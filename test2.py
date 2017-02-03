#!/usr/bin/python

import numpy as np
import rings as ri
import rationals as QQ
import graded_linear_algebra as gla
import a8_categories as cat

'''print("Testing A8Category")'''
K=QQ.QQ()
V=gla.sph(K,2)
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
A=cat.A8Category(K,objects,morphisms,operations)
'''print(A.base)
print(A.objects)
print(A.hom(1,1))
print(A.mu((1,1,1)).gr_map(0)[0,0])
print("Testing A8Modules")'''
M=A.yoneda(1)
'''print(M.mod(1))'''
### Need to test ltimes
#M4=M2.twist(1)
'''for n in M.modules:
    print(n,M.modules[n])'''
M1=M.twist(1)
print(M1)
'''for n in M1.modules:
    print(n,M1.modules[n])'''
#print(M1.operations)
#print(M1.mu((2,)))
#print(M1.cpx(1).cohomology())
#print(M1.cpx(2).cohomology())
#M2=M.twist(2)
#for n in M2.modules:
#    print(n,M2.modules[n])
#M3=M1.twist(2)
#for n in M3.modules:
#    print(n,M3.modules[n])
'''for n in M4.modules:
    print(n,M4.modules[n])'''

    
    
