#!/usr/bin/python

import rings as ri
import finite_fields as FF
import graded_linear_algebra as gla
import numpy as np

K=FF.FF(2)
V=gla.GradedVectorSpace(K)
V.graded_dim={0:1,1:2,2:3,3:4}
F=gla.GradedLinearMap(1,V,V)
z1=K.num(np.array([[1],[0]]))
z2=K.num(np.array([[0,1],[0,0],[0,1]]))
z3=K.num(np.array([[0,1,0],[1,1,1],[1,0,1],[0,1,0]]))
F.graded_map={0:z1,1:z2,2:z3}
for n in F.gr_map(0):
    for m in n:
        print(m)
D=gla.CochainComplex(V,F)
print(D.cohomology())

'''K=QQ.QQ()
Z=gla.GradedVectorSpace(K)
Z.graded_dim={0:1,1:1,2:1}
V=gla.GradedVectorSpace(K)
V.graded_dim={2:1}
W=gla.sph(K,2)
F=gla.GradedLinearMap(0,V,W)
F.graded_map={2:K.num(np.array([[1]]))}
Z.eye().otimes(F).verify()'''


