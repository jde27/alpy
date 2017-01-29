#!/usr/bin/python

import numpy as np
import rings as ri
import rationals as QQ
import graded_linear_algebra as gla
import a8_categories as cat

print("Testing A8Category")
K=QQ.QQ()
V=gla.sph(K,2)
F=gla.sph_op(V)
objects={1}
morphisms={(1,1):V}
operations={(1,1,1):F}
A=cat.A8Category(K,objects,morphisms,operations)
print(A.base)
print(A.objects)
print(A.hom(1,1))
print(A.mu(1,1,1).gr_map(0)[0,0])
print("Testing A8Modules")
M=A.yoneda(1)
print(M.mod(1))
