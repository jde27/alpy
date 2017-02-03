#!/usr/bin/python

import rings as ri
import rationals as QQ
import graded_linear_algebra as gla

K=QQ.QQ()
V=gla.sph(K,2)
F=gla.sph_op(V)
(V.eye().otimes(F)).verify()
print({p+q for p in [0,2] for q in [0,2,4]})
