#!/usr/bin/python

import rings as ri
import rationals as QQ
import graded_linear_algebra as gla
import numpy as np

K=QQ.QQ()
Z=gla.GradedVectorSpace(K)
Z.graded_dim={0:1,1:1,2:1}
V=gla.GradedVectorSpace(K)
V.graded_dim={2:1}
W=gla.sph(K,2)
F=gla.GradedLinearMap(0,V,W)
F.graded_map={2:K.num(np.array([[1]]))}
Z.eye().otimes(F).verify()
