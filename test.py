#!/usr/bin/python

import finite_fields as fi
import rationals as rat
import rings as ri
import numpy as np
import graded_linear_algebra as gla

L=fi.FF(3)
K=rat.QQ()
V=gla.GradedVectorSpace(K)
V.graded_dim={0:1,2:1}
print("V:{0:1,2:1}=",V)
W=gla.GradedVectorSpace(K)
W.graded_dim={0:1}
#print(W)
X=V.oplus(W)
print("X:{0:2,2:1}=",X)
Y=V.otimes(X)
print("Y:{0:2,2:3,4:1}=",Y)
print("X+Y:{0:4,2:4,4:1}=",X.oplus(Y))
print(V.gr_dim(0),V.gr_dim(1))
M=V.eye()
print("1,1=",M.gr_map(0)[0,0],M.gr_map(2)[0,0])
print(type(M.gr_map(0)[0,0].numerator))
R=W.shift()
print("R:{-1:1}=",R,"W:{0:1}=",W)
for i in V.graded_dim:
    print(V.gr_dim(i))
print(V.base)
'''
M=np.array([[1,8],[20,0]])
N=K.num(M)
for i in range(0,2):
    for j in range(0,2):
        print(N[i,j])
for i in range(0,2):
    for j in range(0,2):
        print((N.dot(N).dot(N))[i,j])
print("Check that the original matrix is unchanged.")
for i in range(0,2):
    for j in range(0,2):
        print(N[i,j])
v=ri.Number(L,(3,27))
P=L.num(np.array([[19,v],[2,1]]))
for i in range(0,2):
    for j in range(0,2):
        print((P.dot(P).dot(P))[i,j])
print(v.I())
'''
