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
print("numpy.int32 ",type(M.gr_map(0)[0,0].numerator))
R=W.shift()
print("R:{-1:1}=",R,"W:{0:1}=",W)
for i in V.graded_dim:
    print(V.gr_dim(i))
print(V.base)
U=gla.GradedVectorSpace(K)
U.graded_dim={0:1,1:1,3:2,4:3}
F=gla.GradedLinearMap(1,U,U)
print("F.degree=1,F.source,F.target",F.degree,F.source,F.target)
M=K.num(np.array([[1]]))
x=K.num(1)/K.num(2)
N=K.num(np.array([[1,0],[1,0],[0,x]]))
F.graded_map={0:M,3:N}
print("gr_map",F.gr_map(0)[0,0])
A,B=F.ker_im()
print(A,B)
C=gla.CochainComplex(U,F)
print(C.cohomology())
print("------------------")
G=gla.GradedLinearMap(0,U,V)
G.graded_map={0:M}
H=F.otimes(G)
A,B=H.ker_im()
print(H.source)
print(H.target)
print("deg",H.degree)
for i in range(0,3):
    for j in range(0,2):
        print(H.gr_map(3)[i,j])
print(A,B)
print("F")
F.verify()
print("G")
G.verify()
print("H")
H.verify()
I=gla.GradedLinearMap.block(F,F,F,F)
print("I")
I.verify()
t,s=I.gr_map(0).shape
for i in range(t):
    print("Row")
    for j in range(s):
        print(I.gr_map(0)[i,j])
#J=I.koszulify()
#J=I.rejig_1()
J=G.rejig_3()
print("J")
J.verify()
for m in J.graded_map:
    print(J.gr_map(m)[0,0])
#print(I.source,I.target)
print(G.source,G.target)
print(J.source,J.target)
print("Testing addition")
u=gla.GradedVectorSpace(L)
u.graded_dim={0:1,1:1,3:2,4:3}
f=gla.GradedLinearMap(1,u,u)
m=L.num(np.array([[1]]))
o=L.num(1)/L.num(2)
n=L.num(np.array([[1,0],[1,0],[0,o]]))
f.graded_map={0:m,3:n}
r=f+f+f
A,B=F.ker_im()
print(A,B)
A,B=r.ker_im()
print(A,B)
print("Testing composition")
f.verify()
rr=f*f
rr.verify()
print("Still testing composition")
g=gla.GradedLinearMap(0,u,u)
m=L.num(np.array([[2]]))
g.graded_map={0:m}
g.verify()
(g*g).verify()
print(g.gr_map(0)[0,0])
print((g*g).gr_map(0)[0,0])
(f.oplus(f)).verify()

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
