#!/usr/bin/python

import numpy as np
import graded_linear_algebra as gla
import rings as ri
import finite_fields as ff

'''This should take three graded vector spaces
and return an array whose nth entry
is an array telling us how to shuffle the basis
of (A*B)*C to get the basis of A*(B*C).
'''

def shuffle(n,A,B,C):
    '''Args:
            n [int]
            A,B,C [gla.GradedVectorSpace]

    '''
    m1=min(A.graded_dim.keys())
    m2=min(B.graded_dim.keys())
    m3=min(C.graded_dim.keys())
    r1=max(A.graded_dim.keys())
    r2=max(B.graded_dim.keys())
    r3=max(C.graded_dim.keys())
    pre_shuffle={}
    post_shuffle={}
    ctr=0
    for x in range(m1,r1+1):
        for u in range(0,A.gr_dim(x)):
            for y in range(m2,r2+1):
                for v in range(0,B.gr_dim(y)):
                    for w in range(0,C.gr_dim(n-x-y)):
                        post_shuffle[(x,u,y,v,n-x-y,w)]=ctr
                        ctr+=1
    ctr=0
    for i in range(m1+m2,r1+r2+1):
        for x in range(m1,r1+1):
            for u in range(0,A.gr_dim(x)):
                for v in range(0,B.gr_dim(i-x)):
                    for w in range(0,C.gr_dim(n-i)):
                        pre_shuffle[ctr]=(x,u,i-x,v,n-i,w)
                        ctr+=1
    shuffle_array=[]
    for index in pre_shuffle:
        shuffle_array.append(post_shuffle[pre_shuffle[index]])
    return shuffle_array

def shuffle_map(f,A,B,C):
    '''f is a graded linear map f:A*(B*C)-->W
    and we want to get the map g:(A*B)*C-->W
    obtained by precomposing with the canonical shuffle.
    '''
    d=f.degree
    T=f.target
    g=gla.GradedLinearMap(d,(A.otimes(B)).otimes(C),T)
    for n in f.graded_map:
        shuffle_array=shuffle(n,A,B,C)
        M=f.gr_map(n)
        g.graded_map[n]=M[:,shuffle_array]
    return g
    


K=ff.FF(2)
V=gla.GradedVectorSpace(K)
V.graded_dim={-1:1,0:1,1:2,2:1,3:1}
W=gla.GradedVectorSpace(K)
W.graded_dim={0:1,2:1}
print(shuffle(3,V,W,W))
mult=gla.GradedLinearMap(0,W.otimes(W),W)
mult.graded_map[0]=K.num(np.array([[1]]))
mult.graded_map[2]=K.num(np.array([[1,1]]))
f=V.eye().otimes(mult)
g=shuffle_map(f,V,W,W)
f.display()
g.display()
(f+g).display()
