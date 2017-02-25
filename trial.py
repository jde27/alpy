#!/usr/bin/python

import finite_fields as ff
import rationals as QQ
import a_infinity as ainf

K=ff.FF(2)
#K=QQ.QQ()
A=ainf.BP(2,3,K,2)
M=A.yoneda(1)
N=A.yoneda(2)
# Enter braid word here:
w=[1,2,1,2,1,2]
P=M
Q=N
for n in w:
    R,S=P.twist(n),Q.twist(n)
    P,Q=R,S
print(P.width(),Q.width())
