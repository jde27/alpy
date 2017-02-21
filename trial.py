#!/usr/bin/python

import fields as fi
import finite_fields as FF
import rationals as QQ
from linear_algebra import *
from a_infinity import *

#K=FF.FF(2)
K=QQ.QQ()
A=BP(2,3,K,2)
A.verify()
M=A.yoneda(1)
N=A.yoneda(2)
P=M.twist(2)
Q=P.twist(2)
#R=Q.twist(2)
#P.verify()
#Q.verify()
#R.verify()
#R.cpx(1).verify()
#R.cpx(2).verify()
#print(P.cpx(1).cohomology())
print(Q.cpx(1).cohomology())
#print(R.cpx(1).cohomology())

