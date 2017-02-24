#!/usr/bin/python

import finite_fields as ff
import a_infinity as ainf

K=ff.FF(2)
A=ainf.BP(2,3,K,2)
A.verify()
M=A.yoneda(1)
N=A.yoneda(2)
P=M.twist(2)
#print('DISPLAYING P')
#P.display()
print('*********************************')
Q=P.twist(2)
#print('DISPLAYING Q')
#Q.display()
print('*********************************')
R=Q.twist(2)
#print('DISPLAYING R')
#R.display()
print('*********************************')
S=R.twist(2).twist(2)
T=S.simplify()
print(S.cpx(1).cohomology())
print(T.cpx(1).cohomology())
T.display()
#S.cpx(1).verify()
#R.cpx(2).verify()
#print(P.cpx(1).cohomology())
#print(Q.cpx(1).cohomology())
#print(R.cpx(1).cohomology())

