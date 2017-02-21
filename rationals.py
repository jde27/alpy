#!/usr/bin/python

import fields as fi
import arithmetic as ar

def rat_normal_form(a,b):
    '''Puts the numerator and denominator into normal form:

    i.e. denominator > 0
    gcd(numerator,denominator)=1
    '''
    if b==0:
        print("Tried to divide by zero!")
    elif a==0:
        return 0,1
    else:
        sign=a*b//abs(a*b)
        numerator=sign*abs(a)//ar.gcd(a,b)
        denominator=abs(b)//ar.gcd(a,b)
        return numerator,denominator

def rat_init(self,params):
    '''Creates a rational number object a/b in normal form'''
    a,b=params
    self.numerator,self.denominator = rat_normal_form(a,b)

def rat_add(x,y):
    '''Addition for rational numbers'''
    a=x.numerator*y.denominator+x.denominator*y.numerator
    b=x.denominator*y.denominator
    return a,b

def rat_sub(x,y):
    '''Subtraction for rational numbers'''
    a=x.numerator*y.denominator-x.denominator*y.numerator
    b=x.denominator*y.denominator
    return a,b

def rat_mul(x,y):
    '''Multiplication for rational numbers'''
    a=x.numerator*y.numerator
    b=x.denominator*y.denominator
    return a,b

def rat_div(x,y):
    '''Division for rational numbers'''
    a,b=x.numerator,x.denominator
    c,d=y.numerator,y.denominator
    if b*c!=0:
        return a*d,b*c
    else:
        print("Cannot divide by zero.")

def rat_inv(x):
    '''Inversion of rational numbers'''
    a,b=x.numerator,x.denominator
    if a!=0:
        return b,a
    else:
        print("Cannot invert zero.")

def rat_neg(x):
    '''Negation of rational numbers'''
    a,b=x.numerator,x.denominator
    return -a,b
        
def rat_eq(x,y):
    '''Tests equality of rational numbers'''
    if type(y) is fi.Number:
        if x.numerator==y.numerator and x.denominator==y.denominator:
            # Since rational numbers are stored with their numerator and
            # denominator coprime and the sign in the numerator,
            # equality holds iff equality holds for numerator and
            # denominator separately.
            return True
        else:
            return False
    else:
        # We may be comparing a rational number to an integer:
        if x.numerator==y and x.denominator==1:
            return True
        else:
            return False

def rat_num(n):
    '''Sends an integer n to the rational number n/1.'''
    return n,1

def rat_print(x):
    '''Prints a rational number nicely in form a/b'''
    if x.denominator!=1:
        return str(x.numerator)+"/"+str(x.denominator)
    else:
        return str(x.numerator)

def QQ():
    '''Defines an instance of the field of rational numbers'''
    return fi.Field(rat_init,rat_add,rat_sub,rat_mul,rat_div,
                    rat_inv,rat_neg,rat_eq,rat_num,rat_print,0)
    
'''
#Some random tests
'''
#Q=QQ()
#x=fi.Number(Q,(5,7))
#y=fi.Number(Q,(7,5))
#z=Q.num(9)
#u=Q.num(3)
#v=Q.num(26)
#w=v.I()
#print((y*z*u/v).numerator,"/",(y*z*u/v).denominator)
#print((x*y).denominator)
#print(w)
#g=(x+y-u)*w
#print(g)
#M=np.array([x, y, z])
#N=np.array([[x],[y],[z]])
#B=M.dot(N)
#print(B[0])
#D=np.array([1,2])
#print(D)
#C=np.array([D])
#print(C)
#M=np.array([[1,2,3],[2,3,4],[3,4,5]])
#print(M)
#C=Q.num(M)
#print(C)
#D=C.dot(C)
#for i in range(0,3):
#    for j in range(0,3):
#        print(D[i][j])
#print(M.dot(M))
#x=fi.Number(Q,(5,7))
#y=fi.Number(Q,(4,9))
#if x-y!=0:
#    print("Correct")
#if x-x!=0:
#    print("Wrong!")

