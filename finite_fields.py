#!/usr/bin/python

import fields as fi
import arithmetic as ar

def FF(p):
    '''Creates an instance of the field Z/p.'''
    
    def ff_init(self,params):
        '''Creates a number in Z/p normal form'''
        self.value=params%p
        
    def ff_add(x,y):
        '''Addition for numbers over Z/p'''
        a,b=x.value,y.value
        return (a+b)%p
    
    def ff_sub(x,y):
        '''Subtraction for numbers over Z/p'''
        a,b=x.value,y.value
        return (a-b)%p
    
    def ff_mul(x,y):
        '''Multiplication for numbers over Z/p'''
        a,b=x.value,y.value
        return (a*b)%p
    
    def ff_div(x,y):
        '''Division for numbers over Z/p'''
        z=ff_inv(y)
        a=x.value
        return (a*z)%p
    
    def ff_inv(x):
        '''Inversion of numbers over Z/p'''
        a=x.value
        if a%p!=0:
            (s,t,g)=ar.extendedEuclideanAlgorithm(p,a%p)
            return g*t
        else:
            print("Can't invert zero!")
            return None

    def ff_neg(x):
        '''Negation of numbers over Z/p'''
        a=x.value
        return p-a
    
    def ff_eq(x,y):
        '''Tests equality of numbers over Z/p'''
        a=x.value
        if type(y) is fi.Number:
            b=y.value
        else:
            b=y
        if a%p==b%p:
            return True
        else:
            return False

    def ff_num(n):
        '''Sends an integer n to its reduction mod p.'''
        return n%p
    
    def ff_print(x):
        '''Prints a number over Z/p.'''
        a=x.value
        return str(a%p)

    return fi.Field(ff_init,ff_add,ff_sub,ff_mul,ff_div,
                    ff_inv,ff_neg,ff_eq,ff_num,ff_print,p)
