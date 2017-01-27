#!/usr/bin/python

import numpy as np
import rings as ri
import arithmetic as ar

def FF(p):
    '''Creates an instance of the field Z/p.'''
    
    def ff_init(self,params):
        '''Creates a number in Z/p normal form'''
        self.value=params%p
        
    def ff_add(x,y):
        '''Addition for numbers over Z/p'''
        return (x+y)%p
    
    def ff_sub(x,y):
        '''Subtraction for numbers over Z/p'''
        return (x-y)%p
    
    def ff_mul(x,y):
        '''Multiplication for numbers over Z/p'''
        return (x*y)%p
    
    def ff_div(x,y):
        '''Division for numbers over Z/p'''
        if y%p!=0:
            return x*ff_inv(y)
        else:
            print("Can't divide by zero!")
            return None
        
    def ff_inv(x):
        '''Inversion of numbers over Z/p'''
        if x%p!=0:
            (x,y,g)=ar.extendedEuclideanAlgorithm(p,x%p)
            return g*y
        else:
            print("Can't invert zero!")
            return None
    
    def ff_eq(x,y):
        '''Tests equality of numbers over Z/p'''
        if x%p==y%p:
            return True
        else:
            return False

    def ff_num(n):
        '''Sends an integer n to its reduction mod p.'''
        return n%p
    
    def ff_print(x):
        '''Prints a number over Z/p.'''
        return str(x%p)

    return ri.Field(ff_init,ff_add,ff_sub,ff_mul,ff_div,
                    ff_inv,ff_eq,ff_num,ff_print,p)
