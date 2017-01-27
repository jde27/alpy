#!/usr/bin/python

import numpy as np

def convert(x,func):
    '''A function to convert integers to numbers over a field.'''
    if type(x) is not np.ndarray:
        if type(x) is not number:
            return func(x)
        else:
            return x
    else:
        N=np.empty(shape=x.shape,dtype='object')
        for s,t in enumerate(x):
            N[s]=convert(t,func)
        return N

class ring():
    '''Creates a ring R (algebraic structure which allows
    addition, subtraction, multiplication).
    
    Attributes:
        num_init [function]: Creates numbers over R from input data.
        num_add, num_sub, num_mult,
            [functions]: Used by numbers over R to perform addition,
                subtraction, and multiplication.
        num_eq [function]: Used by numbers over R to test their
            equality with others.
        mod [function]: Turns an integer into a number over R.
        char [int]: The characteristic of R.
    '''
    def __init__(self,num_init,num_add,num_sub,num_mul,\
                 num_eq,num_mod,num_print,char):
        self.num_init=num_init
        self.num_add=num_add
        self.num_sub=num_sub
        self.num_mul=num_mul
        self.num_eq=num_eq
        self.num_mod=num_mod
        self.num_print=num_print
        self.char=char

    def mod(self,n):
        def converter(x):
            return number(self,self.num_mod(x))
        return convert(n,converter)

class field(ring):
    '''Creates a field k (algebraic structure which allows
    addition, multiplication, subtraction, division).
    
    A field is a ring, and inherits all methods and attributes from that
    class. In addition, it has attributes for dividing and inverting numbers.

    Attributes:
        num_div [function]: Used by numbers over k to perform division.
        num_inv [function]: Used by numbers over k to invert themselves over k.
    '''
    def __init__(self,num_init,num_add,num_sub,num_mul,num_div,\
                 num_inv,num_eq,num_mod,num_print,char):
        super().__init__(num_init,num_add,num_sub,num_mul,\
                       num_eq,num_mod,num_print,char)
        self.num_div=num_div
        self.num_inv=num_inv

        
class number():
    '''Instantiates a number over a specified field, k.

    The methods for adding and multiplying numbers over k (amongst
    other things) are defined in the description of k.
    '''
    def __init__(self,k,*args):
        '''Sets the base field and uses the functionality of k
        to create a new number from whatever arguments it's supposed to.
        '''
        self.base=k
        self.base.num_init(self,*args)

    def __add__(self,other):
        k=self.base
        return number(k,k.num_add(self,other))

    def __sub__(self,other):
        k=self.base
        return number(k,k.num_sub(self,other))
    
    def __mul__(self,other):
        k=self.base
        return number(k,k.num_mul(self,other))

    def __truediv__(self,other):
        k=self.base
        return number(k,k.num_div(self,other))

    def I(self):
        k=self.base
        return number(k,k.num_inv(self))

    def __eq__(self,other):
        k=self.base
        return k.num_eq(self,other)

    def __str__(self):
        k=self.base
        return k.num_print(self)

