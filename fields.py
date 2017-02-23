#!/usr/bin/python

class Field():
    '''Creates a field K (algebraic structure which allows
    addition, subtraction, multiplication, division).
    
    Attributes:

        num_init [function]: Creates numbers over K from input data.

        num_add, num_sub, num_mult, num_div, num_inv, num_neg
            [functions]: Used by numbers over K to perform addition,
                subtraction, multiplication, division,
                inversion and negation.

        num_eq [function]: Used by numbers over K to test their
            equality with others.

        num_num [function]: Turns an integer into a number over K.

        char [int]: The characteristic of K.

    Methods:

        K(n):
            Args: n [int]
            Returns [number]: Turns n into an element of K.
    '''
    def __init__(self,num_init,num_add,num_sub,num_mul,
                 num_div,num_inv,num_neg,num_eq,num_num,
                 num_print,char):
        self.num_init=num_init
        self.num_add=num_add
        self.num_sub=num_sub
        self.num_mul=num_mul
        self.num_div=num_div
        self.num_inv=num_inv
        self.num_neg=num_neg
        self.num_eq=num_eq
        self.num_num=num_num
        self.num_print=num_print
        self.char=char

    def __call__(self,n):
        '''This function converts integers into field elements
        using the function specified by num_num.
        '''
        if type(n) is Number:
            return n
        else:
            return Number(self,self.num_num(n))
        
class Number():
    '''Instantiates a number over a specified field, K.

    The methods for adding and multiplying numbers over K etc
    are defined in the description of K.
    '''
    def __init__(self,K,*args):
        '''Sets the base field and uses the functionality of K
        to create a new number from whatever arguments it's supposed to.
        '''
        self.field=K
        self.field.num_init(self,*args)

    def __add__(self,other):
        K=self.field
        return Number(K,K.num_add(self,other))

    def __sub__(self,other):
        K=self.field
        return Number(K,K.num_sub(self,other))
    
    def __mul__(self,other):
        K=self.field
        if type(other) is Number:
            return Number(K,K.num_mul(self,other))
        else:
            return NotImplemented

    def __truediv__(self,other):
        K=self.field
        return Number(K,K.num_div(self,other))

    def __neg__(self):
        K=self.field
        return Number(K,K.num_neg(self))
    
    def I(self):
        K=self.field
        return Number(K,K.num_inv(self))

    def __eq__(self,other):
        K=self.field
        return K.num_eq(self,other)

    def __str__(self):
        K=self.field
        return K.num_print(self)
