#!/usr/bin/python

import numpy as np

def find_nonzero(n):
    '''A function which finds the index of the first nonzero
    entry in an array

    Input: n [array]
    
    Output: [bool], [int], [object of n]
    
    If there is a nonzero entry in position x,
    find_nonzero returns True, x, n[x]
    otherwise, it returns False, 0, 0
    '''
    ctr=0
    while ctr<n.size:
        if n[ctr]!=0:
            return True,ctr,n[ctr]
        else:
            ctr+=1
    return False,0,0
