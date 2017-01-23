#!/usr/bin/python

### These are some routines involving integers that are required
### for defining operations involving more general rings/fields,

def extendedEuclideanAlgorithm(a,b):
    '''An implementation of the extended Euclidean algorithm,
    taken from Jeremy Kun's blog.
    '''
    if abs(b) > abs(a):
        (x,y,d) = extendedEuclideanAlgorithm(b,a)
        return (y,x,d)

    if abs(b) == 0:
        return (1,0,a)

    x1, x2, y1, y2 = 0, 1, 1, 0
    while abs(b)>0:
        q,r = divmod(a,b)
        x=x2-q*x1
        y=y2-q*y1
        a, b, x2, x2, y2, y1 = b, r, x1, x, y1, y

    return (x2, y2, a)

def gcd(a,b):
    '''Compute the greatest common divisor of two integers'''
    p, q, r = extendedEuclideanAlgorithm(a,b)
    return abs(r)

