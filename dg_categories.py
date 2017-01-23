#!/usr/bin/python

import numpy as np
import rings as ri
import rationals as QQ
import vector_spaces as vs
import graded_vector_spaces as gvs

class dg_cat():
    'Differential graded category'
    def __init__(self,obj,mor,diff,prod):
        '''
        Input:
        obj/objects is a set
        mor/morphisms is a dictionary of gvs.graded_vector_space indexed by pairs of objects
        diff is a dictionary of gr_hom indexed by morphisms
        diff[(A,B)] is a gr_hom with:
           source morphism[(A,B)] and target morphism[(A,B)].shift()
        prod is a dictionary of gr_hom indexed by triples of objects
        prod[(A,B,C)] is a gr_hom with:
           source morphism[(B,C)]*morphism[(A,B)]
           target morphism[(A,C)]
        In all cases, we only bother to include nonzero stuff.
        '''#        for A,B in objects,objects
        self.objects=obj
        self.morphisms=mor
        self.diff=diff
        self.prod=prod
    def hom(self,X,Y):
        '''
        If A=self is a dg-category and X, Y are objects of A,
        A.hom(X,Y) returns the morphism space between X and Y [gr_vs],
        provided that space is defined (otherwise it returns
        the zero vector space).
        '''
        A=self
        if {X,Y}<=A.objects:
            if (X,Y) in A.morphisms:
                return A.morphisms[(X,Y)]
            else:
                return gvs.graded_vector_space.zero(K)
        else:
            print("Cannot take homs between ",X," and ",Y," as they do not live in ",A,".")
    def d(self,X,Y):
        '''
        If A=self is a dg-category and X,Y are objects of A,
        A.d(X,Y) returns the differential:

        A.hom(X,Y)-->A.hom(X,Y)[1]  [gr_hom],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 1).
        '''
        A=self
        if {X,Y}<=A.objects:
            if (X,Y) in A.diff:
                return A.diff[(X,Y)]
            else:
                return zero_gr_hom(1,A.hom(X,Y),A.hom(X,Y))###
    def m(self,X,Y,Z):
        '''
        If A=self is a dg-category and X,Y,Z are objects of A,
        A.m(X,Y,Z) returns the multiplication map:

        A.hom(Y,Z)*A.hom(X,Y)-->A.hom(X,Z)  [gr_hom],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 0).
        '''
        A=self
        if {X,Y,Z}<=A.objects:
            if (X,Y,Z) in A.prod:
                return A.prod[(X,Y,Z)]
            else:
                return zero_gr_hom(0,A.hom(Y,Z)*A.hom(X,Y),A.hom(X,Z))
    def yoneda(self,Z):
        '''
        Generates the right Yoneda module associated to the object Z
        Example:
        A=dg_cat({'X'},...)
        S=A.yoneda('X')
        '''
        A=self
        if Z in A.objects:
            M={}
            diff={}
            act={}
            for Y in A.objects:
                if (Y,Z) in A.morphisms:# Remember we only include nonzero spaces
                    M[Y]=A.hom(Y,Z)
                if (Y,Z) in A.diff:        # and nonzero maps
                    diff[Y]=A.d(Y,Z)
                for X in A.objects:
                    if (X,Y,Z) in A.prod:  # same here
                        act[(X,Y)]=A.m(X,Y,Z)
            return dg_mod(A,M,diff,act)
        else:
            print("Yoneda module does not exist: ",X," is not in ",A,".")
    
class dg_mod():
    '(Right-) differential graded modules over a dg-category'
    def __init__(self,A,M,diff,act):
        '''Input:
        A is dg_cat
        M is dictionary of gr_vs indexed by A.objects
        diff is dictionary of gr_hom indexed by A.objects
        diff[X] is a gr_hom with source M[X] and target M[X].shift()
        act is dictionary of gr_hom indexed by pairs of A.objects
        act[(X,Y)] is gr_hom with source M[Y]*A([X,Y]) and target M[X]
        As usual, we only include nonzero spaces and maps in our dictionaries.
        '''
        self.cat=A
        self.module=M
        self.diff=diff
        self.act=act
        '''
        Let e_i be a basis for M and f_i a basis for A.
        It's important that the matrices for act:M*A-->M are written with respect
        to the basis ordered e_1*f_1, e_1*f_2,...
        because otherwise we need a shuffle to implement
        M_1*A+M_2*A=(M_1+M_2)*A in the direct sum of modules.
        This is also compatible with Kronecker product of matrices.
        '''
    def mod(self,X):
        '''
        If M=self is a dg-module over A and X is an object of A,
        M.mod(X) returns the graded vector space (gr_vs) underlying the
        dg-module associated to the object X, provided that vector space
        is defined (otherwise it returns the zero vector space).
        '''
        M=self
        A=M.cat
        if X in A.objects:
            if X in M.module:
                return M.module[X]
            else:
                return gvs.graded_vector_space.zero(K)
        else:
            print("Module does not exist: ",X," is not in ",A,".")
    def d(self,X):
        '''
        If M=self is a dg-module over A and X is an object of A,
        M.d(X) returns the differential:

        M.mod(X)-->M.mod(X)[1]  [gr_hom],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 1).
        '''
        M=self
        A=M.cat
        if X in A.objects:
            if X in M.diff:
                return M.diff[X]
            else:
                return zero_gr_hom(1,M.mod(X),M.mod(X))
    def m(self,X,Y):
        '''
        If M=self is a dg-module over A and X,Y are objects of A,
        M.m(X,Y) returns the action map:

        M.mod(Y)*A.hom(X,Y)-->M.mod(X)  [gr_hom],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 0).
        '''
        M=self
        A=M.cat
        if {X,Y}<=A.objects:
            if (X,Y) in M.act:
                return M.act[(X,Y)]
            else:
                return zero_gr_hom(0,M.mod(Y),A.hom(X,Y))
    def __add__(self,other):
        'Direct sum of two dg-modules'
        A=self.cat
        if A!=other.cat:
            print("Can't sum dg-modules over different dg-categories!")
        else:
            N={}
            diff={}
            act={}
            for X in A.objects:
                k=self.mod(X)+other.mod(X)
                l=(self.d(X)).oplus(other.d(X))
                # In each case, we only bother to store nonempty data
                if k.graded_pieces:
                    N[X]=k
                if l.matrix:
                    diff[X]=l
                for Y in A.objects:
                    q=(self.m(X,Y)).oplus(other.m(X,Y))
                    if q.matrix:
                        act[(X,Y)]=q
            return dg_mod(A,N,diff,act)
    def ext(self,X):
         return self.d(X).cohomology()
    def total_ext(self):
        total=gr_vs({})
        for X in self.cat.objects:
            total+=self.d(X).cohomology()
        return total
    def twist(Q):
        A=self.cat
        M=self
        if Q in A.objects:
            for X in A.objects:
                a=M.mod(X)+((M.mod(Q))*A.hom(X,Q)).shift()
                if n.graded_pieces:
                    N[X]=n
        else:
            print("Cannot twist around ",X," as it is not an object of ",A,".")
    @classmethod
    def ext(cls,M,N):
        '''
        Computes the Ext group Ext*(M,N) in the category of dg-modules:
        Technically, we don't need this method yet, as we only want to
        compute Ext groups between the Yoneda modules associated to objects X
        and another module M; this information is tautologically contained
        in the module M and can be accessed using the method M.ext(X)
        '''
        
        
#### A useful source of dg-categories

class dynkin_graph():
    def __init__(self,V,A):
        '''
        Input a directed graph in the form:
        V - a set
        A - a dictionary of sets
        
        Example:
        V={'A','B','C'}
        A={'A':{'B','C'},'B':{'C'}}
        encodes the graph with three vertices A,B,C and arrows
        from A to B and C and from B to C.
        '''
        self.vertices,self.arrows=G,A
    def categorify(self,N):
        'Produce the associated dg-category of dimension N'
        obj=self.vertices
        mor={}
        diff={}
        prod={}
        sph=gr_vs({0:1,N:1})  # The cohomology of S^N
        m_sph={0:np.matrix('1'),N:np.matrix('1 (-1)**N')}
        #The sign comes from Seidel's convention [Seidel, p.8 Eq (1.3)]
        tween1=gr_vs({1:1})   # Rank 1, degree 1
        tween2=gr_vs({N-1:1}) # Rank 1, degree N-1
        i=np.matrix('1')
        for X in obj:
            mor[(X,X)]=sph
            prod[(X,X,X)]=gr_hom(0,mor[(X,X)]*mor[(X,X)],mor[(X,X)],m_sph)
        for X in obj:
            for Y in self.arrows[X]:
                mor[(X,Y)]=tween1
                mor[(Y,X)]=tween2
                # Need to get signs right! Have not checked them.
                prod[(X,X,Y)]=gr_hom(0,mor[(X,Y)]*mor[(X,X)],mor[(X,Y)],{1:i})
                prod[(X,Y,X)]=gr_hom(0,mor[(Y,X)]*mor[(X,Y)],mor[(X,X)],{N:i})
                prod[(X,X,Y)]=gr_hom(0,mor[(X,Y)]*mor[(X,X)],mor[(X,Y)],{N-1:i})
                prod[(X,Y,Y)]=gr_hom(0,mor[(Y,Y)]*mor[(X,Y)],mor[(X,Y)],{1:i})
                prod[(Y,X,Y)]=gr_hom(0,mor[(X,Y)]*mor[(Y,X)],mor[(Y,Y)],{N:i})
                prod[(Y,Y,X)]=gr_hom(0,mor[(Y,X)]*mor[(Y,Y)],mor[(Y,X)],{N-1:i})
        for X in obj:
            for Y in self.arrows[X]:
                for Z in (self.arrows[X].keys() & self.arrows[Y].keys()):
                    # Confused about the m_2 products from triangles - have asked Ailsa
                    prod[(X,Y,Z)]=gr_hom(0,mor[(Y,Z)]*mor[(X,Y)],mor[(X,Z)],{})
        return dg_cat(obj,mor,diff,prod)


    
'''
# List of objects
obj=['A','B']

# Data-type:
# graded vector space
#
# This is stored as an array, e.g.
# V=[0,1,0,1]
# where V[0] is the minimal grading of a nonzero graded piece
# V[1] is the dimension of V^{V[0]}
# V[2] is the dimension of V^{V[0]+1}
# etc
#
# A graded module over the category is then an array of
# graded vector spaces indexed by the objects of the category
# and the hom spaces of the category is an array of graded vector spaces
# indexed by pairs of objects of the category
hom={'A':{'A':[0,1,0,1],'B':[0,1]},'B':{'A':[2,1],'B':[0,1,0,1]}}
O=numpy.matrix([[0]])
I=numpy.matrix([[1]])
r=numpy.zeros(shape=[0,1])
c=numpy.zeros(shape=[1,0])
z=numpy.zeros(shape=[0,0])
mu={'A':
    {'A':
     {'A':{0:{0:I,1:c,2:I},1:{0:r,1:z,2:r},2:{0:I,1:c,2:O}},
      'B':{0:{0:I,1:c,2:0}}},
     'B':
     {'A':{2:{0:I}},
      'B':{0:{0:I},1:{0:r},2:{0:O}}}},
    'B':
    {'A':
     {'A':{0:{2:I},1:{2:r},2:{2:O}},
      'B':{0:{0:I}}},
     'B':
     {'A':{2:{0:I,1:c,2:O}},
      'B':{0:{0:I,1:c,2:I},1:{0:r,1:z,2:r},2:{0:I,1:c,2:O}}}}}
#
# mu['A']['B']['C'][k][l] gives the matrix for multiplying
# hom^k(B,C) (x) hom^l(A,B) --> hom^{k+l}(A,C)
#
# act['A']['B'][k][l] gives the matrix for the module action
# M^k(B) (x) hom^l(A,B) --> M^{k+l}(A)
#
# d['A'][k] gives the matrix for the differential
# M^k(A)-->M^{k+1}(A)
#
#
# This function returns the Yoneda module associated to the object L
#
def yoneda(L):
    Y={}
    d={}
    act={}
    for A in obj:
        Y[A]=hom[A][L]
    for B in obj:
        d[B]=[]
        for k in range(1,len(Y[B])-1):
            d[B].append(numpy.zeros(shape=[Y[B][k],Y[B][k+1]]))
            for A in obj:
                act[A]={}
                act[A][B]={}
                act[A][B][k]={}
                for l in range(1,len(hom[A][B])-1):
                    act[A][B][k][l]=mu[A][B][L][k][l]
    return Y,d,act


# If X is a graded vector space starting from degree K
# and we ask for the dimension of the degree s-K piece
# this could be out of range, in which case this function returns zero
def gr_dim( X,s ):
    q=s-X[0]
    if 1 <= q <= len(X)-1:
        return X[q]
    else:
        return 0

# This function returns the sum of two graded vector spaces
def gr_sum( X,Y ):
    Z=[min(X[0],Y[0])] # what is the minimal grading in the sum?
    maxgr=max(X[0]+len(X)-1,Y[0]+len(Y)-1) # what is the maximal grading?
    for g in range(1,maxgr-Z[0]+1):
        Z.append(gr_dim(X,g+Z[0])+gr_dim(Y,g+Z[0]))
    return Z

# This function returns the tensor product of two graded vector spaces
def gr_tensor( X,Y ):
    Z=[X[0]+Y[0]]
    for s in range(1,len(X)+len(Y)-2):
        Z.append(0)
    for s in range(1,len(X)):
        for t in range(1,len(Y)):
            Z[s+t-1]=Z[s+t-1]+X[s]*Y[t]
    return Z
    
# This function returns the shift X[1] of a graded vector space X
def gr_shift( X ):
    X[0]=X[0]-1
    return X

# This function takes a dg-module M and an object Q, and returns the
# dg-module T_Q M. Recall that the underlying vector space associated
# to the object A is:
#
# (T_Q M)(A) = M(A) (+) (M(Q) (x) hom(A,Q))[1]
#
def twist( M,Q ):
    N={}
    for A in obj:
        N[A]=gr_sum(M[A],gr_shift(gr_tensor(M[Q],hom[A][Q])))
    return N
'''
