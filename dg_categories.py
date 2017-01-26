#!/usr/bin/python

import numpy as np
import rings as ri
import rationals as QQ
from graded_vector_spaces import *

class dg_cat():
    '''The class of differential graded categories over a field k.
    
    Usage:
        D=dg_cat(k,objects,morphisms,diff,prod)
    
    Attributes:
        D.base [ri.field]: field of definition of D
        D.objects [list]: list of objects of D
        D.morphisms [dictionary]:
            { (X,Y) [pair of objects in D.objects] :
                morphism space X-->Y [gvs] }
        diff [dictionary]:
            { (X,Y) [pair of objects in D.obj] :
                differential on hom(X,Y) [ghomo]
                                                of degree 1 }
        prod [dictionary]:
            { (X,Y,Z) [triple of objects in D.obj] :
                product hom(Y,Z) (x) hom(X,Y) --> hom(X,Z)
                            [ghomo] of degree 0 }

    Methods:
        D.hom(X,Y):
            Args: X, Y in D.objects
            Returns: [gvs] D.morphisms[(X,Y)] or zero if this does not exist
        D.d(X,Y):
            Args: X, Y in D.objects
            Returns: [ghomo] D.diff[(X,Y)] or zero if this does not exist
        D.m(X,Y,Z):
            Args: X, Y, Z in D.objects
            Returns: [ghomo] D.prod[(X,Y,Z)] or zero if this does not exist
        D.yoneda(X):
            Args: X in D.objects
            Returns: [dg_mod] the Yoneda module associated to X
        
    '''
    def __init__(self,k,objects,morphisms,diff,prod):
        self.base=k
        self.objects=objects
        self.morphisms=morphisms
        self.diff=diff
        self.prod=prod
        
    def hom(self,X,Y):
        '''
        If A=self is a dg-category and X, Y are objects of A,
        A.hom(X,Y) returns the morphism space between X and Y [gvs],
        provided that space is defined (otherwise it returns
        the zero vector space).
        '''
        A=self
        k=A.base
        if {X,Y}<=A.objects:
            if (X,Y) in A.morphisms:
                return A.morphisms[(X,Y)]
            else:
                return gvs(k)
        else:
            print("Cannot take homs between ",X,\
                  " and ",Y," as they do not live in ",A,".")
    def d(self,X,Y):
        '''
        If A=self is a dg-category and X,Y are objects of A,
        A.d(X,Y) returns the differential:

        A.hom(X,Y)-->A.hom(X,Y)[1]  [ghomo],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 1).
        '''
        A=self
        if {X,Y}<=A.objects:
            if (X,Y) in A.diff:
                return A.diff[(X,Y)]
            else:
                return ghomo(1,A.hom(X,Y),A.hom(X,Y))
            
    def m(self,X,Y,Z):
        '''
        If A=self is a dg-category and X,Y,Z are objects of A,
        A.m(X,Y,Z) returns the multiplication map:

        A.hom(Y,Z)*A.hom(X,Y)-->A.hom(X,Z)  [ghomo],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 0).
        '''
        A=self
        if {X,Y,Z}<=A.objects:
            if (X,Y,Z) in A.prod:
                return A.prod[(X,Y,Z)]
            else:
                return ghomo(0,A.hom(Y,Z)*A.hom(X,Y),A.hom(X,Z))
            
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
    '''(Right-) differential graded modules over a dg-category

    Usage:
        dg_mod(cat,modules,diff,act)

    Attributes:
        M.base [ri.field]: the base field of definition (inherited from
                the ambient category)

        M.cat [dg_cat]: the ambient dg-category over which M is a module.

        M.module [dict]:
            { X [in M.cat.objects] : [gvs] the module M(X) }

        M.diff [dict]:
            { X [in M.cat.objects] : [ghomo] the differential for M(X) }
            diff[X] : M.mod(X) --> M.mod(X).shift()

        M.act [dict]:
            { (X,Y) [both in M.cat.objects] :
                [ghomo] the action M(Y) (x) hom(X,Y) --> M(X) }
            act[(X,Y)] : M.mod(Y) * A.hom(X,Y) --> M.mod(X)

    Methods:
        M.mod(X):
            Args: X [in M.cat.objects]
            Returns: [gvs] M.module[X] if defined and zero otherwise

        M.d(X):
            Args: X [in M.cat.objects]
            Returns: [ghomo] M.diff[X] if defined and zero otherwise

        M.m(X,Y):
            Args: X, Y [in M.cat.objects]
            Returns: [ghomo] M.act[(X,Y)] if defined and zero otherwise

        M+N:
            Args: M, N [dg_mod]
            Returns: The direct sum of M and N.

        M.ext(X):
            Args: X [in M.cat.objects]
            Returns: [gvs] the cohomology of (M(X),d(X))

        M.total_ext():
            Returns: [gvs] the sum of M.ext(X) over X in M.cat.objects

        M.twist(X):
            Args: X [in M.cat.objects]
            Returns: [dg_mod] the twist of M around X.
    '''
    def __init__(self,cat,module,diff,act):
        self.base=A.base
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
        M.mod(X) returns the graded vector space (gvs) underlying the
        dg-module associated to the object X, provided that vector space
        is defined (otherwise it returns the zero vector space).
        '''
        M=self
        A=M.cat
        k=self.base
        if X in A.objects:
            if X in M.module:
                return M.module[X]
            else:
                return gvs(k)
        else:
            print("Module does not exist: ",X," is not in ",A,".")
            
    def d(self,X):
        '''
        If M=self is a dg-module over A and X is an object of A,
        M.d(X) returns the differential:

        M.mod(X)-->M.mod(X)[1]  [ghomo],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 1).
        '''
        M=self
        A=M.cat
        if X in A.objects:
            if X in M.diff:
                return M.diff[X]
            else:
                return ghomo(1,M.mod(X),M.mod(X))
            
    def m(self,X,Y):
        '''
        If M=self is a dg-module over A and X,Y are objects of A,
        M.m(X,Y) returns the action map:

        M.mod(Y)*A.hom(X,Y)-->M.mod(X)  [ghomo],

        provided that map is defined (otherwise it returns
        the zero map with the correct source and target, of degree 0).
        '''
        M=self
        A=M.cat
        if {X,Y}<=A.objects:
            if (X,Y) in M.act:
                return M.act[(X,Y)]
            else:
                return ghomo(0,M.mod(Y),A.hom(X,Y))
            
    def __add__(self,other):
        '''Direct sum of two dg-modules'''
        A=self.cat
        if A!=other.cat:
            print("Can't sum dg-modules over different dg-categories!")
        else:
            N={}
            diff={}
            act={}
            for X in A.objects:
                sum_of_mod=self.mod(X)+other.mod(X)
                sum_of_diff=(self.d(X)).oplus(other.d(X))
                # In each case, we only bother to store nonempty data
                if sum_of_mod.gps:
                    N[X]=sum_of_mod
                if sum_of_diff.gps:
                    diff[X]=sum_of_diff
                for Y in A.objects:
                    sum_of_act=(self.m(X,Y)).oplus(other.m(X,Y))
                    if sum_of_act.gps:
                        act[(X,Y)]=sum_of_act
            return dg_mod(A,N,diff,act)
        
    def ext(self,X):
         return self.d(X).cohomology()
     
    def total_ext(self):
        total=gvs(self.base)
        for X in self.cat.objects:
            total+=self.d(X).cohomology()
        return total

    def twist(Q):
        A=self.cat
        M=self
        twisted_module={}
        twisted_diff={}
        twisted_act={}
        if Q in A.objects:
            for X in A.objects:
                # First construct the twisted module
                S=M.mod(X)
                T=((M.mod(Q))*A.hom(X,Q)).shift()
                twisted_module[X]=S+T
                # Then construct the twisted differential
                d1=M.diff(X)
                B=M.act[(X,Q)].shift()########### SHIFT?!?!
                Z=ghomo(1,S,T)
                d2=(M.diff(Q)).otimes(A.hom(X,Q).Id())##### SHIFT/SIGNS?!s
                twisted_diff[X]=ghomo.block(d1,B,Z,d2)
                #Finally, construct the twisted action
                    
                ###### still to do twisted_act
                return dg_mod(A,twisted_module,\
                              twisted_diff,twisted_act)
        else:
            print("Cannot twist around ",Q,\
                  " as it is not an object of ",A,".")
        
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
        sph=gvs(self.base)
        sph.setgps({0:1,N:1})  # The cohomology of S^N
        ###################THIS NEEDS MASSIVE EDITING:
        ########## NO LONGER USE MATRICES
        ############# AND have to use setgps
        m_sph={0:np.matrix('1'),N:np.matrix('1 (-1)**N')}
        #The sign comes from Seidel's convention [Seidel, p.8 Eq (1.3)]
        tween1=gvs(self.base)
        tween1.setgps({1:1})   # Rank 1, degree 1
        tween2=gvs(self.base)
        tween2.setgps({N-1:1}) # Rank 1, degree N-1
        i=np.matrix('1')
        for X in obj:
            mor[(X,X)]=sph
            prod[(X,X,X)]=ghomo(0,mor[(X,X)]*mor[(X,X)],mor[(X,X)])
            prod[(X,X,X)].setgps(m_sph)
        for X in obj:
            for Y in self.arrows[X]:
                mor[(X,Y)]=tween1
                mor[(Y,X)]=tween2
                # Need to get signs right! Have not checked them.
                prod[(X,X,Y)]=ghomo(0,mor[(X,Y)]*mor[(X,X)],mor[(X,Y)],{1:i})
                prod[(X,Y,X)]=ghomo(0,mor[(Y,X)]*mor[(X,Y)],mor[(X,X)],{N:i})
                prod[(X,X,Y)]=ghomo(0,mor[(X,Y)]*mor[(X,X)],mor[(X,Y)],{N-1:i})
                prod[(X,Y,Y)]=ghomo(0,mor[(Y,Y)]*mor[(X,Y)],mor[(X,Y)],{1:i})
                prod[(Y,X,Y)]=ghomo(0,mor[(X,Y)]*mor[(Y,X)],mor[(Y,Y)],{N:i})
                prod[(Y,Y,X)]=ghomo(0,mor[(Y,X)]*mor[(Y,Y)],mor[(Y,X)],{N-1:i})
        for X in obj:
            for Y in self.arrows[X]:
                for Z in (self.arrows[X].keys() & self.arrows[Y].keys()):
                    # Confused about the m_2 products from triangles - have asked Ailsa
                    prod[(X,Y,Z)]=ghomo(0,mor[(Y,Z)]*mor[(X,Y)],mor[(X,Z)],{})
        return dg_cat(obj,mor,diff,prod)
