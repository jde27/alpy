#!/usr/bin/python

import rings as ri

class GradedVectorSpace:
    '''The class of graded vector spaces V=\bigosum_n V^n.
    
    Usage: V=GradedVectorSpace(base)

    Mutability: V.graded_dim is mutable. All methods are designed
                not to affect the caller and this will work provided
                the entries of V.graded_dim are immutable.

    Attributes:

        V.base [ri.Field]: Field over which V is defined.

        V.graded_dim [dict]:
            { n [int] : dimension of V^n [int] }

    Methods:

        V.gr_dim(n):
            Args: n [int]
            Returns: V.graded_dim[n], or zero
                        if n is not in V.graded_dim.

        V.oplus(W):
            Args: W [GradedVectorSpace]
            Returns [GradedVectorSpace]: The direct sum V (+) W.

        V.otimes(W):
            Args: W [GradedVectorSpace]
            Returns [GradedVectorSpace]: The tensor product V (x) W.

        V.eye():
            Returns [GradedLinearMap]: The identity map V-->V.
                    The name "eye" is a bit silly, but consistent
                    with NumPy.
    '''

    def __init__(self,base):
        self.base=base

    def gr_dim(self,n):
        '''A safe way of accessing the graded dimensions of V.
        This will return zero if nothing is stored.'''
        if n in self.graded_dim:
            return self.graded_dim[n]
        else:
            return 0

    def eye(self):
        '''Returns the identity map.'''
        k=self.base
        identity=GradedLinearMap(0,self,self)
        for n in self.graded_dim:
            identity.level[n]=k.num(np.eye(self.graded_dim[n]))
        return identity
    
    def shift(self,n=1):
        '''Returns a copy of the graded vector space
        shifted in degree by n (n=1 by default).

        Recall that V[n]^d=V^{n+d}.
        '''
        shifted_space=GradedVectorSpace(self.base)
        for n in self.graded_dim:
            shifted_space.graded_dim[n-1]=self.graded_dim[n]
        return shifted_space
        
class GradedLinearMap:
    '''The class of graded linear maps F^n:V^n-->W^{n+d} of degree d.

    Usage: F=GradedLinearMap(degree,source,target)

    Attributes:

        F.degree [int]: The degree of F.

        F.source, F.target [GradedVectorSpaces]:
            The source and target of F.

        F.level [dict]:
            { n [int] : F^n: V^n-->W^{n+d} [np.ndarray] }
                (The entries of F^n should be numbers in the field
                F.source.base).

    Methods:

        F.gr_map(n):
            Args: n [int]
            Returns: The nth graded level F.level[n], or a zero
                    matrix of appropriate size if n is not
                    in F.level.
    
        F.otimes(G):
            Args: G [GradedLinearMap]
            Returns [GradedLinearMap]:
                The tensor (Kronecker) product of F and G.
    
        F.koszulify(m=1):
            Args: m [int] (m=1 by default)
            Returns [GradedLinearMap]:
                The graded linear map K with K^n = (-1)**(m*(n-1))) F^n
    
        F.rejig_1(m), F.rejig_2(m), F.rejig_3(m):
            Args: m [int]
            Returns:##############

    Class Methods:

        GradedLinearMap.block(A,B,C,D):
            Args: A, B, C, D [GradedLinearMaps]
            Returns [GradedLinearMap]:
                The block map with matrix
                    [ A B ]
                    [ C D ]
    '''
    def __init__(self,degree,source,target):
        self.degree=degree
        self.source=source
        self.target=target
        self.level={}
        
    def gr_map(self,n):
        '''A safe way of accessing the graded levels F^n.'''
        k=self.source.base
        if n in self.level:
            # The graded levels are mutable, so
            # we return a copy - one safety feature!
            return self.level[n].copy()
        else:
            # If the graded level is not stored
            # we return zero - another safety feature!
            d=self.degree
            d1,d2=self.source.gr_dim(n),self.target.gr_dim(n+d)
            return k.num(np.zeros(shape=d2,d1))

    def otimes(self,other):
        ###
        
    def koszulify(m=1):
        new_map=GradedLinearMap(self.degree,self.source,self.target)
        for n in self.level:
            new_map.level[n]=self.base.num((-1)**(m*(n-1)))*self.gr_map(n)
            # We use gr_map to access self.level safely
            ### Hopefully np.arrays can cope with this
        return new_ghomo

    def rejig_1(m=1):
        '''This implements the canonical isomorphism
        hom(M,N)[m] --> hom(M,N[m])
        '''
        d,M,N=self.degree,self.source,self,target
        new_map=GradedLinearMap(d-m,M,N.shift(m))
        for n in self.level:
            new_map.level[n]=self.gr_map(n)
            # We use gr_map to access self.level safely
        return new_map
    
    def rejig_2(m=1):
        '''This implements the canonical isomorphism
        hom(M,N)[d]-->hom(M[m],N[m])[d]
        f--------->((-1)**(m(d-1))) f
        '''
        d,M,N=self.degree,self.source,self,target
        new_map=GradedLinearMap(d,M.shift(n),N.shift(n))
        for n in self.levels:
            new_map.level[n]=M.base.num((-1)**(m*(d-1)))*self.gr_maps(n)
            # We use gr_maps to access self.level safely
            ### again, we need to make sure ndarrays can cope with this
        return new_map

    def rejig_3(m=1):
        '''This implements the canonical isomorphism
        hom(M,N)-->hom(M[k],N)[k]
        Remember there are two different ways to implement this
        which differ by a sign - the other is rejig_2(rejig_1(self,-n),n)
        '''
        return rejig_1(rejig_2(self,m),-m)

    @classmethod
    def block(cls,A,B,C,D):
        ###
    
class CochainComplex:
    '''The class of cochain complexes.

    Usage: Z=CochainComplex(cochains,differential)

    Attributes:

        Z.cochains [GradedVectorSpace]: The space of cochains.

        Z.differential [GradedLinearMap]: The differential on cochains.
            This should have:   .degree=1,
                                .source=Z.cochains
                                .target=Z.cochains
    '''
    def __init__(self,cochains,differential):
        self.cochains=cochains
        self.differential=differential
    
    def cohomology(self):
        ###
