#!/usr/bin/python

import numpy as np
import rings as ri
import rationals as rat

'''In this module, we define the following categories:

    GradedVectorSpace
    GradedLinearMap
    CochainComplex

and the following functions:

    find_nonzero(array) [used by rank_nullity()]
    rank_nullity(array) [used by GradedLinearMap.ker_im()]
'''

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

        V.shift(m=1):
            Args: m [int] (m=1 by default)
            Returns [GradedVectorSpace]:
                The shift of V by m, where V[m]^d = V^{m+d}
    
        print(V): Prints V in the form
                    { n : dim(V^n) }
    '''

    def __init__(self,base):
        self.base=base
        self.graded_dim={}
        
    def gr_dim(self,n):
        '''A safe way of accessing the graded dimensions of V.
        This will return zero if nothing is stored.'''
        if n in self.graded_dim:
            return self.graded_dim[n]
        else:
            return 0

    def oplus(self,other):
        '''Returns the direct sum of two graded vector spaces.'''
        (V,W)=(self,other)
        K=self.base
        D=GradedVectorSpace(K)
        D.graded_dim=V.graded_dim.copy()
        # Use of .copy() avoids shared references
        for n in W.graded_dim:
            if n in D.graded_dim:
                D.graded_dim[n]+=W.gr_dim(n)
            else:
                D.graded_dim[n]=W.gr_dim(n)
        return D

    def __add__(self,other):
        return self.oplus(other)

    def __eq__(self,other):
        gr_dim_test=[self.gr_dim(n)
                     for n in (self.graded_dim.keys()|other.graded_dim.keys())
                     if self.gr_dim(n)!=other.gr_dim(n)]
        if not gr_dim_test:
            return True
        else:
            return False
    
    def otimes(self,other):
        '''Returns the tensor product of two graded vector spaces.'''
        (V,W)=(self,other)
        K=self.base
        T=GradedVectorSpace(K)
        for p in V.graded_dim:
            for q in W.graded_dim:
                if p+q in T.graded_dim:
                    T.graded_dim[p+q]+=(V.gr_dim(p))*(W.gr_dim(q))
                else:
                    T.graded_dim[p+q]=(V.gr_dim(p))*(W.gr_dim(q))
        return T
    
    def eye(self):
        '''Returns the identity map.'''
        K=self.base
        identity=GradedLinearMap(0,self,self)
        for n in self.graded_dim:
            identity.graded_map[n]=K.num(np.eye(self.graded_dim[n],dtype=int))
        return identity
    
    def shift(self,m=1):
        '''Returns a copy of the graded vector space
        shifted in degree by m (m=1 by default).

        Recall that V[m]^d=V^{m+d}.
        '''
        shifted_space=GradedVectorSpace(self.base)
        for n in self.graded_dim:
            shifted_space.graded_dim[n-m]=self.graded_dim[n]
        return shifted_space
        
    def __str__(self):
        return "%s" % str(self.graded_dim)

class GradedLinearMap:
    '''The class of graded linear maps F^n:V^n-->W^{n+d} of degree d.

    Usage: F=GradedLinearMap(degree,source,target)

    Attributes:

        F.degree [int]: The degree of F.

        F.source, F.target [GradedVectorSpaces]:
            The source and target of F.

        F.graded_map [dict]:
            { n [int] : F^n: V^n-->W^{n+d} [np.ndarray] }
                (The entries of F^n should be numbers in the field
                F.source.base).

    Methods:

        F.gr_map(n):
            Args: n [int]
            Returns [GradedLinearMap]: The nth graded graded map
                    F.graded_map[n], or a zero matrix of appropriate size
                    if n is not in F.graded_map.

        F+G:
            Returns: The sum of linear maps V-->W.
    
        F*G:
            Returns: The composition of linear maps.
    
        F.oplus(G):
            Args: G [GradedLinearMap]
            Returns [GradedLinearMap]: The direct sum of F and G.

        F.otimes(G):
            Args: G [GradedLinearMap]
            Returns [GradedLinearMap]:
                The tensor (Kronecker) product of F and G.
    
        F.koszulify(m=1):
            Args: m [int] (m=1 by default)
            Returns [GradedLinearMap]:
                The graded linear map G with G^n = (-1)**(m*(n-1))) F^n
    
        F.rejig_1(m), F.rejig_2(m), F.rejig_3(m):
            Args: m [int]
            Returns [GradedLinearMap]: Implements the
                canonical isomorphisms:
                    hom(A,B)[k]-->hom(A,B[k]) (rejig_1)
                    hom(M,N)-->hom(M[m],N[m]) (rejig_2)
                    hom(M,N)-->hom(M[m],N)[m] (rejig_3)
    
        F.ker_im()
            Returns a tuple (V,W) [GradedVectorSpace]:
                The kernel/image of F.

    Class Methods:

        GradedLinearMap.block(A,B,C,D):
            Args: A, B, C, D [GradedLinearMaps]
            Returns [GradedLinearMap]:
                The block map with matrix
                    [ A B ]
                    [ C D ]

    Static methods:
    
        GradedLinearMap.shuffle_map(A,B,C):
            Args: A, B, C [GradedVectorSpaces]
            Returns [GradedLinearMap]:
                The permutation matrix which realises the isomorphism
                ( A (x) B ) (x) C ----> A (x) ( B (x) C )
            of graded triple tensor products.
    '''
    def __init__(self,degree,source,target):
        self.degree=degree
        self.source=source
        self.target=target
        self.graded_map={}
        
    def gr_map(self,n):
        '''A safe way of accessing the graded maps F^n.'''
        K=self.source.base
        if n in self.graded_map:
            # The graded maps are mutable, so
            # we return a copy - one safety feature!
            return self.graded_map[n].copy()
        else:
            # If the graded map is not stored
            # we return zero - another safety feature!
            d=self.degree
            (d1,d2)=(self.source.gr_dim(n),self.target.gr_dim(n+d))
            return K.num(np.zeros(shape=(d2,d1),dtype=int))

    def __add__(self,other):
        '''Returns the usual sum of two graded linear maps V-->W.'''
        (F,G)=(self,other)
        (d1,d2)=(self.degree,other.degree)
        V1=F.source
        W1=F.target
        V2=G.source
        W2=G.target
        H=GradedLinearMap(d1,V1,W1)
        if V1==V2 and W1==W2 and d1==d2:
            for n in (F.graded_map or G.graded_map):
                H.graded_map[n]=F.gr_map(n)+G.gr_map(n)
            return H
        else:
            print("Cannot add linear maps of different degree "\
                  "or between different spaces.")
            return None
        
    def __mul__(self,other):
        '''Returns the composition of two graded linear maps.

        The ordering is as follows:
        G:U-->V, F:V-->W gives F*G=H:U-->W
        '''
        (F,G)=(self,other)
        (U,W)=(G.source,F.target)
        (d1,d2)=(G.degree,F.degree)
        H=GradedLinearMap(d1+d2,U,W)
        for n in G.graded_map:
            H.graded_map[n]=F.gr_map(n+d1).dot(G.gr_map(n))
        return H
        
    def oplus(self,other):
        '''Returns the direct sum of two graded linear maps.'''
        (V1,V2)=(self.source,other.source)
        (W1,W2)=(self.target,other.target)
        d=self.degree
        B=GradedLinearMap(d,V2,W1)
        C=GradedLinearMap(d,V1,W2)
        return GradedLinearMap.block(self,B,C,other)

    def otimes(self,other):
        '''Returns the tensor (Kronecker) product of linear maps.'''
        K=self.source.base
        (F,G)=(self,other)
        (d1,d2)=(F.degree,G.degree)
        (V1,V2)=(F.source,G.source)
        (W1,W2)=(F.target,G.target)
        H=GradedLinearMap(d1+d2,V1.otimes(V2),W1.otimes(W2))
        gradings_H={p+q for p in F.graded_map for q in G.graded_map}
        for n in gradings_H:
            S1={p for p in V1.graded_dim}
            S2={n-p for p in V2.graded_dim}
            S3={p-d1 for p in W1.graded_dim}
            S4={n+d2-p for p in W2.graded_dim}
            for p in sorted(S1|S2|S3|S4):
                q=n-p
                A=F.gr_map(p)
                B=G.gr_map(q)
                dom=(V1.gr_dim(p))*(V2.gr_dim(q))
                tar=(W1.gr_dim(p+d1))*(W2.gr_dim(q+d2))
                # np.kron can't cope with zero-dimensional
                # vector spaces, but we need them, e.g. to handle
                # (V-->0) (x) (W-->W')
                if dom!=0 and tar!=0:
                    C=K.num(np.asarray(np.kron(A,B)))
                else:
                    C=K.num(np.zeros(shape=(tar,dom),dtype=int))
                if n in H.graded_map:
                    b,a=H.graded_map[p+q].shape
                    d,c=C.shape
                    Z1=K.num(np.zeros(shape=(b,c),dtype=int))
                    Z2=K.num(np.zeros(shape=(d,a),dtype=int))
                    block_diag=np.asarray(np.bmat([[H.gr_map(p+q),Z1],[Z2,C]]))
                    H.graded_map[p+q]=K.num(block_diag)
                else:
                    H.graded_map[p+q]=C
        return H
                
    def koszulify(self,m=1):
        '''Implement a Koszul sign-change.'''
        F=self
        K=F.source.base
        new_map=GradedLinearMap(F.degree,F.source,F.target)
        for n in F.graded_map:
            if (m*(n-1))%2==1:
                sigma=-1
            else:
                sigma=1
            new_map.graded_map[n]=K.num(sigma)*F.gr_map(n)
            # We use gr_map to access self.graded_map safely
        return new_map

    def rejig_1(self,m=1):
        '''This implements the canonical isomorphism
        hom(M,N)[m] --> hom(M,N[m])
        '''
        (d,M,N)=(self.degree,self.source,self.target)
        new_map=GradedLinearMap(d-m,M,N.shift(m))
        for n in self.graded_map:
            new_map.graded_map[n]=self.gr_map(n)
            # We use gr_map to access self.graded_map safely
        return new_map
    
    def rejig_2(self,m=1):
        '''This implements the canonical isomorphism
        hom(M,N)[d]-->hom(M[m],N[m])[d]
        f--------->((-1)**(m(d-1))) f
        '''
        (d,M,N)=(self.degree,self.source,self.target)
        k=M.base
        new_map=GradedLinearMap(d,M.shift(m),N.shift(m))
        for n in self.graded_map:
            if (m*(d-1))%2==1:
                sigma=-1
            else:
                sigma=1
            new_map.graded_map[n-m]=k.num(sigma)*self.gr_map(n)
            # We use gr_maps to access self.graded_map safely
        return new_map

    def rejig_3(self,m=1):
        '''This implements the canonical isomorphism
        hom(M,N)-->hom(M[m],N)[m]
        Remember there are two different ways to implement this
        which differ by a sign - the other is rejig_2(rejig_1(self,-n),n)
        '''
        return (self.rejig_2(m)).rejig_1(-m)

    def ker_im(self):
        '''Find the kernel and image of graded linear map.

        Here, we will store ker(f^n) in kernel.graded_dim[n]
        and im(f^{n-1}) in image.graded_dim[n] (so that
        cohomology in degree n is given by ker[n]/im[n]).
        '''
        F=self
        K=self.source.base
        d=F.degree
        kernel=GradedVectorSpace(K)
        image=GradedVectorSpace(K)
        for n in F.source.graded_dim:
            M=F.gr_map(n)
            rank,nullity=rank_nullity(M)
            kernel.graded_dim[n]=nullity
            image.graded_dim[n+d]=rank
        return kernel,image

    def verify(self):
        '''This is just a debugging tool to check that the map
        is internally consistent (its domain and range agree with
        the shape of its matrix).'''
        F=self
        d=F.degree
        for n in (F.graded_map or F.source.graded_dim):
            print("deg",n,"tar,src",F.gr_map(n).shape,
                  F.target.gr_dim(n+d),F.source.gr_dim(n))

    def __eq__(self,other):
        degree_test=(self.degree==other.degree)
        gr_map_test=[((self.gr_map(n)==other.gr_map(n)).all())
                     for n in (self.graded_map.keys()|other.graded_map.keys())]
        if degree_test and all(gr_map_test):
            return True
        else:
            return False
        
            
    def __str__(self):
        return "%s" % str(self.graded_map)

    def display(self):
        for n in self.graded_map:
            print("Degree ",n,":\n[",end='')
            for x in self.gr_map(n):
                print()
                for y in x:
                    print(y,end=' ')
            print("\n]\n")
            
    @classmethod
    def block(cls,A,B,C,D):
        '''Returns a block matrix from four compatible
        graded homomorphisms.'''
        (d,d_other)=(A.degree,(B.degree,C.degree,D.degree))
        degree_test=(x==d for x in d_other)
        (V1,W1)=(A.source,A.target)
        (V2,W2)=(D.source,D.target)
        (V1_test,W2_test)=(C.source,C.target)
        (V2_test,W1_test)=(B.source,B.target)
        space_test=(V1==V1_test and V2==V2_test and \
                     W1==W1_test and W2==W2_test)
        all_keys=(A.graded_map.keys()|B.graded_map.keys()
                  |C.graded_map.keys()|D.graded_map.keys())
        if degree_test and space_test:
            block_map=GradedLinearMap(d,A.source+D.source,A.target+D.target)
            for n in all_keys:
                H=np.asarray(np.bmat([[A.gr_map(n),B.gr_map(n)],
                                      [C.gr_map(n),D.gr_map(n)]]))
                block_map.graded_map[n]=H
            return block_map
        else:
            if not space_test:
                print("Cannot take direct sum of homomorphisms: sources "\
                      "and targets are incompatible.")
            else:
                print("Cannot take direct sum of homomorphisms: degrees "\
                      "are incompatible.")
            return None

    @staticmethod
    def shuffle_map(A,B,C):
        '''This function returns the natural graded linear map

        ( A (x) B ) (x) C ----> A (x) ( B (x) C )

        The main point of this is that we have implicitly picked a basis
        for all our vector spaces and the basis order is different
        for these two triple tensor products.
        '''
        def F(x,y):
            return (A.gr_dim(m1+x))*(B.gr_dim(m2+y))*(C.gr_dim(m3+N-x-y))

        def insert_submatrix(N1,N2,top,left):
            '''Inserts the submatrix N1 into N2 at position (top,left).'''
            height,width=N1.shape
            N2[top:top+height,left:left+width]=N1
            return N2
            
        K=A.base
        D=(A.otimes(B)).otimes(C)
        E=A.otimes(B.otimes(C))
        m1=min(A.graded_dim.keys())
        m2=min(B.graded_dim.keys())
        m3=min(C.graded_dim.keys())
        shuffler=GradedLinearMap(0,D,E)
        for n in D.graded_dim:
            delta=D.gr_dim(n)
            M=np.zeros(shape=(delta,delta),dtype=int)
            N=n-m1-m2-m3
            x,y,alpha=0,0,0
            Alpha={}
            while y<=N:
                Alpha[(x,y)]=alpha
                alpha+=F(x,y)
                if y>0:
                    x,y=x+1,y-1
                else:
                    x,y=0,x+1
            x,y,beta=0,0,0
            while x<N:
                M=insert_submatrix(np.eye(F(x,y),dtype=int),M,Alpha[(x,y)],beta)
                beta+=F(x,y)
                if y < N-x:
                    x,y=x,y+1
                else:
                    x,y=x+1,0
            M=insert_submatrix(np.eye(F(x,y),dtype=int),M,Alpha[(x,y)],beta)
            shuffler.graded_map[n]=K.num(M)
        return shuffler
    
        
class CochainComplex:
    '''The class of cochain complexes.

    Usage: Z=CochainComplex(cochains,differential)

    Attributes:

        Z.cochains [GradedVectorSpace]: The space of cochains.

        Z.differential [GradedLinearMap]: The differential on cochains.
            This should have:   .degree=1,
                                .source=Z.cochains
                                .target=Z.cochains
    
    Methods:
        Z.cohomology():
            Returns [GradedVectorSpace]:
                The cohomology (ker(d)/im(d)) of Z.
 
    '''
    def __init__(self,cochains,differential):
        self.cochains=cochains
        self.differential=differential
    
    def verify(self):
        d=self.differential
        Z=self.cochains
        zero_map=GradedLinearMap(2,Z,Z)
        if d*d!=zero_map:
            print("Not a cochain complex")
        else:
            print("Yes a cochain complex")
        
    def cohomology(self):
        '''Computes cohomology of a cochain complex.'''
        d=self.differential
        K=self.cochains.base
        kernel,image=d.ker_im()
        coh=GradedVectorSpace(K)
        for n in kernel.graded_dim:
            coh.graded_dim[n]=kernel.gr_dim(n)-image.gr_dim(n)
            if coh.graded_dim[n]<0:
                print("Cohomology seems to be negative dimensional...")
        return coh


def find_nonzero(A):
    '''A function which finds the index of the first nonzero
    entry in an array

    Input: A [array]
    
    Output: [bool], [int], [entry of A]
    
    If there is a nonzero entry in position x,
    find_nonzero returns True, x, A[x]
    otherwise, it returns False, 0, 0
    '''
    ctr=0
    while ctr<A.size:
        if A[ctr]!=0:
            return True,ctr,A[ctr]
        else:
            ctr+=1
    return False,0,0
    
def rank_nullity(A):
    '''Computes the rank and nullity of a 2-d np.array A
    via Gaussian elimination over the field of definition'''
    # This algorithm changes the input matrix drastically,
    # so we start by warding against mutability backwash:
    B=A.copy()
    # We start with a guess of the rank and nullity,
    # setting them to be the biggest they could possibly be.
    rows,cols=B.shape
    rank,nullity=rows,cols
    for n in range(0,rows):
        # We look for the first nonzero element in row n
        found,s,t=find_nonzero(B[n])
        if found:
            # If there is a nonzero element we know the
            # nullity is at least 1 less than our previous guess.
            nullity-=1
            # We now clear column number s by subtracting a suitable
            # multiple of row n from all subsequent rows.
            for x in range(n+1,rows):
                c=B[x,s]/t
                for y in range(0,cols):
                    B[x,y]-=c*B[n,y]
        else:
            # If we find an empty row, we know the rank is
            # at least 1 less than our previous guess.
            rank-=1
    return rank,nullity

def sph(K,N):
    '''Returns the graded vector space H^*(S^N;K).'''
    V=GradedVectorSpace(K)
    V.graded_dim={0:1,N:1}
    return V

def sph_op(V):
    '''Returns the ring operations on V=H^*(S^N;K).'''
    F=GradedLinearMap(0,V.otimes(V),V)
    K=V.base
    N=max(V.graded_dim.keys())
    i=K.num(np.array([[1]]))
    if N%2==1:
        sigma=-1
    else:
        sigma=1
    j=K.num(np.array([[1,sigma]]))
    F.graded_map={0:i,2:j}
    return F

def pt(K,N):
    '''Returns the graded vector space K in degree N.'''
    V=GradedVectorSpace(K)
    V.graded_dim={N:1}
    return V
    
def simple(V,W,N):
    '''Returns the identity V-->W where V and W are both a single copy of K
    living in degree N.'''
    F=GradedLinearMap(0,V,W)
    K=V.base
    i=K.num(np.array([[1]]))
    F.graded_map={N:i}
    return F
    
